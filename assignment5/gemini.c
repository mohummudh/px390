#include <stdlib.h>
#include <stdio.h>
#include <mkl_lapacke.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

/* Structure for band matrix representation */
struct band_mat {
    long ncol;      /* Number of columns in band matrix */
    long nbrows;      /* Number of rows (bands in original matrix) */
    long nbands_up;   /* Number of bands above diagonal */
    long nbands_low;  /* Number of bands below diagonal */
    double *array;    /* Storage for the matrix in banded format */
    /* Internal temporary storage for solving inverse problem */
    long nbrows_inv;  /* Number of rows of inverse matrix */
    double *array_inv; /* Store the inverse if this is generated */
    int *ipiv;        /* Additional inverse information */
};

/* Define a new type band_mat */
typedef struct band_mat band_mat;

/* Initialize a band matrix of a certain size, allocate memory,
   and set the parameters. */
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
    bmat->nbrows = nbands_lower + nbands_upper + 1;
    bmat->ncol = n_columns;
    bmat->nbands_up = nbands_upper;
    bmat->nbands_low = nbands_lower;
    bmat->array = (double *)malloc(sizeof(double) * bmat->nbrows * bmat->ncol);
    bmat->nbrows_inv = bmat->nbands_up * 2 + bmat->nbands_low + 1;
    bmat->array_inv = (double *)malloc(sizeof(double) * (bmat->nbrows + bmat->nbands_low) * bmat->ncol);
    bmat->ipiv = (int *)malloc(sizeof(int) * bmat->ncol);
    if (bmat->array == NULL || bmat->array_inv == NULL) {
        return 0;
    }
    /* Initialize array to zero */
    long i;
    for (i = 0; i < bmat->nbrows * bmat->ncol; i++) {
        bmat->array[i] = 0.0;
    }
    return 1;
}

/* Finalize function: should free memory as required */
void finalise_band_mat(band_mat *bmat) {
    free(bmat->array);
    free(bmat->array_inv);
    free(bmat->ipiv);
}

/* Get a pointer to a location in the band matrix, using
   the row and column indexes of the full matrix. */
double *getp(band_mat *bmat, long row, long column) {
    int bandno = bmat->nbands_up + row - column;
    if (row < 0 || column < 0 || row >= bmat->ncol || column >= bmat->ncol) {
        printf("Indexes out of bounds in getp: %ld %ld %ld \n", row, column, bmat->ncol);
        exit(1);
    }
    return &bmat->array[bmat->nbrows * column + bandno];
}

/* Return the value of a location in the band matrix, using
   the row and column indexes of the full matrix. */
double getv(band_mat *bmat, long row, long column) {
    return *getp(bmat, row, column);
}

/* Set an element of a band matrix to a desired value based on the pointer
   to a location in the band matrix, using the row and column indexes
   of the full matrix. */
double setv(band_mat *bmat, long row, long column, double val) {
    *getp(bmat, row, column) = val;
    return val;
}

/* Solve the equation Ax = b for a matrix a stored in band format
   and x and b real arrays */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
    /* Copy bmat array into the temporary store */
    int i, bandno;
    for (i = 0; i < bmat->ncol; i++) {
        for (bandno = 0; bandno < bmat->nbrows; bandno++) {
            bmat->array_inv[bmat->nbrows_inv * i + (bandno + bmat->nbands_low)] =
                bmat->array[bmat->nbrows * i + bandno];
        }
        x[i] = b[i];
    }

    long nrhs = 1;
    long ldab = bmat->nbands_low * 2 + bmat->nbands_up + 1;
    int info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low,
                             bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
    return info;
}

/* Print the matrix values in band format*/
int printmat(band_mat *bmat) {
    long i, j;
    for (i = 0; i < bmat->ncol; i++) {
        for (j = 0; j < bmat->nbrows; j++) {
            printf("%ld %ld %g \n", i, j, bmat->array[bmat->nbrows * i + j]);
        }
    }
    return 0;
}

/* Print the full matrix values*/
int printfullmat(band_mat *bmat) {
    int i, j;
    for (i = 0; i < bmat->ncol; i++) {
        for (j = 0; j < bmat->ncol; j++) {
            // printf("%ld %ld %g \n",i,j,getv(bmat, i, j));
            printf("%g\t", getv(bmat, i, j));
        }
        printf("\n");
    }
    return 0;
}

/* Swaps two double pointers */
void MemSwap(double **a, double **b) {
    double *temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

/* Prints the current set of values to output file. On success, will return 0, otherwise it will return 1 */
int printVector(FILE *file, double **vector, int Na, int **pos, int Ny, double time) {
    if (file == NULL) {
        printf("Error opening output file\n");
        return 1;
    } else {
        for (int j = 0; j < Na; j++) {
            fprintf(file, "%.6f, %lu, %lu, %.6f\n", time, (*pos)[j] / Ny, (*pos)[j] % Ny,
                    (*vector)[j]);
        }
        fprintf(file, "\n");
        return 0;
    }
}

// Function to calculate reaction term f(u)
double reaction_term(int i, int lambda, int Nx, int Ny, double *u) {
    // Implement the reaction term f(u) = lambda*u - u^3 here
    return lambda * u[i] - pow(u[i], 3.0);
}

// Function to calculate diffusion coefficient D(i, Nx, Ny)
double diffusion_coefficient(int i, int Nx, int Ny) {
    return 1.0;
}

int main() {
    // Iteration variables
    int i, j;

    // The parameters input file, 'input.txt' will contain the following values:
    // • Nx: number of x grid points.
    // • Ny: number of y grid points.
    // • Na: number of active grid cells.
    // • Lx: Length of domain G in x direction.
    // • Ly: Length of domain G in y direction.
    // • tf : final time for time evolution.
    // • lambda: λ parameter in equations.
    // • tD: diagnostic timestep.

    // Initialize and extract values from input.txt
    int Nx, Ny, Na;
    double Lx, Ly, tf, lambda, tD;
    double dt; // Set the value of dt

    FILE *input;
    input = fopen("input.txt", "r");
    if (input == NULL) {
        printf("Error opening input file\n");
        exit(1);
    } else {
        if (fscanf(input, "%d %d %d %lf %lf %lf %lf %lf", &Nx, &Ny, &Na, &Lx, &Ly, &tf, &lambda, &tD) != 8) {
            fprintf(stderr, "Error reading values from input file\n");
            fclose(input);
            exit(2);
        }

        printf("Nx = %d\n", Nx);
        printf("Ny = %d\n", Ny);
        printf("Na = %d\n", Na);
        printf("Lx = %lf\n", Lx);
        printf("Ly = %lf\n", Ly);
        printf("tf = %lf\n", tf);
        printf("lambda = %lf\n", lambda);
        printf("tD = %lf\n", tD);
    }
    fclose(input);

     // Calculate grid parameters
    double dx = Lx / (Nx - 1);
    double dy = Ly / (Ny - 1);

        // Open coefficients file and check for errors
    FILE *coefficients = fopen("coefficients.txt", "r");
    if (coefficients == NULL) {
        fprintf(stderr, "Error opening coefficients file\n");
        exit(1);
    }

    // Initialize variables for reading and max calculation
    int x, y;
    double max = -INFINITY; // Initialize max to negative infinity
    double coeff;
    double *vector = malloc(sizeof(double)*Na);

    //This variable stores the gridIndex (x,y) of the nth active gridPoint in the form y*Ny+x
    int *gridIndex = malloc(sizeof(double)*Na);
    //This variable stores the Vector index or -1 otherwise
    int *vectorIndex = malloc(sizeof(double)*Nx*Ny);


    // Read file, extract values, update indices, and find max
    for (i = 0; i < Na; i++) {
        if (fscanf(coefficients, "%d %d %lf", &x, &y, &coeff) == 3) {
            vector[i] = coeff;
            gridIndex[i] = y + x * Ny;
            vectorIndex[Ny * x + y] = i;
            max = fmax(max, coeff); // Update max if coefficient is larger
        } else {
            fprintf(stderr, "Error reading values from coefficients file\n");
            fclose(coefficients);
            free(vector);
            free(gridIndex);
            free(vectorIndex);
            exit(2);
        }
    }

    fclose(coefficients);

    // Arrays for solution, reaction term, and diffusion term
    double *u = (double *)malloc(Na * sizeof(double));
    double *f = (double *)malloc(Na * sizeof(double));
    double *D = (double *)malloc(Na * sizeof(double));

    // Initialize solution, reaction term, and diffusion term arrays
    for (i = 0; i < Na; i++) {
        u[i] = vector[i];
        f[i] = reaction_term(i, lambda, Nx, Ny, u);
        D[i] = diffusion_coefficient(i, Nx, Ny);
    }

    // Calculate spatial and non-linear term bounds (using updated max)
    double B1 = ((0.5 * dx * dx * dy * dy) / (dx * dx + dy * dy));
    double B2;
    if (max > sqrt(lambda) && max > 0) {
        B2 = 1 / (max * (sqrt(lambda) + max));
    } else {
        B2 = 1 / (2 * lambda);
    }

    double min = fmin(B1, B2);

    printf("Diagnostic timestep  = %lf\n", tD);
    printf("spatial bound     = %lf\n", B1);
    printf("max bound       = %lf\n", 1 / (max * (sqrt(lambda) + max)));
    printf("non-linear term bound = %lf\n", 1 / (2 * lambda));

    // Set dt based on the minimum timestep
    dt = min;

    // Initialize band matrix for storing the coefficients
    band_mat A;
    if (!init_band_mat(&A, 1, 1, Na)) {
        printf("Error initializing band matrix\n");
        exit(1);
    }

    double *b;
    b = malloc(sizeof(double) * Na);

    // Set up the band matrix coefficients based on diffusion, reaction, and boundary conditions
    for (i = 0; i < Na; i++) {
        // Set diagonal element
        setv(&A, i, i, 2 * lambda / dt + D[i] / (dx * dx) + D[i] / (dy * dy));

        // Set upper diagonal element (if active neighbor exists)
        if (i + 1 < Na && ((gridIndex[i + 1] / Ny) == (gridIndex[i] / Ny))) {
            setv(&A, i, i + 1, -lambda / dt - D[i] / (2 * dx * dx));
        }

        // Set lower diagonal element (if active neighbor exists)
        if (i > 0 && ((gridIndex[i - 1] / Ny) == (gridIndex[i] / Ny))) {
            setv(&A, i, i - 1, -lambda / dt - D[i] / (2 * dx * dx));
        }

        // Set right-hand side vector based on reaction term
        b[i] = u[i] + f[i] * dt;
    }

    // Boundary conditions (assuming Neumann for inactive regions)
    for (i = 0; i < Na; i++) {
        // Check for vertical boundaries (constant x)
        if ((gridIndex[i] % Ny) == 0 || (gridIndex[i] % Ny) == (Ny - 1)) {
            setv(&A, i, i, setv(&A, i, i, 0.0)); // Set diagonal to zero for Neumann
        }
    }

    // Open output file in append mode
    FILE *output = fopen("output.txt", "a");

    // Error handling for file opening
    if (output == NULL) {
        fprintf(stderr, "Error opening output file\n");
        finalise_band_mat(&A);
        free(gridIndex);
        free(vector);
        exit(3);
    }

    // Iteration step
    int Iteration = floor(tf / tD);

    // Allocate memory for nextVector
    double *nextVector = malloc(sizeof(double) * Na);
    if (nextVector == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        fclose(output);
        finalise_band_mat(&A);
        free(gridIndex);
        free(vector);
        exit(3);
    }

    // Print initial solution (assuming u is already initialized)
    if (printVector(output, &u, Na, &gridIndex, Ny, 0.0) == 1) {
        fprintf(stderr, "Error printing values to output file\n");
        fclose(output);
        finalise_band_mat(&A);
        free(gridIndex);
        free(vector);
        free(nextVector);
        exit(3);
    }

    // Time loop for solving the equation
    double t = 0.0;
    while (t < tf) {
        // Solve Ax = b for u
        solve_Ax_eq_b(&A, u, b);

        // Update solution and time
        t += dt;
        for (i = 0; i < Na; i++) {
            u[i] += dt * f[i];
        }

        // Print solution only at diagnostic timestep
        if (fabs(t - fmod(t, tD)) < 1e-6) {
            if (printVector(output, &u, Na, &gridIndex, Ny, t) != 0) {
                fprintf(stderr, "Error printing values to output file\n");
                fclose(output);
                finalise_band_mat(&A);
                free(gridIndex);
                free(vector);
                free(nextVector);
                exit(3);
            }
        }
    }

    // Close output file
    fclose(output);

    // Free allocated memory
    free(nextVector);
    free(gridIndex);
    free(vector);
    // Free memory
    free(u);
    free(f);
    free(D);
    free(vector);
    finalise_band_mat(&A);

    return 0;
}