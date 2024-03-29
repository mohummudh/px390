#include <stdlib.h>
#include <stdio.h>
#include <mkl_lapacke.h>
#include <stdbool.h>
#include <math.h>


/* Structure representing a band matrix.*/
struct band_mat{
    long ncol;        /* Number of columns in band matrix            */
    long nbrows;      /* Number of rows (bands in original matrix)   */
    long nbands_up;   /* Number of bands above diagonal           */
    long nbands_low;  /* Number of bands below diagonal           */
    double *array;    /* Storage for the matrix in banded format  */

    /* Internal temporary storage for solving inverse problem */

    long nbrows_inv;  /* Number of rows of inverse matrix   */
    double *array_inv;/* Store the inverse if this is generated */
    int *ipiv;        /* Additional inverse information         */
};

typedef struct band_mat band_mat;

/* Initialises a band matrix structure.
 *
 * @param bmat Pointer to the band matrix structure to be initialised.
 * @param nbands_lower Number of bands below the diagonal (must be non-negative).
 * @param nbands_upper Number of bands above the diagonal (must be non-negative).
 * @param n_columns Number of columns in the matrix (must be positive).
 *
 * @return 1 on success, 0 on failure (memory allocation error).
 */
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {

    // Check for valid input parameters
    if (nbands_lower < 0 || nbands_upper < 0 || n_columns <= 0) {
    return 0; // Indicate error with return value
    }

    // Set band matrix parameters
    bmat->nbrows = nbands_lower + nbands_upper + 1;
    bmat->ncol = n_columns;
    bmat->nbands_up = nbands_upper;
    bmat->nbands_low = nbands_lower;
    bmat->nbrows_inv = bmat->nbands_up * 2 + bmat->nbands_low + 1;

    // Allocate memory for arrays
    bmat->array = (double *)calloc(bmat->nbrows * bmat->ncol, sizeof(double));
    bmat->array_inv = (double *)calloc((bmat->nbrows + bmat->nbands_low) * bmat->ncol, sizeof(double));
    bmat->ipiv = (int *)malloc(sizeof(int) * bmat->ncol);

    // Check for allocation failure
    if (!bmat->array || !bmat->array_inv || !bmat->ipiv) {
        free(bmat->array); // Safe to call free on NULL
        free(bmat->array_inv);
        free(bmat->ipiv);
        return 0;
    }
    return 1;
}

/* Frees the memory allocated for the band matrix structure.
 * @param bmat Pointer to the band matrix structure to be finalised.
 */
void finalise_band_mat(band_mat *bmat) {
    // Free memory for all allocated arrays
    free(bmat->array);
    free(bmat->array_inv);
    free(bmat->ipiv);

    // Set pointers to NULL to avoid dangling pointers
    bmat->array = NULL;
    bmat->array_inv = NULL;
    bmat->ipiv = NULL;
}

/* Gets a pointer to a location in the band matrix based on full matrix indices.
 *
 * @param bmat Pointer to the band matrix structure.
 * @param row Row index in the full matrix (0-based).
 * @param column Column index in the full matrix (0-based).
 *
 * @return Pointer to the element in the band matrix, or NULL on error.
 */
double *getp(band_mat *bmat, long row, long column) {
    int bandno = bmat->nbands_up + row - column;
    if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
        printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
        exit(1);
    }
    return &bmat->array[bmat->nbrows*column + bandno];
}

/* Retrun the value of a location in the band matrix, using the row and column indexes of the full matrix. */
double getv(band_mat *bmat, long row, long column) {
    return *getp(bmat,row,column);
}

/* Sets the value of an element in the band matrix based on full matrix indices. 
* @param val The value to be set at the specified indices.
*/
void setv(band_mat *bmat, long row, long column, double val) {
    *getp(bmat,row,column) = val;
}

/* Solves the linear system Ax = b for a band matrix A stored in the provided structure.
 *
 * @param bmat Pointer to the band matrix structure.
 * @param x Array to store the solution vector (length ncol).
 * @param b Right-hand side vector (length ncol).
 *
 * @return LAPACK info code: 0 indicates successful solution, otherwise error.
 */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
    // Copy bmat array into temporary storage for LAPACK
    int i,bandno;
    for(i=0;i<bmat->ncol;i++) { 
        for (bandno=0;bandno<bmat->nbrows;bandno++) {
        bmat->array_inv[bmat->nbrows_inv*i+(bandno+bmat->nbands_low)] = bmat->array[bmat->nbrows*i+bandno];
        }   
        x[i] = b[i];
    }

    long nrhs = 1;
    long ldab = bmat->nbands_low*2 + bmat->nbands_up + 1;
    int info = LAPACKE_dgbsv( LAPACK_COL_MAJOR, bmat->ncol, bmat->nbands_low, bmat->nbands_up, nrhs, bmat->array_inv, ldab, bmat->ipiv, x, bmat->ncol);
    return info;
}

/* Calculates the 1D element index corresponding to a specific grid point.
 *
 * @param j X-coordinate of the grid point.
 * @param p Y-coordinate of the grid point.
 * @param P Total number of grid points in the y-direction.
 *
 * @return The 1D element index corresponding to the grid point.
 */
long indx(long j, long p, long P) {
    return j * P + p;
}

/* Calculates the 2D grid point coordinates corresponding to a 1D element index.
 *
 * @param indx The 1D element index.
 * @param P The total number of grid points in the y-direction.
 * @param j Output: will contain the x-coordinate of the grid point.
 * @param p Output: will contain the y-coordinate of the grid point.
 * 
 */
void gridp(long indx, long P, long *j, long *p) {
    // Calculate y-coordinate (p) using modulo operation
    *p = indx % P;
    // Calculate x-coordinate (j) using integer division
    *j = (indx - *p) / P;
}

/* Swaps the contents of two double pointers.
 *
 * @param a Pointer to the first double pointer.
 * @param b Pointer to the second double pointer.
 */
void MemSwap(double** a, double** b) {
    // Temporary variable to hold the value during the swap
    double* tmp;
    // Store the value pointed to by a in the temporary variable
    tmp = *a;
    // Assign the value pointed to by b to the location pointed to by a
    *a = *b;
    // Assign the previously stored value (pointed to by a) to the location pointed to by b
    *b = tmp;
}

/* Prints the contents of a vector along with additional information to a file.
 *
 * @param file Pointer to the output file.
 * @param vector Pointer to a double array representing the vector.
 * @param Na Size of the vector.
 * @param pos Pointer to a long integer array containing additional information.
 * @param Ny Number of elements in the y-dimension (related to pos).
 * @param time A double value representing time.
 *
 * @return 0 on success, 1 on error.
 */
int printVector(FILE *file, double **vector, int Na, long int **pos, int Ny, double time) {
    // Check for null file pointer
    if (file == NULL) {
        printf("Error opening output file\n");
        return 1;
    }

    // Loop through each element in the vector
    for (int j = 0; j < Na; j++) {
        // Calculate x and y coordinates from pos
        long int x = (*pos)[j] / Ny;
        long int y = (*pos)[j] % Ny;

        // Print formatted output with time, coordinates, and vector value
        fprintf(file, "%.6f, %ld, %ld, %.6f\n", time, x, y, (*vector)[j]);
    }

    // Print a newline at the end
    fprintf(file, "\n");

    return 0;
}

int main() {

    /* Open and test input file opening */
    FILE *input = fopen("input.txt", "r");
    if (input == NULL) {
        perror("Error opening input file");
        return 0;
    }

    // The parameters input file, 'input.txt' will contain the following values:
    // • Nx: number of x grid points.
    // • Ny: number of y grid points.
    // • Na: number of active grid cells.
    // • Lx: Length of domain G in x direction.
    // • Ly: Length of domain G in y direction.
    // • tf : final time for time evolution.
    // • lambda: λ parameter in equations.
    // • tD: diagnostic timestep.

    int Nx, Ny, Na;
    double Lx, Ly, tf, lambda, tD;

    // Read the parameters from the file
    if (fscanf(input, "%d %d %d %lf %lf %lf %lf %lf ", &Nx, &Ny, &Na, &Lx, &Ly, &tf, &lambda, &tD) != 8) {
        fprintf(stderr, "Error reading parameters from file\n");
        fclose(input);
        return 0;
    }

    // Close the file
    fclose(input);

    // Description of the matrix structure
    // - Block tridiagoNal: one block on either side of the diagoNal blocks
  
    // Declare a band matrix structure named A - 
    band_mat A;
    // Initialise the band matrix A with specified parameters - how come Na?
    init_band_mat(&A, Na, Na, Na); 

    // Allocate memory for the starting vector and next vectors (in Ax = b)
    // These vectors will onLy be used for active grid cells.

    // Allocate memory for the starting vector (x)
    double *x = malloc(Na * sizeof(double));
    if (x == NULL) {
        fprintf(stderr, "Error allocating memory for starting vector\n");
        exit(1); // Exit with error code
    }

    // Allocate memory for the next vector (b)
    double *b = malloc(Na * sizeof(double));
    if (b == NULL) {
        fprintf(stderr, "Error allocating memory for next vector\n");
        free(x); // Free previousLy allocated memory for x in case of error
        exit(1); // Exit with error code
    }

    // Calculate grid spacing in x and y directions
    double dx = Lx / Nx;
    double dy = Ly / Ny;

    // Initialize arrays for storing matrix and equation indices
    long *index = malloc(Na * sizeof(long)); // Allocate memory for index array
    if (index == NULL) {
        fprintf(stderr, "Error allocating memory for index array\n");
        exit(1); // Exit with error code
    }

    long *equation_indx = malloc(Ny * Ny * sizeof(long)); // Allocate memory for equation_indx array
    if (equation_indx == NULL) {
        fprintf(stderr, "Error allocating memory for equation_indx array\n");
        free(index); // Free previously allocated memory for index in case of error
        exit(1); // Exit with error code
    }

    // Initialize equation_indx array to -1, indicating uNassigned equations
    for (int i = 0; i < Ny * Ny; i++) {
        equation_indx[i] = -1;
    }

    // Calculate maximum value of the solution (potential scaling factor)
    double max_u = sqrt(fabs(lambda));

    // Open the coefficients file for reading
    FILE *coeff = fopen("coefficients.txt", "r");
    if (coeff == NULL) {
        perror("Error opening Coefficients file");
        return 0;
    }

    // Loop through each line in the coefficients file
    for (int i = 0; i < Na; i++) {
        // Read three values (j, p, b) from the file
        int j, p;
        if (fscanf(coeff, "%d %d %lf", &j, &p, &b[i]) != 3) {
            fprintf(stderr, "Error reading values from coefficients file\n");
            fclose(coeff);
            // Free previously allocated memory (if aNy) to avoid leaks
            free(index);
            free(equation_indx);
            free(x);
            free(b);
            return 0;
        }

        // Update max_u if the absolute value of the read coefficient is larger
        if (fabs(b[i]) > max_u) {
            max_u = fabs(b[i]);
        }

        // Print the read coefficient for debugging purposes
        printf("Coeff Read     = %lf\n", b[i]);

        // Calculate the index in the matrix using the function indx(j, p, Ny)
        index[i] = indx(j, p, Ny);

        // Assign the equation index to the corresponding element in equation_indx
        equation_indx[indx(j, p, Ny)] = i;
    }

    // Close the coefficients file
    fclose(coeff);

    // Set correction factor based on the sign of lambda
    int cor = lambda < 0 ? -1 : 1;

    /*Each B is an upper bound for the timestep */
    double B1 = ((0.5*dx*dx*dy*dy)/( dx*dx + dy*dy));
    double B2;
    if(max_u > sqrt(lambda*cor) && max_u > 0){
        B2 = 1/(max_u*(sqrt(lambda*cor)+max_u));
    } else {
        B2 = 1/(2*lambda*cor);
    }

    // Find the minimum value between B1 and B2
    double min_bound = fmin(B1, B2);

    // Print diagnostic information with clear labels
    printf("Max_u:                 %lf\n", max_u);
    printf("Diagnostic timestep:   %lf\n", tD);
    printf("Spatial bound:         %lf\n", B1);
    printf("Maximum bound:         %lf\n", 1 / (max_u * (sqrt(lambda) + max_u)));
    printf("Non-linear term bound: %lf\n", 1 / (2 * lambda));


    // Initialize timestep and factor variables
    double dt = tD;
    int factor = 1;

    // Calculate the smallest timestep that is less than or equal to the minimum bound
    while (dt > min_bound) {
        factor++;
        dt = tD / factor;
    }

    // Print the chosen timestep
    printf("Timestep dt          = %lf\n", dt);

    // Iterate through each equation in the system
    for (int i = 0; i < Na; i++) {
        // Set the diagoNal element of the matrix to 1
        setv(&A, i, i, 1);

        // Determine neighboring active grid points based on index values
        bool top = (index[i] % Ny != Ny - 1) && (index[i + 1] - index[i] == 1);
        bool bottom = (index[i] % Ny != 0) && (index[i] - index[i - 1] == 1);
        bool left = (index[i] >= Ny) && (equation_indx[index[i] - Ny] != -1);
        bool right = (index[i] < (Nx - 1) * Ny) && (equation_indx[index[i] + Ny] != -1);

        // Modify matrix elements based on boundary conditions and neighboring grid points
        if(bottom){
            setv(&A, i , i - 1  , -(dt)/(dy*dy));
            setv(&A, i , i  , getv(&A,i,i) +(dt)/(dy*dy));
        } 
        if(top){
            setv(&A, i , i + 1  , -(dt)/(dy*dy));
            setv(&A, i , i  , getv(&A,i,i) +(dt)/(dy*dy));
        }
        if(right){
            setv(&A, i , equation_indx[index[i] + Ny]  , -(dt)/(dx*dx));
            setv(&A, i , i  , getv(&A,i,i) +(dt)/(dx*dx));
        }
        if(left){
            setv(&A, i , equation_indx[index[i] - Ny]  , -(dt)/(dx*dx));
            setv(&A, i , i  , getv(&A,i,i) +(dt)/(dx*dx));
        }
    }

    // Free the memory allocated for equation_indx
    free(equation_indx);
  
    // Open output file for writing
    FILE *out = fopen("output.txt", "w");

    // Print initial solution
    printVector(out, &b, Na, &index, Ny, 0.0);

    // Calculate number of iterations
    int iter = floor(tf / dt);
    printf("Iterations: %d\n", iter);

    // Time stepping loop
    for (int i = 1; i <= iter; i++) {
        // Update solution using explicit scheme
        for (int j = 0; j < Na; j++) {
        b[j] = (1 + dt * lambda) * b[j] - dt * b[j] * b[j] * b[j];
        }

        // Solve linear system using the updated solution as the right-hand side
        solve_Ax_eq_b(&A, x, b);

        // Swap solution and temporary vectors
        MemSwap(&x, &b);
        // Print solution at specific time intervals
        double tolerance = 1e-5; // A small number close to the precision you expect
        double checkTime = dt * factor;
        if (fmod(i * dt, checkTime) < tolerance || fabs(i * dt - checkTime) < tolerance) {
        printVector(out, &b, Na, &index, Ny, i * dt);
        }
    }

    // Close output file
    fclose(out);

    /* Free memory */
    finalise_band_mat(&A);
    free(x);
    free(b);
    free(index);

    return 0;
}