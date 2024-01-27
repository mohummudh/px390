#include <stdlib.h>
#include <stdio.h>
#include <mkl_lapacke.h>
#include <math.h>

/* Define structure that holds band matrix information */
struct band_mat{
  long ncol;        /* Number of columns in band matrix */
  long nbrows;      /* Number of rows (bands in original matrix) */
  long nbands_up;   /* Number of bands above diagonal */
  long nbands_low;  /* Number of bands below diagonal */
  double *array;    /* Storage for the matrix in banded format */
  /* Internal temporary storage for solving inverse problem */
  long nbrows_inv;  /* Number of rows of inverse matrix */
  double *array_inv;/* Store the inverse if this is generated */
  int *ipiv;        /* Additional inverse information */
};
/* Define a new type band_mat */
typedef struct band_mat band_mat;

/* Initialise a band matrix of a certain size, allocate memory,
   and set the parameters.  */ 
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
  bmat->nbrows = nbands_lower + nbands_upper + 1;
  bmat->ncol   = n_columns;
  bmat->nbands_up = nbands_upper;
  bmat->nbands_low= nbands_lower;
  bmat->array      = (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol);
  bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
  bmat->array_inv  = (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol);
  bmat->ipiv       = (int *) malloc(sizeof(int)*bmat->ncol);
  if (bmat->array==NULL||bmat->array_inv==NULL) {
    return 0;
  }  
  /* Initialise array to zero */
  long i;
  for (i=0;i<bmat->nbrows*bmat->ncol;i++) {
    bmat->array[i] = 0.0;
  }
  return 1;
};

/* Finalise function: should free memory as required */
void finalise_band_mat(band_mat *bmat) {
  free(bmat->array);
  free(bmat->array_inv);
  free(bmat->ipiv);
}

/* Get a pointer to a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double *getp(band_mat *bmat, long row, long column) {
  int bandno = bmat->nbands_up + row - column;
  if(row<0 || column<0 || row>=bmat->ncol || column>=bmat->ncol ) {
    printf("Indexes out of bounds in getp: %ld %ld %ld \n",row,column,bmat->ncol);
    exit(1);
  }
  return &bmat->array[bmat->nbrows*column + bandno];
}

/* Retrun the value of a location in the band matrix, using
   the row and column indexes of the full matrix.           */
double getv(band_mat *bmat, long row, long column) {
  return *getp(bmat,row,column);
}

/* Set an element of a band matrix to a desired value based on the pointer
   to a location in the band matrix, using the row and column indexes
   of the full matrix.           */
double setv(band_mat *bmat, long row, long column, double val) {
  *getp(bmat,row,column) = val;
  return val;
}

/* Solve the equation Ax = b for a matrix a stored in band format
   and x and b real arrays                                          */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
  /* Copy bmat array into the temporary store */
  int i, bandno;
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

int printmat(band_mat *bmat) {
  long i,j;
  for(i=0; i<bmat->ncol;i++) {
    for(j=0; j<bmat->nbrows; j++) {
      printf("%ld %ld %g \n",i,j,bmat->array[bmat->nbrows*i + j]);
    }
  }
  return 0;
}

void swap(double** first, double** second){
  double* temporary;
  temporary = *first;
  *first = *second;
  *second = temporary;
} 

int main() {
  band_mat bmat;

  /* constants from input.txt and potential.txt */
  long int N, I;
  double x_L, x_R, a_L, a_R, b_L, b_R, E_0;

  FILE *input;
  input = fopen("input.txt", "r");

  if (fscanf(input, "%lf %lf %ld %lf %lf %lf %lf %lf %ld", &x_L, &x_R, &N, &a_L, &a_R, &b_L, &b_R, &E_0, &I) != 9){
    fprintf(stderr, "Error reading values from the file.\n");
    fclose(input);
    return 1;
  }
  fclose(input);

  FILE *potential = fopen("potential.txt", "r");
  if (potential == NULL) {
      fprintf(stderr, "Error opening potential.txt\n");
      exit(1);
  }

  double *P = malloc(sizeof(double) * N);
  if (P == NULL) {
      fprintf(stderr, "Memory allocation failed for potentials\n");
      fclose(potential);
      exit(1);
  }

  for (int i = 0; i < N; i++) {
    if (fscanf(potential, "%lf", &P[i]) != 1) {
        fprintf(stderr, "Error reading potentials\n");
        free(P);
        fclose(potential);
        exit(1);
     }
  }
  fclose(potential);

  /* making the matrix */

  double dx = (x_R - x_L) / (N - 1); /* step */
  long ncols = N; /* columns in the matrix equal the number of grid points */

  /* We have a three-point stencil (domain of numerical dependence) of
     our finite-difference equations:
     1 point to the left  -> nbands_low = 1
     1       to the right -> nbands_up  = 1
  */

  long nbands_low = 1;  
  long nbands_up  = 1;

  /* initialise the band matrix */
  init_band_mat(&bmat, nbands_low, nbands_up, ncols);

  double *x = malloc(sizeof(double)*ncols);
  double *b = malloc(sizeof(double)*ncols);

  long i, j;
  /* Loop over the equation number and set the matrix
     values equal to the coefficients of the grid values. 
     Note that boundaries are treated with special cases */

  //srand(2);
  for(i=0; i<ncols; i++) {
    if(i>0)       {setv(&bmat,i,i-1,-1.0/(dx*dx));};
    setv(               &bmat,i,i,   ((2 + P[i]*dx*dx) / (dx*dx)) - E_0);
    if(i<ncols-1) {setv(&bmat,i,i+1,-1.0/(dx*dx));};

    /* Initialise guess for Ψ0 */
    
    b[i] = 1.0;   
  }
  
  if (b_L==0) {
    setv(&bmat, 0,0,1); /* Beginning of main diagonal */
    setv(&bmat, 0,1,0); /* Beginning of upper diagonal */
    b[0] = 0;
  }
  else {
    setv(&bmat, 0,0,((-(2*a_L*dx / b_L) + (2 + P[0]*dx*dx)) / (dx*dx)) - E_0); /* Beginning of main diagonal */
    setv(&bmat, 0, 1, -2 / (dx*dx)); /* Beginning of upper diagonal */
  }

  if (b_R==0) {
    setv(&bmat, ncols - 1, ncols - 1, 1); /* End of main diagonal */
    setv(&bmat, ncols - 1, ncols - 2, 0); /* End of upper diagonal */
    b[ncols - 1] = 0;
  }
  else {
    setv(&bmat, ncols - 1, ncols - 2, -2 / (dx*dx)); /* End of lower diagonal */
    setv(&bmat, ncols - 1, ncols - 1, (((2*a_R*dx / b_R) + (2 + P[ncols - 1]*dx*dx)) / (dx*dx)) - E_0); /* End of main diagonal */
  }

  /*  Print coefficient matrix for debugging: */ 
  // printmat(&bmat);          

  FILE *output = fopen("output.txt", "w");
  if (output == NULL) {
    printf("Error opening output.txt\n");
    exit(1);
  }

  for (i = 0; i < ncols; i++) {
    fprintf(output,"%ld, %.6f, %.6f \n", 0, x_L + i*(dx), b[i]);
  }

/* Calculate Ψ_i for 1 <= i <= I*/
  for (i = 1; i <= I; i++) {
    /* Calculate Ψ_i */
    solve_Ax_eq_b(&bmat, x, b); 

    /* Making next step current step */
    swap(&b, &x);

    /* brute force boundary condition on the left*/
    if (b_L == 0){
          b[0] =0; 
        }

    /* Print Ψ_i */
    for (j = 0; j < ncols; j++) {
      fprintf(output, "%ld, %.20f, %.20f \n", i, x_L + j*(dx), b[j]);
    }
    
    fprintf(output, "\n"); 
  }

  fclose(output);

  /* Memory clean up */
  finalise_band_mat(&bmat);
  free(P);
  free(x);
  free(b);
  
  return(0);
}