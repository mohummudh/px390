/*****************************************************************
 * This programme defines functions used to generate band
 * matrix for problems in which a solution at some grid point
 * k is coupled to a solution at a number of adjasont points.
 * The file also includes a main() function with an
 * example of how to solve a stationary 2D heat equation with a
 * constant source, for two different sources.
 * 
 * Prior to compilation execute following lines on nenneke:
 * module purge
 * module load intel impi imkl
 * Then:
 * Compile:  gcc -o w10wrkshop w10wrkshop.c -lmkl -liomp5 -lm 
 * Run: ./w10wrkshop
 * ****************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <mkl_lapacke.h>
#include <stdbool.h>
#include <math.h>

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

/* Initialise a band matrix of a certain size, allocate memory,
   and set the parameters.  */ 
int init_band_mat(band_mat *bmat, long nbands_lower, long nbands_upper, long n_columns) {
  bmat->nbrows = nbands_lower + nbands_upper + 1;
  bmat->ncol   = n_columns;
  bmat->nbands_up = nbands_upper;
  bmat->nbands_low= nbands_lower;
  //bmat->array      = (double *) malloc(sizeof(double)*bmat->nbrows*bmat->ncol);
  bmat->array      = (double *) calloc(bmat->nbrows*bmat->ncol, sizeof(double));
  bmat->nbrows_inv = bmat->nbands_up*2 + bmat->nbands_low + 1;
  //bmat->array_inv  = (double *) malloc(sizeof(double)*(bmat->nbrows+bmat->nbands_low)*bmat->ncol);
  bmat->array_inv  = (double *) calloc((bmat->nbrows+bmat->nbands_low)*bmat->ncol, sizeof(double));
  bmat->ipiv       = (int *) calloc(bmat->ncol,sizeof(int));
  if (!bmat->array||!bmat->array_inv ||!bmat->ipiv) {
    free(bmat->array);
    free(bmat->array_inv);
    free(bmat->ipiv);
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

void setv(band_mat *bmat, long row, long column, double val) {
  *getp(bmat,row,column) = val;
}

/* Solve the equation Ax = b for a matrix a stored in band format
   and x and b real arrays                                          */
int solve_Ax_eq_b(band_mat *bmat, double *x, double *b) {
  /* Copy bmat array into the temporary store */
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


/*Check that a grid point has valid coordinates */
int is_valid(long j, long p, long J, long P) {
  return (j>=0)&&(j<J)&&(p>=0)&&(p<P);
}


/* Return the 1D element index corresponding to a particular grid point.
   We can rewrite this function without changing the rest of the code if
   we want to change the grid numbering scheme!
   Output: long integer with the index of the point
   Input:
   long j:  The X grid point index
   long p:  The Y grid point index
   long P:  The number of Y points.
*/
long indx( long j, long p, long P) {
  return j*P + p;
}

/* Return the 2D point corresponding to a particular 1D grid index */
void gridp(long indx, long P, long *j, long *p) {
  *p = indx%P;
  *j = (indx-(*p))/P;
}
/*Custom Built Functions*/

void MemSwap(double** a, double** b){
  double* tmp;
  tmp = *a;
  *a = *b;
  *b = tmp;
}
int printVector(FILE *file,double **vector,int na,long int **index, int Ny,double t) {
  if(file == NULL){
    printf("Error opening output file\n");
    return 1;
  }else{
    for(int j = 0; j < na; j++) {
      fprintf(file, "%.6f, %lu, %lu, %.6f\n", t, (*index)[j]/Ny, (*index)[j]%Ny , (*vector)[j]);
    }
    fprintf(file, "\n");
    return 0;
  }
}

/* An example of how to use the band matrix routines to solve a PDE:
   The equation solved is related to the steady state solution of the heat 
   diffusion equation.   
*/
int main() {
  /* We have a three-point stencil (domain of numerical dependence) of
     our finite-difference equations:
     1 point to the left  -> nbands_low = 1
     1       to the right -> nbands_up  = 1
  */

  /*We open and test input file opening*/
  FILE *input;
  input = fopen("input.txt", "r");
  if (input == NULL) {
      perror("Error opening input file");
      return 0;
  }

  int nx, ny, na;
  double lx, ly, tf, lam, td;

    // Read the parameters from the file
  if (fscanf(input, "%d %d %d %lf %lf %lf %lf %lf ", &nx, &ny, &na, &lx, &ly, &tf, &lam, &td) != 8) {
      fprintf(stderr, "Error reading parameters from file\n");
      fclose(input);
      return 0;
  }
  fclose(input); // Close the file



  /*Storage for the source term and the solution */
  //Use for only active grid cells
  double *x = calloc(na, sizeof(double)); //My next vector
  double *b = calloc(na, sizeof(double)); //My vector

  /* Initialising the dx and dy. Zero boundary conditions
  accounted for later on. */
  double dx = lx/(nx);
  double dy = ly/(ny);
  /* Loop over the equation number and set the matrix
     values equal to the coefficients of the grid values 
     note boundaries treated with special cases           */
  int j, p;
  long *index = malloc(sizeof(double)*(na)); //Stores the grid index.
  long *equation_indx = malloc(sizeof(double)*ny*ny); //Stoes vector index.

  /*Initialise Equation Index to -1 initially*/
  for(int i=0;i<nx*ny;i++){
    equation_indx[i] = -1;
  }

  /*Correction sign for taking the absolutve value of lamda*/
  int cor;
  if (lam<0){
      cor = -1;
    }else{
      cor = 1;
    }

  double max_u = sqrt(lam*cor); //Intialise
 /*We open and test input file opening*/
  FILE *coeff;
  coeff = fopen("coefficients.txt", "r");
  if (coeff == NULL) {
      perror("Error opening Coefficients file");
      free(x);
      free(b);
      free(index);
      free(equation_indx);
      return 0;
  }
  /*Reading in coefficients into vector b and indexing it*/
  for(int i = 0; i<na; i++){
    if(fscanf(coeff,"%d %d %lf",&j, &p, &b[i]) ==3 ){
      if(fabs(b[i])>max_u){
        max_u = b[i];
      }
      index[i] = indx(j,p,ny); //indx(j,p,P) = j*P + p
      equation_indx[indx(j,p,ny)] = i;
    }else{
      fprintf(stderr, "Error reading values from coefficients file\n");
      fclose(coeff);
      free(index);
      free(equation_indx);
      free(x);
      free(b);
      return 0;
    }
  }
  fclose(coeff);

  /* The matrix is block tridiagonal: one block on either
  side of the blocks of diagonals */ 

  band_mat bmat; //Initialising the Band Matrix
  init_band_mat(&bmat, ny, ny, na);


/*We select an appropriate timestep for stability inspired from CFL and Mathematical Analysis*/
  double bound = ((0.5*dx*dx*dy*dy)/( dx*dx + dy*dy));
  if(max_u == 0 && lam == 0){
    ; //Doing nothing
  }else{
    bound = fmin(bound,1/(max_u*(sqrt(lam*cor)+max_u)));
  }
  bound = fmin(bound,td);
  int factor = 1;
  double dt = td;
  while(dt > bound){
    factor = factor + 1;
    dt = td/factor;
  }
  // printf("CHECK dt = %lf\n", dt);


// Iterate through each equation in the system
 int top;
  int bottom;
  int left;
  int right;

  for(int i=0;i<na;i++){
    
    if(index[i]%ny == ny - 1){
      //If on top row
      top = -1;
    }else{
      top = equation_indx[index[i] + 1];
    }

    if(index[i]%ny == 0){
      //bottom row
      bottom = -1;
    }else{
      bottom = equation_indx[index[i] - 1];
    }
\
    if(index[i] < ny){
      //Left column
      left = -1;
    }else{
      left = equation_indx[index[i] - ny];
    }

    if(index[i] >= (nx-1)*ny  ){
      //right column
      right = -1;
    }else{
      right = equation_indx[index[i] + ny];
    }

    //Use line below to check
    //printf("left = %d, top = %d, right = %d, bottom = %d \n",left,top,right,bottom);

    
    setv(&bmat, i, i, 1+2*(dt)/(dy*dy)+2*(dt)/(dx*dx));
    if(bottom == -1){
      setv(&bmat, i , i  , getv(&bmat,i,i) - (dt)/(dy*dy));
    }else{
      setv(&bmat, i , bottom  , - (dt)/(dy*dy));
    }
    if(top == -1){
      setv(&bmat, i , i  , getv(&bmat,i,i) - (dt)/(dy*dy));
    }else{
      setv(&bmat, i , top  , - (dt)/(dy*dy));
    }
    if(right == -1){
      setv(&bmat, i , i  , getv(&bmat,i,i) - (dt)/(dx*dx));
    }else{
      setv(&bmat, i , right  , - (dt)/(dx*dx));
    }
    if(left == -1){
      setv(&bmat, i , i  , getv(&bmat,i,i) - (dt)/(dx*dx));
    }else{
      setv(&bmat, i , left  , - (dt)/(dx*dx));
    }
  }



    // for (int i = 0; i < na; i++) {
    //     // Set the diagonal element of the matrix to 1
    //     setv(&bmat, i, i, 1);
    //     bool top;
    //     bool bottom;
    //     bool left;
    //     bool right;
    //     if(index[i] % ny != ny - 1){
    //       top = false;
    //     }else if((equation_indx[index[i]+1] == -1)){
    //       top = false;
    //     }else{
    //       top =true;
    //     }

    //     if(index[i] % ny != 0){
    //       bottom = false;
    //     }else if((equation_indx[index[i]+1] == -1)){
    //       bottom = false;
    //     }else{
    //       bottom =true;
    //     }

    //      if(index[i] >= ny){
    //       left = false;
    //     }else if((equation_indx[index[i] - ny] != -1)){
    //       left = false;
    //     }else{
    //       left =true;
    //     }
        
    //      if(index[i] < (nx - 1) * ny){
    //       right = false;
    //     }else if((equation_indx[index[i] + ny] != -1)){
    //       right = false;
    //     }else{
    //       right =true;
    //     }

        
//         // Checking if the neighbouring cells are active
//         bool top = (index[i] % ny != ny - 1) && (i<na-1) && (index[i + 1] - index[i] == 1); //&& (equation_indx[index[i] + 1] !=-1)
//         bool bottom = (index[i] % ny != 0)  && (index[i] - index[i - 1] == 1); // && bottom = (equation_indx[index[i] - 1] != -1)
//         bool left = (index[i] >= ny) && (equation_indx[index[i] - ny] != -1);
//         bool right = (index[i] < (nx - 1) * ny) && (equation_indx[index[i] + ny] != -1);
// // 
//&& (i>0)
        // Modify matrix elements based on boundary conditions and neighboring grid points
    //     if(bottom){
    //         setv(&bmat, i , i - 1  , -(dt)/(dy*dy));
    //         setv(&bmat, i , i  , getv(&bmat,i,i) +(dt)/(dy*dy));
    //     } 
    //     if(top){
    //         setv(&bmat, i , i + 1  , -(dt)/(dy*dy));
    //         setv(&bmat, i , i  , getv(&bmat,i,i) +(dt)/(dy*dy));
    //     }
    //     if(right){
    //         setv(&bmat, i , equation_indx[index[i] + ny]  , -(dt)/(dx*dx));
    //         setv(&bmat, i , i  , getv(&bmat,i,i) +(dt)/(dx*dx));
    //     }
    //     if(left){
    //         setv(&bmat, i , equation_indx[index[i] - ny]  , -(dt)/(dx*dx));
    //         setv(&bmat, i , i  , getv(&bmat,i,i) +(dt)/(dx*dx));
    //     }
    // }

 
  
  /* File Computation output*/
  
  FILE *out = fopen("output.txt","w");
  if(printVector(out, &b, na, &index,ny,0.0)==1){
    perror("Error printing output\n");
    finalise_band_mat(&bmat);
    free(index);
    free(equation_indx);
    free(x);
    free(b);
    return 0;
  }
  int iter = floor(tf/dt);
  for(int i = 1; i<=iter; i++){
    for(int j = 0; j<na; j++){
      b[j] = (1+dt*lam)*b[j]-dt*b[j]*b[j]*b[j];
    }
    solve_Ax_eq_b(&bmat,x,b);
    MemSwap(&x,&b);
    // for(int k = 0; k<na; k++){
    // fprintf(out, "%.6f, %lu, %lu, %lf\n", i*dt, index[k]/ny, index[k]%ny , b[k]);
    // }

    if (i%factor == 0) {
        if(printVector(out, &b, na, &index, ny, i*dt)==1){
          perror("Error printing to output file\n");
          fclose(out);
          finalise_band_mat(&bmat);
          free(index);
          free(equation_indx);
          free(x);
          free(b);
          return 0;
          }
    }
  }
  fclose(out);

  /* Free memory */
  finalise_band_mat(&bmat);
  free(index);
  free(equation_indx);
  free(x);
  free(b);
  return 0;
}
