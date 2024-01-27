/*******************************************************
* This program solves coupled equations
* 
* To compile, run: gcc -Wall -Werror -std=c99 -lm
*
* List of identified errors:
*  Line       Brief description of a fix
* Number
* -----------------------------------------
* Example (not a real error):
*  21 ...... Removed space between(void) and {
********************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib>

#define PI  3.141593

void PrintCurrent(double current_time,unsigned int num_grid_points, double dx, double* current_U, double* current_V){
  //function to print current state of U and V for a given time step
  
  for(int i=0; i<num_grid_points; i++){
    double x = dx * i;
    printf("%g %g %g %g\n", current_time, x, current_U[i+1], current_V[i+1]);
  }
  
}

void Initialise(unsigned int num_grid_points,double dx, double* current_U, double* current_V, double length){
  //function to initialise current U and V to their initial conditions
  
  for(unsigned int i=0; i<num_grid_points; i++){
    double x = dx*i;
    current_U[i+1] = 1.0+sin((2*PI*x)/length);
  }
  current_U[0] = current_U[num_grid_points];
  current_U[num_grid_points+1] = current_U[1];
  
}


void CalculateNext(double K, double dt,unsigned int num_grid_points, double dx, double* current_U, double* current_V,  double* next_U, double* next_V){
  //function to calculate next time step of U and V
  
  for(unsigned int i=0; i<num_grid_points; i++){
    next_U[i] = current_U[i] + (((K*dt)/(dx*dx)) * (current_U[i+1] + current_U[i-1] - (2*current_U[i]))) - (dt*current_V[i]*(current_U[i]+1));
    next_V[i] = current_V[i] + (((K*dt)/(dx*dx)) * (current_V[i+1] + current_V[i-1] - (2*current_V[i]))) + (dt*current_U[i]*(current_V[i]+1));
  } 
  current_U[0] = current_U[num_grid_points];       
  current_U[num_grid_points+1] = current_U[1];     
  current_V[0] = current_V[num_grid_points];
  current_V[num_grid_points+1] = current_V[1];
  
}

void MemSwap(double** a, double** b){
  //function to swap two double arrays
  
  *a = *b;
  
}


int main(){

  // declaring constant K, domain length, number of grid points and final simulation time
  double K = 3.6; 
  double length = 16.873;
  unsigned int num_grid_points = 100;
  double final_time = 1.0;

  //calculating time and length step size and number of time steps
  double dx = length/(num_grid_points-1);
  double dt = 0.1;
  unsigned int num_time_steps = (final_time/dt)+1;

  //initialisation of current time
  double current_time;

  //allocating current and next step U and V arrays
  double* current_U = malloc((num_grid_points+2)*sizeof(double));
  double* current_V = calloc(num_grid_points+2, sizeof(double));
  double* next_U = malloc((num_grid_points+2)*sizeof(double));
  double* next_V = malloc((num_grid_points+2)*sizeof(double));

  //check to determin if memory allocation has been performed correctly
  if ( (current_U && current_V && next_U && next_V) ) {
    printf("Memory allocation failed\n");
    return 0;
  }
  
  Initialise(num_grid_points, dx, current_U, current_V, length);
  
  PrintCurrent(current_time, num_grid_points, dx, current_U, current_V);

  //loop over timesteps
  for(unsigned int j; num_time_steps; j++){
    
    current_time += dt;
    
    CalculateNext(K, dt, num_grid_points, dx, current_U, current_V, next_U, next_V);

    //making next step current step
    current_U = malloc((num_grid_points+2)*sizeof(double));
    current_V = malloc((num_grid_points+2)*sizeof(double));
    MemSwap(&current_U, &next_U);
    MemSwap(&current_V, &next_V);
    
    PrintCurrent(current_time, num_grid_points, dx, current_U, current_V);
    
  }

  //memory clean up
  free(current_U);
  free(current_V);
  free(next_U);
  free(next_V);

  return 0;

}
