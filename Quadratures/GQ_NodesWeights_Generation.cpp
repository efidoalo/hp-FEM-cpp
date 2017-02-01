/*=====================================================;
 *
 *  File: GQ_NodesWeights_Generation.cpp
 *  Content: Definitions of functions required for 
 *  determining the quadrature weights and node points
 *  to be used during evaluation of bilinear forms
 *  over reference elements
 *  Date: 13/02/2016
 *
 ******************************************************/

#include "GQ_NodesWeights_Generation.h"
#include <cmath>

// Function that returns the value of the nth Legendre Polynomial at the point x in R^1  
double legendre(const double& x,const int& n) 
{   
  double Pprev1 = x, Pprev2 = 1.0, k = 1.0, CurrentP;
  if (n==0) 
  return 1.0;     
  if (n==1) 
  return x; 
  for (int i=0; i<n-1; ++i) { // iterative loop to determine the value of the nth legendre polynomial at the point x for n>1
     CurrentP = ((((2.0*k) + 1.0)/(k+1.0))*x*Pprev1) - ((k/(k+1.0))*Pprev2);
     Pprev2 = Pprev1;
     Pprev1 = CurrentP;
     ++k;
  }  
  return CurrentP;        
}

// Function that returns the value of the Derivative of the nth Legendre Polynomial at the point x in R^1
double diff_legendre(const double& x, const int& n)   
{   
  double PDashprev1 = 1.0, PDashprev2 = 0.0, k = 1.0, CurrentPDash;
  if (n==0) 
  return PDashprev2;     
  if (n==1) 
  return PDashprev1;
  for (int i=0; i<n-1; ++i) {   // iterative loop to determine the value of the 1st derivative of the nth legendre polynomial at the point x
    CurrentPDash = ((((2.0*k) + 1.0)/k)*x*PDashprev1) - (((k+1.0)/k)*PDashprev2);
    PDashprev2 = PDashprev1;
    PDashprev1 = CurrentPDash;
    ++k;
  }  
  return CurrentPDash;    
}

// Function that stores the zeros of the nth order Legendre Polynomial at addr zeros_store contiguously
// These  zeros are used as the n nodes for the Gauss-Legendre quadrature.
// The n weights corresponding to these nodes are stored at addr weights_store contiguously 
void gauss_quadrature(double* zeros_store, double* weights_store, const int&  n)   
{
  const double PI = std::acos(-1.0);
  double* zero_prev_estimates = new double[n];  
  double* zero_current_estimates = new double[n];
  double tol = 1e-30, current_error = 1, update_error;

  for (int i=0; i<n; ++i) 
  zero_prev_estimates[i] = -cos((((2.0*(i+1)) - 1)/(2.0*n))*PI);

  while (current_error>=tol)  { 
    // Update the approximation of zero i, for  i = 1, 2, ... , n ( at subscripts 0,1,2,...,n-1)	  
    for (int i=0; i<n; ++i) 
    zero_current_estimates[i] = zero_prev_estimates[i] - ((legendre(zero_prev_estimates[i],n))/(diff_legendre(zero_prev_estimates[i], n)));
 
    // Error Tracking  
    current_error = fabs(zero_current_estimates[0] - zero_prev_estimates[0]);         // initial current error for this iteration    
    for (int j=1; j<n; ++j) {
      update_error = fabs(zero_current_estimates[j] - zero_prev_estimates[j]);  // error at sucscript j is update_error
      if (update_error>current_error) 
      current_error = update_error;  
    }         
    // Update zero_prev_estimates for use in next iteration     
    for (int i=0; i<n; ++i) 
    zero_prev_estimates[i] = zero_current_estimates[i];                             
  }

  // Copy over zeros to contiguous memory where subscript 0 is located at address double* zero_store that was given as function arguement
  for (int k=0; k<n; ++k) 
  zeros_store[k] = zero_current_estimates[k];
  
  // Place computed weights in to contiguous memory where subscript 0 is locate at address double* weights_store that was given as function arguement
  for (int i=0; i<n; ++i) 
  weights_store[i] = 2.0/((1 - pow(zero_current_estimates[i],2))*(pow(diff_legendre(zero_current_estimates[i],n) , 2)));     
  
  // Deallocate 
  delete[] zero_prev_estimates;
  delete[] zero_current_estimates;

}

std::vector<double> GaussianQuadrature::Nodes = { };
std::vector<double> GaussianQuadrature::Weights = { };


