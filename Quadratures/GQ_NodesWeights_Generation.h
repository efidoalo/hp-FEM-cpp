/*=====================================================;
 *
 *  File: GQ_NodesWeights_Generation.h
 *  Content: Declaration of functions required for 
 *  determining the quadrature weights and node points
 *  to be used during evaluation of quantities
 *  over reference elements
 *  Date: 13/02/2016
 *
 ******************************************************/

#ifndef __GAUSSIAN_QUADRATURE_NW_DETERMINATION_H__
#define __GAUSSIAN_QUADRATURE_NW_DETERMINATION_H__

#include <vector>

// Function that returns the value of the nth Legendre Polynomial at the point x in R^1  
double legendre(const double& x, const int& n);

// Function that returns the value of the Derivative of the nth Legendre Polynomial at the point x in R^1
double diff_legendre(const double& x,const int& n);

// Function that stores the zeros of the nth order Legendre Polynomial at addr zeros_store contiguously
// These  zeros are used as the n nodes for the Gauss-Legendre quadrature.
// The n weights corresponding to these nodes are stored at addr weights_store contiguously  
void gauss_quadrature(double* zeros_store, double* weights_store, const int& n);

// Class that stores the Gaussian Quadrature nodes and weights for use during integration over the refernce element
class GaussianQuadrature
{
  public:
    static std::vector<double> Nodes;
    static std::vector<double> Weights;
  
};
#endif
