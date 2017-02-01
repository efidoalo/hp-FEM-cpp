/*======================================================;
 *
 *  File: PDE_Definition.cpp
 *  Contents: Function definitions required to define 
 *  the PDE under consideration. 
 *
 *  PDE Defined by -a(x)*U'' + b(x)*U' + c(x)*U = f on the itnerval (a,b) 
 *  a(x) must have positive constant value.
 *
 *  Date: 19/07/2016
 *  Author: A.Oldham
 *
 ********************************************************/

#include "PDE_Definition.h"
#include <cmath>

// functions that must be defined by the user:
//
//  a( . , . ) MUST RETURN POSITIVE CONSTANT
double a(const double& x)
{
  return 1;// 0.00001;
}

// 
double b(const double& x)
{
  return 0.0;
}

double c(const double& x)
{
  return 0.0; // 1.0;
}

double f( const double& x)
{
  return   cos( 2*acos(-1.0)*x );  //sin( 2*acos(-1.0)*x ); // //   1.0; // sin(((x - 27.6)/(73.6/(2.0*acos(-1.0)))));
}
//


// OPTIONAL
double analytic_solution(const double& x)
{
  //double numer1 = exp(x/sqrt(a(0.0))), numer2 = exp( (-x)/sqrt(a(0.0)) )*exp( 1/sqrt(a(0.0)) );
  //double denom1 = exp(1.0/sqrt(a(0.0)));
  //double ret = (-1.0*(numer1/(denom1+1.0))) - (numer2/(denom1+1.0)) + 1;
  //return ret;

	  
	  
  return (1.0/(4.0*pow(acos(-1.0),2.0)))*cos( 2*acos(-1.0)*x );
                              //(1.0/(4.0*pow(acos(-1.0),2.0)))*sin( 2*acos(-1.0)*x ); // (1.0/((4.0*(pow(acos(-1.0),2.0)))+1.0))*sin(2*acos(-1.0)*x);                                                                  
		                                                               //    (pow((73.6/(2.0*acos(-1.0))),2.0))*sin((x - 27.6)/(73.6/(2.0*acos(-1.0))));    
                                                                               //  (1.0/(4.0*pow(acos(-1.0),2.0)))*sin( 2*acos(-1.0)*x );
									       //  (1.0/(4.0*pow(acos(-1.0),2.0)))*cos( 2*acos(-1.0)*x ); 


}

// OPTIONAL
double analytic_solution_deriv(const double& x)
{
  double C = -1.0/(sqrt(0.00001)*(exp(1.0/sqrt(0.00001)) + 1.0));
  double v1 = exp(x/sqrt(0.00001)), v2=exp( (1-x)/sqrt(a(0.0)) );
  return C*(v1-v2);
}

// Required to be user defined if nonhomogeneous boundary conditions are being considered
// nonhomogenous dirichlet boundary contition [if required], specified as a function evaluated at the (global) point std::vector<double> Point
// currently for 1D problem -u'' = cos(2pix) on [0,1] 
double nonhomogeneousDirichlet(const std::vector<double>& Point)
{
  double x=Point[0];
  return analytic_solution(x);
}
