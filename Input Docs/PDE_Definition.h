/*======================================================;
 *
 *  File: PDE_Definition.h
 *  Contents: Function declarations required to define 
 *  the PDE under consideration.
 *  Date: 19/07/2016
 *  Author: A.Oldham
 *
 ********************************************************/

#ifndef __PDE_DEFINITION_H_INCLUDED__
#define __PDE_DEFINITION_H_INCLUDED__
#include <vector>
// PDE is of the form
//   -a(x)*U'' + b(x)*U' + c(x)*U = f,    on the itnerval (a,b)
//           where a(x) = Epsilon, a constant greater than zero.


// Declarations of functions defining the PDE
double a(const double& x);
double b(const double& x);
double c(const double& x);
double f(const double& x);


// Declaration of function specified when nonhomogeneous bc's are
// present
double nonhomogeneousDirichlet(const std::vector<double>& Point);

// Declaration of the analytic solution defined on the mesh
double analytic_solution(const double& x);
double analytic_solution_deriv(const double& x);

#endif
