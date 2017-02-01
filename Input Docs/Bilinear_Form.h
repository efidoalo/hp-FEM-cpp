/*=====================================================================;
 *
 *  File: Bilinear_Form.h
 *  Content: Declaration and Definition of the bilinear form.
 *  Date: 12/06/2016
 *  Author: Andrew Oldham
 *
 ***********************************************************************/
#ifndef __BILINEAR_FORM_H__
#define __BILINEAR_FORM_H__

#include <cmath>
#include "../General Finite Element Class/Finite_Element.h"
#include "../Operators/Polynomial_Operators.h"
#include "../Mesh/Mesh.h"
#include <iostream>

// VOLUME INTEGRAL DEFINITIONS [if present in the bilinrear form] G1(.,.) = BF_j'*BF_i' for this example
double BilinearForm_VI(Finite_Element& FE, const int& j, const int& i, std::vector<double>& Point);

// H1(.,.) = f*BF_i for this example
double LinearFunctional_VI(Finite_Element& FE, const int& i, std::vector<double>& Point);


// Returns an estimate of the error commited via the solution on the finite element numbered FE_number (measured in H1 Norm).
// This function uses a posteriori errorr estimators.
double local_error_estimator(Mesh& COMPUTATIONALMESH, double* coeffVals, const int& FE_number, const double& initial_h);

// Function defined on element that is used n the energy norm estimate for the test problem (code verification).
double weight_Fnc(const double& x_j_sub1, const double& x_j_sup1, const double& x);

// returns value to (along with local_error_estimator()) decide whether the finite element should be refined or not.
double element_error_criteria(Mesh& COMPUTATIONALMESH, const int& FE_number, const double& initial_h, 
		const double& global_tol, const int& initial_NoOfElements, const double& C_0);

#endif
