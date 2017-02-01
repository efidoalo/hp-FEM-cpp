/*=================================================;
 *
 *  File: Polynomial_Operators.h
 *  Content: Header file for the Functions/Routines
 *  that perform operations of Polynomials, whereby
 *  a polynomial is represented as a vector<double>
 *  Date: 12/06/2016
 *  Author: Andrew Oldham
 *
 **************************************************/

#ifndef __POLYNOMIAL_OPERATORS_H__
#define __POLYNOMIAL_OPERATORS_H__

#include "BaseRep.h"
#include <utility>
#include <vector>

// Finds BaseRep& Key as the first element of a pair in Map and returns the second value of that pair (which is an integer).
// Assumes that no first data member of any two pairs in the vector is identical. (ie all pairs in vector have distinct first members)
int index(const std::vector< std::pair<BaseRep,int> >& Map, const BaseRep& Key);

// Differentiates the multivariate polynomial given by Operand w.r.t the spatial variable given by IndependentVariable recursively Order times.
std::vector<double> Differentiate(const int& IndependentVariable, const int& Order, std::vector<double> Operand);

// Evaluate the multivariate polynomial Poly at the point Point, [Point.size()=Poly[2]].
double Eval(const std::vector<double>& Point, const std::vector<double>& Poly);

// Linear Combination of multivariate polynomials (of identical degree and number of spatial variables)
std::vector<double> LinearComb(std::vector<double> scalars, std::vector< std::vector<double> > OrderedPolynomials);

#endif
