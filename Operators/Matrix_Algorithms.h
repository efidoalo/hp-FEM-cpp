/*=================================================;
 *
 *  File: Matrix_Algorithms.h
 *  Content: Header file for the Functions/Eoutines 
 *  to invert matrices dependent upon their form.
 *  Date: 12/06/2016
 *  Author: Andrew Oldham
 *
 **************************************************/

#ifndef __MATRIX_ALGORIMTHS_H__
#define __MATRIX_ALGORITHMS_H__

#include "CSRMatrix.h"
#include <iostream>
#include <cmath>


// Ax = b SOLVERS, where appropriate we specify const double& tol [the tolerance] to be used.
//
// Conjugate Gradient Method
void ConjugateGradientMethod(CSRMatrix& M, double* b, double* x, const double& tol);

// Thomas Algorithm
void ThomasAlgorithm(CSRMatrix& M, double* b, double* x);

// Gauss Seidel Method 
void GaussSeidelMethod(double* entries, int* col_no, int* row_start, double* b, 
double* result, const int& n, const double& tol);

// Guassian Elimination (worst case scenario) 
void GaussianElimination(CSRMatrix& M, double* b, double* x);


// Output the input matrix (that is representedin CSR format) to cmd terminal
void CSR_MatrixOutput(double* entries, int* col_no, int* row_start, const int& n, const int& m);


// Interchange row r1 with row r2 [ r1 and r2 are 1 starting, r1<r2]
void elementary_row_operation1(CSRMatrix& M, const int& r1, const int& r2);
// Multiply row r1 by scalar
void elementary_row_operation2(CSRMatrix& M, const int& r1, const double& scalar);
// // Add sclaar times row r1 onto row r2
void elementary_row_operation3(CSRMatrix& M, const int& r1, const int& r2, const double& scalar);

// return the inverse ofthe input matrix
CSRMatrix CSR_INVERSION(double* matrix_entries, int* col_no, int* row_sart, const int& n);

// retunr the transpse of the input matrix
CSRMatrix CSR_TRANSPOSE(double* matrix_entries, int* col_no, int* row_start, const int& n);

// M is a CSRMatrix type that includes zeros in its matrix entries array, this function
// ensures that M has no zero entries in its matrix entries array.
// Returns number of zero values found in matrix M
void CSR_REDUCE(CSRMatrix& M);

// returns true if STIFFNESSMATRIX  is diagonaly dominant, false otherwise.
// If true is returned, then the associated system can be solved via the conjugate gradient method (if symmetric) or via at least gauss seidel method (definitely).
bool DIAGDOMINANT(CSRMatrix& STIFFNESSMATRIX);

// Checks if the square matrix M is tridiagonal
bool TRIDIAG(CSRMatrix& M);

// returns the upper left sxs matrix of M [ this function is used to determine if a matrix is symmetric positive definite], s is 1 starting
CSRMatrix submatrix(CSRMatrix& M, const int& s);

// Returns the submatrix obtained from M by deleting row i and column j [ i and  are 1 starting indices].
CSRMatrix submatrix1(CSRMatrix& M, const int& i, const int& j);
void mvProduct(CSRMatrix& M, double *x, double *store, const int& n);

// Calculate the determinant of the square matrix M
double det(CSRMatrix& M); // slow
double det1(CSRMatrix& M); // fast

#endif

