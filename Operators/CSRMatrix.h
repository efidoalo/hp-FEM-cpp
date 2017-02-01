/*=========================================;
 *
 *  File: CSRMatrix.h 
 *  Content: Declarations for CSRMatrix class 
 *  used for sparse matrix
 *  storage
 *  Date: 13/07/2016
 *  Author: A.Oldham
 *
 ******************************************/

#ifndef __CSRMATRIX_H_INCLUDED
#define __CSRMATRIX_H_INCLUDED

// Struct for storing the connectivityarray in CSR format
struct CSRMatrix
{ 
  // Constructor
  CSRMatrix(double* entr, int* col, int* row, int NoR);
  
  // Default constructor
  CSRMatrix();
  
  // Copy Constructor
  CSRMatrix(const CSRMatrix& MATRIX);

  // Assignment Operator
  CSRMatrix& operator=( CSRMatrix CSRMATRIX_ASSIGN);

  void set_vars(double* entr, int* col, int* row, int NoR);

  double *matrix_entries;
  int *col_no;
  int *row_start;
  int NoOfRows;

  // destructor
  ~CSRMatrix();

};

#endif
