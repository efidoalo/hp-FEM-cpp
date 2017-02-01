/*=========================================;
 *
 *  File: CSRMatrix.cpp 
 *  Content: Definitions for CSRMatrix class 
 *  used for sparse matrix
 *  storage
 *  Date: 13/07/2016
 *  Author: A.Oldham
 *
 ******************************************/
#include "CSRMatrix.h"
// Standard Constructor
CSRMatrix::CSRMatrix(double* entr, int* col, int* row, int NoR)
: matrix_entries(entr), col_no(col), row_start(row), NoOfRows(NoR) { }

// Default Destructor
CSRMatrix::CSRMatrix()
: matrix_entries(0), col_no(0), row_start(0), NoOfRows(0) { }

// Copy Constructor
CSRMatrix::CSRMatrix(const CSRMatrix& MATRIX)
{ 
  matrix_entries = new double[MATRIX.row_start[MATRIX.NoOfRows]-1];
  col_no = new int[MATRIX.row_start[MATRIX.NoOfRows]-1];
  row_start = new int[MATRIX.NoOfRows+1];
  NoOfRows = MATRIX.NoOfRows;
  for (int i=0; i<MATRIX.row_start[MATRIX.NoOfRows]-1; ++i) {
    matrix_entries[i] = MATRIX.matrix_entries[i];
    col_no[i] = MATRIX.col_no[i];
  }
  for (int i=0; i<NoOfRows+1; ++i)
  row_start[i] = MATRIX.row_start[i];
}

// Assignment Operator
CSRMatrix& CSRMatrix::operator=( CSRMatrix CSRMATRIX_ASSIGN)
{
  matrix_entries = new double[CSRMATRIX_ASSIGN.row_start[CSRMATRIX_ASSIGN.NoOfRows]-1];
  col_no = new int[CSRMATRIX_ASSIGN.row_start[CSRMATRIX_ASSIGN.NoOfRows]-1];
  row_start = new int[CSRMATRIX_ASSIGN.NoOfRows+1];
  NoOfRows = CSRMATRIX_ASSIGN.NoOfRows;
  for (int i=0; i<CSRMATRIX_ASSIGN.row_start[CSRMATRIX_ASSIGN.NoOfRows]-1; ++i) {
    matrix_entries[i] = CSRMATRIX_ASSIGN.matrix_entries[i];
    col_no[i] = CSRMATRIX_ASSIGN.col_no[i];
  }
  for (int i=0; i<NoOfRows+1; ++i)
  row_start[i] = CSRMATRIX_ASSIGN.row_start[i];
   
  return *this;
}

void CSRMatrix::set_vars(double* entr, int* col, int* row, int NoR)
{ 
  matrix_entries = entr;
  col_no = col;
  row_start = row;
  NoOfRows =  NoR;
}

//Destructor
CSRMatrix::~CSRMatrix()
{
  delete[] matrix_entries;
  delete[] col_no;
  delete[] row_start; 
}
