/*==========================================;
 *
 *  File: Matrix_Algorithms.cpp
 *  Content: Function/Routine definitions to invert 
 *  matrices dependent upon their form.
 *  File includes definitions of non-templated
 *  functions.
 *  Date: 16/06/2016
 *  Author: Andrew Oldham
 *
 *******************************************/

#include "CSRMatrix.h"
#include "Matrix_Algorithms.h"
#include <cmath>
#include <vector>

// Output the input matrix (that is representedin CSR format) to cmd terminal
void CSR_MatrixOutput(double* entries, int* col_no, int* row_start, const int& n, const int& m)
{  
  std::cout<<"\n\n";
  for (int j=0; j<n; ++j) {
    int row_startindex = row_start[j]-1, row_endindex=row_start[j+1]-2;
    int index=j+1;
    while (row_start[index]==0)
    ++index;
    row_endindex = row_start[index]-2;
    if (row_startindex == -1) {
      for (int i=0; i<m; ++i)
      std::cout<<0<<" ";
    }
    else {
      for (int i=0; i<m; ++i) {      
        bool hit = false;
        for (int iter=row_startindex; iter<=row_endindex; ++iter) {
          if (col_no[iter]==i+1) {
            hit = true;
            std::cout<<entries[iter]<<" ";
	  }
        }
        if (hit==false)
        std::cout<<0<<" ";
      }
    }
    std::cout<<"\n"; 
  } 
  std::cout<<"\n";
}

void GaussSeidelMethod(double* entries, int* col_no, int* row_start, double* b, 
double* result, const int& n, const double& tol)
{ 
  double error = 1.0, diagVal = 0.0, rowSum =0.0, currentError=0.0;
  int index = 0;
  for (int i=0; i<n; ++i)
  result[i] = 0.0;
  while (error>=tol) {    
    for (int i=0; i<n; ++i) {
      index = row_start[i]-1;
      while (col_no[index]!=i+1)
      ++index;
      diagVal = entries[index];
      index = row_start[i]-1;
      rowSum=0.0;
      while (index<row_start[i+1]-1) {
        if (col_no[index]!=i+1)
        rowSum+=entries[index]*result[col_no[index]-1];
	++index;
      }
      result[i] = (1.0/diagVal)*(b[i] - rowSum);
    }
    index = 0;
    error = 0.0;
    for (int i=0; i<n; ++i) {
      currentError = 0.0;
      while (index<row_start[i+1]-1) {
        currentError+=entries[index]*result[col_no[index]-1];
        ++index;
      }
      error+=pow(currentError - b[i],2.0);
    } 
    error = pow(error,0.5);
  }
}

// Store the result of the matrix-vector product A*x in the preallocated
// array store M is an nxn square matrix.
void mvProduct(CSRMatrix& M, double *x, double *store, const int& n)
{
  for (int i=0; i<M.NoOfRows; ++i)
  store[i] = 0.0;
  for (int i=0; i<n; ++i) {
    if (M.row_start[i]!=0) {      
      int endofRowIndex = i+1;
      if (M.row_start[endofRowIndex]==0) 
      ++endofRowIndex;
      for (int j=M.row_start[i]-1; j<M.row_start[endofRowIndex]-1; ++j) 
      store[i]+= M.matrix_entries[j]*x[M.col_no[j]-1];       
    }
  }
}

double dotproduct(double* x1, double* x2, const int& n)
{
  double val = 0.0;
  for (int i=0; i<n; ++i)
  val+=x1[i]*x2[i];
  return val;
}

// Conjugate Gradient Method [ x preallocated ]
void ConjugateGradientMethod(CSRMatrix& M, double* b, double* x, const double& tol)
{ 
  double *residual = new double[M.NoOfRows], *p = new double[M.NoOfRows], *product = new double[M.NoOfRows];
  
  double alpha=0.0, beta=0.0, error=0.0, temp=0.0;
  for (int i=0; i<M.NoOfRows; ++i) {
    x[i]=0.0;
    residual[i]=b[i];
    p[i]=b[i];
    error+=pow(b[i],2.0);
  }
  error = pow(error,0.5);
  while (error>=tol) {
    temp=0.0;
    beta=0.0;
    mvProduct(M, p, product, M.NoOfRows);
    temp = dotproduct(p, product, M.NoOfRows);
    alpha=pow(error,2.0)/temp;
    for (int i=0; i<M.NoOfRows; ++i) {
      x[i]+=alpha*p[i];
      residual[i]-= alpha*product[i];
    }
    beta=dotproduct(residual,residual, M.NoOfRows)/pow(error,2.0);
    error=0.0;
    for (int i=0; i<M.NoOfRows; ++i) {
     p[i] = residual[i] + beta*p[i];
     error+=pow(residual[i],2.0);
    }
    error=pow(error,0.5);
  }
  delete[] residual;
  delete[] p;
}

void ThomasAlgorithm(CSRMatrix& M, double* b, double* x)
{
  double* u = new double[M.NoOfRows-1];
  double* d = new double[M.NoOfRows];
  for (int i=0; i<M.NoOfRows-1; ++i) {
    if (i==0) {
      u[i] = M.matrix_entries[1]/M.matrix_entries[0];
      d[i] = b[i]/M.matrix_entries[0];
    }
    else {
      u[i] = M.matrix_entries[M.row_start[i]+1]/(M.matrix_entries[M.row_start[i]] - (M.matrix_entries[M.row_start[i]-1]*u[i-1] ));
      d[i] = (b[i] - (M.matrix_entries[M.row_start[i]-1]*d[i-1]))/(M.matrix_entries[M.row_start[i]] - (M.matrix_entries[M.row_start[i]-1]*u[i-1] ));
    }
  }
  x[M.NoOfRows-1] = b[M.NoOfRows-1];
  for (int i = M.NoOfRows-2; i>=0; --i) 
  x[i] = d[i] - u[i]*x[i+1]; 
}

// Elementary Row Operations
// interchage rows r1 and r2, where r1 < r2
void elementary_row_operation1(CSRMatrix& M, const int& r1, const int& r2)
{  
  int index=0;
  index = r1;
  while (M.row_start[index]==0)
  ++index;
  int length_r1 = M.row_start[index]-M.row_start[r1-1];
  index = r2;
  while (M.row_start[index]==0)
  ++index;
  int length_r2 = M.row_start[index]-M.row_start[r2-1];
  if (M.row_start[r1-1]==0)
  length_r1=0;
  if (M.row_start[r2-1]==0)
  length_r2=0;
  if ( (length_r1==0) && (length_r2==0) )
  return;
  if ( length_r1==0 ) {
    double tempval=0.0;
    int tempcol = 0;
    int index1=r1;
    while (M.row_start[index1]==0)
    ++index1;
    for (int i=0; i<M.row_start[index]-M.row_start[r2-1]; ++i) {
      tempval = M.matrix_entries[M.row_start[index]-2];
      tempcol = M.col_no[M.row_start[index]-2];
      for (int j = M.row_start[index]-3; j>=M.row_start[index1]-1;--j) {
        M.matrix_entries[j+1] = M.matrix_entries[j];
        M.col_no[j+1] = M.col_no[j];
      }
      M.matrix_entries[M.row_start[index1]-1] = tempval;   
      M.col_no[M.row_start[index1]-1] = tempcol;    
    }
    M.row_start[r1-1] = M.row_start[index1];
    for (int i=index1; i<r2-1; ++i) {
      if (M.row_start[i]!=0)
      M.row_start[i]+=length_r2;
    }
    M.row_start[r2-1]=0;
    return;
  }
  if (length_r2==0) {
    double tempval=0.0;
    int tempcol=0;
    int index1=r1;
    while (M.row_start[index1]==0)
    ++index1;
    for (int i=0; i<length_r1; ++i) {
      tempval = M.matrix_entries[M.row_start[r1-1]-1];
      tempcol = M.col_no[M.row_start[r1-1]-1];
      for (int j=M.row_start[r1-1]-1; j<M.row_start[index]-2; ++j) {
        M.matrix_entries[j] = M.matrix_entries[j+1];
	M.col_no[j] = M.col_no[j+1];
      }
      M.matrix_entries[M.row_start[index]-2]=tempval;
      M.col_no[M.row_start[index]-2]=tempcol;
    } 
    M.row_start[r1-1]=0;   
    for (int i=index1; i<r2-1; ++i) {
      if (M.row_start[i]!=0) 
      M.row_start[i]-=length_r1;
    }
    M.row_start[r2-1]=M.row_start[index]-length_r1;
    return;
  }
  if (length_r1<length_r2) {
    int posDelta = length_r2 - length_r1;
    std::vector<double> tempr2val;
    std::vector<int> tempr2col;
    int index1=r1;
    while (M.row_start[index1]==0)
    ++index1;
    for (int i=M.row_start[r2-1]-1; i<M.row_start[index]-1; ++i) {
      tempr2val.push_back(M.matrix_entries[i]);
      tempr2col.push_back(M.col_no[i]);
    }
    for (int i=M.row_start[index1]-2; i>=M.row_start[r1-1]-1; --i) {
      M.matrix_entries[M.row_start[index]-1-length_r1+(i-(M.row_start[r1-1]-1))] = M.matrix_entries[i];
      M.col_no[M.row_start[index]-1-length_r1+(i-(M.row_start[r1-1]-1))] = M.col_no[i];
    }
    for (int i=M.row_start[r2-1]-2; i>=M.row_start[index1]-1; --i) { 
      M.matrix_entries[M.row_start[index]-length_r1-M.row_start[r2-1]+i] = M.matrix_entries[i];
      M.col_no[M.row_start[index]-length_r1-M.row_start[r2-1]+i] = M.col_no[i];
    }
    for (int i=0; i<tempr2val.size(); ++i) {
      M.matrix_entries[M.row_start[r1-1]-1 + i] = tempr2val[i];
      M.col_no[M.row_start[r1-1]-1 + i] = tempr2col[i];
    }
    for (int i=r1; i<r2; ++i) {
      if (M.row_start[i]!=0)	    
      M.row_start[i]+=posDelta;
    }   
  }
  if (length_r1==length_r2) {
    double temp1=0;
    for (int i=0 ;i<length_r1; ++i) {
      temp1 = M.matrix_entries[M.row_start[r2-1]-1+i];
      M.matrix_entries[M.row_start[r2-1]-1+i] = M.matrix_entries[M.row_start[r1-1]-1+i];
      M.matrix_entries[M.row_start[r1-1]-1+i] = temp1;
      temp1 = M.col_no[M.row_start[r2-1]-1+i];
      M.col_no[M.row_start[r2-1]-1+i] = M.col_no[M.row_start[r1-1]-1+i];
      M.col_no[M.row_start[r1-1]-1+i] = temp1;
    }
  }
  if (length_r1>length_r2) {
    int posDelta = length_r2 - length_r1;
    std::vector<double> tempr1val;
    std::vector<int> tempr1col;
    int index1=r1;
    while (M.row_start[index1]==0)
    ++index1;
    for (int i=M.row_start[r1-1]-1; i<M.row_start[index1]-1; ++i) {
      tempr1val.push_back(M.matrix_entries[i]);
      tempr1col.push_back(M.col_no[i]);
    }
    for (int i=0; i<length_r2; ++i) {
      M.matrix_entries[M.row_start[r1-1]-1+i] = M.matrix_entries[M.row_start[r2-1]-1+i];
      M.col_no[M.row_start[r1-1]-1+i] = M.col_no[M.row_start[r2-1]-1+i];
    }
    for (int i = M.row_start[index1]-1; i<=M.row_start[r2-1]-2; ++i) {
      M.matrix_entries[M.row_start[r1-1]-1+length_r2 + i] = M.matrix_entries[i];
      M.col_no[M.row_start[r1-1]-1+length_r2 + i] = M.col_no[i];
    }
    for (int i = 0; i<length_r1; ++i) {
      M.matrix_entries[M.row_start[index]-1-length_r1 + i] = tempr1val[i];
      M.col_no[M.row_start[index]-1-length_r1 + i] = tempr1col[i];  
    }
    for (int i=r1; i<r2; ++i) {
      if (M.row_start[i]!=0)	    
      M.row_start[i]+=posDelta;
    }
  }
}
void elementary_row_operation1(double* b, const int& r1, const int& r2)
{ 
  double temp = b[r1-1];
  b[r1-1] = b[r2-1];
  b[r2-1] = temp;
}

// Multiply row r1 by scalar
void elementary_row_operation2(CSRMatrix& M, const int& r1, const double& scalar)
{
  if (M.row_start[r1-1]==0)
  return;
  int index=r1;
  while (M.row_start[index]==0)
  ++index;
  for (int i=M.row_start[r1-1]-1; i<M.row_start[index]-1; ++i) 
  M.matrix_entries[i]*=scalar; 
}
void elementary_row_operation2(double* b, const int& r1, const double& scalar)
{ 
  b[r1-1]*=scalar;
}

// Add scalar lots of row 1 onto row 2
void elementary_row_operation3(CSRMatrix& M, const int& r1, const int& r2, const double& scalar)
{
  int index=0;
  index = r1;
  while (M.row_start[index]==0)
  ++index;
  int length_r1 = M.row_start[index]-M.row_start[r1-1];
  index = r2;
  while (M.row_start[index]==0)
  ++index;
  int length_r2 = M.row_start[index]-M.row_start[r2-1];
  if (M.row_start[r1-1]==0)
  length_r1=0;
  if (M.row_start[r2-1]==0)
  length_r2=0;
  if (length_r1==0)
  return;
  if ( (length_r2==0) ) {
    double* tempentries = new double[M.row_start[M.NoOfRows]-1 + length_r1];
    int* tempcol_no = new int[M.row_start[M.NoOfRows]-1 + length_r1]; 
    for (int i=0; i<M.row_start[index]-1; ++i) {
      tempentries[i] = M.matrix_entries[i];
      tempcol_no[i] = M.col_no[i];
    }   
    int indexStart = M.row_start[index]-1;
    for (int i=0; i<length_r1; ++i) {
      tempentries[indexStart + i] = scalar*M.matrix_entries[M.row_start[r1-1]-1+i];
      tempcol_no[indexStart + i] = M.col_no[M.row_start[r1-1]-1+i];
    }
    indexStart+=length_r1;
    for (int i=M.row_start[index]-1; i<M.row_start[M.NoOfRows]-1; ++i) {
      tempentries[indexStart + i-(M.row_start[index]-1)] = M.matrix_entries[i];
      tempcol_no[indexStart + i-(M.row_start[index]-1)] = M.col_no[i];
    }
    M.row_start[r2-1] = M.row_start[index];
    for (int i=index; i<=M.NoOfRows; ++i) {
      if (M.row_start[i]!=0)
      M.row_start[i]+=length_r1;      
    }
    delete[] M.matrix_entries;
    delete[] M.col_no;
    M.matrix_entries = tempentries;
    M.col_no = tempcol_no;
    return;
  }
   
  int distinct_col=0;
  std::vector<int> distinctcoltrack;
  std::vector<double> correspondingvals;
  int index1 = r1;
  while (M.row_start[index1]==0)
  ++index1;
  for (int i=M.row_start[r1-1]-1; i<M.row_start[index1]-1; ++i) {
    bool match=false;
    for (int j=M.row_start[r2-1]-1; j<M.row_start[index]-1; ++j) {
      if (M.col_no[i]==M.col_no[j])
      match = true;
    }
    if (match==false) {
      ++distinct_col;
      distinctcoltrack.push_back(M.col_no[i]);
      correspondingvals.push_back(M.matrix_entries[i]);
    }
  }
  if (distinct_col!=0) {
    double* tempentries = new double[M.row_start[M.NoOfRows]-1+distinct_col];
    int * tempcol_no = new int[M.row_start[M.NoOfRows]-1+distinct_col];
    for (int i=0; i<M.row_start[M.NoOfRows]-1+distinct_col; ++i) {
      tempentries[i] = 0.0;
      tempcol_no[i] = 0;
    }
    std::vector<double> r2val;
    std::vector<double> r2col;
    int dcolindex = 0, r2colindex = 0;
    for (int i=M.row_start[r2-1]-1; i<M.row_start[index]-1; ++i) {
      r2val.push_back(M.matrix_entries[i]);
      r2col.push_back(M.col_no[i]);
    }
    for (int i=0; i<M.row_start[r2-1]-1; ++i) {
      tempentries[i] = M.matrix_entries[i];
      tempcol_no[i] = M.col_no[i];
    }    
    bool incrementdcolindex=false, incrementr2colindex=false, setendreached=false, endreached=false, setdcolend = false, dcolend = false;
    for (int i=0; i<M.row_start[index]-M.row_start[r2-1]+distinct_col; ++i) {    
      if ( (distinctcoltrack[dcolindex]<r2col[r2colindex]) && (dcolend==false) ){
        tempentries[M.row_start[r2-1]-1+i] = scalar*correspondingvals[dcolindex];
	tempcol_no[M.row_start[r2-1]-1+i] = distinctcoltrack[dcolindex];
	if (dcolindex<distinctcoltrack.size()-1)
	incrementdcolindex = true;
	else setdcolend=true;
      }
      if ( (distinctcoltrack[dcolindex]<r2col[r2colindex]) && (dcolend==true) ) {
        double additionalVal=0.0;
	for (int j=M.row_start[r1-1]-1; j<M.row_start[index1]-1; ++j) {
          if (M.col_no[j]==r2col[r2colindex])
          additionalVal = scalar*M.matrix_entries[j];
	}
	tempentries[M.row_start[r2-1]-1+i] = r2val[r2colindex] + additionalVal;
	tempcol_no[M.row_start[r2-1]-1+i] = r2col[r2colindex];
        incrementr2colindex=true;
      }
      if ((distinctcoltrack[dcolindex]>r2col[r2colindex]) && (endreached==false)) {
	double additionalVal=0.0;
	for (int j=M.row_start[r1-1]-1; j<M.row_start[index1]-1; ++j) {
          if (M.col_no[j]==r2col[r2colindex])
          additionalVal = scalar*M.matrix_entries[j];
	}
	tempentries[M.row_start[r2-1]-1+i] = r2val[r2colindex] + additionalVal;
	tempcol_no[M.row_start[r2-1]-1+i] = r2col[r2colindex];
	if (r2colindex==(r2col.size()-1)) {
	  setendreached=true; 
	}
        else incrementr2colindex=true;	
      }
      if ((distinctcoltrack[dcolindex]>r2col[r2colindex]) && (endreached==true)) {
        tempentries[M.row_start[r2-1]-1+i] = scalar*correspondingvals[dcolindex];
        tempcol_no[M.row_start[r2-1]-1+i] = distinctcoltrack[dcolindex];
	incrementdcolindex = true;
      }
      if (incrementdcolindex) {
        ++dcolindex;
	incrementdcolindex=false;
      }
      if (incrementr2colindex) {
        ++r2colindex;
	incrementr2colindex=false;
      }
      if (setendreached==true)
      endreached=true;
      if (setdcolend==true)
      dcolend = true;
    }
    for (int i=M.row_start[index]-1; i<M.row_start[M.NoOfRows]-1; ++i) {
      tempentries[distinct_col+i] = M.matrix_entries[i];
      tempcol_no[distinct_col+i] = M.col_no[i];
    }
    for (int i=index; i<=M.NoOfRows; ++i) {
      if (M.row_start[i]!=0)
      M.row_start[i]+=distinct_col;
    }
    delete[] M.matrix_entries;
    delete[] M.col_no;
    M.matrix_entries = tempentries;
    M.col_no = tempcol_no;
  } 
  else {
    for (int i=M.row_start[r2-1]-1; i<M.row_start[index]-1; ++i) {
      double additionalVal = 0.0;
      for (int j=M.row_start[r1-1]-1; j<M.row_start[index1]-1; ++j) {
        if (M.col_no[j]==M.col_no[i])
        additionalVal = scalar*M.matrix_entries[j];
      }
      M.matrix_entries[i]+=additionalVal;
    }
  } 
}
void elementary_row_operation3(double* b, const int& r1, const int& r2, const double& scalar)
{ 
  b[r2-1]+=scalar*b[r1-1];
}

// Redue the matrix M which may include zero entries (and zero rows)  in its matrix_entries and row_start members to CSR format.
void CSR_REDUCE(CSRMatrix& MATRIX)
{ 
  int zeroCount = 0, iter=0, temp=0;
  bool pastRow = false;
  for (int i=0; i<MATRIX.NoOfRows; ++i) {
    int rowstart = iter+1;
    temp=0;
    if (MATRIX.row_start[i]!=0) {
      int index = i+1;
      while (MATRIX.row_start[index]==0)
      ++index;
      if (pastRow) {
        while ( (MATRIX.matrix_entries[MATRIX.row_start[i]-1+temp]==0) || (fabs(MATRIX.matrix_entries[MATRIX.row_start[i]-1+temp])<1e-14) )
        ++temp;   
      }
      pastRow = false;
      for (int j=MATRIX.row_start[i]-1+temp; j<MATRIX.row_start[index]-1; ++j) {
        while ((MATRIX.matrix_entries[j]==0) || (fabs(MATRIX.matrix_entries[j])<1e-14)) {
          ++j;
          ++zeroCount; 
        }
        if (j<MATRIX.row_start[index]-1) {
          MATRIX.matrix_entries[j-zeroCount] = MATRIX.matrix_entries[j]; 
          MATRIX.col_no[j-zeroCount] = MATRIX.col_no[j];
          ++iter;
        }
        else
        pastRow=true;
      }
      MATRIX.row_start[i] = rowstart;
      if ((pastRow==true) && ((iter+1)==rowstart))
      MATRIX.row_start[i] = 0;
    }
  }
  MATRIX.row_start[MATRIX.NoOfRows] = MATRIX.row_start[MATRIX.NoOfRows] - zeroCount; 
}

// Guassian Elimination (worst case scenario) [M IS ASSUMED TO BE A SQUARE MATRIX], x is pre-allocated and initialized to zeros
void GaussianElimination(CSRMatrix& M, double* b, double* x)
{
  int zeroindex = 0;
  for (int i=0; i<M.NoOfRows; ++i) {
    int nonzero_rowindex = -1; 
    for (int j=zeroindex; j<M.NoOfRows; ++j) {
      if (M.row_start[j]!=0) {
        if (M.col_no[M.row_start[j]-1]==i+1)
        nonzero_rowindex=j;
        break;
      }
    }
    if (nonzero_rowindex>=0) {
      elementary_row_operation1(M,zeroindex+1,nonzero_rowindex+1);
      elementary_row_operation1(b,zeroindex+1,nonzero_rowindex+1);
      double scalar = (1.0)/(M.matrix_entries[M.row_start[zeroindex]-1]);
      elementary_row_operation2(M,zeroindex+1,scalar);
      elementary_row_operation2(b,zeroindex+1,scalar);
      for (int j=0; j<zeroindex; ++j) {
	if (M.row_start[j]!=0) {
	  int indexj = j+1;
	  while (M.row_start[indexj]==0)
          ++indexj;
          for (int k=M.row_start[j]-1; k<M.row_start[indexj]-1; ++k) {
            if (M.col_no[k]==i+1) {
	      scalar = (-1.0)*M.matrix_entries[k];
              elementary_row_operation3(M,zeroindex+1,j+1,scalar);
	      elementary_row_operation3(b,zeroindex+1,j+1,scalar);
	      M.matrix_entries[k] = 0;
	      CSR_REDUCE(M);
	    }
	  }
	}
      }
      for (int j=zeroindex+1; j<M.NoOfRows; ++j) {
	if (M.row_start[j]!=0) {
          if (M.col_no[M.row_start[j]-1]==i+1) {
	    scalar=(-1.0)*M.matrix_entries[M.row_start[j]-1];
            elementary_row_operation3(M,zeroindex+1,j+1,scalar);
	    elementary_row_operation3(b,zeroindex+1,j+1,scalar);
	    M.matrix_entries[M.row_start[j]-1]=0;
	    CSR_REDUCE(M);
          }
        }
      }
      ++zeroindex;
    }
  }  
  // Back substitution
  for (int i=0; i<M.NoOfRows; ++i) {
    if (M.row_start[i]!=0) {
    x[M.col_no[M.row_start[i]-1]-1] = b[i];
    }
  }
}

CSRMatrix CSR_INVERSION(double* matrix_entries, int* col_no, int* row_start, const int& n)
{
  if (n==1) {
    CSRMatrix Inverse;
    Inverse.matrix_entries = new double;
    Inverse.col_no = new int;
    Inverse.row_start = new int[2];
    Inverse.NoOfRows = 1;
    Inverse.matrix_entries[0] = 1.0/matrix_entries[0];
    Inverse.col_no[0] = 1;
    Inverse.row_start[0] = 1;
    Inverse.row_start[1] = 2;
    return Inverse;
  }
  if (n==2) {
    CSRMatrix Inverse;
    std::vector<double> entries;
    std::vector<int> col_no_vals;
    std::vector<int> row_startval;
    double arg1=0,arg2=0;
    bool arg1ZERO=false, arg2ZERO=false;
    if ((col_no[row_start[0]-1]==1) && (col_no[row_start[2]-2]==2))
    arg1 = matrix_entries[row_start[0]-1]*matrix_entries[row_start[2]-2];
    if ((col_no[row_start[1]-2]==2) && (col_no[row_start[1]-1]==1) ) 
    arg2 = matrix_entries[row_start[1]-2]*matrix_entries[row_start[1]-1];
    double detMatrix = arg1 - arg2; 
    if (arg1!=0) {
    entries.push_back((1.0/detMatrix)*matrix_entries[row_start[2]-2]);
    col_no_vals.push_back(1);
    row_startval.push_back(1);
    }
    if (arg2!=0) {
    entries.push_back(-1.0*(1.0/detMatrix)*matrix_entries[row_start[1]-2]);
    col_no_vals.push_back(2);
    entries.push_back(-1.0*(1.0/detMatrix)*matrix_entries[row_start[1]-1]);
    col_no_vals.push_back(1);
    row_startval.push_back(entries.size());
    }
    if (arg1!=0) {
    entries.push_back((1.0/detMatrix)*matrix_entries[row_start[0]-1]);
    col_no_vals.push_back(2);
    }
    row_startval.push_back(entries.size()+1);
    Inverse.matrix_entries = new double[entries.size()];
    Inverse.col_no = new int[col_no_vals.size()];
    Inverse.row_start = new int[row_startval.size()];
    Inverse.row_start = new int[row_startval.size()-1];
    for (int i=0; i<entries.size(); ++i) {
      Inverse.matrix_entries[i] = entries[i];
      Inverse.col_no[i] = col_no_vals[i];
    }
    for (int i=0; i<row_startval.size(); ++i)
    Inverse.row_start[i] = row_startval[i];
    Inverse.NoOfRows = n;
    return Inverse;
  }
  // Otherwise find iverse matrix via Gaussian Elimination
  CSRMatrix Inverse;
  Inverse.matrix_entries = new double[n];
  Inverse.col_no = new int[n];
  Inverse.row_start = new int[n+1];
  Inverse.row_start[n] = n+1;
  Inverse.NoOfRows = n;
  CSRMatrix M;
  M.matrix_entries = matrix_entries;
  M.col_no = col_no;
  M.row_start = M.row_start;
  M.NoOfRows = n;
  for (int i=0; i<n; ++i) {
    Inverse.matrix_entries[i] = 1.0;
    Inverse.col_no[i] = i+1;
    Inverse.row_start[i] = i+1;
  }
  int zeroindex = 0;
  for (int i=0; i<n; ++i) {
    int nonzero_rowindex = -1; 
    for (int j=zeroindex; j<M.NoOfRows; ++j) {
      if (M.row_start[j]!=0) {
        if (M.col_no[M.row_start[j]-1]==i+1)
        nonzero_rowindex=j;
        break;
      }
    }
    if (nonzero_rowindex>=0) {
      elementary_row_operation1(M,zeroindex+1,nonzero_rowindex+1);
      elementary_row_operation1(Inverse,zeroindex+1,nonzero_rowindex+1);
      double scalar = (1.0)/(M.matrix_entries[M.row_start[zeroindex]-1]);
      elementary_row_operation2(M,zeroindex+1,scalar);
      elementary_row_operation2(Inverse,zeroindex+1,scalar);
      for (int j=0; j<zeroindex; ++j) {
	if (M.row_start[j]!=0) {
	  int indexj = j+1;
	  while (M.row_start[indexj]==0)
          ++indexj;
          for (int k=M.row_start[j]-1; k<M.row_start[indexj]-1; ++k) {
            if (M.col_no[k]==i+1) {
	      scalar = (-1.0)*M.matrix_entries[k];
              elementary_row_operation3(M,zeroindex+1,j+1,scalar);
	      elementary_row_operation3(Inverse,zeroindex+1,j+1,scalar);
	      M.matrix_entries[k] = 0;
	      CSR_REDUCE(M);
	      CSR_REDUCE(Inverse);
	    }
	  }
	}
      }
      for (int j=zeroindex+1; j<M.NoOfRows; ++j) {
	if (M.row_start[j]!=0) {
          if (M.col_no[M.row_start[j]-1]==i+1) {
	    scalar=(-1.0)*M.matrix_entries[M.row_start[j]-1];
            elementary_row_operation3(M,zeroindex+1,j+1,scalar);
	    elementary_row_operation3(Inverse,zeroindex+1,j+1,scalar);
	    M.matrix_entries[M.row_start[j]-1]=0;
	    CSR_REDUCE(M);
	    CSR_REDUCE(Inverse);
          }
        }
      }
      ++zeroindex;
    }
  }  
  // Backwards step
  double scalar = 0.0;
  for (int i=M.NoOfRows-1; i>=0; --i) {
    for (int j=i-1; j>=0; --j) {
      int postrowindex = j+1;
      while (M.row_start[postrowindex]==0)
      ++postrowindex;
      for (int k=M.row_start[j]-1; k<M.row_start[postrowindex]-1; ++k) {
        if (M.col_no[k]==i+1) {
          scalar = (-1.0)*M.matrix_entries[k];
	  elementary_row_operation3(M,i+1,j+1,scalar);
	  elementary_row_operation3(Inverse,i+1,j+1,scalar);
	  CSR_REDUCE(M);
	  CSR_REDUCE(Inverse);
	}
      }
    }
  }
  return Inverse;
}

// calculates transpose matrices
CSRMatrix CSR_TRANSPOSE(double* matrix_entries, int* col_no, int* row_start, const int& n)
{ 
  if (n==1) { 
    CSRMatrix Copy;
    Copy.matrix_entries = new double[row_start[n]-1];
    Copy.col_no = new int[row_start[n]-1]; 
    Copy.row_start = new int[n+1];
    Copy.NoOfRows = n;
    for (int i=0; i<row_start[n]-1; ++i) {
      Copy.matrix_entries[i] = matrix_entries[i];
      Copy.col_no[i] = col_no[i];
    }
    for (int i=0; i<n+1; ++i)
    Copy.row_start[i] = row_start[i];
    return Copy;
  }
 
  CSRMatrix Transpose;
  Transpose.matrix_entries = new double[row_start[n]-1];
  Transpose.col_no = new int[row_start[n] - 1];
  Transpose.row_start = new int[n+1];
  Transpose.NoOfRows = n;
  int iter=0;
  for (int i=0; i<n; ++i) {
    Transpose.row_start[i] = iter+1;
    for (int j=0; j<n; ++j) {
      for (int k=row_start[j]; k<row_start[j+1]; ++k) {
        if (col_no[k-1]==i+1) {
          Transpose.matrix_entries[iter] = matrix_entries[k-1];
          Transpose.col_no[iter] = j+1;
	  ++iter;
	}
      }
    }
  }
  Transpose.row_start[n] = iter+1;
  return Transpose;
}

// Check the input matrix for diagonal dominance
bool DIAGDOMINANT(CSRMatrix& STIFFNESSMATRIX)
{
  bool returnval = true;
  for (int i=0; i<STIFFNESSMATRIX.NoOfRows; ++i) {
    double disc_Radius = 0.0, DIAGVAL=0.0;
    for (int j=STIFFNESSMATRIX.row_start[i]-1; j<STIFFNESSMATRIX.row_start[i+1]-1; ++j) {
       if (STIFFNESSMATRIX.col_no[j]!=i+1)
       disc_Radius+=fabs(STIFFNESSMATRIX.matrix_entries[j]);
       else DIAGVAL = fabs(STIFFNESSMATRIX.matrix_entries[j]);
    }
    if (DIAGVAL-disc_Radius<0)
    return false;
  }
  return true;
}

// return upper left sxs sub matrix
CSRMatrix submatrix(CSRMatrix& M, const int& s)
{ 
  int length=0,iter=0;
  for (int i=0; i<s; ++i) {
    for (int j=M.row_start[i]-1; j<M.row_start[i+1]-1; ++j) {
      if (M.col_no[j]<=s) 
      ++length;
      else break;      
    }
  }
  CSRMatrix UpperLeft;
  UpperLeft.matrix_entries = new double[length];
  UpperLeft.col_no = new int[length];
  UpperLeft.row_start = new int[s+1];
  UpperLeft.row_start[s] = length+1;
  UpperLeft.NoOfRows = s;
  for (int i=0; i<s; ++i) {
    UpperLeft.row_start[i] = iter+1;
    for (int j=M.row_start[i]-1; j<M.row_start[i+1]-1; ++j) {
      if (M.col_no[j]<=s) {
        UpperLeft.matrix_entries[iter] = M.matrix_entries[j];
        UpperLeft.col_no[iter] = M.col_no[j];
        ++iter;
      }
      else break;      
    }
  }
  return UpperLeft;
}

// returns submatrix resulting from deleting row i and column j from the input matrix
// what if sub matrix is the zero matrix?
CSRMatrix submatrix1(CSRMatrix& M, const int& i, const int& j)
{
  CSRMatrix SUBMATRIX;
  bool subtract = false;
  for (int k=M.row_start[i-1]-1; k<M.row_start[i]-1; ++k) {
    if (M.col_no[k]==j)
    subtract = true;
  }
  int NoOfEntries_jcol = 0;
  for (int k=0; k<M.row_start[M.NoOfRows]-1; ++k) {
    if (M.col_no[k]==j) 
    ++NoOfEntries_jcol;
  }
  int sublength = NoOfEntries_jcol + M.row_start[i]-M.row_start[i-1];
  if (subtract)
  --sublength;

  if ((M.row_start[M.NoOfRows]-1-sublength)==0) {
    // return the zero matrix
    SUBMATRIX.matrix_entries = new double;
    SUBMATRIX.col_no = new int;
    SUBMATRIX.row_start = new int[M.NoOfRows];
    SUBMATRIX.matrix_entries[0] = 0;
    SUBMATRIX.col_no[0]= 1;
    for (int i=0; i<M.NoOfRows-1; ++i)
    SUBMATRIX.row_start[i]= 0;
    SUBMATRIX.row_start[M.NoOfRows-1] = 1;
    SUBMATRIX.NoOfRows = M.NoOfRows-1;
    return SUBMATRIX;
  }
  else {
    SUBMATRIX.matrix_entries = new double[M.row_start[M.NoOfRows]-1-sublength];
    SUBMATRIX.col_no = new int[M.row_start[M.NoOfRows]-1-sublength];
    SUBMATRIX.row_start = new int[M.NoOfRows];
    SUBMATRIX.NoOfRows = M.NoOfRows-1;
    SUBMATRIX.row_start[M.NoOfRows-1]=M.row_start[M.NoOfRows]-sublength;
    int iter =0,rowtrack=0,zerotrack=0;
    for (int k=0; k<M.NoOfRows; ++k) {
      zerotrack = iter+1;
      SUBMATRIX.row_start[rowtrack] = iter+1;
      if (k+1!=i) {
        for (int m=M.row_start[k]-1; m<M.row_start[k+1]-1; ++m) {
          if (M.col_no[m]!=j) {
            SUBMATRIX.matrix_entries[iter] = M.matrix_entries[m];
            if (M.col_no[m]<j)
	    SUBMATRIX.col_no[iter] = M.col_no[m];
	    else SUBMATRIX.col_no[iter] = M.col_no[m]-1;
	    ++iter;
	  }
        }
        if (iter<zerotrack)
        SUBMATRIX.row_start[rowtrack]=0;
        ++rowtrack;
      }
    }
    return SUBMATRIX;
  }
}

// determinant  [slow]
double det(CSRMatrix& M)
{
  for (int i=0; i<M.NoOfRows; ++i) {
    if (M.row_start[i]==0)
    return 0;
  }
  if (M.NoOfRows==1) 
  return M.matrix_entries[0];  
  if (M.NoOfRows>1) {
    double DET=0.0;
    for (int j=0; j<M.row_start[1]-1; ++j) {
     CSRMatrix currSUB = submatrix1(M, 1,  M.col_no[j]); 
     double detcurrSUB = det(currSUB);
     DET+=pow(-1.0,j+2)*M.matrix_entries[j]*detcurrSUB;
    }
    return DET;
  }
}

// determinant [fast]
double det1(CSRMatrix& M)
{ 
  for (int i=0; i<M.NoOfRows; ++i) {
    bool action_required = false;
    for (int j=i; j<M.NoOfRows; ++j) {
      if ((M.row_start[j]!=0) && (M.col_no[M.row_start[j]-1]==i+1)) {
	if (j==i) {
          action_required = true;
          break;
	}
	else { 
          elementary_row_operation3(M,j+1,i+1,1);
	  CSR_REDUCE(M);
	  action_required = true;
	  break;
	}
      }
    }
    if (action_required) {
      double diagVal = M.matrix_entries[M.row_start[i]-1];
      for (int j=i+1; j<M.NoOfRows; ++j) {
        if ((M.row_start[j]!=0) && (M.col_no[M.row_start[j]-1]==i+1)) {
	  double scalar = (-1.0*M.matrix_entries[M.row_start[j]-1])/diagVal;
          elementary_row_operation3(M,i+1,j+1,scalar);
	  M.matrix_entries[M.row_start[j]-1] = 0;
	  CSR_REDUCE(M);	  
        }
      }   
    }
  }
  CSR_REDUCE(M);
  double DET=1.0;
  for (int i=0; i<M.NoOfRows; ++i) {
    if (M.col_no[M.row_start[i]-1]>i+1) 
    return 0;
    else DET*=M.matrix_entries[M.row_start[i]-1];
  }
  return DET;
}


// determine whether input matrix is tridiagonal
bool TRIDIAG(CSRMatrix& M)
{
  for (int i=0; i<M.NoOfRows; ++i) {
    if (i==0) {
      if ( (M.col_no[0]!=1) || (M.col_no[1]!=2) || (M.row_start[1]!=3) )
      return false;      
    }
    if ( (i>0) && (i<M.NoOfRows-1) ) {
      if ( (M.col_no[M.row_start[i]-1]!=i) || (M.col_no[M.row_start[i]]!=i+1) || (M.col_no[M.row_start[i]+1]!=i+2) || ((M.row_start[i+1] - M.row_start[i])!=3) )
      return false;
    }
    if (i==M.NoOfRows-1) {
      if ( (M.col_no[M.row_start[i]-1]!=i) || (M.col_no[M.row_start[i]]!=i+1) || ((M.row_start[i+1] - M.row_start[i])!=2) )
      return false;
    }
  }
  return true;
}










