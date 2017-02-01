/*===========================================================
 *
 *  File: Basis_Functions.cpp
 *  Content: Source file defining the ordered basis functions
 *  for each of the finite elements.
 *  Date: 04/06/2016 
 *  Author: Andrew Oldham
 *
 * **********************************************************/

#include "Basis_Functions.h"
#include "../Operators/CSRMatrix.h"
#include "../Operators/Matrix_Algorithms.h"
#include <vector>
#include <iostream>

std::vector< std::vector< std::vector<double> > > Basis_Functions::BasisFunctions;
std::vector<int> Basis_Functions::OrderedBasisFunction_Degrees;

// Generate and store the set of basis functions [ defined on the reference element [0,1] ] that have degree BFDeg. 
// If set of basis functions corresponding to BFDeg is already supported, function does nothing.
void Basis_Functions::Generate_Basis_Functions(const int& BFDeg)
{ 
  bool new_request = true;
  for (int i=0; i<Basis_Functions::OrderedBasisFunction_Degrees.size(); ++i) {
    if (Basis_Functions::OrderedBasisFunction_Degrees[i]==BFDeg) 
    new_request = false; 
  }
  if (new_request) {
    std::vector< std::vector<double > > newBasisFunctions;
    double NodalInterval = 1.0/((double)(BFDeg)); 
    std::vector<double> Nodes; // excluding left hand node on the refrence element, 0.
    for (int j=1; j<BFDeg; ++j) 
    Nodes.push_back(j*NodalInterval);       
    Nodes.push_back(1.0);
    CSRMatrix A;
    A.row_start = new int[BFDeg+1];
    A.NoOfRows = BFDeg;
    A.matrix_entries = new double[BFDeg*BFDeg];
    A.col_no = new int[BFDeg*BFDeg];
    int iter=0;
    for (int j=0; j<BFDeg; ++j) {
      A.row_start[j] = (j*BFDeg)+1;
      for (int k=0; k<BFDeg; ++k) {
        A.matrix_entries[iter] = pow(Nodes[j],k+1); 
        A.col_no[iter] = k+1;
        ++iter;
      }
    }     
    A.row_start[BFDeg] = BFDeg*BFDeg + 1;     
    for (int i=0; i<BFDeg+1; ++i) {
      std::vector<double> b;
      if (i==0) {
        for (int j=0; j<BFDeg;++j)
        b.push_back(-1.0);
      }
      else {
        for (int j=1; j<BFDeg+1;++j) {
          if (i!=j)
          b.push_back(0);
	  else b.push_back(1.0);
        }
      }
      double* Coeffs = new double[BFDeg];   
      for (int j=0; j<BFDeg; ++j)
      Coeffs[j]=0.0;
      GaussSeidelMethod(A.matrix_entries, A.col_no, A.row_start, &b[0], Coeffs, BFDeg, 1e-8);
      //GaussianElimination(A, &b[0], Coeffs);
      std::vector<double> currBasisFunction;
      currBasisFunction.push_back(BFDeg);
      currBasisFunction.push_back(0.0);
      currBasisFunction.push_back(1.0);
      if (i==0)
      currBasisFunction.push_back(1.0);
      else
      currBasisFunction.push_back(0.0);
      for (int j=0; j<BFDeg; ++j)
      currBasisFunction.push_back(Coeffs[j]); 
      newBasisFunctions.push_back(currBasisFunction);
      delete[] Coeffs;	
    }
    BasisFunctions.push_back(newBasisFunctions);
    OrderedBasisFunction_Degrees.push_back(BFDeg);
    return;    
  }
  else return;
}

// returns the set of basis functions that all have degree const int& Degree [ defined on the reference element [0,1] ]
std::vector< std::vector<double> > Basis_Functions::get_Basis_Functions(const int& Degree)
{
  int BFIndex = -1;
  for (int i=0; i<OrderedBasisFunction_Degrees.size(); ++i) {
    if (OrderedBasisFunction_Degrees[i]==Degree)
    BFIndex=i;
  }
  if (BFIndex>=0) 
  return BasisFunctions[BFIndex];
  else {
    std::vector< std::vector<double> > NULLPoly =  { {0.0} };
    return NULLPoly;
  }
}
                
// returns the maximum degree of all stored basis functions.
int Basis_Functions::maxDegree()
{
  int maxdeg = -1;
  for (int i=0; i<Basis_Functions::OrderedBasisFunction_Degrees.size(); ++i) {
    if (maxdeg<OrderedBasisFunction_Degrees[i])
    maxdeg =OrderedBasisFunction_Degrees[i];
  }
  return maxdeg;
}



								        	

							 
	                                                  
                               



			       
 
 
                                                                                                                 
 
	           					                                                         

                                                                	
	 










            


 
