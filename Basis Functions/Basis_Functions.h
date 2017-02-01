/*===========================================================
 *
 *  File: Basis_Functions.h
 *  Content: Header file declaring class Basis_Functions
 *  that manages basis function generation and storage through
 *  static data members.
 *  Date: 04/06/2016 
 *  Author: Andrew Oldham
 *
 * **********************************************************/

#ifndef __BASIS_FUNCTIONS_H__
#define __BASIS_FUNCTIONS_H__

#include <vector>

//
// This Class is to be used as a container to hold the local basis function
// definitions defined on the reference element. 
// 
struct Basis_Functions
{
  public:
    // constructs the basis functions of degree Degree on the reference element and 
    // stores these polynomials in the private data member (as the first entry in the vector) OrderedBasisFunction_Degrees
    static void Generate_Basis_Functions(const int& Degree);

    std::vector< std::vector<double> > get_Basis_Functions(const int& Degree);

    static int maxDegree(); // returns the maximum degree of polynomial basis functions that exists on the current mesh
  private:
    static std::vector< std::vector< std::vector<double> > > BasisFunctions;
    static std::vector<int> OrderedBasisFunction_Degrees;

};
#endif



								        	

							 
	                                                  
                               



			       
 
 
                                                                                                                 
 
	           					                                                         

                                                                	
	 










            


 
