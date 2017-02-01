/*=====================================================================;
 *
 *  File: Bilinear_Form.cpp
 *  Content: Definitions of the bilinear form.
 *  Date: 12/06/2016
 *  Author: Andrew Oldham
 *
 ***********************************************************************/
#include "Bilinear_Form.h"
#include "PDE_Definition.h"
#include "../General Finite Element Class/Finite_Element.h"
#include "../Operators/Polynomial_Operators.h"
#include "../Quadratures/GQ_NodesWeights_Generation.h"
#include <cmath>
//


std::vector<double> SpatialVars = {1.0};
std::vector<double> D_Orders = {1.0};

// bilinear form evaluated over the finite element FE with ordered arguements basis function j and i (of the finite element)
double BilinearForm_VI(Finite_Element& FE, const int& j, const int& i, std::vector<double>& Point)
{ 
  std::vector<double> BF_j_Derivative = FE.Differentiate_BF(SpatialVars,D_Orders,j);  
  std::vector<double> BF_i_Derivative =  FE.Differentiate_BF(SpatialVars,D_Orders,i);  
  std::vector<double> BF_i = (FE.get_BFncs())[i-1];
  std::vector<double> BF_j = (FE.get_BFncs())[j-1]; 

  double mesh_Point = Eval( Point, (FE.get_AffineMap())[0] );

  double res = (a(mesh_Point)*Eval(Point, BF_j_Derivative)*Eval(Point, BF_i_Derivative)) + (b(mesh_Point)*Eval(Point, BF_j_Derivative)*Eval(Point, BF_i)) + 
	       (c(mesh_Point)*Eval(Point, BF_j)*Eval(Point, BF_i));

  return res;
}

// linear functional evaluated over the finite element FE, with elements basis function number i
// as the arguement for the linear functional.
double LinearFunctional_VI(Finite_Element& FE, const int& i, std::vector<double>& Point)
{     
  return (Eval( Point, (FE.get_BFncs())[i-1]))*f( Eval(Point, FE.get_AffineMap()[0]) );
}

// returns a value such that if the local_error_estimator returns a value greater than this, then the element is marked for refinement.
// initial_h is the initial vertex interval on the mesh, before any refinement has taken place [when the mesh is uniform].
// global tol refers to the total error the solution should obey to, measured in the H1 norm over the entire grid and bounded using
// a posteriori techniques. [ currently implemented for the 1D case]
double element_error_criteria(Mesh& COMPUTATIONALMESH, const int& FE_number, const double& initial_h, 
	         	             const double& global_tol, const int& initial_NoOfElements, const double& C_0)
{
  Finite_Element currFE = COMPUTATIONALMESH.find_FE(FE_number);
  double initNoOfElements = initial_NoOfElements, vertex_interval = (currFE.get_VertexPoints())[1][0] - (currFE.get_VertexPoints())[0][0],
	                   iter=0.0, returnval=0.0;
  bool continue_search = true;
  while (continue_search) {
    ++iter;
    double test_vertexinterval = initial_h/iter;
    if (fabs(test_vertexinterval-vertex_interval)<1e-12)
    continue_search=false; 
  }
  returnval = pow( (global_tol*C_0*2.0), 2.0)/(initNoOfElements*iter);
  return returnval;
}

// local error estimator based on a posteriori resudual estimates, to be used to
// obtain a bound on the error of our approximated solution on the finite element
// numbered (1 tstarting) FE_number in the mesh. Returns estimate of the square of the l2 norm of the local error
double local_error_estimator(Mesh& COMPUTATIONALMESH, double* coeffVals, const int& FE_number, const double& initial_h)
{
  Finite_Element currFE1 = COMPUTATIONALMESH.find_FE(FE_number);
  int bfdegree = currFE1.get_BFDegree();
  double local_ErrEstimate = 0.0;
  double degree = (double)bfdegree;
  double* ConnectivityArray_entries = COMPUTATIONALMESH.ConnArray_entryptr();
  int* ConnectivityArray_row_start = COMPUTATIONALMESH.ConnArray_rowptr();
  std::vector<double> scalars;
  std::vector< std::vector<double> > secondDerivatives;
  std::vector< std::vector<double> > firstDerivatives;
  std::vector<double> spatialvar = {1.0};
  std::vector<double> derivativeOrder = {2.0};
  std::vector<double> localU_h_2nd_deriv;
  std::vector<double> localU_h_1st_deriv;
  std::vector<double> localU_h;
  std::vector<double> localPoint;
  for (int j=ConnectivityArray_row_start[FE_number-1]-1; j<ConnectivityArray_row_start[FE_number]-1; ++j) 
  scalars.push_back(coeffVals[((int)(ConnectivityArray_entries[j]))-1]);

  for (int j=0; j<(bfdegree+1); ++j) { 
    std::vector<double> currsecondderiv = currFE1.Differentiate_BF(spatialvar, derivativeOrder, j+1);
    derivativeOrder = {1.0};
    std::vector<double> currfirstderiv = currFE1.Differentiate_BF(spatialvar, derivativeOrder, j+1);
    secondDerivatives.push_back(currsecondderiv);
    firstDerivatives.push_back(currfirstderiv);
    derivativeOrder = {2.0};
  }
  localU_h = LinearComb(scalars, currFE1.get_BFncs());
  localU_h_1st_deriv = LinearComb(scalars, firstDerivatives);
  localU_h_2nd_deriv = LinearComb(scalars, secondDerivatives);
  for (int j=0; j<GaussianQuadrature::Nodes.size(); ++j)  {
    localPoint = {(0.5*GaussianQuadrature::Nodes[j]) + 0.5};
    local_ErrEstimate+= GaussianQuadrature::Weights[j]*pow( (initial_h/sqrt(bfdegree*(bfdegree+1.0)))*( f( Eval(localPoint, (currFE1.get_AffineMap())[0]) ) - 
			  (-a(0.0)*Eval(localPoint, localU_h_2nd_deriv) +b(Eval(localPoint,currFE1.get_AffineMap()[0]))*Eval(localPoint,localU_h_1st_deriv) 
		           +c(Eval(localPoint,currFE1.get_AffineMap()[0]))*Eval(localPoint,localU_h)) ) , 2.0 );
  }
  local_ErrEstimate*=(0.5)*currFE1.detJ();

  return local_ErrEstimate;
 
}

// Weight Funciton used during test problem.
double weight_Fnc(const double& x_j_sub1, const double& x_j, const double& x)
{
  return (x_j - x)*(x - x_j_sub1);
}





	

