/*============================================================;
 *
 *  File: Mesh_Generator.h
 *  Content: Header File declaring functions for mesh 
 *  generation and refinement.
 *  Stiffness matrix and load vector assembly functuions
 *  are also declared. 
 *  Functions measuring the norm of the error of the
 *  approximated solution (energy and l2 norm) are declared.
 *  Date: 12/06/2016
 *  Author: Andrew Oldham
 *
 *************************************************************/

#ifndef __MESH_GENERATOR_H__
#define __MESH_GENERATOR_H__

#include "../Mesh/Mesh.h"
#include "Matrix.h"
#include <list>

// Generate the initial Mesh (using info from input file ...Msc/Dissertation/Input Docs/Computational Domain.dat )
Mesh GenerateMesh();

// Unique method of determining ordered pairs of integers (used for determining the maximum number of 
// nonzero entries in the global stiffness matrix). Assumes range of values that node1 and node2 
// can take is finite, which is a safe computational assumption. Hash function folds the finite
// discrete square into a 1D straight line.
int hashFunc(const int& node1, const int& node2, const int& unit_dist);

// Ordered Insertion
void orderedInsert(std::list<int>& Container, const int& val);

// Calculatethe bilinear form a(bf_j, bf_i) defined over the reference element FE
// for basis functions (1 starting) j and i  
double a(Finite_Element& FE, const int& j, const int& i);

// returns true if the stiffness matrix (STIFFNESSMATRIX) is symmetric positive definite, false otherwise.
bool SYMM_PD(CSRMatrix& STIFFNESSMATRIX);

// Function that return the stiffness matrix in CSR format
CSRMatrix GenerateStiffnessMatrix(Mesh& COMPUTATIONALMESH);

// Generale load vector
double* GenerateLoadVector(Mesh& COMPUTATIONALMESH);

// ImplemeNt boundary conditions (specified in input file Computational Domain.dat)
void BC_IMPOSITION(Mesh& COMPUTATIONALMESH, CSRMatrix& STIFFNESSMATRIX, double* LOADVEC); 

// Evaluate the solution at the point Point in the global mesh, given the coefficients
// of the solution in coeffs and the finite elements as part of the mesh A.
double EvaluateSolution(Mesh& A, double* coeffs, std::vector<double> Point);


// Returns a sharp upper bound on the energy norm of the error of the aproximated solution.
// Is only used for the test case in order to verify program.
double EnergyNormOfErr_Estimator(Mesh& A, double* coeffs, double* L2_PROJECTION_OF_f);

// Function to determine if mesh refinement was initially specified.
bool MeshRefinement();

// Refines the mesh using h-p refinement 
void RefineMesh(Mesh& COMPUTATIONALMESH, double* coeffVals, const double& initial_h, const int& initial_NoOfElements, 
		              const double& tol,const double& SmoothnesParam, const double& markerParam, const double& C_0, double* L2_PROJECTION_OF_f=0);

// refinementVect input arguement will list the numbered elements in the Mesh that are to be refined, 1 starting finite element numbering. 
void MarkElements(Mesh& COMPUTATIONALMESH, double* coeffVals, std::vector<int>& refinementVect, const double& markParam, double* L2_PROJECTION_OF_f);

// Gives a measure of smoothness of the approximate solution on the finite element FE (numbered element_index+1) in the mesh.
// depending on the return value, p refinement is applied if this returned value is sufficiently close to 1 [approximation is quite smooth on the element].
// h refinement is applied if the returned value is sufficiently close to 0 [ie approximation is not very smooth on the element]
double smoothness_measure(const int& element_index, Mesh& FE, double* coeffVals);

// Estimation of the linfinty Norm of the polynomial Poly on the interval [0,1]. ie supp { |Poly(x)| : x E (0,1) }
double linfinityNorm(std::vector<double>& Poly); 
                                                

// Uses h_refinement to refine the element numbered (i+1) in the mesh
void h_refinement(Mesh& A, const int& i, std::unordered_map<int,Finite_Element>& ElementContainer, Finite_Element& refineFE, int& newNoOfNodes, 
  int& refine_element_track, std::vector<double>& conn_entries, std::vector<int>& conn_col_no, std::vector<int>& conn_row_start, int& currNodeNo, int& currColNo,
   double& h_);

// Uses p refinement to refine the finite element numbered i+1 in the mesh.
void p_refinement(Mesh& A, const int& i, std::unordered_map<int,Finite_Element>& ElementContainer, Finite_Element& refineFE, int& newNoOfNodes, 
	int& refine_element_track, std::vector<double>& conn_entries, std::vector<int>& conn_col_no, int& currNodeNo, int& currColNo, double& h_);

// returns the coefficient values defining the linear combination of basis functuons that represents the L2 projection of f. Only used
// to verify program in the test case problem
double* L2_projection_f(Mesh& COMPUTATIONALMESH);


// Function that returns the actual (analytic) energy norm of the error induced from the FE solution. This function is custom built for the 1 test problem
// that is used to verify the adaptive FEM procedure. Only used in test case problem.
double analytic_energy_norm_Err(Mesh& COMPUTATIONALMESH, double* coeffVals);

// A posteriori error estimate of error in the H1 norm
double H1NormofErr_Est(Mesh& COMPUTATIONALMESH, double* coeffs, const double& C_0);

// returns problem dependent coercivity constant.
double get_Coercivity_const();
#endif
