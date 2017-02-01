/*===========================================================
 *
 *  File: Finite_Element.cpp
 *  Content: Definitions for general finite element class.
 *  Date: 08/06/2016 
 *  Author: Andrew Oldham
 *
 * **********************************************************/


#include "Finite_Element.h"
#include "../Operators/CSRMatrix.h"
#include "../Operators/Matrix_Algorithms.h"
#include "../Operators/Polynomial_Operators.h"
#include "../Basis Functions/Basis_Functions.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>

// Class Constructor
Finite_Element::Finite_Element(const int& FE_TypeSpecifier, 
		const int& BF_DEGREE, 
	        const std::vector< std::vector<double> >& Vertex_Points, 
		const int& Dim)
: FEType(FE_TypeSpecifier), BFDegree(BF_DEGREE), Dimension(Dim)
{ 
  for (int i1=0; i1<Vertex_Points.size();++i1)
  VertexPoints.push_back(Vertex_Points[i1]);

  // Define the appropriate Affine Mapping from the reference element on to this particular finite element
  switch (FE_TypeSpecifier)
  {
    case 0:
      { 
	double point1 = (Vertex_Points[0])[0];
	double point2 = (Vertex_Points[1])[0];
        std::vector<double> AffineMap_1DLine = {1.0, 0.0, 1.0, Vertex_Points[0][0], Vertex_Points[1][0] - Vertex_Points[0][0]};      
	std::vector<double> InverseAffineMap_1DLine = {1.0, 0.0, 1.0, ((-1.0*(Vertex_Points[0][0]))/(Vertex_Points[1][0]-Vertex_Points[0][0])),
		                                                      (1.0/(Vertex_Points[1][0]-Vertex_Points[0][0]))};
	 
	InverseAffineMap = {InverseAffineMap_1DLine};//};
	AffineMap = {AffineMap_1DLine}; 
    
      }
      break;

    case 1:
      {
        std::vector<double> AffineMap_2DTriangle_xComponent = {1.0, 0.0, 2.0, Vertex_Points[0][0], Vertex_Points[2][0] - Vertex_Points[0][0],
		                                                                                   Vertex_Points[1][0] - Vertex_Points[0][0]};
	std::vector<double> AffineMap_2DTriangle_yComponent = {1.0, 0.0, 2.0, Vertex_Points[0][1], Vertex_Points[2][1] - Vertex_Points[0][1],
		                                                                                   Vertex_Points[1][1] - Vertex_Points[0][1]};
	AffineMap = {AffineMap_2DTriangle_xComponent, AffineMap_2DTriangle_yComponent};
      }
      break;

    case 2:
      {
        std::vector<double> AffineMap_2DParralelogram_xComponent = {1.0, 0.0, 2.0, Vertex_Points[0][0], Vertex_Points[3][0] - Vertex_Points[0][0],
		                                                                                        Vertex_Points[1][0] - Vertex_Points[0][0]};
	std::vector<double> AffineMap_2DParralelogram_yComponent = {1.0, 0.0, 2.0, Vertex_Points[0][1], Vertex_Points[3][1] - Vertex_Points[0][1],
		                                                                                        Vertex_Points[1][1] - Vertex_Points[0][1]};
	AffineMap = {AffineMap_2DParralelogram_xComponent, AffineMap_2DParralelogram_yComponent};
      }
      break;

    case 3:
      {
        std::vector<double> AffineMap_3DTetrahedron_xComponent = {1.0, 0.0, 3.0, Vertex_Points[0][0], Vertex_Points[3][0] - Vertex_Points[0][0],
		                                                                                      Vertex_Points[2][0] - Vertex_Points[0][0],
												      Vertex_Points[1][0] - Vertex_Points[0][0]};
	std::vector<double> AffineMap_3DTetrahedron_yComponent = {1.0, 0.0, 3.0, Vertex_Points[0][1], Vertex_Points[3][1] - Vertex_Points[0][1],
		                                                                                      Vertex_Points[2][1] - Vertex_Points[0][1],
												      Vertex_Points[1][1] - Vertex_Points[0][1]}; 
	std::vector<double> AffineMap_3DTetrahedron_zComponent = {1.0, 0.0, 3.0, Vertex_Points[0][2], Vertex_Points[3][2] - Vertex_Points[0][2],
		                                                                                      Vertex_Points[2][2] - Vertex_Points[0][2],
												      Vertex_Points[1][2] - Vertex_Points[0][2]};
	AffineMap = {AffineMap_3DTetrahedron_xComponent, AffineMap_3DTetrahedron_yComponent, AffineMap_3DTetrahedron_zComponent};
      }
      break;

    defult: 
      {
	std::cout<<"\nError, Type Specifier in Constructor of Finite_Element class not supported (must be integer value of 0,1,2, or 3).\n";
	exit(EXIT_FAILURE);
      }
      break;
  }
  
  Basis_Functions BF; 
  BF.Generate_Basis_Functions(BF_DEGREE);
  BasisFunctions = BF.get_Basis_Functions(BF_DEGREE);

}

// Copy Constructor
Finite_Element::Finite_Element(const Finite_Element& COPY)
{
  Dimension = COPY.get_Dim(); // Spatial dimension of the finite elment. 
  FEType = COPY.get_FETYPE();    // 0=Straight Line, 1=2D Triangle, 2=2D Parallelogram, 3=3D Tetrahedron
  BFDegree = COPY.get_BFDegree(); 
  AffineMap = COPY.get_AffineMap();
  InverseAffineMap = COPY.get_InverseAffineMap();
  BasisFunctions = COPY.get_BFncs(); 
  VertexPoints = COPY.get_VertexPoints();
}

// Assignment Operator
Finite_Element& Finite_Element::operator=(Finite_Element FE_ASSIGN)
{
  Dimension = FE_ASSIGN.get_Dim(); // Spatial dimension of the finite elment. 
  FEType = FE_ASSIGN.get_FETYPE();    // 0=Straight Line, 1=2D Triangle, 2=2D Parallelogram, 3=3D Tetrahedron
  BFDegree = FE_ASSIGN.get_BFDegree(); 
  AffineMap = FE_ASSIGN.get_AffineMap();
  InverseAffineMap = FE_ASSIGN.get_InverseAffineMap();
  BasisFunctions = FE_ASSIGN.get_BFncs(); 
  VertexPoints = FE_ASSIGN.get_VertexPoints();
  return *this;
}

// Defualt Constructor
Finite_Element::Finite_Element()
: FEType(-1), BFDegree(-1), Dimension(-1) { }

// compute and return jacobian of reference elem ---> Mesh Affine Transformation.
CSRMatrix Finite_Element::Jacobian()
{
  CSRMatrix JACOBIAN;
  std::vector<double> entries;
  std::vector<int> col_notrack;
  std::vector<int> row_starttrack;
  int entryindex = 0;
  for (int j=0; j<Dimension; ++j) {
    row_starttrack.push_back(entryindex+1);
    for (int i=0; i<Dimension; ++i) {
      int order=1, var=i+1;	    
      std::vector<double> diffVal =  Differentiate(var,order,AffineMap[j]);
      if (diffVal[3]!=0) {
        entries.push_back(diffVal[3]);
        col_notrack.push_back(i+1);
	++entryindex;
      } 
    }
  }
  row_starttrack.push_back(entries.size()+1);
  JACOBIAN.matrix_entries = new double[entries.size()];
  JACOBIAN.col_no = new int[col_notrack.size()];
  JACOBIAN.row_start = new int[Dimension+1];
  JACOBIAN.NoOfRows = Dimension;
  for (int i=0; i<entries.size(); ++i) {
    JACOBIAN.matrix_entries[i] = entries[i];
    JACOBIAN.col_no[i] = col_notrack[i];
  }
  for (int i=0; i<Dimension+1; ++i)
  JACOBIAN.row_start[i] = row_starttrack[i];
  return JACOBIAN;
}

// return determinant of jacobian
double Finite_Element::detJ()
{
  CSRMatrix J = Jacobian();
  return det(J);
}

/* Returns the representation of the differential operator defined by the input arguements
   ie the input vectors
              SpatialVars1       D_Orders
               [1 2 3 1 1]     [2 1 3 1 1]
    index:      4 3 2 1 0       4 3 2 1 0

    would return the representation
                   [0 0 1 2 2 2 0 0]
            index:  7 6 5 4 3 2 1 0
*/
std::vector<int> findRepresentation(const std::vector<double>& SpatialVars1, const std::vector<double>& D_Orders)
{
 std::vector<int> representation;
 int size = 0;
 for (int i=0; i<D_Orders.size(); ++i) {
   for (int j=0; j<D_Orders[i]; ++j)
   representation.push_back(SpatialVars1[i]-1);
 }
 return representation;
}

// Differentiate Basis function BFIdentifier (1starting) w.r.t the global spatial variables defined by SpatialVars and D_Orders.
// The returned vector is this derivative represented in terms of local reference element coordinates as a polynomial.
std::vector<double> Finite_Element::Differentiate_BF(const std::vector<double>& SpatialVars1, const std::vector<double>& D_Orders, const int& BFIdentifier)
{ 
  int derivativeOrder = 0;
  for (int i=0; i<D_Orders.size(); ++i)
  derivativeOrder+=D_Orders[i];

  if (derivativeOrder==1) {
   CSRMatrix JACOBIAN = Jacobian();
   CSRMatrix JTRANSPOSE = CSR_TRANSPOSE(JACOBIAN.matrix_entries, JACOBIAN.col_no, JACOBIAN.row_start, Dimension);
   CSRMatrix JTINVERSE =  CSR_INVERSION(JTRANSPOSE.matrix_entries, JTRANSPOSE.col_no, JTRANSPOSE.row_start, Dimension);
   double spatialVar = SpatialVars1[0]; 
   std::vector<double> Basis_Function = BasisFunctions[BFIdentifier-1];
   std::vector< std::vector<double> > BF_LocalDerivatives;
   std::vector<double> scalars;
   for (int i=0; i<Dimension; ++i) {
     int order=1, var=i+1;    
     std::vector<double> BasisFunctionDerivative = Differentiate(var,order,Basis_Function);
     BF_LocalDerivatives.push_back(BasisFunctionDerivative);
     bool hit=false;
     for (int j=JTINVERSE.row_start[((int)spatialVar) -1]; j<JTINVERSE.row_start[(int)spatialVar]; ++j) {
       if (JTINVERSE.col_no[j-1]==(i+1)) {
         hit = true;
         scalars.push_back(JTINVERSE.matrix_entries[j-1]);
       }
     }
     if (hit==false)
     scalars.push_back(0.0);
   }
   std::vector<double> localCoordinateBFDerivativeRepresentation =  LinearComb(scalars, BF_LocalDerivatives);
   return localCoordinateBFDerivativeRepresentation;
  }
  else {
    CSRMatrix A;
    int A_dim = (int)pow((double)Dimension, (double)derivativeOrder), iter=0, subiter=0;
    A.matrix_entries = new double[A_dim*A_dim];
    A.col_no = new int[A_dim*A_dim];
    A.row_start = new int[A_dim+1];
    A.NoOfRows = A_dim;
    BaseRep localDerivative(Dimension, derivativeOrder);
    for (int i=0; i<A_dim; ++i) {
      A.row_start[i] = iter+1;
      subiter=0;
      BaseRep globalDerivative(Dimension, derivativeOrder);
      while (subiter<A_dim) {
        double currVal = 1.0;
        for (int j=0; j<globalDerivative.Representation.size(); ++j) 
        currVal*= (Differentiate(localDerivative.Representation[j]+1,1,AffineMap[globalDerivative.Representation[j]]))[3]; 
        A.matrix_entries[iter] = currVal;
	A.col_no[iter] = subiter+1;
	++subiter;
	++globalDerivative;
	++iter;
      }
      ++localDerivative;
    }
    A.row_start[A_dim] = (A_dim*A_dim)+1;
    CSRMatrix AInverse = CSR_INVERSION(A.matrix_entries, A.col_no, A.row_start, A.NoOfRows);
    std::vector<int> diffOperator = findRepresentation(SpatialVars1,  D_Orders);
    std::vector< std::vector<double> > localDerivatives;
    BaseRep derivativeTrack(Dimension, derivativeOrder);
    std::vector<double> operand = BasisFunctions[BFIdentifier-1];
    for (int i=0; i<A.NoOfRows; ++i) {
      for (int j=0; j<derivativeTrack.Representation.size(); ++j) 
      operand = Differentiate(derivativeTrack.Representation[j]+1,1,operand);
      localDerivatives.push_back(operand);
      operand = BasisFunctions[BFIdentifier-1];
      ++derivativeTrack; 
    }
    derivativeTrack.zero();
    int row_index = 0;
    bool increment_index=false;
    for (int i=0; i<diffOperator.size(); ++i) {
      if (diffOperator[i]!=0)
      increment_index=true;
    }
    while (increment_index) {
      ++derivativeTrack;
      bool increment_index_switch = true;
      for (int i=0; i<derivativeTrack.Representation.size(); ++i) {
        if (derivativeTrack.Representation[i]!=diffOperator[i])
        increment_index_switch = false;
      }
      if (increment_index_switch)
      increment_index=false;
      ++row_index;
    }
    std::vector<double> scalars;
    for (int i=0; i<A.NoOfRows; ++i)
    scalars.push_back(0.0);
    for (int i=AInverse.row_start[row_index]-1; i<AInverse.row_start[row_index+1]-1; ++i) 
    scalars[AInverse.col_no[i]-1]=AInverse.matrix_entries[i];
    
    std::vector<double> localderivativeRepresentation = LinearComb(scalars, localDerivatives );
    return  localderivativeRepresentation;
  }
}
