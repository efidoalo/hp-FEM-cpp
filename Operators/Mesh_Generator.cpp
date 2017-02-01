/*===============================================;
 *
 *  File: Mesh_Generator.cpp
 *  Content: Source File containing the definition 
 *  of the function that generates the 
 *  computational mesh.
 *  Date: 12/06/2016
 *  Author: Andrew Oldham
 *
 *************************************************/

#include "Mesh_Generator.h"
#include "../Mesh/Mesh.h"
#include "CSRMatrix.h"
#include "Matrix_Algorithms.h"
#include "../Basis Functions/Basis_Functions.h"
#include "../Input Docs/Bilinear_Form.h"
#include "../Quadratures/GQ_NodesWeights_Generation.h"
#include "../Input Docs/PDE_Definition.h"
#include "Polynomial_Operators.h"
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <unordered_map>
#include <iostream>
#include <list>
#include <cmath>

// Generate initiall uniform 1D Mesh.
Mesh GenerateMesh()
{ 
  std::ifstream DataInputStream("Input Docs/Computational Domain.dat",std::ios_base::in);
  if (DataInputStream) {
    // Read input file into block of dynamically allocated memory (buffer)
    DataInputStream.seekg(0,DataInputStream.end);
    int length = DataInputStream.tellg();
    DataInputStream.seekg(0,DataInputStream.beg);
    char* buffer = new char[length];
    DataInputStream.read(buffer,length);

    // Store information from input file. 
    int pos=0;
    char DimText[4] = {'d','i','m','='};
    char DensityText[14] = {'d','e','n','s','i','t','y','_','p','a','r','a','m','='};
    char BasisFunctionText[22] = {'b','a','s','i','s','_','f','u','n','c','t','i','o','n','_','d','e','g','r','e','e','='};
    char QuadratureText[30] = {'Q','u','a','d','r','a','t','u','r','e',' ','P','o','i','n','t','s',' ','p','e','r',' ','E','l','e','m','e','n','t','='};
    char BoundaryText[8] = {'B','O','U','N','D','A','R','Y'};
    char InteriorText[8] = {'I','N','T','E','R','I','O','R'};
    // Get address in buffer where the boundary and interior points are defined as well as the dimension of the domain and the fixed density paramter
    char* DIM_Addr = std::search(buffer,buffer+length*(sizeof(char)), DimText, DimText + 4*sizeof(char)) + 4*sizeof(char);
    char* DP_Addr = std::search(DIM_Addr,buffer+length*(sizeof(char)),DensityText, DensityText + 14*sizeof(char)) + 14*sizeof(char);
    char* BFDeg_Addr = std::search(DP_Addr,buffer+length*(sizeof(char)), BasisFunctionText, BasisFunctionText+ 22*sizeof(char)) + 22*sizeof(char);
    char* QP_Addr = std::search(BFDeg_Addr,buffer+length*(sizeof(char)), QuadratureText, QuadratureText + 30*sizeof(char)) + 30*sizeof(char);
    char* BP_Addr = std::search(QP_Addr,buffer+length*(sizeof(char)),BoundaryText, BoundaryText +8*sizeof(char)) + 9*sizeof(char);
    char* IP_Addr = std::search(BP_Addr, buffer+length*(sizeof(char)),InteriorText, InteriorText + 8*sizeof(char)) + 9*sizeof(char);
    std::string numeralContainer;
    while (*DIM_Addr!='\n') {
      numeralContainer.push_back(*DIM_Addr);
      ++DIM_Addr;
    }
    int Dim = atoi(numeralContainer.c_str());
    numeralContainer.clear();
    while (*DP_Addr!='\n') {
      numeralContainer.push_back(*DP_Addr);
      ++DP_Addr;
    }
    double densityParam = strtod(numeralContainer.c_str(), 0);
    numeralContainer.clear();
    while (*BFDeg_Addr!='\n') {
      numeralContainer.push_back(*BFDeg_Addr);
      ++BFDeg_Addr;
    } 
    int BFDeg = atoi(numeralContainer.c_str()); 
    numeralContainer.clear();
    while (*QP_Addr!='\n') {
      numeralContainer.push_back(*QP_Addr);
      ++QP_Addr;
    }
    int NoOfQuadraturePoints = atoi(numeralContainer.c_str());
    double* zeros_store = new double[NoOfQuadraturePoints];
    double* weights_store = new double[NoOfQuadraturePoints];
    gauss_quadrature(zeros_store, weights_store, NoOfQuadraturePoints);
    for (int i=0; i<NoOfQuadraturePoints; ++i) {
       GaussianQuadrature::Nodes.push_back(zeros_store[i]);
       GaussianQuadrature::Weights.push_back(weights_store[i]);    
    } 
    delete[] zeros_store;
    delete[] weights_store;
    numeralContainer.clear();
    while (*BP_Addr!=' ') {
      numeralContainer.push_back(*BP_Addr);
      ++BP_Addr;
    }
    double LeftBoundaryPoint = strtod(numeralContainer.c_str(),0);
    numeralContainer.clear();
    ++BP_Addr;
    while (*BP_Addr!='\n') {
       numeralContainer.push_back(*BP_Addr);
       ++BP_Addr;
    }
    double RightBoundaryPoint = strtod(numeralContainer.c_str(),0);
    numeralContainer.clear();

    // For the case of a 1D computational Domain [no interior points guaranteed], generate initial uniform grid
    if (Dim==1) {
      Basis_Functions::Generate_Basis_Functions(BFDeg); // generate the initial basis functions on the reference element (degree specified in input file)
      double NoOfElements = ceil((RightBoundaryPoint-LeftBoundaryPoint)/densityParam);
      double uniformVertexInterval = (RightBoundaryPoint-LeftBoundaryPoint)/NoOfElements;
      CSRMatrix ConnectivityMatrix;
      ConnectivityMatrix.matrix_entries = new double[((int)NoOfElements)*(BFDeg+1)];
      ConnectivityMatrix.col_no = new int[((int)NoOfElements)*(BFDeg+1)];
      ConnectivityMatrix.row_start = new int[((int)NoOfElements)+1];
      ConnectivityMatrix.NoOfRows = (int)NoOfElements;
      for (int i=0; i<NoOfElements; ++i) 
      ConnectivityMatrix.row_start[i] = i*(BFDeg+1) + 1;
      ConnectivityMatrix.row_start[(int)NoOfElements] = (((int)NoOfElements)*(BFDeg+1))+1;
      std::unordered_map<int,Finite_Element> FEContainer;
      int iter=0, NodeNumber=1;
      for (int i=0; i<NoOfElements;++i) {
	std::vector< std::vector<double> > Vertex_Points;
	double leftVertex = LeftBoundaryPoint + (i*uniformVertexInterval), rightVertex = LeftBoundaryPoint + ((i+1.0)*uniformVertexInterval);
	if (i==NoOfElements-1)
        rightVertex = RightBoundaryPoint;      
        std::vector<double> leftPoint = {leftVertex};
        std::vector<double> rightPoint = {rightVertex};
        Vertex_Points = {leftPoint, rightPoint};
	int LINESPECIFIER=0; // defines the 1D line finite element type.
        Finite_Element currFE(LINESPECIFIER,BFDeg,Vertex_Points,Dim);
	for (int j=0; j<BFDeg+1; ++j) {
          ConnectivityMatrix.matrix_entries[iter] = NodeNumber;
	  ConnectivityMatrix.col_no[iter]=j+1;
          ++iter;
	  ++NodeNumber;
	}
	--NodeNumber;
        FEContainer.emplace((i+1),currFE);
      }
 
      Mesh GeneratedMesh(NoOfElements, NodeNumber, FEContainer, ConnectivityMatrix,uniformVertexInterval);
      return GeneratedMesh;
    }

  } 
}

// Compute and return the hash value of the position (row = node1, column = node2) in a unit_dist x unit_dist
// square matrix. Hash value returned is unique to the position in the square matrix.
int hashFunc(const int& node1, const int& node2, const int& unit_dist)
{ 
  return (node1-1)*unit_dist + node2;
}

// Inserts val into the std::list container, preserving order, if val is not already contained within it.
// Otherwise Container remains unchanged
void orderedInsert(std::list<int>& Container, const int& val)
{ 
  if (Container.size()==0) {
    Container.push_back(val);
    return;
  }
  std::list<int>::iterator LeftIter = Container.begin();
  std::list<int>::iterator RightIter = Container.end();
  if (Container.size()==1) {   
    if (*LeftIter<val)
    Container.push_back(val);
    if (val<*LeftIter)
    Container.insert(LeftIter, val);
    return;
  }
  --RightIter;
  if (val<*LeftIter) {
    Container.push_front(val);
    return;
  }
  if (*RightIter<val) {
    Container.push_back(val);
    return;
  }
  if ( (val==*LeftIter) || (val==*RightIter) )
  return;
  int distance = 0;
  while (true)
  { 
    if ( *(++LeftIter)==*RightIter ) {
      Container.insert(RightIter, val);
      return;
    }
    --LeftIter;
    std::list<int>::iterator temp = LeftIter;
    distance = 0;
    while (temp!=RightIter) {
      ++temp;
      ++distance;
    }
    temp = LeftIter;
    std::advance(temp, distance/2);
    
    if (*temp==val)
    return;
    if (val<*temp)
    RightIter = temp;
    if (*temp<val)
    LeftIter = temp;  
  }  
}


// Returns the value of the bilinear form ( a( . , . ) ) , evaluated over the finite element FE, 
// whose ordered arguements are the jth and ith basis functions of the finite element FE.
double a(Finite_Element& FE, const int& j, const int& i)
{ 
  double volIntegral = 0.0;
  std::vector< std::vector<double> > nodes;
  std::vector< double > leftBoundary = (FE.get_VertexPoints())[0];
  std::vector< double> rightBoundary = (FE.get_VertexPoints())[1];
  std::vector<double> currNode;
  for (int k=0; k<GaussianQuadrature::Nodes.size(); ++k) {
    currNode = {(0.5*GaussianQuadrature::Nodes[k]) + 0.5}; // quadrature points rescaled for use on the refrence element.
    nodes.push_back(currNode);
  } 

  for (int k=0; k<GaussianQuadrature::Nodes.size(); ++k) 
  volIntegral+=GaussianQuadrature::Weights[k]*BilinearForm_VI(FE, j, i, nodes[k]);
  
  volIntegral*=(0.5*FE.detJ());
  return volIntegral;
}

// Evaluates the linear functional l( . ) that is defined in the weak formulation of the PDE over the finite element FE,
// taking as an arguement the ith basis function of the finite element FE.
double l(Finite_Element& FE, const int& i)
{
  double volIntegral = 0.0;
  std::vector< std::vector<double> > nodes;
  for (int k=0; k<GaussianQuadrature::Nodes.size(); ++k) {
    std::vector<double> currNode = {(0.5*GaussianQuadrature::Nodes[k]) + 0.5}; // quadrature points rescaled for use on the refrence element.
    nodes.push_back(currNode);
  } 
  for (int k=0; k<GaussianQuadrature::Nodes.size(); ++k) 
  volIntegral+= GaussianQuadrature::Weights[k]*LinearFunctional_VI(FE, i, nodes[k]); 
  volIntegral*=(0.5*FE.detJ());
  return volIntegral;
}

// Generates and returns stiffness matrix.
CSRMatrix GenerateStiffnessMatrix(Mesh& COMPUTATIONALMESH)
{ 
  int n = COMPUTATIONALMESH.get_NoOfNodes();
  int NoOfElements = COMPUTATIONALMESH.get_NoOfElements();
  int* ConnArray_row_start = COMPUTATIONALMESH.ConnArray_rowptr();
  double* ConnArray_mentries = COMPUTATIONALMESH.ConnArray_entryptr();
  CSRMatrix STIFFNESSMATRIX; 
  std::list<int> HashedNodePairs;
  for (int i=0; i<NoOfElements; ++i) {
    for (int j=(ConnArray_row_start[i]-1); j<(ConnArray_row_start[i+1]-1); ++j) {
      for (int k=(ConnArray_row_start[i]-1); k<(ConnArray_row_start[i+1]-1); ++k) {
        orderedInsert(HashedNodePairs, hashFunc(ConnArray_mentries[j], ConnArray_mentries[k],n) );          	
      }
    }
  }
  STIFFNESSMATRIX.matrix_entries = new double[HashedNodePairs.size()];
  for (int i=0; i<HashedNodePairs.size(); ++i)
  STIFFNESSMATRIX.matrix_entries[i] = 0.0;
  STIFFNESSMATRIX.col_no = new int[HashedNodePairs.size()];
  STIFFNESSMATRIX.row_start = new int[n+1];
  STIFFNESSMATRIX.NoOfRows = n;
  std::list<int>::iterator it = HashedNodePairs.begin();
  int row_startIndex = 0; 
  for (int i=0; i<n; ++i) {
    int lowerBound = i*n;
    while (*it<=lowerBound) {
      ++it;
      ++row_startIndex;
    }
    STIFFNESSMATRIX.row_start[i] = row_startIndex+1; 
  }
  STIFFNESSMATRIX.row_start[n] = HashedNodePairs.size() +1;
  Finite_Element currFE;
  int ConnectivityRowIndex=-1, currValIndex=0, positionHash=0;
  double currVal=0.0;
  for (int i=0; i<NoOfElements; ++i) {	  
    currFE = COMPUTATIONALMESH.find_FE(i+1);
    ConnectivityRowIndex = ConnArray_row_start[i] - 1;
    for (int j=1; j<=(currFE.get_BFDegree()+1); ++j) {
      for (int k=1; k<=(currFE.get_BFDegree()+1); ++k) {

        currVal = a(currFE, j, k);
	positionHash = hashFunc(ConnArray_mentries[ConnectivityRowIndex+k-1],
                                ConnArray_mentries[ConnectivityRowIndex+j-1],
				n );
	it = HashedNodePairs.begin();
	currValIndex = 0;
	while (*it!=positionHash) {
          ++it;
	  ++currValIndex;
	}
        STIFFNESSMATRIX.matrix_entries[currValIndex] += currVal;
	STIFFNESSMATRIX.col_no[currValIndex] = ConnArray_mentries[ConnectivityRowIndex+j-1];
      }
    }
  }
  return STIFFNESSMATRIX;
}

// Genertes and returns load vector
double* GenerateLoadVector(Mesh& COMPUTATIONALMESH)
{
  double *F = new double[COMPUTATIONALMESH.get_NoOfNodes()];
  for (int i=0; i<COMPUTATIONALMESH.get_NoOfNodes(); ++i) 
  F[i]=0.0;
  double *ConnectivityArray_entriesptr = COMPUTATIONALMESH.ConnArray_entryptr();
  int *ConnectivityArray_rowptr = COMPUTATIONALMESH.ConnArray_rowptr();
  int NodesPerElement =-1, NoOfElements = COMPUTATIONALMESH.get_NoOfElements();
  Finite_Element currFE;
  for (int i=0; i<NoOfElements; ++i) {
    Finite_Element currFE = COMPUTATIONALMESH.find_FE(i+1);
    NodesPerElement = currFE.get_BFDegree() + 1;
    int ConnectivityRowIndex = ConnectivityArray_rowptr[i] - 1;
    for (int j=0; j<NodesPerElement; ++j) {
      F[(int)(ConnectivityArray_entriesptr[ConnectivityRowIndex+j])-1]+=l(currFE,j+1);
    }   
  }
  return F;
}


// Impose boundary conditions on system.
void BC_IMPOSITION(Mesh& COMPUTATIONALMESH, CSRMatrix& STIFFNESSMATRIX, double* LOADVEC)
{ 
  // Assume homogeneous Dirichlet boundary conditions
  bool homDirichlet = true;
  std::ifstream DataInputStream("Input Docs/Computational Domain.dat",std::ios_base::in);
  if (DataInputStream) {
    DataInputStream.seekg(0,DataInputStream.end);
    int length = DataInputStream.tellg();
    DataInputStream.seekg(0,DataInputStream.beg);
    char* buffer = new char[length];
    DataInputStream.read(buffer,length);
    int pos=0;
    char BCTYPE[24] = {'h','o','m','o','g','e','n','e','o','u','s',' ','d','i','r','i','c','h','l','e','t',' ','=',' '};
    // Get address in buffer where the type of boundary conditions are specified
    char* BCTYPE_Addr = std::search(buffer,buffer+length*(sizeof(char)), BCTYPE, BCTYPE + 24*sizeof(char)) + 24*sizeof(char);
    if (*BCTYPE_Addr=='n')
    homDirichlet = false;
    // HOMOGENEOUS DIRICHLET BC'S
    if (homDirichlet) {
      STIFFNESSMATRIX.matrix_entries[0] = 1.0; LOADVEC[0] = 0.0;
      for (int i=1; i<STIFFNESSMATRIX.row_start[1]-1;++i)
      STIFFNESSMATRIX.matrix_entries[i] = 0.0;
      for (int i=1; i<STIFFNESSMATRIX.NoOfRows; ++i) {
        if (STIFFNESSMATRIX.col_no[STIFFNESSMATRIX.row_start[i]-1]==1)
        STIFFNESSMATRIX.matrix_entries[STIFFNESSMATRIX.row_start[i]-1] = 0.0;
	if (STIFFNESSMATRIX.col_no[STIFFNESSMATRIX.row_start[i]-2]==STIFFNESSMATRIX.NoOfRows)
        STIFFNESSMATRIX.matrix_entries[STIFFNESSMATRIX.row_start[i]-2]=0.0;
      }
      int i =STIFFNESSMATRIX.row_start[STIFFNESSMATRIX.NoOfRows - 1] - 1;
      for (i; i<STIFFNESSMATRIX.row_start[STIFFNESSMATRIX.NoOfRows]-2; ++i)
      STIFFNESSMATRIX.matrix_entries[i] = 0.0;
      STIFFNESSMATRIX.matrix_entries[i] = 1.0; LOADVEC[STIFFNESSMATRIX.NoOfRows-1] = 0.0;
      delete[] buffer;
    }
    else {
      // NONHOMOGENEOUS DIRICHLET BC'S
      Finite_Element Element1 = COMPUTATIONALMESH.find_FE(1);
      Finite_Element ElementFinal = COMPUTATIONALMESH.find_FE(COMPUTATIONALMESH.get_NoOfElements());
      std::vector<double> firstNode = (Element1.get_VertexPoints())[0];
      std::vector<double> finalNode = (ElementFinal.get_VertexPoints())[1];
      LOADVEC[0] =  nonhomogeneousDirichlet(firstNode); 
      LOADVEC[STIFFNESSMATRIX.NoOfRows-1] = nonhomogeneousDirichlet(finalNode);
      STIFFNESSMATRIX.matrix_entries[0]=1.0;
      STIFFNESSMATRIX.matrix_entries[STIFFNESSMATRIX.row_start[STIFFNESSMATRIX.NoOfRows]-2]=1.0;
      for (int i=1; i<STIFFNESSMATRIX.row_start[1]-1;++i)
      STIFFNESSMATRIX.matrix_entries[i] = 0.0;
      for (int i=1; i<STIFFNESSMATRIX.NoOfRows-1; ++i) {
	if (STIFFNESSMATRIX.col_no[STIFFNESSMATRIX.row_start[i]-1]==1) {
          LOADVEC[i] = LOADVEC[i] - LOADVEC[0]*STIFFNESSMATRIX.matrix_entries[STIFFNESSMATRIX.row_start[i]-1];
          STIFFNESSMATRIX.matrix_entries[STIFFNESSMATRIX.row_start[i]-1] = 0;
	}
	if (STIFFNESSMATRIX.col_no[STIFFNESSMATRIX.row_start[i+1]-2]==STIFFNESSMATRIX.NoOfRows) {
	  LOADVEC[i] = LOADVEC[i] - LOADVEC[STIFFNESSMATRIX.NoOfRows-1]*STIFFNESSMATRIX.matrix_entries[STIFFNESSMATRIX.row_start[i+1]-2];
          STIFFNESSMATRIX.matrix_entries[STIFFNESSMATRIX.row_start[i+1]-2] = 0.0;
	}        
      }
      int i =STIFFNESSMATRIX.row_start[STIFFNESSMATRIX.NoOfRows-1]-1;
      for ( i; i<STIFFNESSMATRIX.row_start[STIFFNESSMATRIX.NoOfRows]-2;++i)
      STIFFNESSMATRIX.matrix_entries[i] = 0.0; 
      STIFFNESSMATRIX.matrix_entries[i] = 1.0; 
      delete[] buffer;
    }  
  }  
}

// Test to see if the input arguement matrix is symmetrix positive definite.
// returns true if so, false otherwise.
bool SYMM_PD(CSRMatrix& STIFFNESSMATRIX)
{   
   std::ifstream DataInputStream("Input Docs/Computational Domain.dat",std::ios_base::in);
   DataInputStream.seekg(0,DataInputStream.end);
   int length = DataInputStream.tellg();
   DataInputStream.seekg(0,DataInputStream.beg);
   char* buffer = new char[length];
   DataInputStream.read(buffer,length);
   // Check if bilinear form (and thus the stiffness matrix) is symmetric [this information is provided by the user]
   char symmTxt[24] = {'s','y','m','m','e','t','r','i','c','_','b','i','l','i','n','e','a','r','_','f','o','r','m','='};
   char* SYMM_ADDR = std::search(buffer,buffer+length*(sizeof(char)), symmTxt, symmTxt + 24*sizeof(char)) + 24*sizeof(char);
   //Sylvester's criterion for Symmetric positive definiteness
   if (*SYMM_ADDR=='y') {
     for (int i=0; i<STIFFNESSMATRIX.NoOfRows; ++i) {
       double detSUB = 0.0;
       CSRMatrix currSUB = submatrix(STIFFNESSMATRIX, i+1);
      // CSR_MatrixOutput(currSUB.matrix_entries, currSUB.col_no,currSUB.row_start, currSUB.NoOfRows, currSUB.NoOfRows);
       detSUB = det1(currSUB);  //`std::cout<<"\ndet 0f current submatirx : "<<detSUB;
       if (detSUB<=0) 
       return false;       
     }
     return true;
   }
   else return false; 
}

// Evaluate the solution obtained via Fem hp method at the point given by std::vector<double> Point.
// std::vector<double> Point [in] defines the point on the mesh under that we wish to evaluate the solution over.
double EvaluateSolution(Mesh& A, double* coeffs, std::vector<double> Point)
{ 
  int dim = A.find_FE(1).get_Dim();
  int NoOfElements = A.get_NoOfElements();
  int exact_element = -1;
  if (dim==1) {
    bool found_element = false;
    int left_element = 1, right_element = NoOfElements; 
    while (found_element==false) { 
      if ((right_element - left_element)==2) {
        exact_element = left_element+1;
        found_element = true;
      }
      if (((A.find_FE(left_element).get_VertexPoints())[0][0]<=Point[0]) && (Point[0]<=(A.find_FE(left_element).get_VertexPoints())[1][0])) {
        found_element = true;
	exact_element = left_element;
      }
      if (((A.find_FE(right_element).get_VertexPoints())[0][0]<=Point[0]) && (Point[0]<=(A.find_FE(right_element).get_VertexPoints())[1][0])) {
        found_element = true;
	exact_element = right_element;
      }
      if (found_element==false) {
        int temp_element = (left_element + right_element)/2;
	if ((A.find_FE(temp_element).get_VertexPoints())[0][0]<=Point[0])
        left_element = temp_element;
	else right_element = temp_element;
      }
      
    }
    
    Finite_Element RelevantElem = A.find_FE(exact_element);
    double evaluatedVal = 0.0;
    for (int i=0; i<RelevantElem.get_BFDegree()+1; ++i) {
      std::vector<double> inverseMapping = (RelevantElem.get_InverseAffineMap())[0];
      double local_coord = Eval(Point, inverseMapping);
      std::vector<double> localPoint = {local_coord}; 
      evaluatedVal+=coeffs[(int)((A.get_ConnArray()).matrix_entries[(A.get_ConnArray()).row_start[exact_element-1]-1+i])-1]
	                      *Eval(localPoint,RelevantElem.get_BFncs()[i]);
    }
    return evaluatedVal;
  }
}


// Determine if Mesh Refinement is required
bool MeshRefinement()
{ 
  bool hp_refine = false;
  std::ifstream DataInputStream("Input Docs/Computational Domain.dat",std::ios_base::in);
  if (DataInputStream) {
    DataInputStream.seekg(0,DataInputStream.end);
    int length = DataInputStream.tellg();
    DataInputStream.seekg(0,DataInputStream.beg);
    char* buffer = new char[length];
    DataInputStream.read(buffer,length);
    int pos=0;
    char refTxt[19] = {'h','p','_','M','e','s','h','_','R','e','f','i','n','e','m','e','n','t','='};
    // Get address in buffer where the boundary and interior points are defined as well as the dimension of the domain and the fixed density paramter
    char* ref_Addr = std::search(buffer,buffer+length*(sizeof(char)), refTxt, refTxt + 19*sizeof(char)) + 19*sizeof(char);
    if (*ref_Addr=='y')
    hp_refine = true;
    return hp_refine;
  }
  return false;
}

// Mark Elements for refineent. Only used during verification of program on a particular test problem.
void MarkElements(Mesh& COMPUTATIONALMESH, double* coeffVals, std::vector<int>& refinementVect, const double& markParam, double* L2_projection_Of_f)
{
  std::vector<double> localEstimators, localPoint, spatialvars = {1.0}, d_orders = {2.0};
  double currLocalEstimator = 0.0, C_p=0.0, maxEst = -1, leftP=0.0, rightP=0.0;
  double* CArray_entries = COMPUTATIONALMESH.ConnArray_entryptr();
  int* CArray_row_start = COMPUTATIONALMESH.ConnArray_rowptr();
  Finite_Element currFE;
  for (int i=0; i<COMPUTATIONALMESH.get_NoOfElements(); ++i) {
    currFE = COMPUTATIONALMESH.find_FE(i+1);
    leftP = currFE.get_VertexPoints()[0][0];
    rightP = currFE.get_VertexPoints()[1][0];
    C_p = 1.0/(sqrt(a(0.0)*(currFE.get_BFDegree())*(currFE.get_BFDegree()+1.0)));
    std::vector<double> scalars;
    std::vector<double> fproject_scalars;
    std::vector< std::vector<double> > local_2nd_derivatives;
    for (int j=CArray_row_start[i]-1; j<CArray_row_start[i+1]-1; ++j) { 
      scalars.push_back(coeffVals[((int)(CArray_entries[j]))-1]);
      fproject_scalars.push_back( L2_projection_Of_f[((int)(CArray_entries[j]))-1] );
    }
    for (int j=0; j<currFE.get_BFDegree()+1; ++j)
    local_2nd_derivatives.push_back(currFE.Differentiate_BF(spatialvars, d_orders, j+1));
    std::vector<double> localSol = LinearComb(scalars, currFE.get_BFncs());
    std::vector<double> localSol_2nd_deriv = LinearComb(scalars, local_2nd_derivatives);
    std::vector<double> l2_projection_of_f = LinearComb(fproject_scalars, currFE.get_BFncs());
    currLocalEstimator = 0.0;
    for (int j=0; j<GaussianQuadrature::Nodes.size(); ++j) {
      localPoint = {(0.5*GaussianQuadrature::Nodes[j]) + 0.5 };
      currLocalEstimator+= GaussianQuadrature::Weights[j]*(pow( ((Eval(localPoint,l2_projection_of_f) + a(0.0)*Eval(localPoint, localSol_2nd_deriv)-c(0.0)*Eval(localPoint,localSol) )*(sqrt(weight_Fnc(leftP,rightP,Eval(localPoint, currFE.get_AffineMap()[0]))))) , 2.0));
    }
    currLocalEstimator*=(0.5*currFE.detJ());
    currLocalEstimator = sqrt(currLocalEstimator);
    currLocalEstimator*=C_p;
    localEstimators.push_back(currLocalEstimator);
    if (currLocalEstimator>maxEst)
    maxEst = currLocalEstimator;
  }
  for (int i=0; i<COMPUTATIONALMESH.get_NoOfElements(); ++i) {
    if (maxEst<= ((1.0/markParam)*localEstimators[i]) ) 
    refinementVect.push_back(i+1);
  }
}

// Refines the mesh COMPUTATIONALMESH using h-p refinement.
void RefineMesh( Mesh& COMPUTATIONALMESH, double* coeffVals, const double& initial_h, const int& initial_NoOfElements, 
		     const double& tol, const double& SmoothnessParam, const double& markerParam, const double& C_0, double* L2_projection_Of_f)
{ 
  std::vector<int> refinement_vec;
  Finite_Element currFE, refineFE;
  double loclalErrEstimate = 0.0;
  double h=0.0;
  for (int i=0; i<COMPUTATIONALMESH.get_NoOfElements(); ++i) {
    if (local_error_estimator(COMPUTATIONALMESH,coeffVals,i+1,COMPUTATIONALMESH.get_h()) >=
		    element_error_criteria(COMPUTATIONALMESH,i+1,initial_h,tol,initial_NoOfElements, C_0))
    refinement_vec.push_back(i+1);
  }  
 
  //MarkElements(COMPUTATIONALMESH, coeffVals, refinement_vec, markerParam, L2_projection_Of_f); was used for program verification.
  

  std::unordered_map<int,Finite_Element> ElementContainer;
  int refine_element_track = 0, newNoOfElements=0, newNoOfNodes=0, currNodeNo=1, currColNo=1;
  std::vector<double> conn_entries;
  std::vector<int> conn_col_no, conn_row_start;  
  for (int i=0; i<COMPUTATIONALMESH.get_NoOfElements(); ++i) { 
    bool existence = std::binary_search(refinement_vec.begin(), refinement_vec.end(), i+1);
    conn_row_start.push_back(conn_entries.size()+1);
    if (existence) { 
      if (smoothness_measure(i, COMPUTATIONALMESH, coeffVals)>=SmoothnessParam) {
        p_refinement(COMPUTATIONALMESH, i, ElementContainer, refineFE, newNoOfNodes, refine_element_track, conn_entries,  conn_col_no,  currNodeNo,  currColNo, h);
      }
      else {
      h_refinement(COMPUTATIONALMESH, i,ElementContainer, refineFE, newNoOfNodes, refine_element_track, 
		                     conn_entries, conn_col_no, conn_row_start, currNodeNo, currColNo, h);
      }
    }
    else { 
      currFE = COMPUTATIONALMESH.find_FE(i+1);
      double h_element = (currFE.get_VertexPoints())[1][0] - (currFE.get_VertexPoints())[0][0];
      if (h<h_element)
      h=h_element;
      if (i==0)
      newNoOfNodes+=(currFE.get_BFDegree()+1);
      else newNoOfNodes+=(currFE.get_BFDegree());
      ElementContainer[i+1+refine_element_track]=currFE;
      for (int j=0; j<currFE.get_BFDegree()+1; ++j) {
        conn_entries.push_back(currNodeNo);
	conn_col_no.push_back(currColNo);
	++currNodeNo;
	++currColNo;
      }
      currColNo=1;
      --currNodeNo;    
    }
  }
  conn_row_start.push_back(conn_entries.size()+1);
  newNoOfElements = COMPUTATIONALMESH.get_NoOfElements() + refine_element_track;
  CSRMatrix updateConnectivityArray;
  updateConnectivityArray.matrix_entries = new double[conn_entries.size()];
  updateConnectivityArray.col_no = new int[conn_col_no.size()];
  updateConnectivityArray.row_start = new int[conn_row_start.size()];
  for (int j=0; j<conn_entries.size(); ++j) {
    updateConnectivityArray.matrix_entries[j] = conn_entries[j];
    updateConnectivityArray.col_no[j] = conn_col_no[j];
  }
  for (int j=0; j<conn_row_start.size(); ++j) 
  updateConnectivityArray.row_start[j] = conn_row_start[j];
  updateConnectivityArray.NoOfRows = conn_row_start.size()-1;
  COMPUTATIONALMESH.redefineMESH(newNoOfElements, newNoOfNodes, ElementContainer, updateConnectivityArray, h); 
  
}

// h refine element i+1 in the mesh.
void h_refinement(Mesh& A, const int& i, std::unordered_map<int,Finite_Element>& ElementContainer, Finite_Element& refineFE, int& newNoOfNodes, 
	int& refine_element_track, std::vector<double>& conn_entries, std::vector<int>& conn_col_no, std::vector<int>& conn_row_start, int& currNodeNo, int& currColNo,
	double& h_)
{ 
  
  refineFE = A.find_FE(i+1);
  if (i==0)
  newNoOfNodes+=((2*(refineFE.get_BFDegree()))+1);
  else newNoOfNodes+=(2*refineFE.get_BFDegree());
  double leftpoint = (refineFE.get_VertexPoints())[0][0], rightpoint = (refineFE.get_VertexPoints())[1][0];
  double middlepoint = (leftpoint + rightpoint)/2.0;
  if (h_<(middlepoint-leftpoint))
  h_=middlepoint-leftpoint;
  if (h_<(rightpoint-middlepoint))
  h_=(rightpoint-middlepoint);
  std::vector<double> vertexpoint1 = {leftpoint};
  std::vector<double> vertexpoint2 = {middlepoint};
  std::vector<double> vertexpoint3 = {rightpoint};
  std::vector< std::vector<double> > leftFE_VertexPoints = {vertexpoint1, vertexpoint2};
  std::vector< std::vector<double> > rightFE_VertexPoints = {vertexpoint2, vertexpoint3};
  int FETYPE = refineFE.get_FETYPE(), BFDEGREE =refineFE.get_BFDegree(), DIM=refineFE.get_Dim();
  Finite_Element leftElement(FETYPE,  BFDEGREE, leftFE_VertexPoints, DIM);
  Finite_Element rightElement(FETYPE, BFDEGREE, rightFE_VertexPoints, DIM);
  ElementContainer[i+1+refine_element_track]=leftElement;
  ElementContainer[i+2+refine_element_track]=rightElement;
  for (int j=0; j<refineFE.get_BFDegree()+1; ++j) {
    conn_entries.push_back(currNodeNo);
    conn_col_no.push_back(currColNo);
    ++currNodeNo;
    ++currColNo;
  }
  currColNo=1;
  --currNodeNo;
  conn_row_start.push_back(conn_entries.size()+1);
  for (int j=0; j<refineFE.get_BFDegree()+1; ++j) {
    conn_entries.push_back(currNodeNo);
    conn_col_no.push_back(currColNo);
    ++currNodeNo;
    ++currColNo;
  }
  currColNo=1;
  --currNodeNo;
  ++refine_element_track;

  return;
}

// p refine element i+1 in the mesh
void p_refinement(Mesh& A, const int& i, std::unordered_map<int,Finite_Element>& ElementContainer, Finite_Element& refineFE, int& newNoOfNodes, 
	int& refine_element_track, std::vector<double>& conn_entries, std::vector<int>& conn_col_no, int& currNodeNo, int& currColNo, double& h_)
{
	
  refineFE = A.find_FE(i+1);
  if (i==0)
  newNoOfNodes+=refineFE.get_BFDegree()+2;
  else newNoOfNodes+=refineFE.get_BFDegree()+1;
  std::vector<double> leftVertex = (refineFE.get_VertexPoints())[0];
  std::vector<double> rightVertex = (refineFE.get_VertexPoints())[1];
  if (h_<(rightVertex[0]-leftVertex[0])) 
  h_=rightVertex[0] - leftVertex[0];
  std::vector< std::vector<double> > VertexPoints = {leftVertex, rightVertex};
  int refinedBasisFunctionDegree = refineFE.get_BFDegree() + 1;
  Basis_Functions::Generate_Basis_Functions(refinedBasisFunctionDegree);
  Finite_Element refinedElement(0, refinedBasisFunctionDegree, VertexPoints, 1);
  ElementContainer[i+1+refine_element_track] = refinedElement;
  for (int j=0; j<refinedElement.get_BFDegree()+1; ++j) {
    conn_entries.push_back(currNodeNo);
    conn_col_no.push_back(currColNo);
    ++currNodeNo;
    ++currColNo;
  }
  --currNodeNo;
  currColNo=1;
  return;
}

// Estimation of linfinitynorm of the polynomial Poly on the interval (0,1) 
double linfinityNorm(std::vector<double>& Poly)
{
  if ((int)Poly[0]==1) {
    std::vector<double> rightPoint = {1.0};
    std::vector<double> leftPoint = {0.0};
    return std::max(fabs(Eval(leftPoint, Poly)), fabs(Eval(rightPoint, Poly)) );
  }
  else {
    double linfinityapprox = 0.0, x_arg =0, curreval = 0.0;
    std::vector< double> currpoint;
    for (int i=0; i<101; ++i) {
      currpoint = { x_arg };
      curreval = fabs(Eval(currpoint, Poly));
      if (curreval>linfinityapprox)
      linfinityapprox = curreval;
      x_arg+=0.01;
    }
    return linfinityapprox;
  }
}

// returns estimate of the smoothness of the FEM solution on element element_index+1 in the mesh.
// 0 <= return value <= 1. 
// 0 << return value <= 1 suggests smooth solution over the element - thus p refine
// 0 <= return value <<1 suggests irregular/nonsmooth solution over the element - thus h refine
double smoothness_measure(const int& element_index, Mesh& COMPUTATIONALMESH, double* coeffVals)
{
  Finite_Element FE = COMPUTATIONALMESH.find_FE(element_index+1);
  int NoOfBasisFunctions = FE.get_BFDegree()+1;
  int* CArray_row_start = COMPUTATIONALMESH.ConnArray_rowptr();
  double* CArray_entries = COMPUTATIONALMESH.ConnArray_entryptr();
  std::vector< std::vector<double> > localBasisFunctions = (FE.get_BFncs());
  std::vector<double> scalars, svar = {1.0}, dorder={1.0};
  int row_start_index = CArray_row_start[element_index]-1;
  for (int j = row_start_index; j<row_start_index+NoOfBasisFunctions; ++j) 
  scalars.push_back(coeffVals[ ((int)(CArray_entries[j]))-1 ]);
  std::vector<double> approximate_solution = LinearComb(scalars,localBasisFunctions); // locally on th element 
  std::vector< std::vector<double> > localfirst_derivatives;
  for (int j=0; j<(FE.get_BFDegree()+1); ++j)
  localfirst_derivatives.push_back(FE.Differentiate_BF(svar, dorder, j+1));
  std::vector<double> approxSol_derivative = LinearComb(scalars, localfirst_derivatives); 
  double lInfinityNormSquared = pow(linfinityNorm(approximate_solution),2.0);
  double l2NormapproxSol = 0.0, l2NormapproxSolderiv=0.0, h=((FE.get_VertexPoints())[1][0]-(FE.get_VertexPoints())[0][0]);
  std::vector<double> currPoint;
  for (int j=0; j<GaussianQuadrature::Nodes.size(); ++j) { 
    currPoint = {(0.5*GaussianQuadrature::Nodes[j]) + 0.5};
    l2NormapproxSol+= GaussianQuadrature::Weights[j]*pow(Eval(currPoint, approximate_solution),2.0);
    l2NormapproxSolderiv+= GaussianQuadrature::Weights[j]*pow(Eval(currPoint, approxSol_derivative),2.0);
  }
  l2NormapproxSol*=(0.5*FE.detJ());
  l2NormapproxSolderiv*=(0.5*FE.detJ());
  double coth_1_ = (exp(1.0) + exp(-1.0))/(exp(1.0) - exp(-1.0));
  double smoothnessMeasure = lInfinityNormSquared/(coth_1_*( ((1.0/h)*l2NormapproxSol) + (h*l2NormapproxSolderiv) ) );
  return smoothnessMeasure;  
}


// Estimates (returns a sharp upper bound on) the Energy Norm of the error in the calculated approximate/discrete solution.
double EnergyNormOfErr_Estimator(Mesh& A, double* coeffs, double* L2_projection_Of_f)
{
  double eNormErrEst = 0.0, currElement_p=0.0, C_p=0.0, h2_NormOf_Resid=0.0, h2_NormOf_f=0.0, rightP=0.0, leftP=0.0;
  Finite_Element currFE;
  double* CArray_entries = A.ConnArray_entryptr();
  int* CArray_row_start = A.ConnArray_rowptr();
  std::vector<double> scalars, l2project_f_scalars, u_h, l2project_f, u_h_2nd_deriv, ref_elem_point = {0.0}, spatialVariable = {1.0}, derivativeOrder2 = {2.0};
  std::vector< std::vector<double> > second_Derivatives_OnRefElem;
  for (int i=0; i<A.get_NoOfElements(); ++i) {
    scalars.clear();
    l2project_f_scalars.clear();
    u_h.clear();
    l2project_f.clear();
    u_h_2nd_deriv.clear();
    second_Derivatives_OnRefElem.clear();
    currFE = A.find_FE(i+1);
    h2_NormOf_Resid = 0.0;
    h2_NormOf_f = 0.0;
    currElement_p = currFE.get_BFDegree();
    leftP = (currFE.get_VertexPoints()[0])[0];
    rightP = (currFE.get_VertexPoints()[1])[0];
    
    C_p = 1.0/(a(0.0)*currElement_p*(currElement_p+1.0));
    for (int j=CArray_row_start[i]-1; j<CArray_row_start[i+1]-1; ++j) {
      scalars.push_back(coeffs[ (int)(CArray_entries[j]) - 1 ]);
      l2project_f_scalars.push_back(L2_projection_Of_f[(int)(CArray_entries[j]) - 1 ]);
    }
    for (int j=0; j<(currFE.get_BFDegree() +1); ++j)
    second_Derivatives_OnRefElem.push_back(currFE.Differentiate_BF(spatialVariable, derivativeOrder2, j+1));  
    u_h = LinearComb(scalars, currFE.get_BFncs());
    l2project_f = LinearComb(l2project_f_scalars, currFE.get_BFncs());
    u_h_2nd_deriv = LinearComb(scalars, second_Derivatives_OnRefElem); 
    for (int j=0; j<GaussianQuadrature::Nodes.size(); ++j) {
      ref_elem_point = { (0.5*GaussianQuadrature::Nodes[j]) + 0.5 };
      h2_NormOf_Resid+= GaussianQuadrature::Weights[j]*(pow( (Eval(ref_elem_point,l2project_f) +a(0.0)*Eval( ref_elem_point, u_h_2nd_deriv) - c(0.0)*Eval(ref_elem_point, u_h))*( sqrt(weight_Fnc(leftP,rightP,Eval(ref_elem_point,currFE.get_AffineMap()[0]))) ) , 2.0));

      h2_NormOf_f+= GaussianQuadrature::Weights[j]*(pow( ( f(Eval(ref_elem_point,currFE.get_AffineMap()[0]))-Eval(ref_elem_point,l2project_f) )*
			                                     ( sqrt(weight_Fnc(leftP,rightP,Eval(ref_elem_point,currFE.get_AffineMap()[0]))) ) , 2.0));
    } 
    h2_NormOf_Resid*=(0.5*currFE.detJ());
    h2_NormOf_f*=(0.5*currFE.detJ());
    eNormErrEst+= C_p*(h2_NormOf_Resid +  h2_NormOf_f);
  }
  return sqrt(eNormErrEst);
}

// returns the coefficients defining the linear combination that represents the l2 projection of the function f.
// Only used during test problem
double* L2_projection_f(Mesh& COMPUTATIONALMESH)
{
  int *CArray_row_start = COMPUTATIONALMESH.ConnArray_rowptr();
  double *CArray_entries = COMPUTATIONALMESH.ConnArray_entryptr();
  int NoOfElements = COMPUTATIONALMESH.get_NoOfElements(), NoOfNodes = COMPUTATIONALMESH.get_NoOfNodes(); 
  std::list<int> HashedNodePairs;
  for (int i=0; i<NoOfElements; ++i) {
    for (int j=(CArray_row_start[i]-1); j<(CArray_row_start[i+1]-1); ++j) {
      for (int k=(CArray_row_start[i]-1); k<(CArray_row_start[i+1]-1); ++k) {
        orderedInsert(HashedNodePairs, hashFunc(CArray_entries[j], CArray_entries[k],NoOfNodes) );          	
      }
    }
  }
  CSRMatrix M;
  M.matrix_entries = new double[HashedNodePairs.size()];
  M.col_no = new int[HashedNodePairs.size()];
  M.row_start = new int[NoOfNodes+1];
  M.NoOfRows = NoOfNodes;
  double* b = new double[NoOfNodes];
  for (int i=0; i<HashedNodePairs.size(); ++i) {
    M.matrix_entries[i] = 0.0;
    M.col_no[i] = 0;
  }
  int row_startIndex = 0, currValIndex=0; 
  std::list<int>::iterator it = HashedNodePairs.begin();
  for (int i=0; i<NoOfNodes; ++i) {
    int lowerBound = i*NoOfNodes;
    while (*it<=lowerBound) {
      ++it;
      ++row_startIndex;
    }
    M.row_start[i] = row_startIndex+1;
    b[i]=0.0; 
  }
  M.row_start[NoOfNodes] = HashedNodePairs.size() +1;
  Finite_Element currFE;
  for (int i=0; i<NoOfElements; ++i) {
    currFE = COMPUTATIONALMESH.find_FE(i+1);
    for (int j=0; j<currFE.get_BFDegree()+1; ++j) {
      double L2_InnerProduct_f_j=0.0;
      std::vector<double> localPoint;
      for (int g=0; g<GaussianQuadrature::Nodes.size(); ++g) {
        localPoint = { (0.5*GaussianQuadrature::Nodes[g]) + 0.5};
        L2_InnerProduct_f_j+=GaussianQuadrature::Weights[g]*f(Eval(localPoint, currFE.get_AffineMap()[0]))*Eval(localPoint, currFE.get_BFncs()[j]); 
      }
      L2_InnerProduct_f_j*=(0.5*currFE.detJ());
      b[((int)(CArray_entries[CArray_row_start[i]-1+j]))-1]+=L2_InnerProduct_f_j;
      for (int k=0; k<currFE.get_BFDegree()+1; ++k) {
        double L2_InnerProduct_j_k = 0.0;
	std::vector<double> localPoint;
	for (int g=0; g<GaussianQuadrature::Nodes.size(); ++g) {
          localPoint = { (0.5*GaussianQuadrature::Nodes[g]) + 0.5};
          L2_InnerProduct_j_k+=GaussianQuadrature::Weights[g]*(Eval(localPoint, currFE.get_BFncs()[j])*Eval(localPoint, currFE.get_BFncs()[k]));
	}
	L2_InnerProduct_j_k*=(0.5*currFE.detJ());
	int positionHash = hashFunc(CArray_entries[CArray_row_start[i]-1+k],
                                    CArray_entries[CArray_row_start[i]-1+j],
			            NoOfNodes );
	it = HashedNodePairs.begin();
	currValIndex = 0;
	while (*it!=positionHash) {
          ++it;
	  ++currValIndex;
	}
	M.matrix_entries[currValIndex]+=L2_InnerProduct_j_k;
	M.col_no[currValIndex] = CArray_entries[CArray_row_start[i]-1+j];
      }
    }
  }
  double *coefficients = new double[NoOfNodes];
  for (int i=0; i<NoOfNodes; ++i)
  coefficients[i]=0.0;
  ConjugateGradientMethod(M, b, coefficients, 1e-20); 
  return coefficients;
}

// returns analytic energy Norm of the error. Only used for code verification.
double analytic_energy_norm_Err(Mesh& COMPUTATIONALMESH, double* coeffVals)
{ 
  double a_=a(0.0), c_=c(0.0), l2Norm_err =0.0, l2Norm_err_deriv=0.0, analyticEnergy_Err=0.0;
  std::vector<double> scalars, u_h, u_h_deriv, svars={1.0}, dorder={1.0}, localPoint;
  std::vector< std::vector<double> > local_BF_derivatives;
  Finite_Element currFE;
  double *CArray_entries = COMPUTATIONALMESH.ConnArray_entryptr();
  int  *CArray_row_start = COMPUTATIONALMESH.ConnArray_rowptr();
  for (int i=0; i<COMPUTATIONALMESH.get_NoOfElements(); ++i) {
    scalars.clear();
    u_h.clear();
    u_h_deriv.clear();
    local_BF_derivatives.clear();
    currFE = COMPUTATIONALMESH.find_FE(i+1);
    l2Norm_err =0.0;
    l2Norm_err_deriv=0.0;
    for (int j=CArray_row_start[i]-1; j<CArray_row_start[i+1]-1; ++j) 
    scalars.push_back(coeffVals[((int)(CArray_entries[j])) - 1]);
    for (int j=0; j<currFE.get_BFDegree()+1; ++j)
    local_BF_derivatives.push_back(currFE.Differentiate_BF(svars, dorder, j+1));
    u_h = LinearComb(scalars,currFE.get_BFncs());
    u_h_deriv = LinearComb(scalars,local_BF_derivatives);
    for (int j=0; j<GaussianQuadrature::Nodes.size(); ++j) {
      localPoint = {(0.5*GaussianQuadrature::Nodes[j])+0.5};
      l2Norm_err+=GaussianQuadrature::Weights[j]*(pow( analytic_solution(Eval(localPoint,currFE.get_AffineMap()[0])) - Eval(localPoint,u_h) , 2.0));
      l2Norm_err_deriv+=GaussianQuadrature::Weights[j]*(pow( analytic_solution_deriv(Eval(localPoint,currFE.get_AffineMap()[0])) - Eval(localPoint,u_h_deriv), 2.0));
    }
    l2Norm_err*=(0.5*currFE.detJ());
    l2Norm_err_deriv*=(0.5*currFE.detJ());
    analyticEnergy_Err+=((a_*l2Norm_err_deriv) + (c_*l2Norm_err));
  }
  return sqrt(analyticEnergy_Err);
}

// returns a sharp upper bound on the H1Norm of the error in the FEM solution defined by coeffVals and COMP'LMESH. C_0 refers to coercivity const
double H1NormofErr_Est(Mesh& COMPUTATIONALMESH, double* coeffVals, const double& C_0)
{
  
  double C_i_1=0.5, local_p=0.0,local_c=0.0, h=COMPUTATIONALMESH.get_h(), localResidNorm =0.0, H1ErrNorm_Est=0.0;
  std::vector<double> scalars, u_h, u_h_deriv1, u_h_deriv2, svars={1.0}, dorder1={1.0}, dorder2={2.0}, localPoint;
  std::vector< std::vector<double> > local_deriv1, local_deriv2;
  double *CArray_entries = COMPUTATIONALMESH.ConnArray_entryptr();
  int *CArray_row_start = COMPUTATIONALMESH.ConnArray_rowptr(); 
  Finite_Element currFE;
  for (int i=0; i<COMPUTATIONALMESH.get_NoOfElements(); ++i) {
    scalars.clear();
    local_deriv1.clear();
    local_deriv2.clear();
    currFE=COMPUTATIONALMESH.find_FE(i+1);
    for (int j=CArray_row_start[i]-1; j<CArray_row_start[i+1]-1; ++j)
    scalars.push_back(coeffVals[((int)(CArray_entries[j]))-1]);
    for (int j=0; j<currFE.get_BFDegree()+1; ++j) {
      local_deriv1.push_back(currFE.Differentiate_BF(svars, dorder1, j+1));
      local_deriv2.push_back(currFE.Differentiate_BF(svars, dorder2, j+1));
    }
    u_h = LinearComb(scalars,currFE.get_BFncs());
    u_h_deriv1 = LinearComb(scalars, local_deriv1);
    u_h_deriv2 = LinearComb(scalars, local_deriv2);
    localResidNorm=0.0;
    local_p=currFE.get_BFDegree();
    local_c=h/(sqrt(local_p*(local_p+1.0)));
    for (int j=0; j<GaussianQuadrature::Nodes.size(); ++j) {
      localPoint = { (0.5*GaussianQuadrature::Nodes[j]) + 0.5};
      localResidNorm+=GaussianQuadrature::Weights[j]*( pow( local_c*( f(Eval(localPoint, currFE.get_AffineMap()[0])) 
				                                      + a(0.0)*Eval(localPoint,u_h_deriv2)                                                                                                                                    - b(Eval(localPoint,currFE.get_AffineMap()[0]))*Eval(localPoint,u_h_deriv1)                                                                                             - c(Eval(localPoint, currFE.get_AffineMap()[0]))*Eval(localPoint,u_h) ) , 2.0) );
    } 
    localResidNorm*=(0.5*currFE.detJ());
    H1ErrNorm_Est+=localResidNorm;
  }
  H1ErrNorm_Est = sqrt(H1ErrNorm_Est);
  H1ErrNorm_Est*=(C_i_1/C_0);
  return H1ErrNorm_Est;
}

// returns 0 if failed. otherqise returns coercivity constant
double get_Coercivity_const()
{
  std::ifstream DataInputStream("Input Docs/Computational Domain.dat",std::ios_base::in);
  double C_pf=0.0;
  if (DataInputStream) {
    // Read input file into block of dynamically allocated memory (buffer)
    DataInputStream.seekg(0,DataInputStream.end);
    int length = DataInputStream.tellg();
    DataInputStream.seekg(0,DataInputStream.beg);
    char* buffer = new char[length];
    DataInputStream.read(buffer,length);

    // Store information from input file. 
    int pos=0;
    char BoundaryText[8] = {'B','O','U','N','D','A','R','Y'};
    char* BP_Addr = std::search(buffer,buffer+length*(sizeof(char)),BoundaryText, BoundaryText +8*sizeof(char)) + 9*sizeof(char);
    std::string numeralContainer;
    while (*BP_Addr!=' ') {
      numeralContainer.push_back(*BP_Addr);
      ++BP_Addr;
    }
    double LeftBoundaryPoint = strtod(numeralContainer.c_str(),0);
    numeralContainer.clear();
    ++BP_Addr;
    while (*BP_Addr!='\n') {
       numeralContainer.push_back(*BP_Addr);
       ++BP_Addr;
    }
    double RightBoundaryPoint = strtod(numeralContainer.c_str(),0);
    C_pf=0.5*(pow(RightBoundaryPoint-LeftBoundaryPoint,2.0));
    double C_0 =a(0.0)/(1.0+C_pf);
    return C_0; 
  }
  else return 0.0;
}

