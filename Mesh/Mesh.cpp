/*=============================================;
 *
 *  File: Mesh.cpp
 *  Content: Source File defining the Mesh class data structure.
 *  Date: 10/06/2016
 *  Author: A.Oldham
 *
 **********************************************/ 

#include "Mesh.h"
#include "../Operators/CSRMatrix.h"
#include "../General Finite Element Class/Finite_Element.h"
#include <unordered_map>

// Find and return finite element number index in the mesh. (index is 1 starting integer) 
Finite_Element Mesh::find_FE(const int& index) 
{
  return ElementContainer[index];
}

// Copy constructor
Mesh::Mesh(const Mesh& MESHCOPY)
{
  NoOfElements = MESHCOPY.get_NoOfElements();
  NoOfNodes = MESHCOPY.get_NoOfNodes();
  ConnectivityArray = MESHCOPY.get_ConnArray();
  ElementContainer = MESHCOPY.get_FEContainer();
  h = MESHCOPY.get_h();
}

// Assignment Opertor
Mesh& Mesh::operator=(Mesh MESHASSIGN)
{  
  NoOfElements = MESHASSIGN.get_NoOfElements();
  NoOfNodes = MESHASSIGN.get_NoOfNodes();
  ConnectivityArray = MESHASSIGN.get_ConnArray();
  ElementContainer = MESHASSIGN.get_FEContainer();
  h= MESHASSIGN.get_h();
  return *this;
}

// redefine the mesh
void Mesh::redefineMESH(const int& NoOFElements, const int& NoOFNodes, std::unordered_map<int, Finite_Element>& EC, CSRMatrix& CArray, double& h_)
{
  NoOfElements = NoOFElements;
  NoOfNodes = NoOFNodes;
  ConnectivityArray = CArray;
  ElementContainer = EC;
  h=h_;
}
//Destructor 
Mesh::~Mesh()
{
  
}
 
