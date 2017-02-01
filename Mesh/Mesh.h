/*=============================================;
 *
 *  File: Mesh.h
 *  Content: Header File for the Mesh class data structure.
 *  Date: 09/06/2016
 *  Author: A.Oldham
 *
 **********************************************/  

#ifndef __MESH_H_INCLUDED__
#define __MESH_H_INCLUDED__

#include "../Operators/CSRMatrix.h"
#include "../General Finite Element Class/Finite_Element.h"
#include "../Operators/Matrix_Algorithms.h"
#include <unordered_map>

// Mesh Class Declaration
class Mesh
{
  public:
    // Constructor
    Mesh(const int& NumberOfElements, const int& NumberOfNodes,
	 const std::unordered_map<int,Finite_Element>& Container,
	 const CSRMatrix& ConnectivityMatrix, const double& h_)
    : NoOfElements(NumberOfElements), NoOfNodes(NumberOfNodes),
    ElementContainer(Container), ConnectivityArray(ConnectivityMatrix),h(h_) { }

    // Default Constructor (ElementContainer and CSRMatrix initialized using default constructors)
    Mesh()
    : NoOfElements(0), NoOfNodes(0)
    { }
    
    // Copy Constructor
    Mesh(const Mesh& MESHCOPY);

    // Assignment Opertor
    Mesh& operator=(Mesh MESHASSIGN);

    // Find Finite_Element number index in the mesh, where index is a 1 starting integer representing a unique element in the mesh.
    // NOT RANGE CHECKED. 
    Finite_Element find_FE(const int& index);
    
    // get data members
    CSRMatrix get_ConnArray() const { return ConnectivityArray; }
    int get_NoOfNodes() const { return NoOfNodes; }
    int get_NoOfElements() const { return NoOfElements; }
    double get_h() const { return h; }
    std::unordered_map<int,Finite_Element> get_FEContainer() const { return ElementContainer; }
    int* ConnArray_rowptr() { return ConnectivityArray.row_start; } // returns row_start data member of ConnectivityArray
    double* ConnArray_entryptr() { return ConnectivityArray.matrix_entries; } // returns matrix_entries data member of ConnectivityArray

    // Redefines the current mesh
    void redefineMESH(const int& NoOfElements, const int& NoOfNodes, std::unordered_map<int, Finite_Element>& EC, CSRMatrix& CArray, double& h_);

    // Destructor
    ~Mesh();
    
  private:
    int NoOfElements;
    int NoOfNodes;
    std::unordered_map<int,Finite_Element> ElementContainer;
    // Maximum distance between vertex points.
    double h;
    CSRMatrix ConnectivityArray;
};
#endif
