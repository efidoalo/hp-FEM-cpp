/*===========================================================
 *
 *  File: Finite_Element.h
 *  Content: Header file for the general finite element class.
 *  Date: 04/06/2016 
 *  Author: Andrew Oldham
 *
 * **********************************************************/

#ifndef __FE_H__
#define __FE_H__

#include "../Operators/CSRMatrix.h"
#include "../Operators/Matrix_Algorithms.h"
#include <vector>
#include <iostream>

class Finite_Element
{
  public:
    // Standard Constructor
    Finite_Element(const int& FE_TypeSpecifier,                              // Specifies the finite element (ie 1D line, 2D Triangle, 3D Tet, etc) 
		   const int& Basis_Function_Degree,                         // Specifies the set of basis functions (ie, linear, cubic etc)
		   const std::vector< std::vector<double> >& Vertex_Points,  // Vector defining the ordered vertices of the finite element.
		   const int& Dim);                                          // Spatial dimension of the element.
    
    // Copy Constructor
    Finite_Element(const Finite_Element& COPY);

    // Assignment Operator
    Finite_Element& operator=( Finite_Element FE_ASSIGN);

    // Defualt Constructor
    Finite_Element();

    // Differentiate Basis function BFIdentifier (1starting) w.r.t the global spatial variables defined by SpatialVars and D_Orders.
    // The returned vector is this derivative represented in terms of local reference element coordinates as a polynomial.
    std::vector<double> Differentiate_BF(const std::vector<double>& SpatialVars, const std::vector<double>& D_Orders, const int& BFIdentifier);
    
    // Returns a the jacobian of the affine mapping (from reference element to global cooridate space that defines this element
    // in the mesh) in the form of a CSRMatrix struct object.
    CSRMatrix Jacobian();

    // Calculates and returns the determinant of the Jacobian matrix of the affine map that maps local reference element coords to global/Mesh coords
    double detJ();

    // Get Data Members
    std::vector< std::vector<double> > get_BFncs() const { return BasisFunctions; }
    std::vector< std::vector<double> > get_AffineMap()  const { return AffineMap; }
    std::vector< std::vector<double> > get_VertexPoints() const { return VertexPoints; }
    std::vector< std::vector<double> > get_InverseAffineMap() const { return InverseAffineMap; }

    int get_BFDegree() const { return BFDegree; }
    int get_Dim() const { return Dimension; }
    int get_FETYPE() const { return FEType; }


  private:
    int Dimension; // Spatial dimension of the finite elment. 
    int FEType;    // 0=Straight Line, 1=2D Triangle, 2=2D Parallelogram, 3=3D Tetrahedron
    int BFDegree;  // Basis Function Index, together with the FEType data member defines the set of basis functions to be used. 
                   // eg FEType=0, BFIndex =2, would induce BasisFunctions[0][2] to be used which are the lagrangian cubic basis 
		   // functions on the straight line, see "Basis_Functions.h" header file for more details of the BasisFunctions container.
    std::vector< std::vector<double> > AffineMap;    // Defines the Affine Mapping from the reference element to global cooridnates.
    std::vector< std::vector<double> > InverseAffineMap; // Defines the inverse Afifne map, mapping global coordinates to the reference element.
    std::vector< std::vector<double> > VertexPoints; // Defines the ordered vertices of this FE element in the mesh.
    std::vector< std::vector<double> > BasisFunctions;    
};

#endif
