/*================================================
 *
 *  File: Headers.h
 *  Content: Header File that includes all necessary
 *  headers in order for the hp-FEM solver to run
 *  successfully
 *  Date: 31/07/2016
 *  Author: A.Oldham
 *
 * **********************************************/

#ifndef __HEADERS_H_INCLUDED__
#define __HEADERS_H_INCLUDED__

#include <iostream>
#include <vector>
#include <utility>
#include <cstdlib>
#include <list>
#include <fstream>
#include <cmath>
#include "General Finite Element Class/Finite_Element.h"
#include "Operators/BaseRep.h"
#include "Operators/Mesh_Generator.h"
#include "Operators/Matrix_Algorithms.h"
#include "Operators/Polynomial_Operators.h"
#include "Mesh/Mesh.h"
#include "Basis Functions/Basis_Functions.h"
#include "Quadratures/GQ_NodesWeights_Generation.h"
#include "Basis Functions/Basis_Functions.h"
#include "Input Docs/Bilinear_Form.h"
#include "Input Docs/PDE_Definition.h"

#endif
