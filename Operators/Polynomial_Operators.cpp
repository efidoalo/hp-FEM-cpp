/*=================================================;
 *
 *  File: Polynomial_Operators.cpp
 *  Content: Source file defining the Functions/Routines 
 *  that perform operations of Polynomials, whereby
 *  a polynomial is represented as a vector<double>
 *  Date: 12/06/2016
 *  Author: Andrew Oldham
 *
 **************************************************/

#include "Polynomial_Operators.h"
#include "BaseRep.h"
#include <vector>
#include <utility>
#include <iostream>
// Recursive multiplication, modified factorial function. Returns Operand*(Operand-1)*....*(Operand-(Condition-1)) [Operand>=Condition]
int fact1(int Operand, const int& Condition)
{
  int returnVal = 1;
  for (int i=0; i<Condition; ++i) 
  returnVal*=(Operand-i);
  return returnVal;
}	

int index(const std::vector< std::pair<BaseRep,int> >& Map, const BaseRep& Key)
{  
  for (int i=0; i<Map.size(); ++i) {
    if ((Map[i].first)==Key) {
      return Map[i].second; 
    }
  }
}

// Differentiates the multivariate polynomial given by Operand w.r.t the spatial variable given by IndependentVariable recursively Order times.
std::vector<double> Differentiate(const int& IndependentVariable,const int& Order, std::vector<double> Operand)
{                 
  int NoOfSpatialVariables = Operand[2];
  int k = Operand[0];                // Defines the polynomial space, either P_k or Q_k
  std::vector<double> returnPoly;    // Polynomial to be returned. ( differential of Operand)
  int indextrack =0;                 // integer to track the index of the original Polynomial Operand
  std::vector< std::pair<BaseRep,int> > BaseReptoIndex;     
  for (int i=0; i<3; ++i)
  returnPoly.push_back(Operand[i]); 
  int upperbound = (int)(pow(k+1,NoOfSpatialVariables));
  
  // If the polynomial space is P_k (as opposed to Q_k, a tensor product)
  if (Operand[1]==0) {
    BaseRep CoeffTrack(k+1,NoOfSpatialVariables);
    BaseRep Derivative(k+1,NoOfSpatialVariables);
    int currentdegree=0;
    for (int i=0; i<upperbound; ++i) {
      currentdegree=0;
      for (int j=0; j<CoeffTrack.Representation.size(); ++j)
      currentdegree+=CoeffTrack.Representation[j];
      
      if (currentdegree<=k) {     
	std::pair<BaseRep,int> currPair(CoeffTrack,indextrack);
	BaseReptoIndex.push_back(currPair); 
	returnPoly.push_back(0.0);
        ++indextrack;
	if ( ((CoeffTrack.Representation[CoeffTrack.Representation.size()-IndependentVariable])>(Order-1)) && (Operand[ indextrack + 2]!=0) )	{
           int spatialDIndex = CoeffTrack.Representation.size()-IndependentVariable;
	   for (int j=0; j<CoeffTrack.Representation.size(); ++j) {
             if (j!=spatialDIndex) 
	     Derivative.Representation[j] = CoeffTrack.Representation[j];
	     else
	     Derivative.Representation[j] = CoeffTrack.Representation[j] - Order;
	   }
	   int DerivativeIndex = index(BaseReptoIndex, Derivative);
	   returnPoly[DerivativeIndex+ 3] = Operand[ indextrack + 2]*((double)(fact1(CoeffTrack.Representation[spatialDIndex], Order)));
	}
      }
      ++CoeffTrack;
    }
   // std::cout<<"hf "<<returnPoly.size()<<" "<<returnPoly[0]<<returnPoly[1]<<returnPoly[2]<<returnPoly[3]<<returnPoly[4]<<" ";
    return returnPoly;   
  }

  // If the Polynomial space is Q_k (Tensor Product)
  if (Operand[1]==1) {
    BaseRep CoeffTrack(k+1, NoOfSpatialVariables);
    BaseRep Derivative(k+1, NoOfSpatialVariables);
    for (int i=0; i<(int)pow(k+1,NoOfSpatialVariables); ++i) {
      std::pair<BaseRep,int> currPair(CoeffTrack,i);
      BaseReptoIndex.push_back(currPair);
      returnPoly.push_back(0.0);
      int spatialDIndex = CoeffTrack.Representation.size()-IndependentVariable;
      if ( (CoeffTrack.Representation[spatialDIndex]>(Order-1)) &&
             (Operand[ i + 3]!=0) )	{
        for (int k=0; k<CoeffTrack.Representation.size(); ++k) {
          if (k!=spatialDIndex) 
	  Derivative.Representation[k] = CoeffTrack.Representation[k];
	  else
	  Derivative.Representation[k] = CoeffTrack.Representation[k] - Order;
	}
	int DerivativeIndex = index(BaseReptoIndex, Derivative);
        returnPoly[DerivativeIndex + 3] = Operand[ i + 3]*fact1(CoeffTrack.Representation[spatialDIndex], Order);
      }
      ++CoeffTrack;
    }
    return returnPoly;
  }
}

// Evaluate the multivariate polynomial Poly at the point Point, [Point.size()=Poly[2]].
double Eval(const std::vector<double>& Point, const std::vector<double>& Poly)
{
  int k = Poly[0]; // Poly is in the polynomial space P_k or Q_k
  int indextrack = 0, localdegree=0;
  double eval=0.0, localval=1.0;
  BaseRep CoeffTrack(k+1, Poly[2]);
  if (Poly[0]==0) 
  return Poly[3];
  
  if (Poly[1]==0) {
    for (int i=0; i<(int)pow(k+1,Poly[2]); ++i) {
      localdegree=0;
      for (int k1=0; k1<CoeffTrack.Representation.size(); ++k1)
      localdegree+=CoeffTrack.Representation[k1];
      if (localdegree<k+1) {
	localval=1.0;
        for (int j=0; j<Point.size(); ++j) 
        localval*=pow(Point[j],CoeffTrack.Representation[CoeffTrack.Representation.size()-j-1]);
        eval+=localval*Poly[3+indextrack];
	++indextrack;
      }
      ++CoeffTrack;
    }
    return eval;
  }
  if (Poly[1]==1) { 
    for (int i=0; i<(int)pow(k+1,Poly[2]); ++i) {
      localval=1.0;
      for (int j=0; j<Point.size(); ++j) 
      localval*=pow(Point[j],CoeffTrack.Representation[CoeffTrack.Representation.size()-j-1]);
      eval+=localval*Poly[3+indextrack];
      ++indextrack; 
      ++CoeffTrack;
    }
    return eval;
  }
}

// Linear Combination of multivariate polynomials (of identical degree and numbe of spatial variables)
// scalars.size() must equal OrderedPolynomials.size()
std::vector<double> LinearComb(std::vector<double> scalars, std::vector< std::vector<double> > OrderedPolynomials)
{
  std::vector<double> LinearCombination;
  for (int i=0; i<3; ++i) 
  LinearCombination.push_back((OrderedPolynomials[0])[i]); // Polynomials all of identical degree and number of spatial variables

  for (int i=0; i<OrderedPolynomials[0].size()-3; ++i) {	  
    LinearCombination.push_back(0.0);
    for (int j=0; j<OrderedPolynomials.size(); ++j) 
    LinearCombination[3+i] += scalars[j]*OrderedPolynomials[j][3+i]; 
  }
  return LinearCombination;

}
