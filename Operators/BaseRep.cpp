/*===========================================================================
 *
 *  File: BaseRep.cpp
 *  Content: Source .cpp file defining class to be used for representation 
 *  of arbitrary base arithmetic
 *  Date: 05/06/2016
 *  Author: Andrew Oldham
 *
 * **************************************************************************/

#include "BaseRep.h"
#include <iostream>
// Constructor for the BaseRep Class
BaseRep::BaseRep(const int& base, const int& length)
: Base(base) 
{
  // Initialize Representation to zero.
  for (int i=0; i<length; ++i)
  Representation.push_back(0);

}

// Default Constructor
BaseRep::BaseRep()
{
  Base = 0;
}

// Copy COnstructor
BaseRep::BaseRep(const BaseRep& COPY)
{
  Base = COPY.get_Base();
  for (int i=0; i<COPY.Representation.size(); ++i)
  Representation.push_back(COPY.Representation[i]);
}

// Assignnment Operator
BaseRep& BaseRep::operator=(BaseRep ASSIGN)
{
  Base = ASSIGN.get_Base();
  for (int i=0; i<ASSIGN.Representation.size(); ++i) 
  Representation.push_back(ASSIGN.Representation[i]);
  return *this;
}

// Increment Representation data member by 1. Undefined behaviour for incrementing past base^(length) [ inclusive of base^length ]
BaseRep& BaseRep::operator++() 
{ 
  bool carry=false;
  if (Representation[0]<(Base-1)) {
    Representation[0] = Representation[0]+1;
  }
  else {
    Representation[0] = 0;
    carry = true;
    int RepIndex=1;
    while (carry) {
      if (RepIndex>=Representation.size())
      return *this; 
      if (Representation[RepIndex]<(Base-1)) {
        Representation[RepIndex] = Representation[RepIndex]+1;
	carry=false;
      }
      else {
        Representation[RepIndex]=0;
        ++RepIndex;
      }	
    }
  }
  return *this;
}

// Comparison Operators
bool operator==(const BaseRep& BR1, const BaseRep& BR2)
{
  bool val = true;
  for (int i=0; i<BR1.Representation.size(); ++i) {
    if (BR1.Representation[i]!=BR2.Representation[i])
    val=false;
  }
  return val;
}
bool operator!=(const BaseRep& BR1, const BaseRep& BR2)
{
  return !(BR1==BR2);
}

// sets Representation to 0.
void BaseRep::zero()
{
  for (int i=0; i<Representation.size(); ++i)
  Representation[i] = 0;
}	

// Destructor
BaseRep::~BaseRep()
{
}







