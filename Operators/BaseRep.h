/*===========================================================================
 *
 *  File: BaseRep.h
 *  Content: Header File declaring class to be used for representation of 
 *  arbitrary base arithmetic
 *  Date: 05/06/2016
 *  Author: Andrew Oldham
 *
 * **************************************************************************/

#ifndef __BASEREP_H_INCLUDED__
#define __BASEREP_H_INCLUDED__

#include <vector> 
#include <cmath>

// Class to be used for representation of arithmetic in any integer base.
class BaseRep
{
  public:
    // Constructor
    BaseRep(const int& base,     // Defines the base used when representing in this class instance 
	    const int& length);  // Defines the fixed number of digits used to represent integers in the given base.
    
    // Default Constructor
    BaseRep(); 
    
    // Copy COnstructor
    BaseRep(const BaseRep& COPY); 
    
    // Assignment Operator
    BaseRep& operator=(BaseRep ASSIGN); 

    // Overloaded Increment Operator
    BaseRep& operator++();        
    
    // sets Representation to 0.
    void zero();
    int get_Base() const { return Base; }

    std::vector<int> Representation;   // Vector that defines the current Integer stored in this class instance, 
                                       // Representation[0] is the least significant bit/trit/quit/.. etc

    // Destructor
    ~BaseRep();  

  private:
    int Base;
};

// Comparison Operators
bool operator==(const BaseRep& BR1, const BaseRep& BR2);
bool operator!=(const BaseRep& BR1, const BaseRep& BR2);

#endif
