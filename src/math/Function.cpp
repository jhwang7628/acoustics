//////////////////////////////////////////////////////////////////////
// Function.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "Function.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
Function::Function()
  : _derivatives( 0 )
{
}

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
Function::Function( const RealFunction &function,
                    const vector<RealFunction> &derivatives )
  : _function( function ),
    _derivatives( derivatives )
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
Function::~Function()
{
}

