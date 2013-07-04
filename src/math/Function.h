//////////////////////////////////////////////////////////////////////
// Function.h: Interface for the function class
//
//////////////////////////////////////////////////////////////////////

#ifndef FUNCTION_H
#define FUNCTION_H

#include <TYPES.h>

#include <utils/Evaluator.h>

#include <vector>

//////////////////////////////////////////////////////////////////////
// Function class
//
// Simple class for storing a function pointer and a pointer to
// function derivatives
//////////////////////////////////////////////////////////////////////
class Function {
	public:
		Function();

        Function( const RealFunction &function,
                  const std::vector<RealFunction> &derivatives );

        // Destructor
        virtual ~Function();

        RealFunction                 &function()
                                      {
                                        return _function;
                                      }

        RealFunction                 &derivative( int order )
                                      {
                                        return _derivatives[ order - 1 ];
                                      }

        std::vector<RealFunction>    &derivatives()
                                      {
                                        return _derivatives;
                                      }

        int                          &supportSamples()
                                      {
                                        return _supportSamples;
                                      }

        REAL                         &support()
                                      {
                                        return _support;
                                      }

	protected:

	private:
        RealFunction                 _function;
        std::vector<RealFunction>    _derivatives;

        int                          _supportSamples;
        REAL                         _support;

};

#endif
