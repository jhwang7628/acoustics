////////////////////////////////////////////////////////////////////////////////
// File: type_trait.hpp
// Copyright (c) 2007 by Changxi Zheng
////////////////////////////////////////////////////////////////////////////////

#ifndef CARBINE_TYPE_TRAIT_HPP
#   define CARBINE_TYPE_TRAIT_HPP

#include <complex>
#include "null_type.hpp"

namespace carbine
{
    template<typename T>
    struct TrivialType
    { typedef NullType    type; };

    template<>
    struct TrivialType<double>
    { typedef double      type; };

    template<>
    struct TrivialType<float>
    { typedef float       type; };

    template<>
    struct TrivialType< std::complex<double> >
    { typedef double      type; };

    template<>
    struct TrivialType< std::complex<float> >
    { typedef float       type; };
}

#endif
