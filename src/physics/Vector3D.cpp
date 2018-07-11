/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#include "Vector3D.h"

//------------------------------------------------------------------------------

//! Specific version for the parameter double since it throws a compiler warning
//! if double values are compared directly.
//! \param vector Second operand for equality.
//! \return True if values of vectors are equal, false otherwise
bool operator== (const Vector3D<double> &vector1,
		 const Vector3D<double> &vector2)
{

    for (int i=0; i<3; i++){
	if (! (FP_ZERO == std::fpclassify(vector1[i] - vector2[i])) )
	    return false;
    }

    return true;
}

//------------------------------------------------------------------------------

//! \return True if all values of vector are equal to parameter.
bool operator==(const Vector3D<double> &vector1, const double value)
{
    for(int i=0; i<3; i++){
	if (! (FP_ZERO == std::fpclassify(vector1[i] -  value) ) )
	    return false;
    }

    return true;
}

//------------------------------------------------------------------------------

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
