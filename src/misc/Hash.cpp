/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "Hash.h"

//------------------------------------------------------------------------------

//! Compute a unique hash value for a Vector3D object. Even for permutations of
//! the same indices a different value has to be delivered.
//! \param vector 3D Vector whose hash shall be computed.
//! \return Calculated hash value.
size_t Hash::operator()(const Vector3D<indexType> &vector) const
{
    return ( ( std::hash<int>()(vector[0])
	       ^ (std::hash<int>()(vector[1]) << 1) ) >> 1)
	^ (std::hash<int>()(vector[2]) << 1 );
}

//------------------------------------------------------------------------------

//! Compute a unique hash value for a Vector3D object. Even for permutations of
//! the same indices a different value has to be delivered.
//! \param vector 3D Vector whose hash shall be computed.
//! \return Calculated hash value.
size_t Hash::operator()(const std::tuple<std::string, elementType> &entry) const
{

    return ( std::hash<std::string>{}( std::get<0>(entry) + "__" +
                                      std::to_string(std::get<1>(entry))
                                      )
             );
}

//------------------------------------------------------------------------------

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
