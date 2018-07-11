/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


/**
Copyright © 2017 Hermann Detz

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
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
