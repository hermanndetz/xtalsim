/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#ifndef __HASH_H__
#define __HASH_H__

#include <functional>
#include <tuple>
#include <string>

#include <physics/Vector3D.h>
#include <physics/Element.h>

//! \brief Provides methods to calculate hashes.

class Hash{
public:
    //! Generates hash for 3D vector.
    size_t operator()(const Vector3D<indexType> &vector) const;
    size_t operator()(const std::tuple<std::string, elementType> &entry) const;
    size_t operator()(const std::tuple<elementType,elementType> &entry) const;
};

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
