/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __RANGE_3D_H__
#define __RANGE_3D_H__

#include <string>
#include <sstream>

#include <physics/Vector3D.h>

//! Value ranges in three dimensions and configurable data type.

//! Defines for three dimensions a lower and upper limit as well as a boolean
//! variable which defines if the range shall be applied.
template <class T>
class Range3D{
 public:

    //! If true apply defined limit in that dimension.
    Vector3D<bool> apply;
    //! Lower limit in all three dimensions.
    Vector3D<T> start;
    //! Upper limit in all three dimensions.
    Vector3D<T> stop;

    //! Constructor
    Range3D();
    //! Constructor
    Range3D(const Vector3D<T> &_start, const Vector3D<T> &_stop,
	    const Vector3D<bool> &_apply);
    //! Constructor
    Range3D(const T _start[], const T _stop[], const bool _apply[]);
    //! Copy Constructor
    Range3D(const Range3D<T> &range);
    //! Destructor
    ~Range3D();

    //! Adapt vector such that it lies within given limit.
    inline void fitInBox(const Vector3D<T> &limit);
    //! Adapt vector only in one dimension to lie within given limit.
    inline void fitInBox(const Vector3D<T> &limit, const uint8_t dimension);
    
    //! Return values as string.
    std::string str(void) const;
};

#include "Range3D.inl"

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
