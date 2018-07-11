/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __VECTOR_3D_H__
#define __VECTOR_3D_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <typeinfo>

typedef int indexType;
typedef double spaceType;

template <class T>
class Range3D;

//! Coordinates in different data types.

//! Stores x, y and z coordinate of a specific point. Since different data types
//! are required this class is defined as template
template <class T>
class Vector3D{

private:
    //! Coordinates in a three dimensional space.
    T values_[3];
    
public:
    //! Constructor
    Vector3D(const T _x=0, const T _y=0, const T _z=0);
    //! Constructor
    Vector3D(const T *_values);
    //! Copy Constructor
    Vector3D(const Vector3D<T> &vector);
    //! Destructor
    ~Vector3D();

    //! Return squared length of vector
    inline double squaredLength(void) const;
    //! Return angle between two vectors.
    inline double angle(const Vector3D<T> &vector) const;
    //! Return in product of two vectors.
    inline double inProduct(const Vector3D<T> &vector) const;
    //! Checks if the vector is inside a specified range.
    inline bool isInRange(const Range3D<T> &range) const;
    //! Adapt vector such that it lies within given limit.
    inline void fitInBox(const Vector3D<T> &limit);
    //! Adapt vector only in one dimension to lie within given limit.
    inline void fitInBox(const Vector3D<T> &limit, const uint8_t dimension);
    //! Compute smallest difference between two vectors.
    Vector3D<T> getMinimalDistance(const Vector3D<T> &vector,
			     const Vector3D<T> &limit= Vector3D<T>()) const;
    
    //! Check equality of two vectors.
    //    bool operator==(const Vector3D<T> &vector) const;
    //! Check if all entries have the given value.
    //    bool operator==(const T value) const;
    //! Check if all entries are smaller than counterparts in second vector.
    bool operator<(const Vector3D<T> &vector) const;
    //! Check if all entries are smaller than a specific value
    bool operator<(const T value) const;
    //! Check if all entries are smaller or equal their counterparts in second
    //! vector. 
    bool operator<=(const Vector3D<T> &vector) const;
    //! Check if all entries are smaller or equal a specific value
    bool operator<=(const T value) const;
    //! Check if all entries are bigger than counterparts in second vector.
    bool operator>(const Vector3D<T> &vector) const;
    //! Check if all entries are bigger than a specific value
    bool operator>(const T value) const;
    //! Set values according to parameter vector.
    void operator=(const Vector3D<T> &vector);
    //! Set values according to parameter array.
    void operator=(const T *_values);
    //! Return value specified by value
    T & operator[](uint8_t index);
    //! Return value specified by value
    T const & operator[](uint8_t index) const;
    
    //! Multiply by constant.
    void operator*=(const double factor);
    //! Add vector to current one.
    void operator+=(const Vector3D<T> &vector);
    //! Add numerical value to each direction.
    void operator+=(const T value);
    //! Add vector to current one and return result.
    Vector3D<T> operator+(const Vector3D<T> &vector) const;
    //! Add numerical value to each direction and return result.
    Vector3D<T> operator+(const T value) const;
    //! Subtract vector from current one.
    void operator-=(const Vector3D<T> &vector);
    //! Subtract numerical value from each direction.
    void operator-=(const T value);
    //! Subtract vector from current one and return result.
    Vector3D<T> operator-(const Vector3D<T> &vector) const;
    //! Subtract numerical value from each direction and return result.
    Vector3D<T> operator-(const T value) const;
    
    //! Return values as string "(x,y,z)"
    const std::string str(void) const;
    //! Return values as string "x <tabulator> y <tabulator> z"
    const std::string strTab(void) const;
};

// Forward declaration. Only with this constellation it was possible to specify
// specific double implementations which very necessary to prevent <double>
// comparison warnings.
bool operator==(const Vector3D<double> &vector1, const double value);
bool operator==(const Vector3D<double> &vector1,
		const Vector3D<double> &vector2);

// Unfortunately it is necessary to include the code of template member
// functions in the header file.
// Other solutions are documented at
// http://www.codeproject.com/Articles/48575/How-to-define-a-template-class-in-a-h-file-and-imp 
#include "Vector3D.inl"

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
