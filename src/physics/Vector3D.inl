/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __VECTOR_3D_INL__
#define __VECTOR_3D_INL__

template <class T>
Vector3D<T>::Vector3D(const T x, const T y, const T z)
{
    values_[0]=x;
    values_[1]=y;
    values_[2]=z;
}

//------------------------------------------------------------------------------

template <class T>
Vector3D<T>::Vector3D(const T *values)
{
    values_[0]=values[0];
    values_[1]=values[1];
    values_[2]=values[2];
}

//------------------------------------------------------------------------------

template <class T>
Vector3D<T>::Vector3D(const Vector3D<T> &vector)
{
    for (int i=0; i<3; i++){
	this->values_[i] = vector[i];
    }
}

//------------------------------------------------------------------------------

template <class T>
Vector3D<T>::~Vector3D() {}
    
//------------------------------------------------------------------------------

//! Calculates the squared length of the vector, i.e. length^2. This prevents
//! the usage of sqrt() which is very costly function and helps to save some
//! time.
//! \return (Length of vector)^2
template <class T>
inline double Vector3D<T>::squaredLength(void) const
{
    return(values_[0]*values_[0]+values_[1]*values_[1]+values_[2]*values_[2]);
}

//------------------------------------------------------------------------------

//! Calculates the angle between the current vector and a specified one. 
//! \param vector Vector to which the angle shall be computed.
//! \return Angle between two vectors in radians.
template <class T>
inline double Vector3D<T>::angle(const Vector3D<T> &vector) const
{

    double inProduct;
    
    if ( ( (*this) == (T)0 ) || (vector == (T)0) )
	return 0;

    inProduct = this->inProduct(vector);
    return acos( (inProduct) /
		 (sqrt(this->squaredLength()) * sqrt(vector.squaredLength())) );
}

//------------------------------------------------------------------------------

//! Calculates the in-product of the current vector with a specified one.
//! \param vector Second vector taken for in-product.
//! \return In product of the two vectors.
template <class T>
inline double Vector3D<T>::inProduct(const Vector3D<T> &vector) const
{
    double result=0;

    for(int i=0; i<3; i++){
	result += this->values_[i]*vector[i];
    }

    return result;
}

//------------------------------------------------------------------------------

//! \param range The range defining the search space.
//! \return True if coordinates inside given range, false otherwise.
template <class T>
inline bool Vector3D<T>::isInRange(const Range3D<T> &range) const
{
    
    for (int i=0; i<3; i++){
	if (! range.apply[i])
	    continue;

	if (range.start[i] < range.stop[i]){
	    if ((values_[i] < range.start[i]) || (values_[i] > range.stop[i]))
		return false;
	}
	else{
	    if ((values_[i] < range.start[i]) && (values_[i] > range.stop[i]))
		return false;
	}
    }

    return true;
    
}

//------------------------------------------------------------------------------

//! Checks if the values of the vector is bigger than the specified limit. If
//! this is the case the corresponding value of the limit is subtracted until
//! the value is within the limits.
//! \param limit The range defining the maximal size.
template <class T>
inline void Vector3D<T>::fitInBox(const Vector3D<T> &limit)
{

    if (limit.squaredLength() > 0) {
    for (int i=0; i<3; i++){
	while( values_[i] >= limit[i] )
	    values_[i] -= limit[i];

	while( values_[i] < 0 )
	    values_[i] += limit[i];
    }
    }
    
}

//------------------------------------------------------------------------------

//! \param limit The range defining the maximal size.
//! \param dimension The dimension that shall be checked.
template <class T>
inline void Vector3D<T>::fitInBox(const Vector3D<T> &limit,
				  const uint8_t dimension)
{
    if (limit.squaredLength() > 0) {
    while( values_[dimension] >= limit[dimension] )
	values_[dimension] -= limit[dimension];

    while( values_[dimension] < 0 )
	values_[dimension] += limit[dimension];
    }
}

//------------------------------------------------------------------------------

//! In this case a periodic box in all directions is assumed. Therefore the
//! plain difference between the vectors may not result to the minimal distance.
//! \param vector Second vector.
//! \param limit Size of periodic box.
//! \return 3D Vector determining the minimial distance.
template <class T>
Vector3D<T> Vector3D<T>::getMinimalDistance(const Vector3D<T> &vector,
					      const Vector3D<T> &limit) const
{
    Vector3D<T> result;

    // If no box specified do not use periodic boundaries.
    if (limit==0)
	return (*this - vector);

    for (int i=0; i<3; i++){
	T tmp;
	
	result[i] = vector[i] - values_[i];

	tmp = vector[i] - limit[i] - values_[i];
	if (std::abs(result[i]) > std::abs(tmp))
	    result[i] = tmp;

	tmp = vector[i] + limit[i] - values_[i];
	if (std::abs(result[i]) > std::abs(tmp))
	    result[i] = tmp;
    }

    return (result);

}

//------------------------------------------------------------------------------

//! \param vector Second operand for equality.
//! \return True if values of vectors are equal, false otherwise
template <class T>
bool operator==(const Vector3D<T> &vector1, const Vector3D<T> &vector2)
{
    
    for (int i=0; i<3; i++){
	// if (typeid(vector[0]).name() == "f"){
	//     if (! (FP_ZERO == std::fpclassify(this->values_[i] - vector[i])) )
	// 	return false;
	// }
	// else{
	if (vector1[i] != vector2[i])
	    return false;
	    //}
    }

    return true;
}

//------------------------------------------------------------------------------

//! \return True if all values of vector are equal to parameter.
template <class T>
bool operator==(const Vector3D<T> &vector1, const T value)
//bool Vector3D<T>::operator==(const T value) const
{
    for(int i=0; i<3; i++){
	if (vector1[i] != value)
	    return false;
    }

    return true;
}

//------------------------------------------------------------------------------

//! \param vector Vector that shall be assigned to current one.
template <class T>
void Vector3D<T>::operator=(const Vector3D<T> &vector)
{
    for(int i=0; i<3; i++){
	this->values_[i] = vector[i];
    }
}

//------------------------------------------------------------------------------

//! \warning The parameter has to contain at least three elements. Otherwise the
//! border of the array is passed and arbitrary values are read.
//! \param values Array containing at least three values.
template <class T>
void Vector3D<T>::operator=(const T *values)
{
    for(int i=0; i<3; i++){
	this->values_[i] = values[i];
    }
}
  
//------------------------------------------------------------------------------

//! \param vector Vector containing reference values.
//! \return True if all values of vector are smaller than corresponding values
//! of the second vector.
template <class T>
bool Vector3D<T>::operator<(const Vector3D<T> &vector) const
{

    for(int i=0; i<3; i++){
	if (this->values_[i] >= vector[i])
	    return false;
    }

    return true;
}

//------------------------------------------------------------------------------

//! \param value Reference value.
//! \return True if all values of vector are smaller than given value.
template <class T>
bool Vector3D<T>::operator<(const T value) const
{

    for(int i=0; i<3; i++){
	if (this->values_[i] >= value)
	    return false;
    }

    return true;
}
    
//------------------------------------------------------------------------------

//! \param vector Vector containing reference values.
//! \return True if all values of vector are smaller than corresponding values
//! of the second vector.
template <class T>
bool Vector3D<T>::operator<=(const Vector3D<T> &vector) const
{

    for(int i=0; i<3; i++){
	if (this->values_[i] > vector[i])
	    return false;
    }

    return true;
}

//------------------------------------------------------------------------------

//! \param value Reference value.
//! \return True if all values of vector are smaller than given value
template <class T>
bool Vector3D<T>::operator<=(const T value) const
{

    for(int i=0; i<3; i++){
	if (this->values_[i] > value)
	    return false;
    }

    return true;
}

//------------------------------------------------------------------------------

//! \param vector Vector containing reference values.
//! \return True if all values of vector are smaller than corresponding values
//! of the second vector.
template <class T>
bool Vector3D<T>::operator>(const Vector3D<T> &vector) const
{

    for(int i=0; i<3; i++){
	if (this->values_[i] <= vector[i])
	    return false;
    }

    return true;
}

//------------------------------------------------------------------------------

//! \param value Reference value.
//! \return True if all values of vector are bigger than given value
template <class T>
bool Vector3D<T>::operator>(const T value) const
{

    for(int i=0; i<3; i++){
	if (this->values_[i] <= value)
	    return false;
    }

    return true;
}

//------------------------------------------------------------------------------

//! Prints values of the vector to a stream.
template <class T>
const std::string Vector3D<T>::str(void) const
{
    std::ostringstream msg;

    msg << '(' << values_[0] << ',' << values_[1] << ',' << values_[2] << ')';
    return msg.str();
}

//------------------------------------------------------------------------------

//! Prints values of the vector to a stream separated by tabs.
template <class T>
const std::string Vector3D<T>::strTab(void) const
{
    std::ostringstream msg;

    msg << values_[0] << '\t' << values_[1] << '\t' << values_[2];
    return msg.str();
}

//------------------------------------------------------------------------------

//! Return value of Vector as reference, making it possible to assign a value to
//! it directly.
//! \param index Value that shall be returned.
template <class T>
T & Vector3D<T>::operator[](uint8_t index)
{
    return values_[index];
}

//------------------------------------------------------------------------------

//! Return value of Vector as constant reference.
//! \param index Value that shall be returned.
template <class T>
T const & Vector3D<T>::operator[](uint8_t index) const
{
    return values_[index];
}

//------------------------------------------------------------------------------

//! Multiply whole vector with a constant factor
//! \param factor Constant factor used for multiplication.
template <class T>
void Vector3D<T>::operator*=(const double factor)
{
    for(int i=0; i<3; i++){
	values_[i] = (T)(values_[i]*factor);
    }
}

//------------------------------------------------------------------------------

//! Adds a vector to the current one and stores the result in the current
//! vector.
//! \param vector Vector that shall be added.
template <class T>
void Vector3D<T>::operator+=(const Vector3D<T> &vector)
{
    for(int i=0; i<3; i++){
	this->values_[i] += vector[i];
    }
}

//------------------------------------------------------------------------------

//! Adds a value to each direction of the vector.
//! \param value Value that shall be added
template <class T>
void Vector3D<T>::operator+=(const T value)
{
    for(int i=0; i<3; i++){
	this->values_[i] += value;
    }
}

//------------------------------------------------------------------------------

//! Adds a vector to the current one and returns result;
//! \param vector Vector that shall be added.
template <class T>
Vector3D<T> Vector3D<T>::operator+(const Vector3D<T> &vector) const
{
    Vector3D<T> result;
    
    for(int i=0; i<3; i++){
	result[i] = this->values_[i] + vector[i];
    }

    return result;
}

//------------------------------------------------------------------------------

//! Adds a value to each direction of the vector and returns result.
//! \param value Value that shall be added
template <class T>
Vector3D<T> Vector3D<T>::operator+(const T value) const
{
    Vector3D<T> result;
    
    for(int i=0; i<3; i++){
	result[i] = this->values_[i] + value;
    }

    return result;
}

//------------------------------------------------------------------------------

//! Subtracts a vector to the current one and stores the result in the current
//! vector.
//! \param vector Vector that shall be subtracted.
template <class T>
void Vector3D<T>::operator-=(const Vector3D<T> &vector)
{
    for(int i=0; i<3; i++){
	this->values_[i] -= vector[i];
    }
}

//------------------------------------------------------------------------------

//! Subtracts a value to each direction of the vector.
//! \param value Value that shall be subtracted.
template <class T>
void Vector3D<T>::operator-=(const T value)
{
    for(int i=0; i<3; i++){
	this->values_[i] -= value;
    }
}

//------------------------------------------------------------------------------

//! Subtracts a vector to the current one and returns result;
//! \param vector Vector that shall be subtracted.
template <class T>
Vector3D<T> Vector3D<T>::operator-(const Vector3D<T> &vector) const
{
    Vector3D<T> result;
    
    for(int i=0; i<3; i++){
	result[i] = this->values_[i] - vector[i];
    }

    return result;
}

//------------------------------------------------------------------------------

//! Subtracts a value to each direction of the vector and returns result.
//! \param value Value that shall be subtracted.
template <class T>
Vector3D<T> Vector3D<T>::operator-(const T value) const
{
    Vector3D<T> result;
    
    for(int i=0; i<3; i++){
	result[i] = this->values_[i] - value;
    }
    
    return result;
}

//------------------------------------------------------------------------------

//! Prints values of the vector to a stream.
//! \param vector The fector to be printed.
template <class T>
std::ostream& operator<<(std::ostream& os, Vector3D<T> const & vector)
{

    os << "(" << vector[0] << "," << vector[1] << "," << vector[2] << ")";

    return os;
}

//------------------------------------------------------------------------------

//! Prints values of the vector to a file stream.
//! \param vector The fector to be printed.
template <class T>
std::ofstream& operator<<(std::ofstream& os, Vector3D<T> const & vector)
{

    os << "(" << vector[0] << "," << vector[1] << "," << vector[2] << ")";

    return os;
}

//##############################################################################

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
