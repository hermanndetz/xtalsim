/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#ifndef __RANGE3D_INL__
#define __RANGE3D_INL__

//! Set all limits to zero and to not apply the limits.
template <class T>
Range3D<T>::Range3D():apply(false, false, false), start((T) 0, (T) 0, (T) 0),
		   stop((T) 0, (T) 0, (T) 0)
{
    // nothing to be done
}

//------------------------------------------------------------------------------

template <class T>
Range3D<T>::Range3D(const Vector3D<T> &_start, const Vector3D<T> &_stop,
	    const Vector3D<bool> &_apply):
    apply(_apply), start(_start), stop(_stop)
{
    // nothing to be done
}

//------------------------------------------------------------------------------

template <class T>
Range3D<T>::Range3D(const T _start[], const T _stop[], const bool _apply[]):
    apply(_apply), start(_start), stop(_stop)
{
    // nothing to be done
}

//------------------------------------------------------------------------------

template <class T>
Range3D<T>::Range3D(const Range3D<T> &range):
    apply(range.apply), start(range.start), stop(range.stop)
{
    // nothing to be done
}
    
//------------------------------------------------------------------------------

template <class T>
Range3D<T>::~Range3D() { }

//------------------------------------------------------------------------------

//! Prints values of the range to a stream.
template <class T>
std::string Range3D<T>::str(void) const
{
    std::ostringstream msg;

    msg << "start: " << start << ", stop: " << stop << ", apply: " << apply;
    return msg.str();
}

//------------------------------------------------------------------------------

//! Checks if the values of the vectors are bigger than the specified limit. If
//! this is the case the corresponding value of the limit is subtracted until
//! the value is within the limits.
//! \param limit The range defining the maximal size.
template <class T>
inline void Range3D<T>::fitInBox(const Vector3D<T> &limit)
{
    start.fitInBox(limit);
    stop.fitInBox(limit);
}

//------------------------------------------------------------------------------

//! \param limit The range defining the maximal size.
//! \param dimension The dimension that shall be checked.
template <class T>
inline void Range3D<T>::fitInBox(const Vector3D<T> &limit,
				  const uint8_t dimension)
{
    start.fitInBox(limit, dimension);
    stop.fitInBox(limit, dimension);
}

//##############################################################################

//! Prints values of the vector to a stream.
template <class T>
std::ostream& operator<<(std::ostream& os, Range3D<T> const & range)
{

    os << "start: " << range.start << ", stop: " << range.stop <<
	", apply: " << range.apply;

    return os;
}

//------------------------------------------------------------------------------

//! Prints values of the vector to a file stream.
template <class T>
std::ofstream& operator<<(std::ofstream& os, Range3D<T> const & range)
{

    os << "start: " << range.start << ", stop: " << range.stop <<
	", apply: " << range.apply;

    return os;
}

//------------------------------------------------------------------------------

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
