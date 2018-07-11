/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __TUPLES_H_
#define __TUPLES_H_

#include <cstdint>
#include <tuple>

//! Tuple to be used for journals etc. 
//! Uses long long unsigned int to be compatible with pugixml.
typedef std::tuple<long long unsigned int, double> UDTuple;

#endif


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
