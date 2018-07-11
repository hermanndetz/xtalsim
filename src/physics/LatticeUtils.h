/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __LATTICE_UTILS_H__
#define __LATTICE_UTILS_H__

#include <cmath>

namespace Zincblende{

    //! Converts squared bond length to lattice constant in zincblende lattice.
    double squaredBondLengthToLatticeConstant (double squaredBondLength);
    //! Converts lattice constant to bond length in zincblende lattice.
    double latticeConstantToBondLength (double latticeConstant);
}

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
