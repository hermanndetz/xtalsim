/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "LatticeUtils.h"

//------------------------------------------------------------------------------

//! \param squaredBondLength (bond length) * (bond length)
double Zincblende::squaredBondLengthToLatticeConstant (
						 double squaredBondLength)
{

    return std::sqrt(squaredBondLength) / std::sqrt(3.0 * 0.25 * 0.25);
    
}   

//------------------------------------------------------------------------------

//! \param latticeConstant Lattice constant that shall be transformed.
double Zincblende::latticeConstantToBondLength (double latticeConstant)
{

    return latticeConstant * std::sqrt(3.0 * 0.25 * 0.25);
    
}   

//------------------------------------------------------------------------------

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
