/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "ArrheniusFactor.h"

double ArrheniusFactor::Boltzmann(const double energyDifference,
				  const double temperature)
{
    if (temperature < temperatureLimit)
	return 0;
    
    return scaling*exp(-1 * energyDifference / (kB * temperature));
}

//------------------------------------------------------------------------------

double ArrheniusFactor::BoseEinstein(const double energyDifference,
				  const double temperature)
{
    if (temperature < temperatureLimit)
	return 0;
    
    return scaling*(1.0 / (exp(energyDifference / (kB * temperature)) -1 ) );
}

//------------------------------------------------------------------------------

double ArrheniusFactor::FermiDirac(const double energyDifference,
				  const double temperature)
{
    if (temperature < temperatureLimit)
	return 0;
    
    return scaling*(1.0 / (1.0 + exp(energyDifference / (kB * temperature)) ) );
}

//------------------------------------------------------------------------------

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
