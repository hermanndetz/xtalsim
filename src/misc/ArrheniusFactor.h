/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __ARRHENIUS_FACTOR_H__
#define __ARRHENIUS_FACTOR_H__

#include <math.h>

//! \brief Provides functions to compute temperature based Arrhenius factors.

//! This namespace defines different methods to calculate an arrhenius factor
//! which is used to accept a certain degree of energetic not optimal
//! configurations. The value is computed based on 
//! the current lattice temperature. The higher the latter the higher the
//! probability that non optimal results are accepted.

namespace ArrheniusFactor{
    //! Boltzmann constant
    const double kB = 8.6173324e-5; // eV/K
    //! \brief Temperature below this limit is assumed 0.
    
    //! Is required to prevent arbitrary small values and division by zero
    //! exceptions.
    const double temperatureLimit = 1e-6; // K
    //! Scaling factor for results.
    const double scaling=0.25;

    //! Determine value according to Boltzmann statistic.
    double Boltzmann(const double energyDifference, const double temperature);
    //! Determine value according to Bose-Einstein statistic.
    double BoseEinstein(const double energyDifference,
			const double temperature);
    //! Determine value according to Fermi-Dirac statistic.
    double FermiDirac(const double energyDifference, const double temperature);
}

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
