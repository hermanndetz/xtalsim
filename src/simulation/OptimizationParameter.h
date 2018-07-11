/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __OPTIMIZATION_PARAMETER_H__
#define __OPTIMIZATION_PARAMETER_H__

#include <limits>

#include <misc/Configuration.h>

//! \brief Stores parameter for the optimization steps.

class OptimizationParameter{
public:

    //! Maximal displacement in MMC relaxation.
    double mmcMaxDisplacement;
    //! Maximal scaling swing around 1 for atomic layer scaling.
    double maxScaling;
    //! Overall optimization steps.
    int runCount;
    //! MMC runs that shall be carried out with one call of function.
    int mmcRunCount;

    // parameters for dynamic optimization
    //! Minimal displacement in MMC relaxation.
    double mmcMinDisplacement;
    //! Maximal scaling swing around 1 for atomic layer scaling.
    double minScaling;
    //! Optimization steps carried out before energy gain is checked.
    int checkCount;
    //! Factor the energy has to drop before parameters are reduced.
    double energyDropFactor;
    //! Factor used to reduce the parameters.
    double reductionFactor;
    //! Minimal energy gain that has to be achieved in first checkCount runs.
    double minEnergy;
   //! Probability with which an anion is added for suitable bond partners.
    double anionPassivationProbability; 

    //! Maximum number of concurrent threads.
    uint16_t maxThreadCount;
    
    //! Constructor
    OptimizationParameter();

    //! Destructor
    ~OptimizationParameter();

    //! Fill variables of class from configuration object.
    void set(const Configuration &config);
    
};

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
