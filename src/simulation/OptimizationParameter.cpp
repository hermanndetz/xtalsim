/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#include "OptimizationParameter.h"

//------------------------------------------------------------------------------

OptimizationParameter::OptimizationParameter():
    mmcMaxDisplacement(0), maxScaling(0), runCount(0), mmcRunCount(1),
    mmcMinDisplacement(0), minScaling(0), checkCount(1),
    energyDropFactor(2), reductionFactor(2),
    minEnergy(0), anionPassivationProbability(1.0),	maxThreadCount(1)
{
    //nothing to be done
}

//------------------------------------------------------------------------------

OptimizationParameter::~OptimizationParameter() {}

//------------------------------------------------------------------------------

//! \param config Configuration holding desired values.
void OptimizationParameter::set(const Configuration &config)
{
   
    this->mmcMinDisplacement = config.minDisplacement;
    this->mmcMaxDisplacement = config.maxDisplacement;
    this->mmcRunCount = config.mmcRunCount;
    
    this->minScaling = config.minScaling;
    this->maxScaling = config.maxScaling;

    this->runCount = config.runCount;
    this->checkCount = config.checkCount;

    this->reductionFactor = config.reductionFactor;
    this->energyDropFactor = config.energyDropFactor;
    this->minEnergy = config.minEnergy;
    this->anionPassivationProbability = config.anionPassivationProbability;

    this->maxThreadCount = config.maxThreadCount;
    
}

//------------------------------------------------------------------------------

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
