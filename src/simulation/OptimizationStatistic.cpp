/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#include "OptimizationStatistic.h"

//------------------------------------------------------------------------------

OptimizationStatistic::OptimizationStatistic(const std::string _name):
    name(_name), acceptCount(0), rejectCount(0)
{
    //nothing to be done
}

//------------------------------------------------------------------------------

OptimizationStatistic::~OptimizationStatistic() {}

//------------------------------------------------------------------------------

void OptimizationStatistic::merge (const OptimizationStatistic &statistic)
{
    this->acceptCount += statistic.acceptCount;
    this->rejectCount += statistic.rejectCount;
}

//------------------------------------------------------------------------------

void OptimizationStatistic::clear (void)
{
    this->acceptCount = 0;
    this->rejectCount = 0;
}

//------------------------------------------------------------------------------

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
