/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __OPTIMIZATION_STATISTIC_H__
#define __OPTIMIZATION_STATISTIC_H__

#include <string>

//! \brief Stores statistical parameters for a single optimization action.

//! This class solely stores the statistics of an action and no other properties
//! such as the function to call or the probability. This was done to increase
//! the speed as these actions are called a lot, meaning that the access to the
//! function to call has to be as fast as possible.

class OptimizationStatistic{
public:

    //! Name of Action.
    std::string name;
    //! Number of runs leading to a better overall energy.
    unsigned long acceptCount;
    //! Number of runs leading to a worse overall energy.
    unsigned long rejectCount;
   
    //! Constructor
    OptimizationStatistic(const std::string _name="undef");

    //! Destructor
    ~OptimizationStatistic();

    //! Add values of other statistic to current one.
    void merge(const OptimizationStatistic &statistic);

    //! Clear current statistic.
    void clear(void);
    
};

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
