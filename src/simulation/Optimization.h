/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __OPTIMIZATION_H__
#define __OPTIMIZATION_H__

#include <vector>
#include <memory>
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>
#include <string>
#include <iostream>
#include <fstream>

#include <easyloggingcpp/easylogging++.h>

#include <misc/XmlHandler.h>
#include <misc/Journal.h>
#include <physics/Range3D.h>
#include <simulation/TersoffPotential.h>
#include <simulation/OptimizationParameter.h>
#include <simulation/OptimizationStatistic.h>

class SimulationBox;

typedef void (SimulationBox::*OptimizationFunction) (
			       const TersoffPotential&,
			       const OptimizationParameter &parameter,
			       OptimizationStatistic &statistic,
			       Range3D<indexType>);

//typedef std::function<void(Range3D<indexType>,
//			   const TersoffPotential&,
//			   const double)> OptimizationFunction;


//! \brief Carries out optimizations of simulation box.

//! It is possible to register single actions together with the probability it
//! shall be called. The actions are stored as function pointers with a single
//! paramter, a Range3D object, definine the range on which the function shall
//! be applied.

class Optimization{

private:

    //! Available actions.
    //std::vector<void (SimulationBox::*) (Range3D<indexType>&)> actions_;
    std::vector<OptimizationFunction> actions_;
    //! Statistics for each single action.
    std::vector<OptimizationStatistic> statistics_;
    //! Probabilites with which corresponding actions shall be triggered.
    std::vector<double> probabilities_;
    //! Sum of probabilites.
    std::vector<double> probabilitiesSum_;
    
    //! Logger name.
    const std::string logName_;
    
    //! Random number generator.
    std::mt19937 random_;

    //! Update probability vectors.
    void checkProbabilities(void);

    //! Returns index of random action.
    inline uint16_t getRandomAction (void);
    
public:

    //! Constructor
    Optimization(const char *logName="Optimization");

    //! Desctructor
    ~Optimization();
    
    //! Adds an optimization action.
//    void registerAction(void (SimulationBox::*fp) (Range3D<indexType>&),
    void registerAction(const std::string name,
			const OptimizationFunction fp,const double probability);

    //! Start the static optimization.
    void runStatic(SimulationBox &simbox,
	     const Range3D<indexType> &range,
	     const TersoffPotential &tersoff,
	     const OptimizationParameter &parameter,
	     const std::string journalPreamble = "");

    //! Start the dynamic optimization.
    void runDynamic(SimulationBox &simbox,
		      const Range3D<indexType> &range,
		      const TersoffPotential &tersoff,
		      OptimizationParameter parameter,
		      const std::string journalPreamble = "");
    
    //! Print statistic of Optimization actions so far.
    void printStatistic(void) const;

    //! Resets statistic to its initial state.
    void clearStatistic(void);
};

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
