/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __INTERFACE_TOOL_H__
#define __INTERFACE_TOOL_H__

#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <limits>
#include <tuple>
#include <algorithm>
#include <random>
#include <unordered_set>

#include <easyloggingcpp/easylogging++.h>

#include <misc/BondInfo.h>
#include <misc/Journal.h>
#include <misc/XmlHandler.h>
#include <misc/JsonHandler.h>
#include <misc/Tuples.h>
#include <misc/Hash.h>
#include <physics/Range3D.h>
#include <physics/Vector3D.h>
#include <physics/Lattice.h>
#include <physics/Atom.h>
#include <physics/Material.h>
#include <physics/SimulationBox.h>
#include <physics/ExchangeReaction.h>
#include <simulation/Optimization.h>
#include <simulation/OptimizationParameter.h>
#include <simulation/TersoffPotential.h>


//******************************************************************************

const std::vector<std::tuple<int,int>> cationPositionsInPlane = {
    std::make_tuple(2,2),
    std::make_tuple(-2,2),
    std::make_tuple(2,2),
    std::make_tuple(2,-2),
    std::make_tuple(-2,-2),
    std::make_tuple(2,-2),
    std::make_tuple(-2,-2),
    std::make_tuple(-2,2)};

const std::vector<std::tuple<int,int>> cationPositionsOutOfPlane = {
    std::make_tuple(0,2),
    std::make_tuple(2,0),
    std::make_tuple(0,-2),
    std::make_tuple(-2,0)};

class ExchangeReaction;

//! Defines possible metrics to determine optimum position of atoms at interface.
typedef enum {
    ITM_Energy = 0,
    ITM_BondStrain,
    ITM_GrowthStrain,
    ITM_DistanceStrain
} InterfaceToolMetric;

//! \brief Handles interface related tasks in simulation box.

//! Handles all interface related tasks like creation in different shapes and
//! forms. These method was chosen to reduce the code size of the SimulationBox
//! class.

class InterfaceTool{
private:

    //! The simulation box whose interfaces are handled.
    SimulationBox *simbox_;

    //! Tersoff potentials used to optimize interfaces.
    const TersoffPotential *tersoff_;

    //! Exchange reaction for anion exchange
    ExchangeReaction *anionExchange_;

    //! Journal logging all tried positions.
    Journal<std::string> posJournal_;
    //! Journal logging all tried positions.
    Journal<std::string> choiceJournal_;
    //! Journal logging total energy or strain of interface layers.
    Journal<UDTuple> metricJournal_;
    
    //! random number generator
    //std::default_random_engine random_;
    std::mt19937 random_;

    //! Stores all indices already tried for a single atom.
    std::unordered_set<Vector3D<indexType>, Hash> checkedIndices_;

    //! Logger name.
    const char *logName_;

    //! Defines after how many tries search for a suitable position is aborted.
    const int unsuitedPositionsLimit_=20;
    //! Counts how many unsuitable postions have been encountered.
    int unsuitedPositionsCount_;

    //! stores index to assign modification order to atoms
    uint32_t modificationIndex_{0};
    
    //! Find the next suitable location for atom.
    Vector3D<indexType> findSuitablePosition(const Vector3D<indexType> position,
				   const Range3D<indexType> interfaceRange,
				   const int direction,
				   const Material *material);

    //! Move atom to next position along a given direction.
    Vector3D<indexType> moveAtom(const Vector3D<indexType> position,
				 const Range3D<indexType> interfaceRange,
				 const int direction,
				 const Material *material);
    
    //! Add anions if possible
    void addAnions(const Vector3D<indexType> position,const Material *material,
                   double passivationProbability, bool requireMatchingCation=true);
    
public:

    //! Constructor
    InterfaceTool(SimulationBox &simbox, const TersoffPotential &tersoff,
		  const char *logName="InterfaceTool");

    //! Constructor
    InterfaceTool(SimulationBox &simbox, const TersoffPotential &tersoff,
          ExchangeReaction &anionExchange, const char *logName="InterfaceTool");

    //! Destructor
    ~InterfaceTool();
    
    //! Search for an energetical optimized interface roughness.
    void createOptimizedRoughness(Range3D<indexType> interfaceRange,
				  const unsigned int atomsToDeployCount,
				  const unsigned int positionsToCheckCount,
				  const Material &material,
				  const OptimizationParameter &parameter,
                  const MaterialCollection &collection,
                  const InterfaceToolMetric metric=ITM_Energy,
				  const std::string &journalPreamble="");
};

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
