/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "BondTool.h"

//------------------------------------------------------------------------------

//! The constructor registers initializes the logger.
//! \param logName Optional parameter defining the logger's name.
//! \param simbox Simulation box the tool operates on.
//! is supposed to grow.
BondTool::BondTool(const SimulationBox &simbox,
		       const char *logName):
    simbox_(&simbox), logName_(logName)
{
    // nothing to be done
}

//------------------------------------------------------------------------------

BondTool::~BondTool() {}

//------------------------------------------------------------------------------

//! Analyzes a particular layer and returns bond information.
//! \param layerId Layer index of which bonds shall be analyzed.
//! \return Bond information about requrested layer.
LayerBondInfo BondTool::AnalyzeLayer(indexType layerId) const {

    LayerBondInfo result;

    for (auto atom: simbox_->getLattice().getAtomsInLayer(layerId,
                                                          simbox_->getOutOfPlaneDimension()) ) {

        NeighborList neighbors = atom->getNeighbors();

        for (auto neighborIndex: neighbors) {
            auto neighbor = simbox_->getLattice()(neighborIndex);

            auto entry = std::make_tuple(atom->getElementId(), neighbor->getElementId());
            
            if (result.find(entry) == result.end())
                result[entry] = 1;
            else
                result.at(entry) ++;
        }
    }

    return result;
}

//------------------------------------------------------------------------------

//! Analyzes the complete Lattice of a SimulationBox and returns 
//! bond information.
//! \return Composition information about SimulationBox.
BondInfo BondTool::AnalyzeStructure (void) const {
    BondInfo result;

    Vector3D<indexType> latticeSize = simbox_->getLattice().getSize();

    for (indexType i = 0; i < latticeSize[simbox_->getOutOfPlaneDimension()]; i++) {
        result << AnalyzeLayer(i);
    }

    return result;
}

//------------------------------------------------------------------------------

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:

