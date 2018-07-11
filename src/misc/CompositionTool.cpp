/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "CompositionTool.h"

//------------------------------------------------------------------------------

//! The constructor registers initializes the logger.
//! \param logName Optional parameter defining the logger's name.
//! \param simbox Simulation box the tool operates on.
//! is supposed to grow.
CompositionTool::CompositionTool(const SimulationBox &simbox,
		       const char *logName):
    simbox_(&simbox), logName_(logName)
{
    // nothing to be done
}

//------------------------------------------------------------------------------

CompositionTool::~CompositionTool() {}

//------------------------------------------------------------------------------

//! Analyzes a particular layer and returns composition information.
//! \param layerId Layer index of which composition shall be analyzed.
//! \return Composition information about requrested layer.
LayerCompositionInfo CompositionTool::AnalyzeLayer(indexType layerId) const {

    LayerCompositionInfo result;

    for (auto atom: simbox_->getLattice().getAtomsInLayer(layerId,
                                                          simbox_->getOutOfPlaneDimension()) ) {

        auto entry = std::make_tuple(atom->getMaterial()->getName(), atom->getElementId());
        
        if (result.find(entry) == result.end())
            result[entry] = 1;
        else
            result.at(entry) ++;
    }

    return result;
}

//------------------------------------------------------------------------------

//! Analyzes the complete Lattice of a SimulationBox and returns compositoin information.
//! \return Composition information about SimulationBox.
CompositionInfo CompositionTool::AnalyzeStructure (void) const {
    CompositionInfo result;

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
