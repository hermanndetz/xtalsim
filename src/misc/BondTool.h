/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __BOND_TOOL_H__
#define __BOND_TOOL_H__

#include <iostream>

#include <misc/BondInfo.h>
#include <physics/SimulationBox.h>


class BondTool {
private:
    //! Lattice the strain operations work on.
    SimulationBox const * const simbox_;

    //! Logger name.
    const std::string logName_;

    //! Analyzes a particular layer and returns composition information.
    LayerBondInfo    AnalyzeLayer(indexType layerId) const;

public:
    //! Constructor
    BondTool (const SimulationBox &simbox, const char *logName="BondTool");

    //! Destructor
    ~BondTool ();

    //! Analyzes the complete Lattice of a SimulationBox and returns compositoin information.
    BondInfo AnalyzeStructure (void) const;
};

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:

