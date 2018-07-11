/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __STRAIN_TOOL_H__
#define __STRAIN_TOOL_H__

#include <vector>
#include <tuple>
#include <sstream>

#include <physics/Lattice.h>
#include <physics/Atom.h>
#include <physics/Field3D.h>
#include <physics/Vector3D.h>
#include <physics/MaterialCollection.h>
#include <physics/Material.h>
#include <misc/LayerStrainInfo.h>

//! Vector containing strain infos for each single bond.
typedef std::vector<std::tuple<Vector3D<indexType>,
			       Vector3D<indexType>, double>> StrainInfo;

//! Defines methods for strain calculation
typedef enum {
    SCM_BondStrain,
    SCM_GrowthStrain,
    SCM_DistanceStrain
} StrainCalculationModes;

//! This include must preceed typedefs due to circular dependencies.
#include <physics/SimulationBox.h>

//! Number of digits after the decimal point in strain output files.
const uint8_t strainOutputPrecision=10;

//! \brief Implements several ways to calculate strain on SimulationBox.

//! This class was created to separate the strain calculation from the code in
//! the SimulationBox. Several ways to calculate the strain are implemented. The
//! results can be written to files or returned by the functions.

class StrainTool{
private:

    //! Lattice the strain operations work on.
    SimulationBox * simbox_;

    //! Logger name.
    const std::string logName_;

    //! Temporary storage of strain data
    std::shared_ptr<Field3D<double>> strainField_{};

    //! true, if strainField_ contains valid data
    bool strainFieldReady_;
    
    //! Calculate strain of a single Layer
    void calculateStrainLayer(const MaterialCollection &materials,
			      const std::string &outputFileName,
			      const int layer,
			      const double refDistance,			      
			      double &bondStrain,
			      double &growthStrain,
			      double &distanceStrain,
                  bool writeToFile=true) const;
    
public:

    //! Constructor
    StrainTool(SimulationBox &simbox,
	   const char *logName="StrainTool");

    //! Constructor
    StrainTool(std::shared_ptr<SimulationBox> simbox,
	   const char *logName="StrainTool");

    //! Desctructor
    ~StrainTool();

    //! Calculate strain and write it to file.
    LayerStrainInfo calculateStrain(const MaterialCollection &materials,
			 const std::string &outputFileName,
			 double refDistance=-1.0, bool writeToFile=true) const;

    //! Calculate strain for single atom
    double calculateStrainAtom(const MaterialCollection &materialsInSimbox, 
            const MaterialCollection &otherMaterials,
            const Atom *atom, StrainCalculationModes mode,
            double refDistance=-1.0) const;

    //! Calculate strain for each bond and return info.
    StrainInfo getStrainInfo (const MaterialCollection &materials,
			      double refDistanceParam) const;

    //! Calculate strain state of each atom and return as scalar field.
    std::shared_ptr<Field3D<double>> getStrainField (const MaterialCollection &materials,
                double refDistance);
};

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
