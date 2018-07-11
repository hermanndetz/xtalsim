/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <vector>
#include <tuple>
#include <string>
#include <sstream>

#include <physics/AtomState.h>
#include <physics/Vector3D.h>

//! \brief Store content of a configuration file.

class Configuration{
public:

    //! Cations of desired material.
    std::vector< std::tuple<std::string,double> > cations;
    //! Anions of desired material.
    std::vector< std::tuple<std::string,double> > anions;
    //! Material name.
    std::string materialName;
    //! File containing material definition.
    std::string materialFile;
    //! Material that shall be searched in file.
    std::string materialSearchName;
    //! Periodic Table file name.
    std::string periodicTableFileName;
    //! XYZ file name.
    std::string xyzFileName;
    //! elastic constant.
    double c11;
    //! elastic constant.
    double c12;
    //! elastic constant.
    double c44;
    //! Lattice constant.
    double latticeConstant;
    //! Type of lattice.
    std::string latticeType;
    //! Size of add structure.
    Vector3D<indexType> size;
    //! Dimension the crystal is allowed to grow
    uint8_t growthDimension;
    //! The radius used to generate neighbor list.
    double neighborRadius;
    //! The atomic layer count used to generate neighbor list.
    indexType neighborLayers;
    //! Lattice temperature in Kelvin.
    double latticeTemperature;
    
    //! Start index.
    int startIndex;
    //! Stop index.
    int stopIndex;    
    //! Probability to trigger MMC relaxation.
    double mmcProbability;
    //! How often mmc relaxation shall be carried out per function call.
    int mmcRunCount;
    //! Minimal displacement in MMC relaxation.
    double minDisplacement;
    //! Maximal displacement for MMC relaxation.
    double maxDisplacement;
    //! Probability to trigger Scaling of atomic layers.
    double scalingProbability;
    //! Maximal scaling swing around 1 for atomic layer scaling.
    double minScaling;
    //! Maximal swing around 1 of scaling parameter.
    double maxScaling;
    //! How many optimization steps shall be carried out.
    int runCount;
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
    //! If false, passivating anions can also be introduced if no matching cation is present.
    bool requireMatchingCation;
    //! Maximal number of concurrent threads.
    uint16_t maxThreadCount;
    
    //! Number of cations introduced when creating interface.
    int interfaceAtoms;
    //! Number of positions tested for each cation.
    int interfacePositions;
  
    //! Input file name.
    std::string inputFileName;
    //! Output file name.
    std::string outputFileName;
    //! Log file name.
    std::string logFileName;
    //! Name of file storing all Tersoff parameters.
    std::string tersoffFileName;
    //! Preamble of output files.
    std::string outputPreamble;
    //! Preamble of journal files.
    std::string journalPreamble;
    
    //! Extend lattice.
    bool extend;
    //! Calculate energy.
    bool calculateEnergy;
    //! Print additional information.
    bool print;
    //! Print debug messages.
    bool verbose;
    //! Disable screen output.
    bool quiet;
    
    //! Render bonds with nearest neighbor atoms
    bool renderBonds;
    //! Visualize strain of individual bonds.
    bool visualizeBondStrain;
    //! Reference lattice constant for strain calculationn
    double refLatticeConstant;
    //! Render only atoms and bonds, which were modified
    bool modifiedOnly;
    //! Filter atoms by modification state
    std::vector<std::string> filterModificationState;
    //! Render layers containing modified atoms, but filter modified atoms
    bool modifiedNegative;
    //! Render only atoms up to a certain modification index
    uint32_t modificationIndexMax;
    //! Initial focal point of camera in Visualization.
    Vector3D<double> cameraFocalPoint;
    //! Initial position of camera in Visualization.
    Vector3D<double> cameraDirection;
    //! Initial zoom in Visualization.
    double cameraZoom;
    //! If true animation of modifications only shall be executed.
    bool modificationAnimation;
    //! Defines how many succeeding modifications shall be skipped when one was displayed.
    uint32_t animationStep;
    //! Prefix for animation files.
    std::string animationFilePrefix;
    //! Filename where screenshot is stored.
    std::string screenshotFileName;
    //! Number of rotational pictures made.
    uint16_t rotationCount;
    
    //! Use dynamic optimization.
    bool dynamicOptimization;
    //! Use static optimization.
    bool staticOptimization;
    
    //! Constructor
    Configuration();

    //! Copy constructor
    Configuration(const Configuration &config);

    //! Destructor
    ~Configuration();

    //! Assignment
    void operator=(const Configuration &config);
    
    //! Returns filter definition for atom states
    std::bitset<AtomState::StateCount> getModificationState(void) const;

    //! Return liste of all parameters with description as string.
    const std::string str(void) const;
    
};

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
