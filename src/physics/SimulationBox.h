/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __SIMULATION_BOX_H__
#define __SIMULATION_BOX_H__

#include <vector>
#include <cmath>
#include <memory>
#include <algorithm>
#include <random>
#include <chrono>
#include <iostream>
#include <fstream>
#include <tuple>
#include <thread>

#include <easyloggingcpp/easylogging++.h>

#include <physics/Atom.h>
#include <physics/PeriodicTable.h>
#include <physics/MaterialCollection.h>
#include <physics/Material.h>
#include <physics/Vector3D.h>
#include <physics/Range3D.h>
#include <physics/Lattice.h>
#include <physics/LatticeUtils.h>
#include <simulation/TersoffPotential.h>
#include <simulation/TersoffParameter.h>
#include <simulation/OptimizationParameter.h>
#include <simulation/OptimizationStatistic.h>
#include <misc/CsvHandler.h>
#include <misc/InterfaceTool.h>
#include <misc/LayerStrainInfo.h>
#include <misc/StrainTool.h>
#include <misc/ArrheniusFactor.h>
#include <misc/BondTool.h>
#include <misc/CompositionTool.h>

//! \brief Representation of the simulation box.

//! Contains the lattice as well as additional information. At the moment it is
//! only possible to grow the simulation box in one direction and with one in
//! plane lattice constant. The latter is fixed when the atomic layers are
//! created in the box.

//! \todo overall: improveme by using functions from library "algorithm", e.g.,
//! for_each, ...

class SimulationBox{
private:

    //! Description of the simulation Box.
    std::string description_;
    //! atomic layer representation of atoms
    Lattice lattice_;

    //! dimension the lattice is allowed to grow
    uint8_t outOfPlaneDimension_;
    //! in plane lattice constant
    double inPlaneLatticeConstant_;
    //! size in real space
    Vector3D<spaceType> size_;
    //! Materials used in the simulation box.
    MaterialCollection materials_;
    
    //! random number generator
    //std::default_random_engine random_;
    std::mt19937 random_;

    //! If true an update of the neighbor list is necessary.
    bool updateNeighborList_;

    //! layer count used to generate neighbor list
    indexType neighborListLayers_;
    //! real space radius used to generate neighbor list
    double neighborListRadius_;
    //! periodic boundary condition used to generate neighbor list
    bool neighborListPeriodicBoundaries_;
    //! if true neighbor list is created in parallel in separate thread
    //    bool neighborListThreadRunning_;
    
    //! Logger name.
    const std::string logName_;
    
    //! Get surrounding range of lattice position.
    Range3D<indexType> getSurrounding(const Vector3D<indexType>&positionLattice) const;

    //! Pick random interval in specified range in specified dimension.
    void getRandomInterval(Range3D<indexType> &range,const indexType dimension);

    void createLattice(const Vector3D<indexType> layerCounts,
                       const Material &material,
                       const Lattice::Type latticeType);

    //! Generate neighbor list periodically.
    void generateNeighborsConcurrent(void);
    //! Update neighborlists in each atom.
    void updateNeighbors(void);

    //! Returns energy for given range.
    void getEnergyByReference(const TersoffPotential &tersoff, double &energy,
                   Range3D<indexType> range) const;
    
public:
    
    //! Constructor
    SimulationBox(const uint8_t outOfPlaneDimension,
                  const Vector3D<spaceType> size={0,0,0},
                  const double inPlaneLatticeConstant=-1,
                  const double latticeTemperature=0,
                  const char *logName="SimulationBox");

    //! Destructor
    ~SimulationBox();

    //! Returns description of Simulation Box.
    std::string getDescription(void) const;
    //! Sets description of Simulation Box.
    void setDescription(const std::string description);
    
    //! Returns lattice temperature.
    double getLatticeTemperature(void) const;
    //! Sets lattice temperature.
    void setLatticeTemperature(const double temperature);
    
    //! Return out of plane dimension.
    uint8_t getOutOfPlaneDimension(void) const;
    //! Return in plane lattice constant.
    double getInPlaneLatticeConstant(void) const;
    //! Return size of simulation box in real space.
    Vector3D<spaceType> getSize(void) const;
    //! Return atomic lattice.
    Lattice & getLattice(void);
    //! Return atomic lattice.
    Lattice const & getLattice(void) const;
    //! Return stored materials.
    MaterialCollection & getMaterials(void);
    //! Return stored materials.
    MaterialCollection const & getMaterials(void) const;
    //! Adds a material to the collection stored in simbox
    void addMaterial(const Material &material);

    //! Returns the number of valid atoms in a range.
    uint32_t getAtomCount(const Range3D<indexType>&range=Range3D<indexType>())
	const;
    
    //! Resets modification attribute of all atoms in range.
    void clearModification(const Range3D<indexType>&range=Range3D<indexType>());
    
    //! Introduce Zincblende Lattice structure.
    void createZincblende(const Vector3D<indexType> layerCounts,
                          const Material &material);
    //! Introduce Zincblende Lattice structure.
    void createHalfHeusler(const Vector3D<indexType> layerCounts,
			  const Material &material);
    //! Introduce Zincblende Lattice structure.
    void createFullHeusler(const Vector3D<indexType> layerCounts,
			  const Material &material);

    //! Start generation of neighbor list in a separate thread.
    void startConcurrentNeighborListGeneration(void);
    //! Create neighbor lists in atoms using atoms inside radius.
    void generateNeighbors(indexType layers=1, const double radius=0,
			   bool periodicBoundaries=true);
    //! Create neighbor lists.
    void generateNeighbors(bool force=false);
    
    double getEnergyParallel(const TersoffPotential &tersoff,
                           const Range3D<indexType> range = Range3D<indexType>(),
                             const int maxThreadCount=1) const;
    //! Calculate energy of a specified range.
    double getEnergy(const TersoffPotential &tersoff,
                     Range3D<indexType> range = Range3D<indexType>()) const;
    //! Calculate energy of a single atom.
    double getEnergySingleAtom(const Atom *atom,
			      const TersoffPotential &tersoff) const;

    //! Calculate strain and write it to file.
    LayerStrainInfo calculateStrain(const MaterialCollection &materials,
			 const std::string &outputFileName,
			 double refDistance=-1, bool writeToFile=true);

    //! Calculate strain of a single atom.
    double calculateStrainSingleAtom(const Atom *atom, StrainCalculationModes mode, const MaterialCollection &collection);

    //! Calculate strain for each bond
    StrainInfo getStrainInfo (const MaterialCollection &materials,
				    const double refDistanceParam);

    //! Obtain composition info and write it to file.
    CompositionInfo analyzeComposition(const std::string &outputFileName) const;

    //! Generate bond statistics and write them to file.
    BondInfo analyzeBonds(const std::string &outputFileName) const;

    //! Create energetic optimized interface roughness.
    void createOptimizedRoughness(Range3D<indexType> interfaceRange,
				  const unsigned int atomsToDeployCount,
				  const unsigned int positionsToCheckCount,
				  const Material &material,
				  const OptimizationParameter &parameter,
				  const TersoffPotential &tersoff,
                  const MaterialCollection &collection,
				  const std::string &journalPreamble="" );

    //! Relax atoms using Metropolis Monte Carlo approach in separate threads.
    void mmcRelaxParallel(const TersoffPotential &tersoff,
		  const OptimizationParameter &parameter,
		  OptimizationStatistic &statistic,
                  Range3D<indexType> range = Range3D<indexType>());
    //! Relax atoms using Metropolis Monte Carlo approach.
    void mmcRelax(const TersoffPotential &tersoff,
		  const OptimizationParameter &parameter,
		  OptimizationStatistic &statistic,
		  Range3D<indexType> range = Range3D<indexType>());
    //! Scale simulation box.
    void scale(const TersoffPotential &tersoff,
	       const OptimizationParameter &parameter,
	       OptimizationStatistic &statistic,
	       Range3D<indexType> range = Range3D<indexType>());

    //! Write atoms to .xyz File.
    void writeToXYZ(const std::string fileName,
		    const Range3D<indexType> &range = Range3D<indexType>(),
		    const bool checkModification = false,
		    const std::string materialName = "")
	const;
};

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
