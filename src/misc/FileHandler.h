/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __FILE_HANDLER_H__
#define __FILE_HANDLER_H__

#include <string>
#include <exception>
#include <sstream>
#include <vector>
#include <memory>
#include <unordered_map>
#include <map>
#include <limits>

#include <physics/SimulationBox.h>
#include <physics/Lattice.h>
#include <physics/MaterialCollection.h>
#include <physics/Material.h>
#include <physics/PeriodicTable.h>
#include <physics/Element.h>
#include <physics/Atom.h>
#include <physics/Vector3D.h>
#include <simulation/TersoffPotential.h>
#include <simulation/TersoffParameter.h>
#include <misc/Configuration.h>
#include <misc/Journal.h>

//! \brief Defines interfaces for different file handlers.

class FileHandler{
    
 public:

    //! Clear DOM structure.
    virtual void clear(void) =0;
    //! Load XML file into DOM structure.
    virtual void load(const std::string &fileName) =0;
    //! Export DOM structure to XML file.
    virtual void save(const std::string &fileName) const =0;

    //! Write Elements to DOM structure.
    virtual void set(const PeriodicTable &periodicTable) =0;
    //! Fill PeriodicTable with Elements read from XML file.
    virtual void get(PeriodicTable &periodicTable) const =0;

    //! Write Tersoff Potential to DOM structure.
    virtual void set(const TersoffPotential &potential) =0;
    //! Fill Tersoff Potential from DOM structure.
    virtual void get(TersoffPotential &potential) const =0;

    //! Write Material Collection to DOM structure.
    virtual void set(const MaterialCollection &collection) =0;
    //! Fill Material Collectoin from DOM structure.
    virtual void get(MaterialCollection &collection) const =0;
    
    //! Write lattice to DOM structure.
    virtual void set(const Lattice &lattice) =0;
    //! Fill lattice structure with atoms read from XML file.
    virtual void get(Lattice &lattice,
	     const MaterialCollection &materials) const =0;

    //! Write Simulation Box to DOM structure.
    virtual void set(const SimulationBox &simbox) =0;
    //! Fill Simulationg Box with atoms read from XML file.
    virtual void get(std::shared_ptr<SimulationBox> &simbox) const = 0;
    
    //! Fill Configuration from DOM structure.
    virtual void get(Configuration &conf) const =0;
    
};

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
