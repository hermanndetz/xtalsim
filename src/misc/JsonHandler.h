/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __JSON_HANDLER_H__
#define __JSON_HANDLER_H__

#include <string>
#include <iostream>
#include <fstream>
#include <tuple>

#include <jsoncpp/json.h>
#include <easyloggingcpp/easylogging++.h>

#include <physics/SimulationBox.h>
#include <physics/Lattice.h>
#include <physics/MaterialCollection.h>
#include <physics/Material.h>
#include <physics/PeriodicTable.h>
#include <physics/Element.h>
#include <simulation/TersoffPotential.h>
#include <simulation/TersoffParameter.h>
#include <misc/Configuration.h>
#include <misc/Journal.h>
#include <misc/FileHandler.h>

//! \brief JSON file handler

//! Provides methods to read and write JSON files. The data are internally
//! represented using a document object model (DOM) structure.
//! This software is based on JsonCpp
//! (https://github.com/open-source-parsers/jsoncpp)

class JsonHandler : public FileHandler {

private:

    //! Logger name.
    const std::string logName_;
    
    //! Root of DOM structure.
    Json::Value jsonDoc_;

    //! Read Journal entries as Integer.
    int readEntries(Journal<int> &journal, Json::Value entries) const;
    //! Read Journal entries as Double.
    int readEntries(Journal<double> &journal, Json::Value entries) const;
    //! Read Journal entries as std::string.
    int readEntries(Journal<std::string> &journal, Json::Value entries) const;

    //! Write Journal Entries.
    template <class T>
    int writeEntries(const Journal<T> &journal, Json::Value &entryNode);
    
public:

    //! Constructor
    JsonHandler(const char *logName="JsonHandler");
    //! Destructor
    ~JsonHandler();

    //! Clear DOM structure.
    void clear(void);
    //! Load JSON file into DOM structure.
    void load(const std::string &fileName);
    //! Export DOM structure to JSON file.
    void save(const std::string &fileName) const;

    //! Write Elements to DOM structure.
    void set(const PeriodicTable &periodicTable);
    //! Fill PeriodicTable with Elements read from XML file.
    void get(PeriodicTable &periodicTable) const;

    //! Write Tersoff Potential to DOM structure.
    void set(const TersoffPotential &potential);
    //! Fill Tersoff Potential from DOM structure.
    void get(TersoffPotential &potential) const;

    //! Write Material Collection to DOM structure.
    void set(const MaterialCollection &collection);
    //! Fill Material Collectoin from DOM structure.
    void get(MaterialCollection &collection) const;

    //! Write lattice to DOM structure.
    void set(const Lattice &lattice);
    //! Fill lattice structure with atoms read from XML file.
    void get(Lattice &lattice, const MaterialCollection &materials) const;

    //! Write Simulation Box to DOM structure.
    void set(const SimulationBox &simbox);
    //! Fill lattice structure with atoms read from XML file.
    void get(std::shared_ptr<SimulationBox> & simbox) const;

    //! Write CompositionInfo to DOM structure.
    void set(const CompositionInfo &composition);
    //! Fill vector of CompositionInfos from XML file.
    void get(std::vector<CompositionInfo> &compositions) const;
    
    //! Write Journal to DOM structure.
    template <class T>
    void set(const Journal<T> &journal);
    //! Fill Journal from DOM structure.
    template <class T>
    void get(Journal<T> &journal, const std::string name) const;
    
    //! Fill Configuration from DOM structure.
    void get(Configuration &conf) const;
    
};

#include "JsonHandler.inl"

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
