/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __XML_HANDLER_H__
#define __XML_HANDLER_H__

#include <string>
#include <exception>
#include <sstream>
#include <vector>
#include <memory>
#include <unordered_map>
#include <map>
#include <limits>
#include <tuple>

#include <pugixml/pugixml.hpp>
#include <easyloggingcpp/easylogging++.h>

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
#include <misc/CompositionInfo.h>
#include <misc/Configuration.h>
#include <misc/Journal.h>
#include <misc/FileHandler.h>
#include <misc/Tuples.h>

//! \brief XML file handler

//! Provides methods to read and write XML files. The data are internally
//! represented using a document object model (DOM) structure.
//! This software is based on pugixml library (http://pugixml.org). pugixml is
//! Copyright (C) 2006-2015 Arseny Kapoulkine. 

class XmlHandler : public FileHandler{
 private:

    //! Root of DOM tree structure.
    pugi::xml_document xmlDoc_;

    //! Logger name.
    const std::string logName_;

    //! Read Journal entries as Integer.
    int readEntries(Journal<int> &journal, pugi::xml_node entry) const;
    //! Read Journal entries as Double.
    int readEntries(Journal<double> &journal, pugi::xml_node entry) const;
    //! Read Journal entries as std::string.
    int readEntries(Journal<std::string> &journal, pugi::xml_node entry) const;
    //! Read Journal entries as UDTuple.
    int readEntries(Journal<UDTuple> &journal, pugi::xml_node entry) const;

    //! Write Journal Entries.
    template <class T>
    int writeEntries(const Journal<T> &journal,
				  pugi::xml_node &journalNode);

    //! Write std::string Journal Entries.
    int writeEntries(const Journal<std::string> &journal,
				  pugi::xml_node &journalNode);
    //! Write UDTuple Journal Entries.
    int writeEntries(const Journal<UDTuple> &journal,
                pugi::xml_node &journalNode);
    
 public:
   
    //! Constructor
    XmlHandler(const char *logName="XmlHandler");
    //! Destructor
    ~XmlHandler();

    //! Clear DOM structure.
    void clear(void);
    //! Load XML file into DOM structure.
    void load(const std::string &fileName);
    //! Export DOM structure to XML file.
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
    void get(std::shared_ptr<SimulationBox> &simbox) const;

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


//------------------------------------------------------------------------------

//! \brief XML exception

//! The used XML parser (pugixml) indicates errors during operation by return
//! codes, this is 
//! changed at the interface to exception-based notification.
class XmlException : public std::exception {
    
 public:

    //! Identifies source of the Exception.
    enum class Id{
	FileNotFound,
	OutOfMemory,
	Tag,
	Declaration,
	PcData,
	StartElement,
	Attribute,
	EndElement,
	Mismatch,
	NoElement,
	Save,
	Journal,
	Unknown
    } ;
    
    //! Constructor
    XmlException(const std::string &message="",
		 const std::string &fileName="",
		 const XmlException::Id id=XmlException::Id::Unknown);

    //! Return description of execution.
    const char * what () const throw ();
    //! Return more detailled error message.
    const std::string getMessage() const noexcept;
    //! Return name of file which caused the exception.
    const std::string getFileName () const noexcept;
    //! Return Exception identification number.
    const XmlException::Id getId() const noexcept;
    
 private:

    std::string message_; //!< error message
    std::string fileName_; //!< name of file which caused the exception
    XmlException::Id id_; //! < unique ID identifying the exception source
    
};

#include "XmlHandler.inl"

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
