/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __PERIODIC_TABLE_H__
#define __PERIODIC_TABLE_H__

#include <vector>
#include <exception>

#include <easyloggingcpp/easylogging++.h>

#include <physics/Element.h>

//! \brief Stores available elements.

//! Elements with additional information are stored in this class. They can be
//! read from an XML file but also can be later extended. Currently a method to
//! delete elements is not planned.
//! This class was implemented as singleton as only exactly one Periodic Table
//! is valid and also realeases us from providing it to each function as
//! parameter. Furthermore this class is supposed to be filled once at the
//! startup of the application.

class PeriodicTable{
private:
    //! Constructor
    PeriodicTable();

    //! all elements
    std::vector<Element> elements_;

    //! Logger name.
    const std::string logName_;
    
public:

    //! Destructor
    ~PeriodicTable();

    //! Returns unique object of PeriodicTable.
    static PeriodicTable& getInstance(void)
    {
        static PeriodicTable pt;
        return pt;
    }
    
    //! Delete Copy Constructor.
    PeriodicTable(PeriodicTable const&)  = delete;
    //! Delete assignment operator.
    void operator=(PeriodicTable const&) = delete;

    //! Add Element to the periodic table.
    void add (const Element &element);

    //! Get all Elements.
    std::vector<Element> const & get(void) const;
    //! Select element by Element ID.
    const Element & getById(const elementType id) const;
    //! Select element by Element Name.
    const Element & getByName(const std::string &name) const;
    //! Select element by Element short name.
    const Element & getBySymbol(const std::string &symbol) const;
    //! Select element by proton and neutron count
    const Element & getByProtonNeutron(const uint8_t protonCount,
					  const uint8_t neutronCount) const;
    
    //! Number of elements.
    uint16_t size() const;
};


//------------------------------------------------------------------------------

//! \brief Periodic Table exception

//! Exceptions thrown by PeriodicTable class
class PeriodicTableException : public std::exception {
   
 public:
   
    //! Identifies source of the Exception.
    enum class Id{
    	ElementId,
    	Name,
    	Symbol,
    	ProtonNeutron,
    	Unknown
    } ;
    
    //! Constructor
    PeriodicTableException(const std::string &message="",
			   const PeriodicTableException::Id id =
			   PeriodicTableException::Id::Unknown ); 

    //! Return description of execution.
    const char * what () const throw ();
    //! Return more detailled error message.
    const std::string getMessage() const noexcept;
    //! Return Exception identification number.
    const PeriodicTableException::Id getId() const noexcept;
    
 private:

    std::string message_; //!< error message
    PeriodicTableException::Id id_; //!< unique ID identifying source
    
};


#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
