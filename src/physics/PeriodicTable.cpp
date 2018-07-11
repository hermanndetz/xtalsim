/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "PeriodicTable.h"

//------------------------------------------------------------------------------

//! The constructor registers initializes the logger.
//! \param logName Optional parameter defining the logger's name.
PeriodicTable::PeriodicTable(void): logName_("PeriodicTable")
{
    
    // reserve room for every known element as we store pointers to the object
    // and vector may move its members if its capacity is reached, which would
    // invalidate the stored pointers.
    elements_.reserve(118);
}

//------------------------------------------------------------------------------

PeriodicTable::~PeriodicTable() {}

//------------------------------------------------------------------------------

//! \param element Element to add.
void PeriodicTable::add (const Element &element)
{
    CLOG(TRACE, logName_) << "adding element " << element.name << " ("
			  << element.symbol << ") to table";
    elements_.push_back(element);

}

//------------------------------------------------------------------------------

//! \param id Unique Element ID.
//! \throws PeriodicTableException
const Element & PeriodicTable::getById(const elementType id) const
{

    CLOG(TRACE, logName_) << "Searching for Element ID '" << id << "'";
    
    for (auto const &el: elements_)
	if (el.id == id) return el;

    CLOG(ERROR, logName_) << "Element ID '" << id <<
	"' not found in periodic table";
    throw PeriodicTableException("Element not found by element ID",
				 PeriodicTableException::Id::ElementId);    
}

//------------------------------------------------------------------------------

//! \warning This method is only assured to
//! deliver the right element if element names are unique.
//! \param name Name of the desired element.
//! \throws PeriodicTableException
const Element & PeriodicTable::getByName(const std::string &name) const
{
    //    for (int i=0; i<elements_.size(); i++){
    //	if (elements_[i].name == name) return elements_[i];
    //}

    CLOG(TRACE, logName_) << "Searching for Element by name '" << name << "'";
    
    for (auto const &el: elements_)
	if (el.name == name) return el;

    CLOG(ERROR, logName_) << "Element name '" << name <<
	"' not found in periodic table";
    throw PeriodicTableException("Element not found by name",
				 PeriodicTableException::Id::Name);
    
}

//------------------------------------------------------------------------------

//! \warning This method is only assured to
//! deliver the right element if element symbols are unique.
//! \param symbol Symbol of the desired element.
//! \throws PeriodicTableException
const Element & PeriodicTable::getBySymbol(const std::string &symbol) const
{
    CLOG(TRACE, logName_) << "Searching for Element by symbol '" << symbol
                          << "'";
    
    for (auto const &el: elements_)
	if (el.symbol == symbol) return el;

        CLOG(ERROR, logName_) << "Element symbol '" << symbol <<
	"' not found in periodic table";
    throw PeriodicTableException("Element not found by symbol",
				 PeriodicTableException::Id::Symbol);    
}

//------------------------------------------------------------------------------

//! \warning This method is only assured to
//! deliver the right element if only one element with supplied combination of
//! proton and neutron count exists.
//! \param protonCount Amounts of protons in the core.
//! \param neutronCount Amounts of neutrons in the core.
//! \throws PeriodicTableException
const Element & PeriodicTable::getByProtonNeutron(const uint8_t protonCount,
				 const uint8_t neutronCount) const
{

    CLOG(TRACE, logName_) << "Searching for Element by proton count " <<
        protonCount << " and neutron count " << neutronCount;

    
    for (auto const &el: elements_){
	if ((el.protonCount == protonCount) &&
	    (el.neutronCount == neutronCount))
	    return el;
    }

    CLOG(ERROR, logName_) << "Element not found by proton count " << protonCount
			  << " and neutron count " << neutronCount;
    throw PeriodicTableException("Element not found by proton/neutron count",
				 PeriodicTableException::Id::ProtonNeutron);
}

//------------------------------------------------------------------------------

//! Return a constant reference to all available Elements. Please note that it
//! is not possible to edit the elements in that way.
//! \return Vector of all Elements currently stored.
std::vector<Element> const & PeriodicTable::get(void) const
{
    return elements_;
}

//------------------------------------------------------------------------------

//! Return the number of elements, which are loaded in this Periodic Table.
//! \return Number of all Elements currently stored.
uint16_t PeriodicTable::size() const
{
    return elements_.size();
}

//##############################################################################

PeriodicTableException::PeriodicTableException(const std::string &message,
				      const PeriodicTableException::Id id)
    :message_(message), id_(id)
{
    // nothing to be done
}

//------------------------------------------------------------------------------

//! Returns a short description of the exception. This is the implementation of
//! the virtual function inherited from std::exception
//! \remark This function does not throw exceptions!
const char * PeriodicTableException::what () const throw ()
{
    static std::string text;
    text = "Periodic Table Exception -- " + message_;
    return text.c_str();
}

//------------------------------------------------------------------------------

//! Returns a more detailled description of the exception, making it possible to
//! use a single exception class for multiple errors.
//! \remark This function does not throw exceptions!
const std::string PeriodicTableException::getMessage () const noexcept
{
    return message_;
}

//------------------------------------------------------------------------------

//! Returns the ID of the exception source, making possible to easily
//! distinguish different errors.
//! \remark This function does not throw exceptions!
const PeriodicTableException::Id PeriodicTableException::getId () const noexcept
{
    return id_;
}

//##############################################################################

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
