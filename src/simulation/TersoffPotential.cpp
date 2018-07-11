/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "TersoffPotential.h"

//------------------------------------------------------------------------------

//! The constructor registers initializes the logger.
//! \param logName Optional parameter defining the logger's name.
TersoffPotential::TersoffPotential(const char *logName): logName_(logName)
{
    // nothing to be done
}

//------------------------------------------------------------------------------

TersoffPotential::~TersoffPotential() {}

//------------------------------------------------------------------------------

//! \param parameter Tersoff Parameter to add.
void TersoffPotential::add (const TersoffParameter &parameter)
{
    parameters_.emplace (getId_(parameter.elementIdLow,parameter.elementIdHigh),
			 parameter);
}

//------------------------------------------------------------------------------

TersoffParameterMap const & TersoffPotential::get(void) const
{
    return parameters_;
}

//------------------------------------------------------------------------------

//! \param elementId1 Unique ID of first element.
//! \param elementId2 Unique ID of second element.
TersoffParameter const & TersoffPotential::get(const elementType elementId1,
				 const elementType elementId2) const
{
    CLOG(TRACE, logName_) << "trying to find potential for elements " <<
	elementId1 << " and " << elementId2;
    try{
        return parameters_.at(getId_(elementId1, elementId2));
    }
    catch(std::exception &e) {
        std::ostringstream text;
        text << "no potential found for ID " << elementId1 <<
            " and " << elementId2;
        CLOG(ERROR, logName_) << text.str();
        
        throw TersoffPotentialException(text.str(),
			   TersoffPotentialException::Id::NoPotential);
    }
}

//------------------------------------------------------------------------------

//! \param elementId1 Unique element ID of first involved element.
//! \param elementId1 Unique element ID of second involved element.
//! \todo improve speed by removing if, doubles memory requirements but should
//! improve speed, not introduced yet because more care has to be taken when
//! exporting parameters then
inline uint32_t TersoffPotential::getId_(const elementType elementId1,
				 const elementType elementId2) const
{

    if (elementId1 < elementId2)
	return ( (elementId1 << 16) + elementId2);
    else
	return ( (elementId2 << 16) + elementId1);

}

//##############################################################################

TersoffPotentialException::TersoffPotentialException(const std::string &message,
			   const TersoffPotentialException::Id id):
    message_(message), id_(id) {}

//------------------------------------------------------------------------------

//! Returns a short description of the exception. This is the implementation of
//! the virtual function inherited from std::exception
//! \remark This function does not throw exceptions!
const char * TersoffPotentialException::what () const throw ()
{
    static std::string text;
    text = "Tersoff Potential Exception: " +message_;
    return text.c_str();
}

//------------------------------------------------------------------------------

//! Returns a more detailled description of the exception, making it possible to
//! use a single exception class for multiple errors.
//! \remark This function does not throw exceptions!
const std::string TersoffPotentialException::getMessage () const noexcept
{
    return message_;
}

//------------------------------------------------------------------------------

//! Returns the ID of the exception source, making possible to easily
//! distinguish different errors.
//! \remark This function does not throw exceptions!
const TersoffPotentialException::Id TersoffPotentialException::getId () const
    noexcept
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
