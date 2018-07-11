/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __TERSOFF_POTENTIAL_H__
#define __TERSOFF_POTENTIAL_H__

#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <cmath>

#include <easyloggingcpp/easylogging++.h>

#include <physics/Atom.h>
#include <physics/Vector3D.h>
#include <simulation/TersoffParameter.h>

typedef std::unordered_map<uint32_t,TersoffParameter> TersoffParameterMap;

class SimulationBox;

//! \brief Stores Tersoff Parameters.

//! Different Tersoff Parameter objects are stored together in this class. It
//! provides functions to easily find the desired parameters. A major
//! requirement for this class is speed, i.e., it should work as fast as
//! possible. 

class TersoffPotential{
 private:
    //! all Tersoff Parameters
    TersoffParameterMap parameters_;

    //! Logger name.
    const std::string logName_;

    //! Generates unique ID for specified element IDs.
    inline uint32_t getId_(const elementType elementId1,
			   const elementType elementId2) const;
    
 public:

    //! Constructor
    TersoffPotential(const char *logName="TersoffPotential");
    //! Destructor
    ~TersoffPotential();

    //! Add Tersoff Parameter to the Tersoff Potential.
    void add (const TersoffParameter &parameter);
    //! Return all stored Tersoff Parameters.
    TersoffParameterMap const & get(void) const;
    //! Select Tersoff parameter by element IDs.
    TersoffParameter const & get(const elementType elementId1,
				 const elementType elementId2) const;
    
};

//------------------------------------------------------------------------------

//! \brief XML exception

//! The used XML parser (pugixml) indicates errors during operation by return
//! codes, this is 
//! changed at the interface to exception-based notification.
class TersoffPotentialException : public std::exception {
    
 public:

    //! Identifies source of the Exception.
    enum class Id{
    NoPotential,
	Unknown
    } ;
    
    //! Constructor
    TersoffPotentialException(const std::string &message="",
		 const TersoffPotentialException::Id id=TersoffPotentialException::Id::Unknown);

    //! Return description of execution.
    const char * what () const throw ();
    //! Return more detailled error message.
    const std::string getMessage() const noexcept;
    //! Return Exception identification number.
    const TersoffPotentialException::Id getId() const noexcept;
    
 private:

    std::string message_; //!< error message
    TersoffPotentialException::Id id_; //! < unique ID identifying the exception source
    
};

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
