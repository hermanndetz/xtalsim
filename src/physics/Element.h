/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __ELEMENT_H__
#define __ELEMENT_H__

#include <stdint.h>
#include <string>

#include <misc/Color.h>

typedef uint16_t elementType;

//! \brief Stores information regarding a single element.

//! Specific information for a single element are stored here. Please note that
//! the same atomic number can appear due the usage of isotopes. Therefore a
//! unique identification number is generated using both the proton and neutron
//! count of the element.

class Element{
    
 public:

    elementType id; //!< unique element ID
    uint8_t protonCount; //!< number of protons in core
    uint8_t neutronCount; //!< number of neutrons in core

    std::string name; //!< name of the element
    std::string symbol; //!< chemical symbol
    
    uint8_t period; //!< period in periodic table
    uint8_t group; //!< group in periodic table

    double mass; //!< mass of element
    double weight; //!< weight of element

    Color color; //!< color of element in visualization
    
    //! Constructor
    Element(std::string _name="undef", std::string _symbol="undef",
	    uint8_t _protonCount=0, uint8_t _neutronCount=0,
	    uint8_t _period=0, uint8_t _group=0,
	    double _mass=0.0, double _weight=0.0, std::string _color="#000000");

    //! Destructor
    ~Element();

    //! Generate unique Element Id.
    static elementType getId(uint8_t count);
    //! Generate unique Element Id.
    static elementType getId(uint8_t proton, uint8_t neutron);
    
};

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
