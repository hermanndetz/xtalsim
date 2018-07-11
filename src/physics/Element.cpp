/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "Element.h"

//------------------------------------------------------------------------------

//! Automatically initializes the internal id with the given Proton and Neutron
//! count.
Element::Element(std::string _name, std::string _symbol,
		 uint8_t _protonCount, uint8_t _neutronCount,
		 uint8_t _period, uint8_t _group,
		 double _mass, double _weight, std::string _color)
    :protonCount(_protonCount), neutronCount(_neutronCount), name(_name),
     symbol(_symbol), period(_period), group(_group), mass(_mass),
     weight(_weight), color(_color)
{
    id = getId(protonCount, neutronCount);
}

//------------------------------------------------------------------------------

Element::~Element() {}

//------------------------------------------------------------------------------

//! The generated ID is used to uniquely identify an element, which is only
//! possible using both the neutron and proton count, as isotopes share the name
//! and the abbreviation.
//! \param count Proton Count of Element.
//! \return Unique Element ID.
elementType Element::getId(uint8_t count)
{
    return (elementType)((count << 8) + count);
}

//------------------------------------------------------------------------------

//! The generated ID is used to uniquely identify an element, which is only
//! possible using both the neutron and proton count, as isotopes share the name
//! and the abbreviation.
//! \param protonCount Proton count of Element.
//! \param neutronCount Neutron count of Element.
//! \return Unique Element ID.
elementType Element::getId(uint8_t protonCount, uint8_t neutronCount)
{
    return (elementType)((protonCount << 8) + neutronCount);
}

//------------------------------------------------------------------------------

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
