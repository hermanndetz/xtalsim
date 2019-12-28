/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#ifndef __COLOR_H__
#define __COLOR_H__

#include <string>
#include <stdint.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdlib.h>

#include <physics/Vector3D.h>

typedef uint8_t colorType;

//! Default color is grey. red = green = blue = defaultColorValue
const colorType defaultColorValue=211;

//! \brief Stores color information in RGB format.

class Color {
private:
    //! color components [0]: red, [1]: green, [2]: blue
    Vector3D<colorType> rgb_; 

public:

    //! Index of component red in RGB vector.
    static const uint8_t Red = 0;
    //! Index of component green in RGB vector.
    static const uint8_t Green = 1;
    //! Index of component blue in RGB vector.
    static const uint8_t Blue = 2;    
        
    //! Constructor - default initialization for 'gray'
    Color (const colorType red=defaultColorValue,
	   const colorType green=defaultColorValue,
	   const colorType blue=defaultColorValue);
    Color (const std::string webcolor);

    //! Copy Constructor
    Color (const Color &color);

    //! Destructor
    ~Color ();

    //! Assignment operator overloading;
    void operator=(const Color &col);

    //! Convert color into string.
    std::string str (void) const;
    //! Returns all rgb values of color.
    Vector3D<colorType> get (void) const;
    //! Select single component. 
    colorType get (const uint8_t index) const;
    //! Select single component converted to double [0:1]
    double get_double (const uint8_t index) const;
};

#endif


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
