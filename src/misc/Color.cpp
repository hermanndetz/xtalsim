/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "Color.h"

//------------------------------------------------------------------------------

//! Automatically initializes color components, if given as parameters.
Color::Color (const colorType red, const colorType green,
	      const colorType blue) : rgb_(red,green,blue)
{
    // nothing to be done
}

//------------------------------------------------------------------------------

//! Automatically initializes color components, if given as parameters.
//! Default value is RGB (211, 211, 211) = \#D3D3D3 = gray.
//! \param webcolor Hex string representation of desired color.
Color::Color (const std::string webcolor)
{
    // auto getComponents = [this, webcolor] (uint8_t offset) {
    //     char tmpChar[3]{};

    //     for (int i : {0, 1, 2}) {
    //         tmpChar[0] = webcolor[offset+i*2];
    //         tmpChar[1] = webcolor[offset+i*2+1];
    //         this->rgb_[i] = strtoul((const char *)&tmpChar, NULL, 16);
    //     }
    // };

    // if (webcolor.length() == 6) {
    //     getComponents(0);
    // } else if ((webcolor.length() == 7) && (webcolor[0] == '#')) {
    //     getComponents(1);
    // }

    rgb_[0] = defaultColorValue;
    rgb_[1] = defaultColorValue;
    rgb_[2] = defaultColorValue;

    unsigned int offset = 0;
    
    if (webcolor[0] == '#')
        offset=1;

    if (webcolor.length() == (6 + offset)) {
        for (int i=0; i< 3; i++)
            rgb_[i] = (colorType) std::stoi(webcolor.substr(offset+2*i, 2), 0, 16);
    }
}

//------------------------------------------------------------------------------

Color::Color (const Color &color)
{
    rgb_ = color.rgb_;
}

//------------------------------------------------------------------------------

//! Destructor
Color::~Color () {

}

//------------------------------------------------------------------------------

void Color::operator=(const Color &color)
{
    rgb_ = color.rgb_;
}

//------------------------------------------------------------------------------

//! Returns the color description as web-compatible string '\#rrggbb'.
//! \return Color code as string.
std::string Color::str (void) const
{
    std::stringstream stream;

    stream << "#";

    for (int i : {0, 1, 2}) {
        //! Attention! Cast to int necessary, otherwise uint8_t is cast
        //! to std::string directly.
        stream << std::setfill('0') << std::setw(2) << std::hex << (int)rgb_[i];
    }

    return stream.str();
}

//------------------------------------------------------------------------------

Vector3D<colorType> Color::get (void) const
{
    return rgb_;
}

//------------------------------------------------------------------------------

colorType Color::get (const uint8_t index) const
{
    return rgb_[index];
}

//------------------------------------------------------------------------------

double Color::get_double (const uint8_t index) const
{
    return (double)(rgb_[index])/255.0;
}

//##############################################################################


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
