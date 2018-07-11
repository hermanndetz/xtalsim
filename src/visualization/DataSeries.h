/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#ifndef __DATA_SERIES__
#define __DATA_SERIES__

#include <misc/Color.h>

//! This is a container class to store information about a
//! single data series in a plot. It stores the information,
//! which column is plotted against which column of a table.
//!
//! It is implemented as a class because it will be extended
//! with more elaborate contructors, which will automatically
//! assign color codes.
class DataSeries {
private:
    uint32_t    tableIndex;
    uint32_t    xIndex;
    uint32_t    yIndex;
    Color       colorCode;

public:
    //! Constructor
    DataSeries (uint32_t x, uint32_t y, Color &color);
    DataSeries (uint32_t table, uint32_t x, uint32_t y, Color &color);
    DataSeries (const DataSeries &ds);

    //! Destructor
    ~DataSeries ();

    //! Returns color of object.
    Color       GetColor ();
    //! Returns table index of object.
    uint32_t    GetTableIndex ();
    //! Returns index in x direction of object.
    uint32_t    GetXIndex ();
    //! Returns index in y direction of object.
    uint32_t    GetYIndex ();
    //! Returns content  of the object as string.
    std::string str(void) const;
};

#endif


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
