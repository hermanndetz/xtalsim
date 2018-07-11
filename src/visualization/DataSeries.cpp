/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "DataSeries.h"

#include <string>

//! Constructor
DataSeries::DataSeries (uint32_t x, uint32_t y, Color &color) :
    tableIndex{}, xIndex(x), yIndex(y), colorCode(color) {

}

//------------------------------------------------------------------------------
//! Constructor
DataSeries::DataSeries (uint32_t table, uint32_t x, uint32_t y, Color &color) :
    tableIndex(table), xIndex(x), yIndex(y), colorCode(color) {

}

//------------------------------------------------------------------------------

//! Constructor
DataSeries::DataSeries (const DataSeries &ds) {
    tableIndex = ds.tableIndex;
    xIndex = ds.xIndex;
    yIndex = ds.yIndex;
    colorCode = ds.colorCode;
}

//------------------------------------------------------------------------------

//! Destructor
DataSeries::~DataSeries () {

}

//------------------------------------------------------------------------------

//! \return color
Color DataSeries::GetColor () {
    return colorCode;
}

//------------------------------------------------------------------------------

//! \return table index
uint32_t DataSeries::GetTableIndex () {
    return tableIndex;
}

//------------------------------------------------------------------------------

//! \return index in x direction
uint32_t DataSeries::GetXIndex () {
    return xIndex;
}

//------------------------------------------------------------------------------

//! \return index in y direction
uint32_t DataSeries::GetYIndex () {
    return yIndex;
}

//------------------------------------------------------------------------------

//! \return content of object in string format
std::string DataSeries::str(void) const {
    std::stringstream stream;

    stream << tableIndex << ":" << xIndex << ":" << yIndex << " - " <<
        colorCode.str();

    return stream.str();
}

//##############################################################################


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
