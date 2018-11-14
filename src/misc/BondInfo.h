/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/



#ifndef __BONDINFO_H__
#define __BONDINFO_H__

#include <unordered_map>
#include <vector>

#include <misc/Hash.h>
#include <physics/Element.h>
#include <physics/PeriodicTable.h>
#include <visualization/DataSeries.h>

#include "projectConfigure.h"

#ifdef __VTK__
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkLongLongArray.h>
#include <vtkStringArray.h>
#endif

typedef std::unordered_map<std::tuple<elementType,elementType>, uint32_t, Hash> LayerBondInfo;

class BondInfo {
private:
    std::vector<LayerBondInfo> data;
    std::vector<DataSeries> dataSeries;

public:
    //! Constructor
    BondInfo();
    //! Constructor
    BondInfo(const BondInfo &src);
    //! Destructor
    ~BondInfo();

    //! Appends bond info of one layer to internal data.
    BondInfo & operator << (LayerBondInfo layer);

    //! Returns number of layers stored in internal data structure.
    uint32_t getLayerCount (void) const;
    //! Returns the maximum number of bonds per element pair per layer.
    uint32_t getMaxBondCount (void) const;
    //! Returns bond info of one given layer.
    LayerBondInfo getLayer (uint32_t layerIndex) const;
    //! Prints internal data to std::cout.
    void dump (void) const;

    //! Returns data series description for XYPlot.
    std::vector<DataSeries> get(void) const;

#ifdef __VTK__
    vtkSmartPointer<vtkTable> getCompressed(void);
    vtkSmartPointer<vtkTable> getSparse(void);
    vtkSmartPointer<vtkTable> getSparse(std::vector<std::tuple<elementType,elementType>> filterBonds, bool interpolate);
#endif
};

#endif


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:

