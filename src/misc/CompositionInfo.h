/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __COMPOSITION_INFO_H__
#define __COMPOSITION_INFO_H__

#include <iostream>
#include <unordered_map>
#include <vector>

#include <physics/PeriodicTable.h>
#include <physics/Element.h>
#include <visualization/DataSeries.h>
#include <misc/Hash.h>

#include "projectConfigure.h"

#ifdef __VTK__
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkLongLongArray.h>
#include <vtkStringArray.h>
#endif

// Key: <Material name, elementID>
typedef std::unordered_map<std::tuple<std::string, elementType>, uint32_t, Hash> LayerCompositionInfo;

class CompositionInfo {
private:
    std::vector<LayerCompositionInfo> data;
    std::vector<DataSeries> dataSeries;

public:
    //! Constructor
    CompositionInfo();
    //! Constructor
    CompositionInfo(const CompositionInfo &src);
    //! Destructor
    ~CompositionInfo();

    //! Appends composition info of one layer to internal data.
    CompositionInfo & operator << (LayerCompositionInfo layer);

    //! Returns number of layers stored in internal data structure.
    uint32_t getLayerCount (void) const;
    //! Returns the maximum number of atoms per element per layer.
    uint32_t getMaxAtomCount (void) const;
    //! Returns composition info of one given layer.
    LayerCompositionInfo getLayer (uint32_t layerIndex) const;
    //! Prints internal data to std::cout.
    void dump (void) const;

    //! Returns data series description for XYPlot.
    std::vector<DataSeries> get(void) const;


#ifdef __VTK__
    vtkSmartPointer<vtkTable> getCompressed(void);
    vtkSmartPointer<vtkTable> getSparse(void);
    vtkSmartPointer<vtkTable> getSparse(std::vector<elementType> filterElements, bool interpolate);
#endif
};

#endif


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
