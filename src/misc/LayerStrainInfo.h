/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __LAYER_STRAIN_INO__
#define __LAYER_STRAIN_INO__

#include <tuple>
#include <vector>

#include "projectConfigure.h"

#include <physics/Vector3D.h>
#include <visualization/DataSeries.h>

#ifdef __VTK__
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkLongLongArray.h>
#include <vtkDoubleArray.h>
#endif

//! Vector containing strain info for each layer.
typedef std::tuple<indexType, double, double, double> LayerStrainInfoData;

class LayerStrainInfo {
private:
    std::vector<LayerStrainInfoData> data;
    std::vector<DataSeries> dataSeries;

public:
    //! Constructor
    LayerStrainInfo();
    //! Constructor
    LayerStrainInfo(const LayerStrainInfo &src);
    //! Destructor
    ~LayerStrainInfo();

    //! Appends strain info of one layer to internal data.
    LayerStrainInfo & operator << (LayerStrainInfoData layer);

    //! Returns number of layers stored in internal data structure.
    uint32_t getLayerCount (void) const;
    //! Returns composition info of one given layer.
    LayerStrainInfoData getLayer (uint32_t layerIndex) const;

    //! Returns data series description for XYPlot.
    std::vector<DataSeries> get(void) const;

    //! Returns maximum srain value (independent of column)
    double getMaximumStrain(void) const;

    //! Returns minimum srain value (independent of column)
    double getMinimumStrain(void) const;

#ifdef __VTK__
    vtkSmartPointer<vtkTable> getTable();
#endif
};


#endif


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
