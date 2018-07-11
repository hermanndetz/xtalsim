/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "LayerStrainInfo.h"

//------------------------------------------------------------------------------

//! Constructor
LayerStrainInfo::LayerStrainInfo() {

}

//------------------------------------------------------------------------------

//! Constructor
LayerStrainInfo::LayerStrainInfo(const LayerStrainInfo &src) {
    data = src.data;
    dataSeries = src.dataSeries;
}

//------------------------------------------------------------------------------

//! Destructor
LayerStrainInfo::~LayerStrainInfo() {

}

//------------------------------------------------------------------------------

//! Appends strain info of one layer to internal data.
//! \param layer Strain info of a layer.
LayerStrainInfo & LayerStrainInfo::operator << (LayerStrainInfoData layer) {
    data.push_back(layer);
}

//------------------------------------------------------------------------------

//! Returns number of layers stored in internal data structure.
//! \return Number of layers stored in CompositionInfo object.
uint32_t LayerStrainInfo::getLayerCount (void) const {
    return data.size();
}

//------------------------------------------------------------------------------

//! Returns composition info of one given layer.
//! \param layerIndex Index of which the strain information is requested.
//! \return Strain info about one layer.
LayerStrainInfoData LayerStrainInfo::getLayer (uint32_t layerIndex) const {
    return data[layerIndex];
}

//------------------------------------------------------------------------------

std::vector<DataSeries> LayerStrainInfo::get() const {
    return dataSeries;
}

//------------------------------------------------------------------------------

//! Returns maximum srain value (independent of column)
//! \return Maximum strain value.
double LayerStrainInfo::getMaximumStrain(void) const {
    double result{};

    for (auto layer : data) {
        if (std::get<1>(layer) > result)
            result = std::get<1>(layer);
        if (std::get<2>(layer) > result)
            result = std::get<2>(layer);
        if (std::get<3>(layer) > result)
            result = std::get<3>(layer);
    }

    return result;
}

//------------------------------------------------------------------------------

//! Returns minimum srain value (independent of column)
//! \return Minimum strain value.
double LayerStrainInfo::getMinimumStrain(void) const {
    double result{};

    for (auto layer : data) {
        if (std::get<1>(layer) < result)
            result = std::get<1>(layer);
        if (std::get<2>(layer) < result)
            result = std::get<2>(layer);
        if (std::get<3>(layer) < result)
            result = std::get<3>(layer);
    }

    return result;
}

//------------------------------------------------------------------------------

#ifdef __VTK__

vtkSmartPointer<vtkTable> LayerStrainInfo::getTable() {
    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

    vtkSmartPointer<vtkLongLongArray> colLayer = vtkSmartPointer<vtkLongLongArray>::New();
    vtkSmartPointer<vtkDoubleArray> colBondStrain = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> colGrowthStrain = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> colRefDistanceStrain = vtkSmartPointer<vtkDoubleArray>::New();

    colLayer->SetName("layer");
    colBondStrain->SetName("bond strain");
    colGrowthStrain->SetName("growth strain");
    colRefDistanceStrain->SetName("ref dist strain");

    for (auto entry: data) {
        colLayer->InsertNextValue(std::get<0>(entry));
        colBondStrain->InsertNextValue(std::get<1>(entry));
        colGrowthStrain->InsertNextValue(std::get<2>(entry));
        colRefDistanceStrain->InsertNextValue(std::get<3>(entry));
    }

    table->AddColumn(colLayer);
    table->AddColumn(colBondStrain);
    table->AddColumn(colGrowthStrain);
    table->AddColumn(colRefDistanceStrain);

    Color colBond = Color("#D73027");
    DataSeries dsBond = DataSeries(0, 0, 1, colBond);
    Color colGrowth = Color("#4575B4");
    DataSeries dsGrowth = DataSeries(0, 0, 2, colGrowth);
    Color colRefDistance = Color("#FEE090");
    DataSeries dsRefDistance = DataSeries(0, 0, 3, colRefDistance);

    dataSeries.push_back(dsBond);
    dataSeries.push_back(dsGrowth);
    dataSeries.push_back(dsRefDistance);

    return table;
}

#endif

//------------------------------------------------------------------------------


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
