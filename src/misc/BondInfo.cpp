/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "BondInfo.h"

//------------------------------------------------------------------------------

//! Constructor
BondInfo::BondInfo() {

}

//------------------------------------------------------------------------------

//! Constructor
BondInfo::BondInfo(const BondInfo &src) {
    data = src.data;
}

//------------------------------------------------------------------------------

//! Destructor
BondInfo::~BondInfo() {

}
//
//------------------------------------------------------------------------------

//! Appends composition info of one layer to internal data.
//! \param layer Composition info of a layer.
//! \todo return type does not match body
//! \todo maybe use other operator? This one is used in general when printing
//! information.
BondInfo & BondInfo::operator << (LayerBondInfo layer) {
    data.push_back(layer);
}

//------------------------------------------------------------------------------

//! Returns number of layers stored in internal data structure.
//! \return Number of layers stored in BondInfo object.
uint32_t BondInfo::getLayerCount (void) const {
    return data.size();
}

//------------------------------------------------------------------------------

//! Returns the maximum number of bonds per element pair per layer.
//! \return Maximum number of bonds per layer per element pair.
uint32_t BondInfo::getMaxBondCount(void) const {
    uint32_t result{};

    for (auto layer: data) {
        for (auto element: layer) {
            if (element.second > result)
                result = element.second;
        }
    }

    return result;
}

//------------------------------------------------------------------------------

//! Returns bond info of one given layer.
//! \param layerIndex Index of which the bond information is requested.
//! \return Bond info about one layer.
LayerBondInfo BondInfo::getLayer (uint32_t layerIndex) const {

    return data[layerIndex];
}

//------------------------------------------------------------------------------

//! Prints internal data to std::cout.
void BondInfo::dump (void) const {
    uint32_t i = 0;

    for (auto layer: data) {
        std::cout << "Layer " << i << std::endl;
        for (auto entry: layer) {
            auto key= std::get<0>(entry);
            std::cout << "Element 1: " << std::get<0>(key)
                      << "\tElement 2: " << std::get<1>(key)
                      << "\tCount: " << std::get<1>(entry) << std::endl;
        }

        i++;
    }
}

//------------------------------------------------------------------------------

std::vector<DataSeries> BondInfo::get(void) const {
    return dataSeries;
}

//------------------------------------------------------------------------------

#ifdef __VTK__

vtkSmartPointer<vtkTable> BondInfo::getCompressed(void) {
    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
    vtkSmartPointer<vtkLongLongArray> colLayer = vtkSmartPointer<vtkLongLongArray>::New();
    vtkSmartPointer<vtkStringArray> colElement = vtkSmartPointer<vtkStringArray>::New();
    vtkSmartPointer<vtkLongLongArray> colCount = vtkSmartPointer<vtkLongLongArray>::New();

    const PeriodicTable &pt = PeriodicTable::getInstance();
    
    colLayer->SetName("layer");
    colElement->SetName("elements");
    colCount->SetName("count");

    uint32_t i = 0;

    for (auto layer: data) {
        for (auto entry: layer) {
            auto key = std::get<0>(entry);
            colLayer->InsertNextValue(i);
            colElement->InsertNextValue(pt.getById(std::get<0>(key)).symbol + "-" + pt.getById(std::get<1>(key)).symbol);
            colCount->InsertNextValue(std::get<1>(entry));
        }

        i++;
    }

    table->AddColumn(colLayer);
    table->AddColumn(colElement);
    table->AddColumn(colCount);

    return table;
}

#endif

//------------------------------------------------------------------------------

#ifdef __VTK__

vtkSmartPointer<vtkTable> BondInfo::getSparse(void) {
    return getSparse(std::vector<std::tuple<elementType,elementType>>{}, false);
}

#endif

//------------------------------------------------------------------------------

#ifdef __VTK__

vtkSmartPointer<vtkTable> BondInfo::getSparse(
                std::vector<std::tuple<elementType,elementType>> filterBonds,
                bool interpolate) {

    const PeriodicTable &pt = PeriodicTable::getInstance();
    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

    // first element contains x-values, all other y-values
    std::vector<vtkSmartPointer<vtkLongLongArray>> colValues;

    // map contains elementId and column index of result table
    LayerBondInfo translation;

    uint32_t layerIndex = 0;

    colValues.push_back(vtkSmartPointer<vtkLongLongArray>::New());
    colValues[0]->SetName("layer");
    
    for (auto layer: data) {
        colValues[0]->InsertNextValue(layerIndex);

        for (auto entry: layer) {
            auto key = std::get<0>(entry);
            // The elementID is taken from the first element here.
            // It is later only used to define the color of the
            // data series.
            elementType elementID = std::get<0>(key);
            
            if (std::find(filterBonds.begin(), filterBonds.end(), key) != filterBonds.end())
                continue;

            // add element for specific material if not available
            if (translation.find(key) == translation.end()) {
                auto columnIndex = colValues.size();
                
                translation[key] = columnIndex;
                colValues.push_back(vtkSmartPointer<vtkLongLongArray>::New());

                std::string name = pt.getById(std::get<0>(key)).symbol + "-" + pt.getById(std::get<1>(key)).symbol;
                colValues.back()->SetName(name.c_str() );

                // as not element of each material appears in each layer
                // initialize complete vector with zeros.
                for (auto i=0; i<data.size(); i++)
                    colValues.back()->InsertNextValue(0);
                
                Color col = Color(pt.getById(elementID).color);
                DataSeries tmp = DataSeries(0, 0, columnIndex, col);
                dataSeries.push_back(tmp);
            }
            
            // add values for material/element combination in current layer
            auto search = translation.find(key);

            colValues[search->second]->SetValue(layerIndex, std::get<1>(entry));

            if (interpolate == true) {
                switch (layerIndex%2) {
                case 0:
                    colValues[search->second]->SetValue(layerIndex+1, std::get<1>(entry));
                    break;
                case 1:
                    colValues[search->second]->SetValue(layerIndex-1, std::get<1>(entry));
                    break;
                }
            }

        }  //for (auto entry: layer)

        layerIndex++;
        
    } // for (auto layer: data)
   
    for (auto i = 0; i < colValues.size(); i++)
        table->AddColumn(colValues[i]);

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

