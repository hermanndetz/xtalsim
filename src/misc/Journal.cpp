/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "Journal.h"

//------------------------------------------------------------------------------

#ifdef __VTK__
template <>
void Journal<UDTuple>::get (vtkSmartPointer<vtkTable> table, std::vector<std::string>columnTitles) {
    table->Initialize();
    vtkSmartPointer<vtkUnsignedLongLongArray> colX = vtkSmartPointer<vtkUnsignedLongLongArray>::New();
    vtkSmartPointer<vtkDoubleArray> colY = vtkSmartPointer<vtkDoubleArray>::New();
    colX->Initialize();
    colY->Initialize();

    if (columnTitles.size() >= 2) {
        colX->SetName(columnTitles[0].c_str());
        colY->SetName(columnTitles[1].c_str());
    } else {
        colX->SetName("X");
        colY->SetName("Y");
    }

    for (auto const &entry: entries_) {
        colX->InsertNextValue((long long unsigned int)std::get<0>(entry));
        colY->InsertNextValue((double)std::get<1>(entry));
    }

    table->AddColumn(colX);
    table->AddColumn(colY);
}
#endif

//------------------------------------------------------------------------------


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
