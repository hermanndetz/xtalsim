/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __JOURNAL_INL__
#define __JOURNAL_INL__

template <class T>
Journal<T>::Journal(const std::string &name,
			 const std::string &description):
    name_(name), description_(description)
{
    // nothing to be done    
}

//------------------------------------------------------------------------------

template <class T>
Journal<T>::~Journal()
{
    // nothing to be done
}

//------------------------------------------------------------------------------

template <class T>
std::string Journal<T>::getDescription(void) const
{
    return description_;
}

//------------------------------------------------------------------------------

template <class T>
void Journal<T>::setDescription(const std::string description)
{
    description_ = description;
}

//------------------------------------------------------------------------------

template <class T>
std::string Journal<T>::getName(void) const
{
    return name_;
}

//------------------------------------------------------------------------------

template <class T>
void Journal<T>::setName(const std::string name)
{
    name_ = name;
}

//------------------------------------------------------------------------------

template <class T>
void Journal<T>::clear(void)
{
    entries_.clear();
}

//------------------------------------------------------------------------------

// \param entry Entry that shall be added to the Journal.
template <class T>
void Journal<T>::add(T entry)
{
    entries_.push_back(entry);
}

//------------------------------------------------------------------------------

template <class T>
const std::vector<T>& Journal<T>::getEntries(void) const
{
    return entries_;
}

//------------------------------------------------------------------------------

#ifdef __VTK__
template <class T> void Journal<T>::get (vtkSmartPointer<vtkTable> table, std::vector<std::string>columnTitles)
{
    table->Initialize();

    // Ugly code duplication, necessary in lack of better solution!
    // The specialized function is not called by the compiler because
    // parameters and return types do not differ.
    // 
    // vtkChart needs a table with vtkUnsignedLongLongArray or
    // vtkDoubleArray. If given a vtkVariantArray, it fails with the
    // not so meaningful error "vtkPlotLine ...: No X column is set (index 0)."
    //
    // \todo This works as long as the only compound type is an UDTuple!
    if (std::is_compound<T>::value) {
        vtkSmartPointer<vtkUnsignedLongLongArray>colX = vtkSmartPointer<vtkUnsignedLongLongArray>::New();
        vtkSmartPointer<vtkDoubleArray>colY = vtkSmartPointer<vtkDoubleArray>::New();

        colX->Initialize();
        colY->Initialize();

        if (columnTitles.size() > 0) {
            colX->SetName(columnTitles[0].c_str());

            if ((columnTitles.size() > 1) && (std::is_compound<T>::value))
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
    } else {
        vtkSmartPointer<vtkVariantArray> colX = vtkSmartPointer<vtkVariantArray>::New();
        colX->Initialize();

        if (columnTitles.size() > 0)
            colX->SetName(columnTitles[0].c_str());
        else
            colX->SetName("X");

        for (auto const &entry: entries_)
            colX->InsertNextValue(std::get<0>(entry));

        table->AddColumn(colX);
    }
}
#endif

//------------------------------------------------------------------------------


#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
