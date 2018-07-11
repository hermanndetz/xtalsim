/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __JOURNAL_H__
#define __JOURNAL_H__

#include <stdio.h>

#include <string>
#include <vector>
#include <type_traits>

#include <misc/Tuples.h>

#include "projectConfigure.h"

#ifdef __VTK__
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkTable.h>
#include <vtkUnsignedLongLongArray.h>
#include <vtkVariant.h>
#include <vtkVariantArray.h>
#endif

//! \brief Stores intermediate results of an internal procedure.

//! Mainly used for debug purposes but also to retrace the behaviour of the
//! probabilistic procedures. The results are supposed to be stored later on in
//! an XML file to be analyzed.
template <class T>
class Journal{
 private:

    //! Name of Journal.
    std::string name_;
    //! Optional description of the stored data.
    std::string description_;

    //! Stored messages.
    std::vector<T> entries_;
    
 public:
    //! Constructor
    Journal(const std::string &name="undef",
		const std::string &description="");
    //! Desctructor
    ~Journal();

    //! Return description.
    std::string getDescription(void) const;
    //! Return name.
    std::string getName(void) const;

    //! Set description of Journal.
    void setDescription(const std::string description);
    //! Set name of Journal.
    void setName(const std::string name);
    
    //! Add entry to log.
    void add(T entry);
    //! Return all stored entries.
    const std::vector<T>& getEntries(void) const;

    //! Clear all stored entries.
    void clear(void);

#ifdef __VTK__
    void get (vtkSmartPointer<vtkTable> table, std::vector<std::string>columnTitles=std::vector<std::string>());
#endif
};

//******************************************************************************

#include "Journal.inl"

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
