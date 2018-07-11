/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __MATERIAL_COLLECTION_H__
#define __MATERIAL_COLLECTION_H__

#include <list>
#include <algorithm>
#include <sstream>
#include <string>

#include <easyloggingcpp/easylogging++.h>

#include <physics/Material.h>
#include <physics/PeriodicTable.h>

//! \brief Stores available materials.

//! Available Materials are stored in this function. It is used to pass them to
//! a function and provides an interface for searching and selecting specific
//! materials.

class MaterialCollection{
private:
    //! All available Materials. This has to be a list and not a std::vector
    //! since only a pointer is stored in the atom class. Therefore the address
    //! of each single element must not be altered.
    std::list<Material> materials_;

    //! Logger name.
    const char* logName_;
    
public:

    //! Constructor
    MaterialCollection(const char *logName="MaterialCollection");
    //! Destructor
    ~MaterialCollection();

    //! Add Material to the collection.
    void add (const Material &material);

    //! Get all Materials.
    std::list<Material> const & get(void) const;
    //! Select Material by name.
    const Material & getByName(const std::string &name) const;
    //! Select Material by name or adds it if not available.
    const Material & getByNameOrAdd(const Material &material);
    //! Select Material by cations and anions.
    const Material & getByCationsAnions(const std::vector<elementType> &cations,
				 const std::vector<elementType> &anions) const;
    //! Select Material by element ID.
    const Material & getByElementID(const elementType elementA,
				    const elementType elementB) const;

    //! Get information about stored materials.
    const std::string str(void) const;
};

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
