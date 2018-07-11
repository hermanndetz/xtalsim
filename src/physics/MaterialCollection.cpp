/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "MaterialCollection.h"

//------------------------------------------------------------------------------

//! The constructor registers initializes the logger.
//! \param logName Optional parameter defining the logger's name.
MaterialCollection::MaterialCollection(const char *logName):
    logName_(logName)
{
    //nothing to be done
}

//------------------------------------------------------------------------------

MaterialCollection::~MaterialCollection() {}

//------------------------------------------------------------------------------

//! \param material Material to add.
void MaterialCollection::add (const Material &material)
{
    CLOG(TRACE, logName_) << "adding material " << material.str() <<
	" to collection";
    materials_.push_back(material);
}

//------------------------------------------------------------------------------

//! \return All materials stored in the collection.
std::list<Material> const & MaterialCollection::get(void) const
{
    return materials_;
}

//------------------------------------------------------------------------------

//! \warning This method is only assured to
//! deliver the right material if the names are unique.
//! \param name Name of the desired material.
//! \throws MaterialException if material not found.
//! \ŗeturn Reference to the found entry inside the collection.
const Material & MaterialCollection::getByName(const std::string &name) const
{

    CLOG(TRACE, logName_) << "Searching for Material by name '" << name << "'";
    
    auto mat = find(materials_.begin(), materials_.end(), name);

    if (mat == materials_.end()) {
	CLOG(ERROR, logName_) << "Material " << name
			      << "' not found in collection";
        throw MaterialException("Material not found by name",
				MaterialException::Id::NameCollection);
    }

    return (*mat);
}

//------------------------------------------------------------------------------

//! Checks if a material is part of the collection. If not the material is
//! added.
//! \param material Material that shall be found or added.
//! \ŗeturn Reference to the found or new entry inside the collection.
const Material & MaterialCollection::getByNameOrAdd(const Material &material)
{

    CLOG(TRACE, logName_) << "Searching for Material by name "
			  << material.getName() << "' or add it";

    auto mat = find(materials_.begin(), materials_.end(), material.getName());

    if (mat == materials_.end()) {
	CLOG(DEBUG, logName_) << "Material '" << material.getName()
			    << "' not found in collection.";
	CLOG(DEBUG, logName_) << "Going to add it now!";
	    
	add(material);
	mat = find(materials_.begin(), materials_.end(), material.getName());
    }

    return (*mat);

}

//------------------------------------------------------------------------------

//! Searches for a material having the specified cations and anions. They have
//! to fit perfectly, i.e., the amount as well as the type has to be the same.
//! \param cations Desired cations.
//! \param anions Desired anions.
//! \return Found Material.
//! \throws MaterialException if no fitting material found.
const Material & MaterialCollection::getByCationsAnions(
		       const std::vector<elementType> &cations,
		       const std::vector<elementType> &anions) const
{
    std::string cationString = "[";
    std::string anionString = "[";
    bool allFound;
    const PeriodicTable &pt = PeriodicTable::getInstance();
    
    for (auto i: cations) {
        cationString += std::to_string(i);
        cationString += ("(" + pt.getById(i).symbol + ")");
        cationString += ",";
    }

    cationString.back()= ']';
    
    for (auto i: anions) {
        anionString += std::to_string(i);
        anionString += ("(" + pt.getById(i).symbol + ")");
        anionString += ",";
    }

    anionString.back()= ']';

    CLOG(TRACE, logName_) << "Trying to find material by cations " <<
	cationString << " and anions " << anionString;

    for(auto it=materials_.begin(); it != materials_.end(); ++it){
        const Material &material = *it;

        if (material.getCations().size() != cations.size())
            continue;

        if (material.getAnions().size() != anions.size())
            continue;

        allFound = true;
        for (auto elementID : cations){
            if ( ! material.hasCation(elementID)){
            allFound = false;
            break;
            }
        }

        if ( ! allFound) continue;

        allFound = true;
        for (auto elementID : anions){
            if ( ! material.hasAnion(elementID)){
            allFound = false;
            break;
            }
        }
        
        if ( ! allFound) continue;

        CLOG(TRACE, logName_) << "Found suitable material with name " <<
            (*it).getName() << "!";
        return (*it);
    }
    
    // Disabled log output here.
    // This function is called with two element IDs
    // and tries both combinations as cations and anions (AB and BA)
    // Therefore, this error is expected to happen frequently.
    //CLOG(INFO, logName_) << "No material with cations " <<
	//cationString << " and anions " << anionString << " found!";
    throw MaterialException("Material not found by cations " + 
            cationString + " and anions " + anionString,
			    MaterialException::Id::CationAnionCollection);
}

//------------------------------------------------------------------------------

const Material & MaterialCollection::getByElementID(const elementType elementA, 
        const elementType elementB) const {

    for (auto material = materials_.begin(); material != materials_.end();
	 ++material) {
        // pretend elementA is the cation and elementB is the anion
        if ((material->isOnlyCation(elementA) == true) &&
            (material->isOnlyAnion(elementB) == true)) {
            return *material;
        }

        // pretend elementA is the anion and elementB is the cation
        if ((material->isOnlyCation(elementB) == true) &&
            (material->isOnlyAnion(elementA) == true)) {
            return *material;
        }
    }

    CLOG(ERROR, logName_) << "No material consisting of " <<
	elementA << " and " << elementB << "found!";
    throw MaterialException("Material not found by IDs",
			    MaterialException::Id::CationAnionCollection);

}

//------------------------------------------------------------------------------

//! Get information about stored materials.
const std::string MaterialCollection::str(void) const {
    std::stringstream result{};

    result << "MaterialCollection contains " << materials_.size() << " materials" << std::endl;

    for (auto material: materials_) {
        result << "Material: " << material.str(true) << std::endl;
    }

    return result.str();
}

//------------------------------------------------------------------------------


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
