/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#include "ExchangeReaction.h"

//******************************************************************************

//! Constructor
ExchangeReaction::ExchangeReaction(const char *logName) : logName_(logName)
{
    random_.seed(std::random_device()());
}

//------------------------------------------------------------------------------

//! Destructor
ExchangeReaction::~ExchangeReaction() {

}

//------------------------------------------------------------------------------

//! \param elementId1 Unique element ID of first involved element.
//! \param elementId1 Unique element ID of second involved element.
inline std::string ExchangeReaction::getId_(const std::string material1,
                                            const std::string material2) const
{
    
	return ( material1 + '_' + material2);

}

//------------------------------------------------------------------------------

//! Extracts element ID from unique ID for reaction
//! \param id Unique reaction ID
//! \param first returns first element, if true, otherwise second
//! \return element id of first or second element
inline std::string ExchangeReaction::getMaterial_(const std::string id, bool first) {

    std::size_t found = id.find("_");
    
    if (first)
        return id.substr(0,found);
    else
        return id.substr(found+1);

}

//------------------------------------------------------------------------------

//! Adds an entry for an exchange reaction with a given probability
//! \param elements Tuple with element IDs of form (orignal, replacement)
//! \param probability [0:1]
void ExchangeReaction::addReaction (const std::string material1, const std::string material2, 
                                    const double probability) 
{
    CLOG(DEBUG, logName_) << "adding reaction " << material1 <<
        "->" << material2 << ": " << probability;

    reactions_[getId_(material1, material2)] = probability;
}

//------------------------------------------------------------------------------

//! Adds an entry for an exchange reaction with a given probability
//! \param reactionDefinition Reaction definition in the form
//! Material1_Material2_probability.
void ExchangeReaction::addReaction (const std::string reactionDefinition)
{

    CLOG(DEBUG, logName_) << "got reaction definition " << reactionDefinition;
    
    std::size_t found = reactionDefinition.find("_");
    std::string originalMaterial = reactionDefinition.substr(0,found);

    std::size_t found2 = reactionDefinition.find("_",found+1);
    std::string replacementMaterial = reactionDefinition.substr(found+1, found2-found-1);
    double probability = std::stod(reactionDefinition.substr(found2+1));

    addReaction(originalMaterial,replacementMaterial, probability);

}

//------------------------------------------------------------------------------

//! Replaces anions in the specified layer according to the defined reactions
//! If the reaction is not defined (no probability is set for a given species
//! pair), this function does not do anything.
//! \param layerID layer index, where pattern is generated.
//! \param dimension normal to the layer plane. (0 = (100), 1 = (010), 2 = (001))
//! \param original identifies the material to be replaced
//! \param replacement identifies the material that replaces original anions
void ExchangeReaction::replaceAnions (SimulationBox *simbox, const indexType layerID,
                                      const indexType dimension, 
                                      const std::string original, const std::string replacement) 
{
    Lattice &lattice = simbox->getLattice();

    CLOG(TRACE, logName_) << "Requested exchange reaction from material " <<
        original << " to replacement " << replacement;

    auto reaction = reactions_.find(getId_(original, replacement));

    if (reaction != reactions_.end()) {
        // reaction is defined -> do the work
        CLOG(DEBUG, logName_) << "Performing exchange reaction at probability " <<
            reactions_[getId_(original, replacement)];

        std::vector<Atom *> atoms = lattice.getAtomsInLayer(layerID, dimension);

        CLOG(DEBUG, logName_) << "Atoms in layer: " << atoms.size();

        std::uniform_real_distribution<double> distribution(0.0,1.0);

        uint32_t count = 0;
        std::string id = getId_(original, replacement);
        const Material &material = simbox->getMaterials().getByName(replacement);

        for (auto atom: atoms) {
            if (atom->getMaterial()->getName() == original) {

                if (distribution(random_) <= reactions_[id]) {
                    atom->setElement(material.getRandomAnion(random_),
                                     & material, AtomState::ModifiedExchangeReaction);
                    count ++;
                }
            }
        }

	    CLOG(DEBUG, logName_) << "Exchanged " << count << " atoms in layer " << layerID;
    }
}

//------------------------------------------------------------------------------

//! Returns information about stored reactions
//! \return string containing information about stored reactions
std::string ExchangeReaction::str () {
    std::stringstream result;

    for (auto reaction: reactions_) {
        result << getMaterial_(reaction.first, true) << " -> " <<
            getMaterial_(reaction.first, false) << ": " <<
            reaction.second;
    }

    return result.str();
}

//##############################################################################


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
