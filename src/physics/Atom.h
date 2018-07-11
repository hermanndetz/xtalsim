/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __ATOM_H__
#define __ATOM_H__

#include <string>
#include <random>
#include <vector>
#include <unordered_map>
#include <mutex>

#include <misc/Hash.h>
#include <physics/AtomState.h>
#include <physics/Vector3D.h>
#include <physics/Material.h>

//! \todo change this to Atom *?? This way safer however the Atom objects should
//! not be moved anyway.
//typedef std::unordered_map< Vector3D<indexType>, double, Hash > NeighborList;
typedef std::vector< Vector3D<indexType> > NeighborList;

//! This class stores all the information of a single atom.

//! \todo optimization of neighbor list possible, e.g., by cleverer update
//! function, storing and updating the distance, ....

class Atom{
 private:

    elementType elementId_; //!< unique element id
    Vector3D<spaceType> position_; //!< position in space
    Vector3D<spaceType> previousPosition_; //!< previous position in space
    Vector3D<indexType> index_; //!< indices of atomic layer

    bool wasModified_; //!< indicates if element of atoms was changed
    bool wasMoved_; //!< indicates if atom was moved in real space
    std::bitset<AtomState::StateCount> state_;

    uint32_t modificationIndex_; //!< stores the order in which the atoms were modified
    
    const Material *material_; //!< material the atom belongs to

    //! Random number generator.
    std::mt19937 random_;

    //! List of atoms within a specific radius around the atom.
    NeighborList neighbors_;

    //! controls concurrent access to neighbor list.
    std::mutex neighborMutex_;

 public:

    //! Constructor
    Atom(const elementType elementId, const Vector3D<spaceType> position,
         const Vector3D<indexType> index, const Material *material,
         const bool wasModified=false, const bool wasMoved=false, const uint32_t modificatonIndex=0, 
         const std::bitset<AtomState::StateCount> state=AtomState_Default);

    //! Copy Constructor
    Atom(const Atom &atom);
    
    //! Desctructor
    ~Atom();

    //! Displace atom randomly in space.
    void displace(const double distance);
    
    //! Return Element ID of the atoms element.
    elementType getElementId(void) const;
    //! Return position of the atom in space.
    Vector3D<spaceType> getPosition(void) const;
    //! Set position of atom in space.
    void setPosition(const Vector3D<spaceType>& position,
		     const bool savePosition=false);
    //! Set position of atom in space to previous location.
    void setPreviousPosition(void);
    //! Return position of the atom in the lattice.
    Vector3D<indexType> getIndex(void) const;
    //! Return pointer to material the atom belongs to.
    const Material *getMaterial(void) const;
    //! Set element ID and material.
    void setElement(const elementType elementId, const Material *material, const AtomState state=AtomState::ModifiedUnknown);
    
    //! Indicates if atom element was changed.
    bool wasModified(void) const;
    //! Returns modification state.
    std::bitset<AtomState::StateCount> getState(void) const;
    //! Add modification state.
    void addModificationState(const AtomState state);
    //! Reset modification attribute and clear modification index.
    void clearModification(void);
    //! Reset modification state.
    void clearState(void);
    //! Set modification order.
    void setModificationOrder(uint32_t index);
    //! Get modification order.
    uint32_t getModificationOrder(void) const;
    //! Indicates if atom was moved in real space.
    bool wasMoved(void) const;
    //! Reset modification attribute.
    void clearMovement(void);
    
    //! Add entry to neighbor list.
    void addNeighbor(const Vector3D<indexType> &position);
    //! Clears neighbor list.
    void clearNeighbors(void);
    //! Clears neighbor list.
    void updateNeighbors(NeighborList neighbors);
    //! Returns neighbor list.
    NeighborList const & getNeighbors(void) const;

    //! Return info about object as string.
    std::string str(void) const;

    //! Assignment operator
    Atom & operator=(const Atom &atom);
};

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
