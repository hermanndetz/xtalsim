/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "Atom.h"

//------------------------------------------------------------------------------

//! Initialize the random number genrator.
//! \param elementId Unique ID of the atoms Element.
//! \param position Position of atom in real space.
//! \param index Position of atom inside lattice.
//! \param material Pointer to the material the atom belongs to.
//! \param wasModified Marker indicating if atoms was modified recently.
//! \param state Defines, if atom was modified through interface generator, exchange reaction, ...
Atom::Atom(const elementType elementId, const Vector3D<spaceType> position,
	   const Vector3D<indexType> index, const Material *material,
           const bool wasModified, const bool wasMoved, const uint32_t modificatonIndex, 
           std::bitset<AtomState::StateCount> state):
    elementId_(elementId), position_(position),
    previousPosition_(position), index_(index),
    wasModified_(wasModified), wasMoved_(wasMoved), material_(material),
    modificationIndex_(modificatonIndex), state_{state}
{
    random_.seed(std::random_device()());
}

//------------------------------------------------------------------------------

Atom::Atom(const Atom &atom):
    elementId_(atom.elementId_), position_(atom.position_),
    previousPosition_(atom.position_), index_(atom.index_),
    wasModified_(atom.wasModified_), material_(atom.material_),
    neighbors_(atom.neighbors_), modificationIndex_(atom.modificationIndex_),
    state_(atom.state_)
{
    random_.seed(std::random_device()());    
}
    
//------------------------------------------------------------------------------

Atom::~Atom() {}

//------------------------------------------------------------------------------

//! \param distance Maximal amount of displacement.
void Atom::displace(const double distance)
{

    std::uniform_real_distribution<double> vectorDistribution(-1.0,1.0);
    std::uniform_real_distribution<double> lengthDistribution(0,1.0);
    
    Vector3D<spaceType> displacement = {
	(double)vectorDistribution(random_),
	(double)vectorDistribution(random_),
	(double)vectorDistribution(random_)};
    
    // set length of vector to distance, which is the maximum
    displacement *= ( (distance*distance) / (displacement.squaredLength()) );
    // vary it arbitrarly beteween 0 and distance
    displacement *= lengthDistribution(random_);

    previousPosition_ = position_;
    position_ += displacement;
    
}

//------------------------------------------------------------------------------

//! \return Unique Element ID of the atoms element.
elementType Atom::getElementId(void) const
{
    return elementId_;
}

//------------------------------------------------------------------------------

//! \return Position of atom in real space.
Vector3D<spaceType> Atom::getPosition(void) const
{
    return position_;
}

//------------------------------------------------------------------------------

//! \param position New position in real space.
//! \param savePosition If true the old position is backuped before it is
//! overwritten.
void Atom::setPosition(const Vector3D<spaceType>& position,
		       const bool savePosition)
{
    //if ( ! (position == position_))
    wasMoved_ = true;

    if (savePosition)
        previousPosition_ = position_;
    position_ = position;
}

//------------------------------------------------------------------------------

//! Set position of atom to previous location. This becomes necessary when
//! optimization moves did not result in an energy gain.
void Atom::setPreviousPosition(void)
{
    position_ = previousPosition_;
}

//------------------------------------------------------------------------------

Vector3D<indexType> Atom::getIndex(void) const
{
    return index_;
}

//------------------------------------------------------------------------------

//! \param position Position in real space of added neighbor.
void Atom::addNeighbor(const Vector3D<indexType> &position)
{
    std::lock_guard<std::mutex> guard(neighborMutex_);
    neighbors_.push_back(position);
}

//------------------------------------------------------------------------------

void Atom::clearNeighbors(void)
{
    std::lock_guard<std::mutex> guard(neighborMutex_);
    neighbors_.clear();
}

//------------------------------------------------------------------------------

void Atom::updateNeighbors(NeighborList neighbors)
{
    std::lock_guard<std::mutex> guard(neighborMutex_);
    neighbors_.swap(neighbors);
}

//------------------------------------------------------------------------------

//! \return List of all neighbors of the atom. The list may not be edited.
NeighborList const & Atom::getNeighbors(void) const
{
    return neighbors_;
}

//------------------------------------------------------------------------------

const Material * Atom::getMaterial(void) const
{
    return material_;
}

//------------------------------------------------------------------------------

//! Sets the element together with the material pointer.
//! \note The modification state is only changed, if old and new element differ.
//! \param elementId New element ID.
//! \param material New material pointer.
//! \param state Sets modification state (default ModifiedUnknown)
void Atom::setElement(const elementType elementId, const Material *material, const AtomState state)
{
    if (elementId_ != elementId) {
        state_.set(state);
    }

    elementId_=elementId;
    material_=material;
    wasModified_=true;
}

//------------------------------------------------------------------------------

//! \return True if atom was modified lately, false otherwise.
bool Atom::wasModified(void) const
{
    return wasModified_;
}

//------------------------------------------------------------------------------

//! resets the modification flag and the modification index
void Atom::clearModification(void)
{
    wasModified_ = false;
    modificationIndex_ = 0;
}

//------------------------------------------------------------------------------

//! Set modification order.
//! \param index the modification index to be stored
void Atom::setModificationOrder(uint32_t index) {
    modificationIndex_ = index;
}

//------------------------------------------------------------------------------

//! Get modification order.
//! \return the modification index
uint32_t Atom::getModificationOrder(void) const {
    return modificationIndex_;
}

//------------------------------------------------------------------------------

//! \return True if atom was modified lately, false otherwise.
bool Atom::wasMoved(void) const
{
    return wasMoved_;
}

//------------------------------------------------------------------------------

//! Returns modification state.
std::bitset<AtomState::StateCount> Atom::getState(void) const {
    return state_;
}

//------------------------------------------------------------------------------

//! Add modification state.
void Atom::addModificationState(const AtomState state) {
    state_.set(state);
}

//------------------------------------------------------------------------------

void Atom::clearMovement(void)
{
    wasMoved_ = false;
}

//------------------------------------------------------------------------------

//! Reset modification state.
void Atom::clearState(void) {
    state_ = AtomState_Default;
}

//------------------------------------------------------------------------------

//! Return info about object as string.
std::string Atom::str(void) const {
    std::stringstream result{};

    result << "ElementID: " << elementId_ << std::endl;
    result << "Position: " << position_.str() << std::endl;
    result << "Prev. Pos.: " << previousPosition_.str() << std::endl;
    result << "Index: " << index_.str() << std::endl;
    result << "Modified: " << wasModified_;

    if (modificationIndex_ > 0)
        result << " (" << modificationIndex_ << ")";
        
    result << std::endl;
    result << "Moved: " << wasMoved_ << std::endl;
    result << "State: " << state_.to_string() << std::endl;
    result << "Material: " << material_->getName() << std::endl;
    result << "Neighbors: ";

    for (auto item: neighbors_) {
        result << item.str() << " ";
    }

    result << std::endl;

    return result.str();
}

//------------------------------------------------------------------------------

//! \param atom Atom whose properties shall be copied.
Atom & Atom::operator=(const Atom &atom)
{
    if (this== &atom)
	return *this;
    
    elementId_ = atom.elementId_;
    position_ = atom.position_;
    previousPosition_ = atom.previousPosition_;
    index_ = atom.index_;
    material_ = atom.material_;
    wasModified_ = atom.wasModified_;
    modificationIndex_ = atom.modificationIndex_;
    state_ = atom.state_;

    neighbors_.clear();
    neighbors_ = atom.neighbors_;

    return *this;
}

//------------------------------------------------------------------------------

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
