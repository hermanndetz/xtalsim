/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#include "Lattice.h"

//------------------------------------------------------------------------------

//! The constructor initializes the logger.
//! \param temperature Lattice Temperature.
//! \param logName Optional parameter defining the logger's name.
Lattice::Lattice(const double temperature,
		 const char *logName):
    size_(0,0,0), temperature_(temperature), logName_(logName)
{
    random_.seed(std::random_device()());
}

//------------------------------------------------------------------------------

//! The atoms inside the lattice get created with 'new' so they have to be
//! deleted if the object is destroyed.
Lattice::~Lattice() {

    for (int i=0; i<size_[0]; i++){
	for (int j=0; j<size_[1]; j++){
	    for (int k=0; k<size_[2]; k++){

		if (lattice_[i][j][k] != nullptr )
		    delete(lattice_[i][j][k]);
	    }
	}
    }
}

//##############################################################################

//! Create a lattice if the internal lattice is still empty.
//! \throws LatticeException If the internal lattice is not empty or the given
//! size has at least in one dimension a zero.
//! \param x,y,z Defines the lattice size in the first, second and third
//! dimension.
//! \param dimension Specifies the growth dimension.
void Lattice::generate(const indexType x, const indexType y, const indexType z,
		       const uint8_t dimension)
{
    generate(Vector3D<indexType>(x,y,z), dimension);
}

//------------------------------------------------------------------------------

//! Create a lattice if the internal lattice is still empty.
//! \throws LatticeException If the internal lattice is not empty or the given
//! size has at least in one dimension a zero.
//! \param size Defines the lattice size. Each dimension has to have a value >=
//! 1
//! \param dimension Specifies the growth dimension.
void Lattice::generate(const Vector3D<indexType> &size, const uint8_t dimension)
{

    CLOG(TRACE, logName_) << "creating lattice of size " << size.str();
    
    // Lattice is not empty so extend it.
    if ( size_ > 0){
	extend(size[dimension], dimension);
	return;
    }

    // parameter contains at least in one dimension a zero
    if (! (size > 0)){
	CLOG(ERROR, logName_) << "Lattice with size zero can not be created!";
	throw LatticeException("Lattice with dimension zero can not be"	\
			       "created!", LatticeException::Id::SizeZero);
    }

    lattice_.resize(size[0]);

    for (unsigned int i=0; i<lattice_.size(); i++){
	lattice_[i].resize(size[1]);

	for (unsigned int j=0; j<lattice_[i].size(); j++){
	    lattice_[i][j].resize(size[2],nullptr);
	}
    }

    size_=size;
    
    CLOG(TRACE, logName_) << "lattice created successfully";
}

//------------------------------------------------------------------------------

//! The lattice is extended homogeniously in one dimension. With this method it
//! is not possible to grow only certain parts of the lattice.
//! \param dimension The desired growth dimension. Must be either 0, 1 or 2.
//! \param amount Number of layers that shall be added.
//! \throws LatticeException if a dimension >2 is provided and if lattice shall
//! be created and size in at least one dimension is zero
void Lattice::extend(indexType amount, uint8_t dimension)
{

    CLOG(TRACE, logName_) << "Extending lattice in dimension " << dimension
			  << " by " << amount;

    // Lattice is a three dimensional vector so if first dimension is zero then
    // complete lattice is empty.
    if (size_[0] == 0){

	CLOG(ERROR, logName_) << "Empty lattice can not be extended.";
	throw LatticeException("Empty lattice can not be extended.",
			       LatticeException::Id::NotExists);
	    
    }

    // extend lattice
    switch(dimension){
    case(0):
	lattice_.resize(size_[0]+amount);
	for (int i=0; i<amount; i++){
	    lattice_[size_[0]+i].resize(size_[1]);

	    for (int j=0; j<size_[2]; j++)
		lattice_[size_[0]+i][j].resize(size_[2], nullptr);
	}
	size_[0]+=amount;	
	break;
    case(1):
	for (int i=0; i<size_[0]; i++){
	    lattice_[i].resize(size_[1]+amount);

	    for (int j=0; j<amount; j++)
		lattice_[i][size_[1]+j].resize(size_[2], nullptr);
	}
	size_[1]+=amount;
	break;
    case(2):
	for (int i=0; i<size_[0]; i++){
	    for (int j=0; j<size_[1]; j++)
		lattice_[i][j].resize(size_[2]+amount);
	}
	size_[2]+=amount;
	break;
    default:
	CLOG(ERROR, logName_) << "Wrong or no dimension provided.";
	throw LatticeException("Wrong or no dimension provided.",
			       LatticeException::Id::WrongDimension);
    }

    CLOG(TRACE, logName_) << "new size of lattice is" << size_.str();
    
}

//------------------------------------------------------------------------------

// Introduces either a single cation or anion layer (depending on the provided
// layer type) at the given layer number in growth direction according to the
// zincblende lattice structure.
// \param i Atom layer index of the newly introduced layer in growth direction.
// \param position Constains the position in space of the layer to be
// introduced.
// \param material Pointer to the material that shall be used.
// \param layerType Defines which type of layer (cation or anion, regular or
// shifted) shall be created. Only supports CationRegular, CationShifted,
// AnionRegular and AnionShifted.
void Lattice::createZincblendeLayer(const indexType i,
                                    Vector3D<spaceType> position,
                                    const Material *material,
                                    const Lattice::Layer layerType,
                                    uint8_t outOfPlaneDimension,
                                    double inPlaneLatticeConstant)
{

    uint8_t inPlaneDimension1 = (outOfPlaneDimension+1)%3;
    uint8_t inPlaneDimension2 = (outOfPlaneDimension+2)%3;
    double atomPitchInPlane = inPlaneLatticeConstant /4;
    elementType elementId;
    indexType j, k;
    Vector3D<indexType> index;

    bool isCationLayer = (layerType == Lattice::Layer::CationRegular);
    isCationLayer = isCationLayer || (layerType == Lattice::Layer::CationShifted);

    bool isShifted = (layerType == Lattice::Layer::CationShifted);
    isShifted = isShifted || (layerType == Lattice::Layer::AnionShifted);
    
    CLOG(DEBUG, logName_) << "introducing zincblende structure at layer  " << i
                          << " with type " << layerType;

    // create zincblende lattice
    //   x   i   x     o ... 1st layer
    // o   n   o       x ... 2nd layer
    //   i   x   i     n ... 3rd layer
    // n   o   n       i ... 4th layer
    //   x   i   x
    // o   n   o
   
    j = isCationLayer?0:1;
    index[outOfPlaneDimension]=i;
    
    for(; j<size_[inPlaneDimension1]; j+=2) {

        position[inPlaneDimension1] = j*atomPitchInPlane;
        index[inPlaneDimension1]=j;

        k= isShifted? (2+j)%4 : j%4 ;
           
        for(; k<size_[inPlaneDimension2]; k+=4) {

            position[inPlaneDimension2] = k*atomPitchInPlane;
            index[inPlaneDimension2]=k;

            if (isCationLayer)
                elementId = material->getRandomCation(random_);
            else
                elementId = material->getRandomAnion(random_);
                    
            this->operator()(index) = new Atom(elementId, position, index,
                                       material);
            CLOG(TRACE, logName_) << "added element " << elementId <<
                " at " << index.str() << " position " << position.str();

            // should not happen at all, BUT:
            // checking pos[2] first because this is the most likely case
            if ((std::isnan(position[2])) || (std::isnan(position[1])) || (std::isnan(position[0]))){
                throw LatticeException("Coordinate in lattice is NAN",
                               LatticeException::Id::PositionNAN);
            }


        }
    }

    CLOG(DEBUG, logName_) << "zincblende structure at layer " << i
                         << " successfully finished.";
}

//------------------------------------------------------------------------------

// Introduces a rocksalt lattice at the given layer number in growth direction.
// If only one element is specified of cations respectively anions the complete
// layer is created from this single element. If two or more were specified the
// first two are used to create an ordered one.
// \param i Atom layer index of the newly introduced layer in growth direction.
// \param position Constains the position in space of the layer to be
// introduced.
// \param material Pointer to the material that shall be used.
// \param layerType Cation* respectively Anion* define which elements shall be
// used. *Regular* means that lattice is started at (0,0), *Shifted* at
// (1,1). Finally *Switched means that the element at (0,0) respectively (1,1)
// is of the type that was defined second (if available), not first as in the
// other cases.
void Lattice::createRocksaltLayer(const indexType i,
                                  Vector3D<spaceType> position,
                                  const Material *material,
                                  const Lattice::Layer layerType,
                                  uint8_t outOfPlaneDimension,
                                  double inPlaneLatticeConstant)
{

    uint8_t inPlaneDimension1 = (outOfPlaneDimension+1)%3;
    uint8_t inPlaneDimension2 = (outOfPlaneDimension+2)%3;
    double atomPitchInPlane = inPlaneLatticeConstant /4;
    elementType elementId;
    indexType j, k;
    Vector3D<indexType> index;
    const MaterialComponents *component;
    int startElementIndex, elementIndex;
    
    bool isCationLayer = (layerType == Lattice::Layer::CationRegular)
        || (layerType == Lattice::Layer::CationRegularSwitched)
        || (layerType == Lattice::Layer::CationShifted)
        || (layerType == Lattice::Layer::CationShiftedSwitched);

    bool isShifted = (layerType == Lattice::Layer::CationShifted)
        || (layerType == Lattice::Layer::CationShiftedSwitched)
        || (layerType == Lattice::Layer::AnionShifted)
        || (layerType == Lattice::Layer::AnionShiftedSwitched);

    bool isSwitched = (layerType == Lattice::Layer::CationRegularSwitched)
        || (layerType == Lattice::Layer::CationShiftedSwitched)
        || (layerType == Lattice::Layer::AnionRegularSwitched)
        || (layerType == Lattice::Layer::AnionShiftedSwitched);
        
    if(isCationLayer)
        component = &(material->getCations());
    else
        component = &(material->getAnions());        
    
    CLOG(DEBUG, logName_) << "introducing rocksalt structure at layer  " << i
                          << " with type " << layerType;

    // create rocksalt lattice
    //                 o ... 1st element
    // o   x   o       x ... 2nd element
    //                 
    // x   o   x      
    //             
    // o   x   o
   
    j = isShifted? 1: 0;
    startElementIndex = isSwitched? 1 : 0;
    index[outOfPlaneDimension]=i;
    
    for(; j<size_[inPlaneDimension1]; j+=2) {

        position[inPlaneDimension1] = j*atomPitchInPlane;
        index[inPlaneDimension1]=j;

        k=isShifted ? 1 : 0;

        elementIndex =(startElementIndex + j/2)%2;
        for(; k<size_[inPlaneDimension2]; k+=2) {

            position[inPlaneDimension2] = k*atomPitchInPlane;
            index[inPlaneDimension2]=k;

            if (component->size()==1)
                elementId = std::get<0>((*component)[0]);
            else
                elementId = std::get<0>((*component)[elementIndex]);
                    
            this->operator()(index) = new Atom(elementId, position, index,
                                       material);
            CLOG(TRACE, logName_) << "added element " << elementId <<
                " at " << index.str();

            elementIndex = (elementIndex+1)%2;
        }
    }

    CLOG(DEBUG, logName_) << "rocksalt structure at layer " << i
                         << " successfully finished.";
}

//##############################################################################

//! Picks an atom in the specified layer. The atom is searched near the
//! coordinates (0,0). The first valid location found is delivered.
//! \param layerNumber Layer index that shall be searched.
//! \param outOfPlaneDimension Dimension the layer is fixed.
//! \return Pointer to picked atom.
const Atom* Lattice::getFirstAtomInLayer(const indexType layerNumber,
				    const uint8_t outOfPlaneDimension) const
{
    uint8_t inPlaneDimension1 = (outOfPlaneDimension+1)%3;
    uint8_t inPlaneDimension2 = (outOfPlaneDimension+2)%3;
    Vector3D<indexType> vector;

    CLOG(TRACE, logName_) << "searching for first atom in layer " << layerNumber
                          << " in dimension " << (int) outOfPlaneDimension;
    
    vector[outOfPlaneDimension] = layerNumber;
    
    for (int i=0; i< size_[inPlaneDimension1]; i++){
	vector[inPlaneDimension1] = i;
	for (int j=0; j< size_[inPlaneDimension2]; j++){
	    vector[inPlaneDimension2] = j;
	    auto atom = this->operator()(vector);
	    if (atom != nullptr)
		return atom;
	}
    }

    return ( nullptr );
}

//------------------------------------------------------------------------------

//! \throws LatticeException if defined position exceeds internal lattice.
//! \param x,y,z Position in first, second and third dimension of desired
//! lattice element.
//! \return Reference to specified lattice location.
Atom* & Lattice::operator()(indexType x, indexType y, indexType z)
{
    return operator()(Vector3D<indexType>(x,y,z));
}

//------------------------------------------------------------------------------

//! \throws LatticeException if defined position exceeds internal lattice.
//! \param position Position in  of desired lattice element.
//! \return Reference to specified lattice location.
Atom* & Lattice::operator()(Vector3D<indexType> position)
{
    if (! (position <= size_)){
	CLOG(ERROR, logName_) << "requested position " << position.str() <<
	    " bigger than lattice " << size_.str();
	throw LatticeException("searched position bigger than lattice",
			       LatticeException::Id::PositionTooBig);
    }

    return lattice_[position[0]][position[1]][position[2]];

    // use const version of function to avoid code duplication
    //return const_cast<Atom* &>( static_cast<const Lattice&>
    //			  (*this).operator()(position) );
}

//------------------------------------------------------------------------------

//! \throws LatticeException if defined position exceeds internal lattice.
//! \param x,y,z Position in first, second and third dimension of desired
//! lattice element.
//! \return Reference to specified lattice location.
const Atom * Lattice::operator()(indexType x, indexType y, indexType z) const
{
    return operator()(Vector3D<indexType>(x,y,z));
}

//------------------------------------------------------------------------------

//! \throws LatticeException if defined position exceeds internal lattice.
//! \param position Position in  of desired lattice element.
//! \return Reference to specified lattice location.
const Atom * Lattice::operator()(Vector3D<indexType> position) const
{
    
    if (! (position <= size_)){
	CLOG(ERROR, logName_) << "requested position " << position.str() <<
	    " bigger than lattice " << size_.str();
	throw LatticeException("searched position bigger than lattice",
			       LatticeException::Id::PositionTooBig);
    }

    return lattice_[position[0]][position[1]][position[2]];
}

//------------------------------------------------------------------------------

//! \return Amount of atomic layers in each spatial direction.
Vector3D<indexType> Lattice::getSize(void) const
{
    return size_;
}

//------------------------------------------------------------------------------

//! \return Maximum coordinate of an atom in each spatial direction.
Vector3D<spaceType> Lattice::getMaxCoordinates(void) const
{
    Vector3D<spaceType> result = {0.0, 0.0, 0.0};
    Vector3D<spaceType> tmp;
   
    for (auto& plane: lattice_) {
        for (auto& row : plane) {
            for (auto& atom : row) {
                // safety check
                if (atom != nullptr )
                {
                    tmp = atom->getPosition();

                    for (auto i: {0, 1, 2}) {
                        if (tmp[i] >= result[i])
                            result[i] = tmp[i];
                    }
                } 
            }
        }
    }

    return result;
}

//------------------------------------------------------------------------------

//! \return Temperature of lattice.
double Lattice::getTemperature(void) const
{
    return temperature_;
}

//------------------------------------------------------------------------------

//! \param temperature New lattice temperature
void Lattice::setTemperature(const double temperature)
{
    temperature_ = temperature;
}

//------------------------------------------------------------------------------
					 
//! Generates a list of all atoms inside the specified range. This is necessary
//! for example in the relaxation algorithms. The order how the atoms are added
//! is not relevant in these cases.
//! \param range Search space for which valid atoms shall be delivered.
//! \param checkModification If true only atoms are added that have recently
//! been modified.
//! \param materialName Only atoms belonging to this material are selected. Is
//! ignored if contains "".
//! \return List of atoms inside specified range.
//! \throws LatticeException if start position exceeds internal lattice.
//! \sa SimulationBox
std::vector<Atom *> Lattice::getAtomList(Range3D<indexType> range,
					 const bool checkModification,
					 const std::string &materialName) const
{

    std::vector<Atom *> atoms;
    indexType i, j, k;
    bool addAtom = false;

    CLOG(TRACE, logName_) << "starting generation of atom list in range " <<
	range.start.str() << " to " << range.stop.str() << " with " <<
	range.apply.str();
    
    for (int i=0; i<3; i++){
	// if no range specified take whole range
	if (range.apply[i] == false){
	    range.start[i] = 0;
	    range.stop[i] = size_[i]-1;
	}
	else{
	    if (range.start[i] >= size_[i])
		throw LatticeException("Start of range bigger than lattice.",
				       LatticeException::Id::RangeTooBig);

	    if (range.stop[i] >= size_[i]){
		CLOG(WARNING, logName_) << "getAtomList: Stop of range bigger than lattice! "
					<< "Setting to maximum!";
		range.stop[i] = (size_[i]-1);
	    }
	}
    }

    for(i = range.start[0]; true; checkOverflow(++i, size_[0]) ){

	for (j = range.start[1]; true; checkOverflow(++j,size_[1])){

	    for (k = range.start[2]; true; checkOverflow(++k, size_[2])){
	
		if (lattice_[i][j][k] != nullptr ){
		    addAtom = true;
		    
		    if (checkModification &&
			! (lattice_[i][j][k])->wasModified()){
			addAtom=false;
		    }

		    if ( (materialName != "") &&
			 (lattice_[i][j][k])->getMaterial()->getName() !=
			 materialName){
			addAtom=false;
		    }

		    if (addAtom) atoms.push_back(lattice_[i][j][k]);
		}

		if (k == range.stop[2]) break;
	    }
	    
	    if (j == range.stop[1]) break;
	}

	if (i == range.stop[0])  break;
    }

    return atoms;
}

//------------------------------------------------------------------------------

//! Get the list of atoms in a given layer.
//! \param layerID index of the layer to be returned
//! \param dimension normal to the layer plane. (0 = (100), 1 = (010), 2 = (001))
//! \return List of atoms in specified layer
std::vector<Atom *> Lattice::getAtomsInLayer(
                indexType layerID, indexType dimension) const
{
    Vector3D<indexType> start{};
    Vector3D<indexType> end{};
    Vector3D<bool> apply{false};

    if (dimension < 3) {
        start[dimension] = layerID;
        end[dimension] = layerID;
        apply[dimension] = true;
    }

    Range3D<indexType> layerDef = Range3D<indexType>(start, end, apply);

    return getAtomList(layerDef, false, "");
}

//------------------------------------------------------------------------------

//! Count atoms according to their state
//! \param metric AtomState to check for
//! \return Number of atoms meeting condition
uint32_t Lattice::countAtomsByState(AtomState metric) {
    uint32_t result{};
    std::bitset<AtomState::StateCount> tmpState;
    tmpState.set(metric);

    for (auto atom: getAtomList()) {

        if ((tmpState & atom->getState()).any())
            result ++;
    }

    return result;
}

//------------------------------------------------------------------------------

//! Count atoms within a layer according to their state
//! \param layerID index of the layer to be returned
//! \param dimension normal to the layer plane. (0 = (100), 1 = (010), 2 = (001))
//! \param metric AtomState to check for
//! \return Number of atoms meeting condition in specified layer
uint32_t Lattice::countAtomsInLayerByState(indexType layerID, indexType dimension,
                AtomState metric) {
    uint32_t result{};

    std::vector<Atom *>atomList = getAtomsInLayer(layerID, dimension);
    std::bitset<AtomState::StateCount> tmpState;
    tmpState.set(metric);

    for (auto atom: atomList) {
        tmpState.reset();
        tmpState.set(metric);
        tmpState &= atom->getState();

        if (tmpState[metric])
            result ++;
    }

    return result;
}

//------------------------------------------------------------------------------

//! Check if value reached limit and if it did reset it to zero.
inline void Lattice::checkOverflow(indexType &value, const indexType limit)
    const
{
    if (value == limit)
	value=0;
}

//------------------------------------------------------------------------------

//! Copy all atoms of a specific range into a vector and return the result.
//! \param range Specific range whose atoms shall be backuped.
//! \return Vector of saved atoms.
std::vector<Atom> Lattice::backupAtoms(const Range3D<indexType> range) const
{
    std::vector<Atom> result;

    CLOG(TRACE, logName_) << "Starting to backup atoms in range "<< range.str();

    for (auto atom: getAtomList(range))
	result.push_back(*atom);
    
    CLOG(TRACE, logName_) << "Successfully backuped " << result.size()
			 << " atoms";
    
    return result;
}

//------------------------------------------------------------------------------

//! Restores Atoms from the given vector into the lattice.
//! \param atoms Vector of atoms that shall be restored.
void Lattice::restoreAtoms(const std::vector<Atom> atoms)
{

    CLOG(TRACE, logName_) << "Starting to restore saved atoms.";
    
    for (auto atom: atoms){
	*(this->operator()(atom.getIndex())) = atom;
    }

    CLOG(TRACE, logName_) << "Successfully restored " << atoms.size()
			 << " atoms.";
    
}

//------------------------------------------------------------------------------

//! Returns the number of valid atoms in the lattice.
//! \remark Calls Lattice::getAtomList
//! \return Number of atoms.
uint16_t Lattice::getAtomCount(Range3D<indexType> range) const
{
    return getAtomList(range).size();
}

//------------------------------------------------------------------------------

//! Get average coordinate
//! \param layer index of layer to be analyzed
//! \param dimension defines orientation of layer plane
//! \param checkModification if true, only modified atoms are considered
//! \return Average coordinate in given dimension
double Lattice::getAvgCoordinate (indexType layer, uint8_t dimension, bool checkModification) const {
    double result = 0.0;

    if (dimension < 3) {
        // code duplication for easier understanding
        // allows modification for loose layer definitions (e.g. cation + anion 
        // layer) in future
        Vector3D<indexType> rangeStart;
        Vector3D<indexType> rangeStop;

        rangeStart[dimension] = layer;
        rangeStop[dimension] = layer;

        Vector3D<bool> apply = Vector3D<bool>(false, false, false);
        apply[dimension] = true;

        Range3D<indexType> layerDef = Range3D<indexType>(rangeStart, rangeStop, apply);

        std::vector<Atom *> atoms = getAtomList(layerDef, false, "");

        uint32_t processedAtoms = 0;

        for (auto atom: atoms) {
            if ((checkModification == true) && (atom->wasModified() == false))
                continue;

            result += atom->getPosition()[dimension];
            processedAtoms ++;
        }

        // dividing by counter rather than atoms.size in case
        // unmodified atoms were filtered out
        if (processedAtoms > 0)
            result /= processedAtoms;
    } else {
        throw LatticeException("Wrong or no dimension provided.",
                       LatticeException::Id::WrongDimension);
    }

    return result;
}

//------------------------------------------------------------------------------

//! Returns the maximum modification index of all atoms in object
//! \return maximum modification index
uint32_t Lattice::getMaxModificationIndex(void) const {
    uint32_t result = 0;

    std::vector<Atom *> atomList = getAtomList();

    for (auto atom: atomList) {
        if (atom->getModificationOrder() > result)
            result = atom->getModificationOrder();
    }

    return result;
}

//##############################################################################

LatticeException::LatticeException(const std::string &message,
				   const LatticeException::Id id)
    :message_(message), id_(id) {}

//------------------------------------------------------------------------------

//! Returns a short description of the exception. This is the implementation of
//! the virtual function inherited from std::exception
//! \remark This function does not throw exceptions!
const char * LatticeException::what () const throw ()
{
    static std::string text;
    text = "Lattice Exception -- " + message_;
    return text.c_str();
}

//------------------------------------------------------------------------------

//! Returns a more detailled description of the exception, making it possible to
//! use a single exception class for multiple errors.
//! \remark This function does not throw exceptions!
const std::string LatticeException::getMessage () const noexcept
{
    return message_;
}

//------------------------------------------------------------------------------

//! Returns the ID of the exception source, making possible to easily
//! distinguish different errors.
//! \remark This function does not throw exceptions!
const LatticeException::Id LatticeException::getId () const noexcept
{
    return id_;
}

//##############################################################################

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
