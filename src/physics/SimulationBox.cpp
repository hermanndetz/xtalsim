/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#include "SimulationBox.h"

//------------------------------------------------------------------------------

//! The constructor registers initializes the logger.
//! \param logName Optional parameter defining the logger's name.
//! \param outOfPlaneDimension One of the three spatial directions the crystal
//! is supposed to grow.
//! \param size Size of the Simulation Box in space.
//! \param latticeTemperature Temperature of the internal lattice.
SimulationBox::SimulationBox(const uint8_t outOfPlaneDimension,
                             const Vector3D<spaceType> size,
                             const double inPlaneLatticeConstant,
                             const double latticeTemperature,
                             const char *logName):
    description_(""),
    lattice_(latticeTemperature),
    outOfPlaneDimension_(outOfPlaneDimension),
    inPlaneLatticeConstant_(inPlaneLatticeConstant), size_(size),
    updateNeighborList_(false), neighborListLayers_(1),
    neighborListRadius_(0), neighborListPeriodicBoundaries_(true),
    //neighborListThreadRunning_(false),
    logName_(logName)
{
    //    random_.seed(std::chrono::system_clock::now().time_since_epoch().count());
    random_.seed(std::random_device()());
    startConcurrentNeighborListGeneration();
}

//------------------------------------------------------------------------------

SimulationBox::~SimulationBox() {}

//------------------------------------------------------------------------------

//! \return Description of Simulation Box.
std::string SimulationBox::getDescription(void) const
{
    return description_;
}

//------------------------------------------------------------------------------

//! \param description New Simulation Box description
void SimulationBox::setDescription(const std::string description)
{
    description_ = description;
}

//------------------------------------------------------------------------------

//! \return Temperature of lattice.
double SimulationBox::getLatticeTemperature(void) const
{
    return lattice_.getTemperature();
}

//------------------------------------------------------------------------------

//! \param temperature New lattice temperature
void SimulationBox::setLatticeTemperature(const double temperature)
{
    lattice_.setTemperature(temperature);
}

//------------------------------------------------------------------------------

uint8_t SimulationBox::getOutOfPlaneDimension(void) const
{
    return outOfPlaneDimension_;
}

//------------------------------------------------------------------------------

double SimulationBox::getInPlaneLatticeConstant(void) const
{
    return inPlaneLatticeConstant_;
}

//------------------------------------------------------------------------------

Vector3D<spaceType> SimulationBox::getSize(void) const
{
    return size_;
}

//------------------------------------------------------------------------------

Lattice & SimulationBox::getLattice(void)
{
    return lattice_;
}

//------------------------------------------------------------------------------

Lattice const & SimulationBox::getLattice(void) const
{
    return lattice_;
}

//------------------------------------------------------------------------------

MaterialCollection & SimulationBox::getMaterials(void)
{
    return materials_;
}

//------------------------------------------------------------------------------

MaterialCollection const & SimulationBox::getMaterials(void) const
{
    return materials_;
}

//------------------------------------------------------------------------------

//! Adds a material to the collection stored in simbox
//! Needed because LatticeGenerator only stores materials that are generated,
//! BUT the InterfaceGenerator can create additional combinations.
//! \param material Material that shall be added to collection.
void SimulationBox::addMaterial(const Material &material) {
    materials_.add(material);
}

//------------------------------------------------------------------------------

//! Returns the number of valid atoms.
//! \return Number of atoms.
uint32_t SimulationBox::getAtomCount(const Range3D<indexType>&range) const
{
    return lattice_.getAtomCount(range);
}

//------------------------------------------------------------------------------

//! \param range Range inside the lattice whose modification attribute shall be
//! cleared.
void SimulationBox::clearModification(const Range3D<indexType> &range)
{
    std::vector<Atom *> atoms = lattice_.getAtomList(range);

    for (auto atom: atoms){
	atom->clearModification();
    }
}

//##############################################################################

//! Create lattice in defined size and introduce atoms, such that the form a
//! zincblende lattice. Please note that the in plane parts are only considered
//! when the lattice is created, later only the out of plane component is
//! considered. Wrapper function for \sa createLattice.
//! \param layerCounts Desired size.
//! \param material Defines the Material that is introduced.
void SimulationBox::createZincblende(const Vector3D<indexType> layerCounts,
				     const Material &material)
{

    CLOG(TRACE, logName_) << "starting creation of " << layerCounts.str()
			  << "layers of zincblende lattice for material "
			  << material.getName();
    
    createLattice(layerCounts, material, Lattice::Type::Zincblende);

    CLOG(TRACE, logName_) << "creation of Zincblende lattice successfully "
                              << "finished.";
}

//------------------------------------------------------------------------------

//! Create half Heusler lattice in defined size. The first cation element is
//! considered to form the ZnS type bond with the anion element.
//! Please note that the in plane parts are only considered
//! when the lattice is created, later only the out of plane component is
//! considered. Wrapper function for \sa createLattice.
//! \param layerCounts Desired size.
//! \param material Defines the Material that is introduced. At least two cation
//! and one anion element have to be defined.
void SimulationBox::createHalfHeusler(const Vector3D<indexType> layerCounts,
				     const Material &material)
{

    CLOG(TRACE, logName_) << "starting creation of " << layerCounts.str()
                          << "layers of half Heusler lattice for material "
                          << material.getName();

    if (material.getCations().size() < 2) {
        CLOG(ERROR, logName_) << "At least two cations have to be "
                              << "specified! ABORTING";
        return;
    }
    
    createLattice(layerCounts, material, Lattice::Type::HalfHeusler);

    CLOG(TRACE, logName_) << "creation of Half Heusler lattice successfully "
                          << "finished.";
    
}

//------------------------------------------------------------------------------

//! Create full Heusler lattice in defined size. The first anion element is
//! considered to sit at (1,1,1) when starting at (0,0,0).
//! Please note that the in plane parts are only considered
//! when the lattice is created, later only the out of plane component is
//! considered.
//! \param layerCounts Desired size.
//! \param material Defines the Material that is introduced. At least two cation
//! and one anion element have to be defined.
void SimulationBox::createFullHeusler(const Vector3D<indexType> layerCounts,
				     const Material &material)
{

    CLOG(TRACE, logName_) << "starting creation of " << layerCounts.str()
			  << "layers of half Heusler lattice for material "
			  << material.getName();

    if (material.getAnions().size() < 2) {
        CLOG(ERROR, logName_) << "At least two cations have to be "
                              << "specified! ABORTING";
        return;
    }
    
    createLattice(layerCounts, material, Lattice::Type::FullHeusler);

    CLOG(TRACE, logName_) << "creation of Full Heusler lattice successfully "
                          << "finished.";
}
    
//------------------------------------------------------------------------------

//! Create lattice in defined size and introduce atoms, such that the form a
//! zincblende lattice. Please note that the in plane parts are only considered
//! when the lattice is created, later only the out of plane component is
//! considered.
//! \param layerCounts Desired size.
//! \param material Defines the Material that is introduced.
void SimulationBox::createLattice(const Vector3D<indexType> layerCounts,
                                  const Material &material,
                                  const Lattice::Type latticeType)
{
    uint8_t inPlaneDimension1 = (outOfPlaneDimension_+1)%3;
    uint8_t inPlaneDimension2 = (outOfPlaneDimension_+2)%3;

    int i;
    double atomPitchInPlane, atomPitchOutOfPlane;
	
    Vector3D<indexType> size;
    Vector3D<spaceType> position;
    const Material *materialPointer;

    std::vector<Lattice::Layer> layerTypes;
    
    switch (latticeType) {
    case Lattice::Type::Zincblende:
        layerTypes.push_back(Lattice::Layer::CationRegular);
        layerTypes.push_back(Lattice::Layer::AnionRegular);
        layerTypes.push_back(Lattice::Layer::CationShifted);
        layerTypes.push_back(Lattice::Layer::AnionShifted);
        break;
    case Lattice::Type::HalfHeusler:       
        layerTypes.push_back(Lattice::Layer::CationRegular);
        layerTypes.push_back(Lattice::Layer::CationRegularSwitched);
        break;
    case Lattice::Type::FullHeusler:
        layerTypes.push_back(Lattice::Layer::CationRegular);
        layerTypes.push_back(Lattice::Layer::AnionShifted);
        layerTypes.push_back(Lattice::Layer::CationRegular);
        layerTypes.push_back(Lattice::Layer::AnionShiftedSwitched);
        break;
    }

    if (inPlaneLatticeConstant_ < 0.0) {
        inPlaneLatticeConstant_ = material.getLatticeConstant();
        atomPitchOutOfPlane = inPlaneLatticeConstant_ /4.0;
    } else {
        atomPitchOutOfPlane =
            material.getLatticeConstant(inPlaneLatticeConstant_)/4.0;
    }
	    
    atomPitchInPlane = inPlaneLatticeConstant_ /4.0;

    CLOG(DEBUG, logName_) << "atom Pitch: " << atomPitchOutOfPlane;

    lattice_.generate(layerCounts, outOfPlaneDimension_);
    size = lattice_.getSize();

    // If lattice was empty before creation.
    if ( !((size_[inPlaneDimension1]+size_[inPlaneDimension2]) > 0)) {
	
        size_[inPlaneDimension1] = size[inPlaneDimension1]* atomPitchInPlane;
        size_[inPlaneDimension2] = size[inPlaneDimension2]* atomPitchInPlane;
    }
       
    i=size[outOfPlaneDimension_]- layerCounts[outOfPlaneDimension_];
    position[outOfPlaneDimension_] = size_[outOfPlaneDimension_];

    CLOG(DEBUG, logName_) << "starting to introduce atoms on overall lattice "
			  << "size " << size.str() << " beginning at layer "
                          << i << " in growth dimension "<< (int)outOfPlaneDimension_;

    // To initialize the Atom object a pointer to the current material is
    // required which has to stay valid until the SimulationBox object is
    // destroyed. Therefore the material is added to the internal
    // MaterialCollection which gives the material a constant address.
    materialPointer = & materials_.getByNameOrAdd(material);
    
    int layerTypesIndex =0;

    for(;i<size[outOfPlaneDimension_]; i+=1) {

        switch(latticeType){
        case Lattice::Type::Zincblende:
            lattice_.createZincblendeLayer(i, position, materialPointer,
                                           layerTypes[layerTypesIndex],
                                           outOfPlaneDimension_,
                                           inPlaneLatticeConstant_);
            layerTypesIndex = (layerTypesIndex+1)%layerTypes.size();
        break;
        case Lattice::Type::HalfHeusler:
            if ((i%2) == 0){
                // rocksalt part
                lattice_.createRocksaltLayer(i, position, materialPointer,
                                             layerTypes[layerTypesIndex],
                                             outOfPlaneDimension_,
                                             inPlaneLatticeConstant_);
                layerTypesIndex = (layerTypesIndex+1)%layerTypes.size();
            }
            else{
                
                // zincblende part
                lattice_.createZincblendeLayer(i, position, materialPointer,
                                               Lattice::Layer::AnionRegular,
                                               outOfPlaneDimension_,
                                               inPlaneLatticeConstant_);
            }
            break;
        case Lattice::Type::FullHeusler:
            lattice_.createRocksaltLayer(i, position, materialPointer,
                                         layerTypes[layerTypesIndex],
                                         outOfPlaneDimension_,
                                         inPlaneLatticeConstant_);
            layerTypesIndex = (layerTypesIndex+1)%layerTypes.size();
            break;
        }
            
        position[outOfPlaneDimension_] += atomPitchOutOfPlane;
    }
    
    size_[outOfPlaneDimension_] = position[outOfPlaneDimension_];
    updateNeighborList_ = true;
    updateNeighbors();

    //    if (neighborListThreadRunning_ == false)
    //        startConcurrentNeighborListGeneration();
        
    CLOG(DEBUG, logName_) << layerCounts[outOfPlaneDimension_] <<
        " atomic layers have been created.";
    
}

//##############################################################################

void SimulationBox::startConcurrentNeighborListGeneration(void){

    //    neighborListThreadRunning_ = true;
    
    // initial neighbor list creation
    //updateNeighbors();
    
    auto thread = std::thread(&SimulationBox::generateNeighborsConcurrent, this);
    thread.detach();
}

//------------------------------------------------------------------------------

//! This method is supposed to be started in a separate thread to continuously
//! update the neighbor lists.
void SimulationBox::generateNeighborsConcurrent(void)
{
    
    while (true) {
        //CLOG(ERROR, logName_) << "concurrent Neighbor list";
        std::this_thread::sleep_for(std::chrono::milliseconds(400));

        // Update of list only necessary when radius used. Otherwise only
        // layer distance is important which does not change over time
        if (neighborListRadius_ > 0)
            updateNeighbors();
    }
    
}

//------------------------------------------------------------------------------
    
//! \param force If true neighbor list is created from scratch.
void SimulationBox::generateNeighbors(bool force)
{
    if (force)
        updateNeighborList_ = true;

    updateNeighbors();
}
 
//------------------------------------------------------------------------------
    
//! Adds atoms that are inside specified layers and radius around one atom to
//! its neighbor list. The list is only recreated/updated if atoms have been
//! moved/added. Otherwise this function does nothing.
//! \param radius Maximal distance between atoms such that they are still
//! counted as neighbors. Ignored if equal to 0.
//! \param layers Maximal allowed index distance such that they are still
//! counted as neighbors. Ignored if equal to 0.
//! \param periodicBoundaries Enables or disables periodic boundaries.
void SimulationBox::generateNeighbors(indexType layers,
				      const double radius, bool periodicBoundaries)

{

    bool parametersChanged = false;
    
    if (layers != neighborListLayers_) {
        neighborListLayers_ = layers;
        parametersChanged = true;
    }

    if (radius != neighborListRadius_) {
        neighborListRadius_ = radius;
        parametersChanged = true;
    }

    if (periodicBoundaries != neighborListPeriodicBoundaries_) {
        neighborListPeriodicBoundaries_ = periodicBoundaries;
        parametersChanged = true;
    }

    // Update not necessary if parameters did not change and only atomic layer
    // distance used.
    if ( (! parametersChanged) && (! (neighborListRadius_ > 0)) )
        return;

    updateNeighborList_ = true;
    updateNeighbors();

}

//------------------------------------------------------------------------------

//! Generates an updated neighbor list for each atom by creating an updated list
//! in parallel and afterwards swapping it in one sweep with the old one.
//! The list is only recreated/updated if atoms have been
//! moved/added. Otherwise this function does nothing.
//! \todo Periodic boundaries can only be switched off, if cut-off is given 
//! in layers but not as radius. JM: can not confirm this, revisit this!!
void SimulationBox::updateNeighbors(void)
{

    // radius is squared to prevent costly sqrt()
    double checkValue = neighborListRadius_*neighborListRadius_;
    std::vector<Atom *> atoms, neighbors;
    
    if (! updateNeighborList_) {
        CLOG(DEBUG, logName_) << "no update of neighbor list required "
                             << "because nothing changed";
        return;
    }
        
    CLOG(TRACE, logName_) << "Creating neighbor list using layer count " <<
	neighborListLayers_ << " and radius " << neighborListRadius_;

    // if layers == x then lattice has at least to be of size 2*x+1
    if ( lattice_.getSize() < (2*neighborListLayers_+1) ){
	CLOG(INFO, logName_) << "Layer count " << neighborListLayers_ << " too big for "
	    "lattice of size " << lattice_.getSize().str();
	neighborListLayers_ = 0;
    }

    // if both smaller or equal to 0
    if ( ! ( (neighborListLayers_ > 0) ||  (neighborListRadius_ > 0) ) ) return;
    
    atoms = lattice_.getAtomList();

    for (auto atom: atoms) {

        int counter=0;
        Range3D<indexType> range = getSurrounding(atom->getIndex());
        NeighborList neighborList;
        
        CLOG(DEBUG, logName_) << "Range to look for neighbors: " <<
            range.str();
    
        neighbors = lattice_.getAtomList(range);

        for (auto neighbor: neighbors) {
            if (neighbor->getIndex() == atom->getIndex()) continue;

            if (neighborListRadius_ > 0){
                Vector3D<spaceType> tmpVector;
                if (neighborListPeriodicBoundaries_){
                    tmpVector =	atom->getPosition().getMinimalDistance(
                                                                       neighbor->getPosition(), size_);
                }
                else{
                    tmpVector =	atom->getPosition().getMinimalDistance(
                                                                       neighbor->getPosition(), Vector3D<spaceType>());
                }

                if ( tmpVector.squaredLength() > checkValue) continue;
            }

            neighborList.push_back(neighbor->getIndex());
            counter ++;
        }

        atom->updateNeighbors(neighborList);
        CLOG(DEBUG, logName_) << "new neighbor list for " <<
            atom->getIndex().str() << " has " << counter << " elements";

    }

    updateNeighborList_ = false;
    CLOG(TRACE, logName_) << "finished update of neighbor lists";
}
    
//------------------------------------------------------------------------------

//! Determines for a given position in the lattice the surrounding atoms within
//! a given radius.
//! \warning If no atoms resides on the determined position nothing is done.
//! \param positionLattice Coordinates inside the lattice that shall be checked.
//! \param radius Maximal distance that shall be covered in the range.
//! \return Range on the lattice containing all positions within the radius
//! around the given location.
Range3D<indexType> SimulationBox::getSurrounding(
				const Vector3D<indexType> &positionLattice) const
{
   
    // radius is squared to prevent costly sqrt()
    double checkValue = neighborListRadius_*neighborListRadius_*1.5;
    Range3D<indexType> range;
    Vector3D<indexType> vector;
    Vector3D<indexType> latticeSize = lattice_.getSize();
    Vector3D<spaceType> positionSpace =(lattice_(positionLattice))->getPosition();
    Vector3D<spaceType> tmpVector;

    CLOG(TRACE, logName_) << "generating surrounding for position " <<
		   positionLattice.str();
    CLOG(TRACE, logName_) << "lattice size " << latticeSize.str();
    
    if ((lattice_(positionLattice)) == nullptr) return range;

	if (neighborListLayers_ > 0) {
	    range.start = positionLattice - neighborListLayers_;
	    range.stop = positionLattice + neighborListLayers_;
	    range.apply = {true, true, true};

	    if (neighborListPeriodicBoundaries_ == true) {
            range.fitInBox(lattice_.getSize());
	    } else {
            for (auto i: {0, 1, 2}) {
                if (range.start[i] < 0)
                    range.start[i] = 0;

                if (range.stop[i] >= latticeSize[i])
                    range.stop[i] = latticeSize[i]-1;
            }
            
	    }
        return range;
    }
    
    // determine lower limits for the search space by going in each direction
    // and searching for the first atom outside the radius
    for (int i=0; i<3; i++) {
        vector=positionLattice;

        while(true) {

            vector[i] -= 1;
            if (neighborListPeriodicBoundaries_ == true)
                vector.fitInBox(latticeSize, i);
            else if (vector[i] < 0){
                vector[i] = 0;
                break;
            }

            if ( vector[i] ==  positionLattice[i]) {
                // if same lattice position reached again take whole lattice
                vector[i] = 0;
                break;
            }
	    
            if (lattice_(vector) != nullptr){
                tmpVector = positionSpace.getMinimalDistance(
                         (lattice_(vector))->getPosition(), size_);

                if ( tmpVector.squaredLength() > checkValue)
                    break;
            }
        }

        range.start[i] = vector[i];
    }

    CLOG(TRACE, logName_) << "finished lower limit";
    
    // determine the upper limit in the same fashion
    for (int i=0; i<3; i++) {
        vector=positionLattice;

        while(true) {

            vector[i]+=1;
            if (neighborListPeriodicBoundaries_ == true)
                vector.fitInBox(latticeSize, i);
            else if (vector[i] >= latticeSize[i]){
                vector[i] = latticeSize[i]-1;
                break;
            }
            
            // when ranges collide
            if (vector[i] == range.start[i]) {
                vector[i] -= 1;
                vector.fitInBox(latticeSize, i);
                
                break;
            }
	    
            if (lattice_(vector) != nullptr) {
                tmpVector = positionSpace.getMinimalDistance(
                         (lattice_(vector))->getPosition(), size_);

                if ( tmpVector.squaredLength() > checkValue)
                    break;
            }
            
        }

        range.stop[i] = vector[i];
    }

    CLOG(TRACE, logName_) << "achieved surrounding from " << range.start.str()
			  << " to " << range.stop.str();
    
    range.apply={true,true,true};
    return range;
    
}

//##############################################################################

//! Runs several relaxation steps in parallel on the simulation box. The number
//! of single pieces that shall be worked on has to be even. Since energy
//! calculation needs also atoms in neighboring slices the single pieces are
//! optimized in an interleaved fashion, i.e., in the first run the 1st, 3rd,
//! 5th .. are processed. In the second step the 2nd, 4th, 6th ... For a
//! description of the parameters \sa mmcRelax
double SimulationBox::getEnergyParallel(const TersoffPotential &tersoff,
                                      const Range3D<indexType> range,
                                      const int maxThreadCount) const
{

    int threadCount = maxThreadCount;
    Vector3D<indexType> size = lattice_.getSize();
    indexType stepWidth = (size[outOfPlaneDimension_]+1)/threadCount;
    
    while (stepWidth < 1){
        threadCount --;
        stepWidth = (size[outOfPlaneDimension_]+1)/threadCount;
    }

    std::thread threads[threadCount];
    double localEnergy[threadCount];

    Range3D<indexType> partialRange;
    partialRange.apply[outOfPlaneDimension_] = true;

    CLOG(DEBUG, logName_) << "Each thread serves " << stepWidth << " layers";

    indexType startLayer = 0;
    for (int threadID=0; threadID<threadCount; threadID++){

        partialRange.start[outOfPlaneDimension_] = startLayer;
        partialRange.stop[outOfPlaneDimension_] = startLayer+stepWidth-1;
        CLOG(TRACE, logName_) << "Starting thread number " << threadID;

        threads[threadID] = std::thread(&SimulationBox::getEnergyByReference, this,
                                        std::ref(tersoff),
                                        std::ref(localEnergy[threadID]), partialRange);
        startLayer+=stepWidth;
    }

    double energy = 0;
    for (int threadID=0; threadID<threadCount; threadID++){
        threads[threadID].join();
        energy += localEnergy[threadID];
    }
        
    CLOG(TRACE, logName_) << "All threads joined";

    return energy;

}

//------------------------------------------------------------------------------

//! Calculates the overall energy of all atoms inside the specified range.
//! \param range Range inside the lattice whose energy shall be computed.
//! \param tersoff Parameters necessary to compute the energy.
//! \return Energy of atoms in specified range.
double SimulationBox::getEnergy(const TersoffPotential &tersoff,
				const Range3D<indexType> range) const
{
    double energy;

    getEnergyByReference(tersoff, energy, range);
    return energy;
}

//------------------------------------------------------------------------------

//! Calculates the overall energy of all atoms inside the specified range.
//! \param range Range inside the lattice whose energy shall be computed.
//! \param tersoff Parameters necessary to compute the energy.
//! \return Energy of atoms in specified range.
void SimulationBox::getEnergyByReference(const TersoffPotential &tersoff,
                              double &energy,
				const Range3D<indexType> range) const
{

    CLOG(TRACE, logName_) << "starting to calculate energy in Range " << range;
    
    energy=0;
    std::vector<Atom *> atoms = lattice_.getAtomList(range);
    
    for (auto atom : atoms){
        energy += getEnergySingleAtom(atom, tersoff);
    }

    CLOG(TRACE, logName_) << "finished energy calculations";

}

//------------------------------------------------------------------------------

//! Calculate the energy of a single atom.
//! \param atom Pointer to the Atom whose energy shall be computed.
//! \param tersoff Parameters necessary to compute the energy.
//! \return Energy of specified atom.
double SimulationBox::getEnergySingleAtom(const Atom *atom,
					 const TersoffPotential &tersoff) const
{

    NeighborList const &neighborList = atom->getNeighbors();
    double lengthVector1{0.0}, xi{0.0}, energy{0.0};
    Vector3D<spaceType> vector1, vector2;
    Vector3D<spaceType> atomPosition (atom->getPosition());
    Vector3D<indexType> atomIndex (atom->getIndex());

    CLOG(TRACE, logName_) << "Start energy calculation for " << atomIndex.str();
    CLOG(TRACE, logName_) << "Element: " << atom->getElementId();

    if (neighborList.size() == 0) {
        CLOG(WARNING, logName_) << "Neighborlist of atom " << atomIndex.str() <<
            " is empty!";
        return 0.0;
    }
    
    for (auto neighbor1: neighborList) {
        if (neighbor1 == atomIndex) continue;

	CLOG(TRACE, logName_) << "processing neighbor at position " <<
	    neighbor1.str();
        
        const TersoffParameter &param = tersoff.get(atom->getElementId(),
                       (lattice_(neighbor1))->getElementId());
        
        vector1 = atomPosition.getMinimalDistance(
                       (lattice_(neighbor1))->getPosition(), size_);
        lengthVector1 = sqrt(vector1.squaredLength());

        if (param.isAboveCutoff(lengthVector1)){
            CLOG(DEBUG, logName_) << "distance between " << atomIndex.str() <<
            " and " << neighbor1.str() << " outside cutoff radius";
                continue;
        }

        xi = 0.0;
        
        // calculate angular term
        for (auto neighbor2: neighborList){
            if (neighbor1 == neighbor2) continue;
            if (neighbor2 == atomIndex) continue;

            vector2 = atomPosition.getMinimalDistance(
                      (lattice_(neighbor2))->getPosition(), size_);
            xi += param.getAngularCoefficient(vector1, lengthVector1, vector2);
        }

        energy += param.getEnergy(lengthVector1, xi);
	
    }

    CLOG(TRACE, logName_) << "Achieved final energy of " << energy << " for "
			 << atomIndex.str();
    CLOG(TRACE, logName_) << "Energy calculation finished successfully";
    return energy;
}

//------------------------------------------------------------------------------

//! \sa StrainTool::calculateStrain
LayerStrainInfo SimulationBox::calculateStrain( const MaterialCollection &materials,
				    const std::string &outputPreamble,
				    const double refDistance, bool writeToFile )
{

    CLOG(TRACE, logName_) << "Start Strain calculations.";
    
    StrainTool strainTool(*this);
    return strainTool.calculateStrain(materials, outputPreamble,
				      refDistance, writeToFile);
}

//------------------------------------------------------------------------------
//! Calculate strain of a single atom.
double SimulationBox::calculateStrainSingleAtom(const Atom *atom, 
        StrainCalculationModes mode, const MaterialCollection &collection) {
    StrainTool strainTool(*this);
    return strainTool.calculateStrainAtom(materials_, collection, atom, mode);
}

//------------------------------------------------------------------------------

//! \sa StrainTool::getStrainInfo
StrainInfo SimulationBox::getStrainInfo(const MaterialCollection &materials,
					double refDistanceParam)
{

    CLOG(TRACE, logName_) << "Start extraction of Strain Infos.";
    
    StrainTool strainTool(*this);
    return strainTool.getStrainInfo(materials, refDistanceParam);
}

//------------------------------------------------------------------------------

//! Obtain composition info and write it to file.
//! \param outputFileName path, where output file is written (if != "")
//! \return CompositionInfo returns the composition information
CompositionInfo SimulationBox::analyzeComposition(const std::string &outputFileName) const
{
    CompositionInfo compInfo;
    CompositionTool compTool = CompositionTool(*this);

    compInfo = compTool.AnalyzeStructure();

    if (outputFileName != "") {
        CsvHandler compInfoCSV = CsvHandler("CSV File");
        compInfoCSV = compInfo;
        compInfoCSV.save(outputFileName, 1, "\t", "Layer\tElement\tCount");
    }

    return compInfo;
}

//------------------------------------------------------------------------------

//! Generate bond statistics and write them to file.
BondInfo SimulationBox::analyzeBonds(const std::string &outputFileName) const
{
    BondInfo bondInfo;
    BondTool bondTool = BondTool(*this);

    bondInfo = bondTool.AnalyzeStructure();

    if (outputFileName != "") {
        CsvHandler bondInfoCSV = CsvHandler("CSV File");
        bondInfoCSV = bondInfo;
        bondInfoCSV.save(outputFileName, std::vector<int32_t>{1, 2}, "\t", "Layer\tElement1\tElement2\tCount");
    }

    return bondInfo;
}

//------------------------------------------------------------------------------

//! \sa InterfaceTool::createOptimizedRoughness
//! \todo Check, if this is needed at all. Not yet prepared for different metrics!
void SimulationBox::createOptimizedRoughness(
				  Range3D<indexType> interfaceRange,
				  const unsigned int atomsToDeployCount,
				  const unsigned int positionsToCheckCount,
				  const Material &material,
				  const OptimizationParameter &parameter,
				  const TersoffPotential &tersoff,
                  const MaterialCollection &collection,
				  const std::string &journalPreamble)
{

    CLOG(TRACE, logName_) << "Starting creation of optimized interface " <<
	" roughness.";
    
    InterfaceTool ifTool(*this, tersoff);
    ifTool.createOptimizedRoughness(interfaceRange,
				    atomsToDeployCount,
				    positionsToCheckCount,
				    material,
				    parameter,
                    collection,
                    ITM_Energy,
				    journalPreamble);

}

//##############################################################################

//! Runs several relaxation steps in parallel on the simulation box. The number
//! of single pieces that shall be worked on has to be even. Since  energy
//! calculation needs also atoms in neighboring slices the single pieces are
//! optimized in an interleaved fashion, i.e., in the first run the 1st, 3rd,
//! 5th .. are processed. In the second step the 2nd, 4th, 6th ... For a
//! description of the parameters \sa mmcRelax
void SimulationBox::mmcRelaxParallel(const TersoffPotential &tersoff,
			     const OptimizationParameter &parameter,
			     OptimizationStatistic &statistic,
			     Range3D<indexType> range)
{
    int partCount = parameter.maxThreadCount*2;

    Vector3D<indexType> size = lattice_.getSize();
    indexType stepWidth = (size[outOfPlaneDimension_]+1)/partCount;

    if (stepWidth < 2) {
        partCount = int((size[outOfPlaneDimension_]+1)/2);

        if ( (partCount%2) == 1 )
            partCount--;

        CLOG(INFO, logName_) <<  "Width of " << stepWidth << " atomic layers is too little for "
                                << " parallel computation. Decreasing parallel worker count to "
                                << partCount/2;
        stepWidth = (size[outOfPlaneDimension_]+1)/partCount;
    }

    int threadCount= partCount/2;
    std::thread threads[threadCount];
    // separate statistic for each thread, prevents error due to parallel writing
    OptimizationStatistic localStatistics[threadCount];

    Range3D<indexType> partialRange;
    partialRange.apply[outOfPlaneDimension_] = true;

    CLOG(DEBUG, logName_) << "Each thread serves " << stepWidth << " layers";
   
    for (int k=0; k<2; k++) {
        indexType startLayer = k*stepWidth;

        for (int threadID=0; threadID<threadCount; threadID++) {
            localStatistics[threadID].clear();
            partialRange.start[outOfPlaneDimension_] = startLayer;
            partialRange.stop[outOfPlaneDimension_] = startLayer+stepWidth-1;
            CLOG(TRACE, logName_) << "Starting thread number " << threadID;

            threads[threadID] = std::thread(&SimulationBox::mmcRelax, this, std::ref(tersoff),
                                            std::ref(parameter), std::ref(localStatistics[threadID]),
                                            partialRange);
            startLayer+=2*stepWidth;
        }

        for (int threadID=0; threadID<threadCount; threadID++)
            threads[threadID].join();

        for (auto localStatistic: localStatistics)
            statistic.merge(localStatistic);
        
        CLOG(TRACE, logName_) << "All threads joined";
    }
    
}

//------------------------------------------------------------------------------

//! Carries out a single round of Metropolis Monte Carlo Simulation. In this
//! procedure one atom after the other is displaced by a random amount in a
//! random direction. After each displacement it is checked if the energy of the
//! atom got better. If this is the case the displacement is fixed, otherwise it
//! is reverted.
//! \param tersoff Tersoff arameters necessary to compute the energy.
//! \param parameter Optimization parameters such as maximal displacement.
//! \param statistic Optimization statistics.
//! \param range Range of the lattice the relaxation is carried out on.
void SimulationBox::mmcRelax(const TersoffPotential &tersoff,
			     const OptimizationParameter &parameter,
			     OptimizationStatistic &statistic,
			     Range3D<indexType> range)
{

    double energyBefore, energyAfter;
    Vector3D<spaceType> position;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    CLOG(TRACE, logName_) << "Starting " << parameter.mmcRunCount <<
        " Metropolis Monte Carlo relaxation runs.";
    
    std::vector<Atom *> atoms = lattice_.getAtomList(range);

    for (int i=0; i<parameter.mmcRunCount; i++) {
        CLOG(TRACE, logName_) << "Starting run number " << i;
        
        int counter=0;
        shuffle (atoms.begin(), atoms.end(), random_);
        
        for (auto atom: atoms){

            CLOG(TRACE, logName_) << "Starting atom number " << counter;
        
            position = atom->getPosition();
            energyBefore=getEnergySingleAtom(atom, tersoff);
            atom->displace(parameter.mmcMaxDisplacement);
            energyAfter=getEnergySingleAtom(atom, tersoff);

            if (energyAfter > energyBefore) {
                //Casually accept moves which worse energy to leave local minima
                if (distribution(random_) >
                    ArrheniusFactor::Boltzmann(
                            std::abs(energyAfter-energyBefore),
                            lattice_.getTemperature()
                                   )) {
                    statistic.rejectCount ++;
                    atom->setPreviousPosition();
                } else
                    statistic.acceptCount ++;
            } else
                statistic.acceptCount ++;	    

            counter++;
        }
    }
    
    updateNeighborList_ = true;
}

//------------------------------------------------------------------------------

//! The chosen interval is stored directly in the given range.
//! \param range Range specifying the maximal interval.
//! \param dimenstion The dimension that shall be processed.
void SimulationBox::getRandomInterval(Range3D<indexType> &range,
				      const indexType dimension)
{
    Vector3D<indexType> latticeSize=lattice_.getSize();
    Range3D<indexType> energyRange;
    bool isCompleteLattice = ( (range.start[dimension] <= 0) &&
	      (range.stop[dimension] >= latticeSize[dimension]) );
    
    if ( (! range.apply[dimension]) || isCompleteLattice ){
        std::uniform_int_distribution<indexType>
            rangeDistr(0, latticeSize[dimension]-1);

        range.apply[dimension] = true;
        range.start[dimension] = rangeDistr(random_);
        range.stop[dimension] = rangeDistr(random_);
    } else if(range.start[dimension] <=
	    range.stop[dimension]) {

        std::uniform_int_distribution<indexType>
            rangeDistrStart(0, range.stop[dimension]-
                   range.start[dimension]);
        
        range.start[dimension] += rangeDistrStart(random_);

        // ensure this way that stop value >= start value
        std::uniform_int_distribution<indexType>
            rangeDistrStop(0, range.stop[dimension]-
                   range.start[dimension]);

        range.stop[dimension] = range.start[dimension] +
            rangeDistrStop(random_);
    } else {
        indexType randomRange = latticeSize[dimension] -
            range.start[dimension]+range.stop[dimension];
        std::uniform_int_distribution<indexType> rangeDistrStart(0,randomRange);

        indexType randomChoice = rangeDistrStart(random_);
        range.start[dimension] += randomChoice;
        range.start[dimension] %= latticeSize[dimension];
        
        // ensure this way that stop value after start value in specified
        // interval
        std::uniform_int_distribution<indexType>
            rangeDistrStop(0, randomRange-randomChoice);

        randomChoice = rangeDistrStop(random_);
        range.stop[dimension] = range.start[dimension] +
            rangeDistrStop(random_);
        range.stop[dimension] %= latticeSize[dimension];	
    }
}

//------------------------------------------------------------------------------

//! Scales certain parts of the simulation box by multiplying the space between
//! the single atomic layers in the out-of-plane direction with a factor. If an
//! energetically worse configuration was achieved the changes are reversed.
//! \param tersoff Parameters necessary to compute the energy.
//! \param parameter Optimization parameters such as maximal displacement.
//! \param statistic Optimization statistics.
//! \param range The are the scaling is applied.
void SimulationBox::scale(const TersoffPotential &tersoff,
			  const OptimizationParameter &parameter,
			  OptimizationStatistic &statistic,
			  Range3D<indexType> range)
{
    std::vector<Atom *> shiftAtoms, scaleAtoms;
    double reference, extremeShift = 0, shift, energyBefore;
    Vector3D<spaceType> position;
    Vector3D<indexType> latticeSize=lattice_.getSize();
    uint8_t refLayerIndex;

    std::uniform_real_distribution<double>
	distribution(0.0, parameter.maxScaling);
   
    // scaling uniformly distributed around 1 (= no change)
    double scaling = (1-parameter.maxScaling/2) + distribution(random_);


    CLOG(TRACE, logName_) << "Start scaling of atomic layers with factor " <<
	scaling;

    getRandomInterval(range, outOfPlaneDimension_);
    
    CLOG(TRACE, logName_) << "Using out-of-plane interval [ " <<
	range.start[outOfPlaneDimension_] << "," <<
	range.stop[outOfPlaneDimension_] << "]";

    
    if (range.start[outOfPlaneDimension_] > 0)
        refLayerIndex = range.start[outOfPlaneDimension_]-1;
    else
        refLayerIndex = latticeSize[outOfPlaneDimension_]-1;
    
    // Pick reference atom which serves as reference point. Its distance to the
    // single atoms is scaled.
    reference = lattice_.getFirstAtomInLayer(refLayerIndex,outOfPlaneDimension_)
        ->getPosition()[outOfPlaneDimension_];

    energyBefore = getEnergy(tersoff);
    scaleAtoms = lattice_.getAtomList(range);

    //--------------------------------------------------------------------------
    // scale
    
    for (auto atom: scaleAtoms) {
        // Scale distance between atoms and reference atom.
        
        double localReference;
        position = atom->getPosition();

        if (position[outOfPlaneDimension_] > reference)
            localReference = reference;
        else
            localReference = reference-size_[outOfPlaneDimension_];
        
        // newpos = (oldpos - reference) * scaling + reference
        // shift = newpos - oldpos
        shift = (scaling -1) * position[outOfPlaneDimension_] +
            (1- scaling) * localReference;
        position[outOfPlaneDimension_] += shift;
        atom->setPosition(position, true);
        
        if (std::abs(shift) > std::abs(extremeShift))
            extremeShift = shift;
        
    }

    //--------------------------------------------------------------------------
    // shift
    
    bool smallerThanSize = (range.stop[outOfPlaneDimension_] <
			    (latticeSize[outOfPlaneDimension_]-1) );
    
    if (range.apply[outOfPlaneDimension_] && smallerThanSize) {
        // Shift all atoms above the specified range to preserve the distance
        // between last scaled layer and following one.
        Range3D<indexType> shiftRange = range;

        shiftRange.start[outOfPlaneDimension_]=
            range.stop[outOfPlaneDimension_]+1;
        shiftRange.stop[outOfPlaneDimension_]=
            latticeSize[outOfPlaneDimension_]-1;

        // CLOG(INFO, logName_) << "Using out-of-plane interval [ " <<
        //     shiftRange.start[outOfPlaneDimension_] << "," <<
        //     shiftRange.stop[outOfPlaneDimension_] << "]";
            
        shiftAtoms = lattice_.getAtomList(shiftRange);

        for(auto atom: shiftAtoms) {
            position = atom->getPosition();
            position[outOfPlaneDimension_] += extremeShift;

            // Some atoms get moved twice, once during the scalig procedure and
            // once in this shift operation. This happens when the start layer >
            // stop layer. For these the second move must not be saved because
            // then returning to the original position is not possible any more!
            if (atom->getIndex().isInRange(range))
                atom->setPosition(position);
            else
                atom->setPosition(position, true);
        }
    }
    
    size_[outOfPlaneDimension_] += extremeShift;

    //--------------------------------------------------------------------------
    // check & revert
    
    double energyAfter = getEnergy(tersoff);
    CLOG(DEBUG, logName_) << "before: " << energyBefore << ", after: " <<
	energyAfter;
    
    if (energyAfter > energyBefore) {
        CLOG(DEBUG, logName_) << "energy did not decrease, scaling is reverted";
        statistic.rejectCount ++;
        size_[outOfPlaneDimension_] -= extremeShift;
        
        for (auto atom: scaleAtoms)
            atom->setPreviousPosition();
            
        for (auto atom: shiftAtoms)
            atom->setPreviousPosition();

    } else
        statistic.acceptCount ++;

    updateNeighborList_ = true;
    CLOG(TRACE, logName_) << "finished scaling of simulation box";
}

//##############################################################################

//! Write all atoms of the simulation box to an .xyz file which can be read by
//! several tools
//! \param fileName Name of the file the data are written to.
//! \param range range of atoms that shall be written.
//! \param checkModification If true only recently modified atoms are
//! considered.
//! \param materialName If specified only atoms belonging to this material are
//! considered.
void SimulationBox::writeToXYZ(const std::string fileName,
			       const Range3D<indexType> &range,
			       const bool checkModification,
			       const std::string materialName) const
{

    const PeriodicTable &pt = PeriodicTable::getInstance();
    std::ofstream xyzFile (fileName, std::ios::out);
    std::vector<Atom *> atoms = lattice_.getAtomList(range, checkModification,
						     materialName);

    CLOG(TRACE, logName_) << "Writing atoms to XYZ file " << fileName;
    
    xyzFile << atoms.size() << std::endl;
    xyzFile << "name\tx\ty\tz\ti\tj\tk" << std::endl;

    for (auto atom : atoms){
        xyzFile << pt.getById(atom->getElementId()).symbol << "\t";
        xyzFile << atom->getPosition().strTab() << "\t";
        xyzFile << atom->getIndex().strTab() << std::endl;
    }
    
    xyzFile.close();

    CLOG(TRACE, logName_) << atoms.size() << " Atoms successfully written.";
}

//------------------------------------------------------------------------------

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
