/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#include "InterfaceTool.h"

//******************************************************************************

//! The constructor registers initializes the logger.
//! \param logName Optional parameter defining the logger's name.
//! \param simbox Simulation box the tool operates on.
//! \param tersoff Tersoff potentials used to calculate energies.
//! is supposed to grow.
InterfaceTool::InterfaceTool(SimulationBox &simbox,
			     const TersoffPotential &tersoff,
			     const char *logName):
    simbox_(&simbox), tersoff_(&tersoff),
    posJournal_("IFTries","Locations tried during Interface"
		"Creation procedure"),
    choiceJournal_("IFChoices","Locations chosen to place atoms on"),
    metricJournal_("IFEnergy","Energy/Atom within interfacial layers"),
    logName_(logName), unsuitedPositionsCount_(0)
{
    random_.seed(std::random_device()());
}

//------------------------------------------------------------------------------

//! The constructor registers initializes the logger.
//! \param logName Optional parameter defining the logger's name.
//! \param simbox Simulation box the tool operates on.
//! \param tersoff Tersoff potentials used to calculate energies.
//! \param anionExchange Exchange reaction defining anion exchange.
//! is supposed to grow.
InterfaceTool::InterfaceTool(SimulationBox &simbox,
			     const TersoffPotential &tersoff,
                 ExchangeReaction &anionExchange,
			     const char *logName):
    simbox_(&simbox), tersoff_(&tersoff), 
    anionExchange_(&anionExchange),
    posJournal_("IFTries","Locations tried during Interface"
		"Creation procedure"),
    choiceJournal_("IFChoices","Locations chosen to place atoms on"),
    metricJournal_("IFEnergy","Energy/Atom within interfacial layers"),
    logName_(logName)
{
    random_.seed(std::random_device()());
}

//------------------------------------------------------------------------------

InterfaceTool::~InterfaceTool() {}

//------------------------------------------------------------------------------

//! Creates energetical optimized interface roughness by introducing one cation
//! after the other, whereat each cation is allowed to try several positions and
//! in addition local relaxations are carried out around the chosen location.
//! Usage:
//! \code
//! InterfaceTool ifTool(...);
//! ifTool.createOptimizedRoughness(...);\endcode
//! \warning Works only properly if start and stop values of the given range
//! point to a cation layer in the lattice.
//! \param interfaceRange Range containing upper and lower limit in the
//! out-of-plane direction. the other entries are ignored.
//! \param atomsToDeployCount Amount of introduced cations.
//! \param positionsToCheckCount Amount of different positions each atom tries.
//! \param material Material the interface roughness is constructed of.
//! \param parameter Optimization parameters used in the relaxation step.
//! \param collecction Material collectoin that can be used to look up additional materials.
//! \param metric Defines metric to determine optimum positions.
//! \param journalPreamble Path preamble that shall be used to store generated
//! journals.
void InterfaceTool::createOptimizedRoughness(
					     Range3D<indexType> interfaceRange,
					     const unsigned int atomsToDeployCount,
					     const unsigned int positionsToCheckCount,
					     const Material &material,
					     const OptimizationParameter &parameter,
                         const MaterialCollection &collection,
                         const InterfaceToolMetric metric,
					     const std::string &journalPreamble)
{

    Lattice &lattice = simbox_->getLattice();
    Vector3D<indexType> latticeSize = lattice.getSize();
    uint8_t outOfPlaneDimension = simbox_->getOutOfPlaneDimension();
    unsigned int atomCount;
    Vector3D<indexType> index;
    Atom *atom;
    Atom *modAtomBestConfiguration;
    std::vector<Atom> bestConfiguration, originalConfiguration;
    Range3D<indexType> range;
    Optimization optimization("localRelaxation");
    const Material *materialPointer;

    std::stringstream journalName{};
    std::stringstream journalDesc{};

    switch (metric) {
    case ITM_Energy:
        journalName << "IF Energy Journal L";
        journalDesc << "Energy/Atom within interfacial layers";
        break;
    case ITM_BondStrain:
        journalName << "IF Bond Strain Journal L";
        journalDesc << "Bond Strain within interfacial layers";
        break;
    case ITM_GrowthStrain:
        journalName << "IF Growth Strain Journal L";
        journalDesc << "Growth Strain within interfacial layers";
        break;
    case ITM_DistanceStrain:
        journalName << "IF Distance Strain Journal L";
        journalDesc << "Distance Strain within interfacial layers";
        break;
    }

    journalName << (indexType)((interfaceRange.stop[outOfPlaneDimension] -
                                interfaceRange.start[outOfPlaneDimension]) / 2 +
                               interfaceRange.start[outOfPlaneDimension]);

    metricJournal_.setName(journalName.str().c_str());
    metricJournal_.setDescription(journalDesc.str().c_str());
    
    uint8_t inPlaneDimension1 = (outOfPlaneDimension+1)%3;
    uint8_t inPlaneDimension2 = (outOfPlaneDimension+2)%3;
    
    std::uniform_int_distribution<int> inPlaneDistribution1
	(0,latticeSize[inPlaneDimension1]-1);
    std::uniform_int_distribution<int> inPlaneDistribution2
	(0,latticeSize[inPlaneDimension2]-1);
    std::uniform_int_distribution<int> directionDistribution(0,3);

    
    CLOG(INFO, logName_) <<"Starting creation of Interface roughness in region "
			  << interfaceRange.str();
    CLOG(DEBUG, logName_) << "Using material " << material.getName();
    
    bool startBigger =(latticeSize[outOfPlaneDimension] <=
		       interfaceRange.start[outOfPlaneDimension]);
    bool stopBigger = (latticeSize[outOfPlaneDimension] <=
		       interfaceRange.stop[outOfPlaneDimension]);
    
    if (startBigger || stopBigger){
	CLOG(WARNING, logName_) << "Layer range outside the specified lattice!"
	    " Aborting operation!";
        return;
    }
  
    for (auto i: {interfaceRange.start[outOfPlaneDimension]-2,
                interfaceRange.start[outOfPlaneDimension]-1,
                interfaceRange.start[outOfPlaneDimension],
                interfaceRange.start[outOfPlaneDimension]+1,
                interfaceRange.start[outOfPlaneDimension]+2}) {
        const Atom * atom = lattice.getFirstAtomInLayer(i, outOfPlaneDimension);
        CLOG(DEBUG, logName_) << "Material in layer " << i << ": "
                             << atom->getMaterial()->getName();
        CLOG(DEBUG, logName_) << "Anion in layer " << i << ": "
                             << std::get<0>(atom->getMaterial()->getAnions()[0]);

    }

    //    std::vector<Atom *> atomsBelow = lattice.getAtomsInLayer(interfaceRange.start[2]-1, outOfPlaneDimension);
    //    std::vector<Atom *> atomsAbove = lattice.getAtomsInLayer(interfaceRange.start[2], outOfPlaneDimension);

    const Atom * atomBelow = lattice.getFirstAtomInLayer(interfaceRange.start[2]-1, outOfPlaneDimension);
    const Atom * atomAbove = lattice.getFirstAtomInLayer(interfaceRange.start[2], outOfPlaneDimension);
  
    // exchange anions at the layer BELOW the interface!
    // ATTENTION! Experimental code!
    CLOG(DEBUG, logName_) << "Exchanging atoms in layer " << interfaceRange.start[2]-1;
    anionExchange_->replaceAnions(simbox_, interfaceRange.start[2]-1, 2,
                                  atomBelow->getMaterial()->getName(),
                                  atomAbove->getMaterial()->getName());
    //                                  std::get<0>(atomsBelow[0]->getMaterial()->getAnions()[0]),
    //                                  std::get<0>(atomsAbove[0]->getMaterial()->getAnions()[0]));

    XmlHandler xmlFile;
    simbox_->writeToXYZ("aaa.xyz");
    
    XmlHandler xmlEJournal(journalName.str().c_str());

    try {
        // if load throws an exception, get is not called
        xmlEJournal.load(journalPreamble + "ifEJournal.xml");
        xmlEJournal.get(metricJournal_, journalName.str().c_str());
    } catch (XmlException e) {
        CLOG(ERROR, logName_) << "Could not load journal file: " << journalPreamble + "ifEJournal.xml";
    }
    
    materialPointer = & simbox_->getMaterials().getByNameOrAdd(material);
    
    interfaceRange.apply[inPlaneDimension1] = false;
    interfaceRange.apply[inPlaneDimension2] = false;
    interfaceRange.apply[outOfPlaneDimension] = true;

    optimization.registerAction(std::string("MMC"), &SimulationBox::mmcRelax,1);
    
    modificationIndex_ = 1;
    atomCount=0;
    range.apply={true,true,true};
    uint32_t atomCheckStep = atomsToDeployCount/10;
    uint32_t atomCheckCount = atomCheckStep;
    while (atomCount < atomsToDeployCount) {

        Vector3D<indexType> bestIndex;
        elementType cationId;
        
        // Start at top of allowed region an descend afterwards until a valid
        // location is found. In this way the trace of an actual atom
        // approaching the surface is modeled accurately.
        index[outOfPlaneDimension] =interfaceRange.stop[outOfPlaneDimension];
        index[inPlaneDimension1] = inPlaneDistribution1(random_);
        index[inPlaneDimension2] = inPlaneDistribution2(random_);

        // maybe integrate this into findSuitablePosition
        atom = lattice(index);
        if (atom == nullptr) continue;	

        CLOG(INFO, logName_) << "processing atom number " << atomCount;
        CLOG(DEBUG, logName_) << "chose index ", index.str();
        
        checkedIndices_.clear();
        cationId = materialPointer->getRandomCation(random_);

        double bestMetric = std::numeric_limits<double>::max();
        int direction = directionDistribution(random_);
        
        posJournal_.add("placing new atom");
        choiceJournal_.add("placing new atom");
       
        unsuitedPositionsCount_=0;
        for (unsigned int i=0; i < positionsToCheckCount; i++) {
            //! \attention Debug output temporarily commented out!
            index = findSuitablePosition(index, interfaceRange,
                         direction, materialPointer);

            if (unsuitedPositionsCount_ > unsuitedPositionsLimit_){
                CLOG(WARNING, logName_) << "Too many unsuited positions encountered."
                                        << " Aborting at iteration # " << i << "!";
                break;
            }
                
            // remember which indices have already been processed
            checkedIndices_.insert(index);
            
            range.start=(index-4);
            range.stop=(index+4);
            range.fitInBox(latticeSize);

            originalConfiguration = lattice.backupAtoms(range);

            CLOG(DEBUG, logName_) << "found suitable location # " << i << " at "
                     << index.str();
            posJournal_.add("suitable location at " + index.str());
            choiceJournal_.add(index.str());
            
            atom = lattice(index);
            atom->setElement(cationId, materialPointer, AtomState::ModifiedInterface);

            addAnions(index, materialPointer, parameter.anionPassivationProbability);
            CLOG(TRACE, logName_) << "starting optimization";
            optimization.runStatic(*simbox_, range, *tersoff_, parameter);

            double tmpMetric{};

            switch (metric) {
            case ITM_Energy:
                tmpMetric = simbox_->getEnergy(*tersoff_, interfaceRange);
                break;
            case ITM_BondStrain:
                tmpMetric = std::abs(simbox_->calculateStrainSingleAtom(atom, SCM_BondStrain, collection));
                break;
            case ITM_GrowthStrain:
                tmpMetric = std::abs(simbox_->calculateStrainSingleAtom(atom, SCM_GrowthStrain, collection));
                break;
            case ITM_DistanceStrain:
                tmpMetric = std::abs(simbox_->calculateStrainSingleAtom(atom, SCM_DistanceStrain, collection));
                break;
            }

        
            if (tmpMetric < bestMetric) {
                bestConfiguration = lattice.backupAtoms(range);
                bestMetric = tmpMetric;
                bestIndex = index;
                modAtomBestConfiguration = atom;
            }

            //! \attention Debug output temporarily commented out!
            lattice.restoreAtoms(originalConfiguration);

            unsuitedPositionsCount_=0;
            //move from current position at least one step
            index = moveAtom(index, interfaceRange, direction, materialPointer);
            //! \attention Debug output temporarily commented out!
            posJournal_.add("moveAtom: chose " + index.str());
        }

        CLOG(DEBUG, logName_) << "Found best position at " << bestIndex.str();
        posJournal_.add("Best position at " + bestIndex.str());
        metricJournal_.add(std::make_tuple(atomCount,bestMetric/simbox_->getAtomCount(interfaceRange)));

        lattice.restoreAtoms(bestConfiguration);
        atomCount++;

        modAtomBestConfiguration->setModificationOrder(modificationIndex_);
        modificationIndex_ += 1;
        
        if (atomCount < atomCheckCount)
            continue;
                
        atomCheckCount += atomCheckStep;
        CLOG(INFO, logName_) << "\e[1m" << std::fixed << std::setprecision(0)
                             << 100.0*((double)atomCount / (double)atomsToDeployCount)
                             << "% of atoms processed" << "\e[0m";
    }

    
    XmlHandler handler("Journal Writer");
    handler.set(posJournal_);
    handler.set(choiceJournal_);
    handler.save(journalPreamble + "posJournal.xml");

    xmlEJournal.set(metricJournal_);
    xmlEJournal.save(journalPreamble + "ifEJournal.xml");

    JsonHandler jHandler;
    jHandler.set(posJournal_);
    jHandler.set(choiceJournal_);
    jHandler.save(journalPreamble + "posJournal.json");    

    CLOG(TRACE, logName_) << "Creation of Interface roughness finished.";
}


//------------------------------------------------------------------------------

//! Start from an initial position and check for a suitable location to place a
//! cation. A position is suitable if all cations two layers below belong to the
//! same material and the current cation belongs to another material.
//! \param position Starting location in the lattice.
//! \param interfaceRange Range containing upper and lower limit in the
//! out-of-plane direction. the other entries are ignored.
//! \param direction One of four directions the atom is allowed to move over the
//! surface.
//! \param material Material the interface roughness is constructed of.
//! \return A suitable position.
Vector3D<indexType> InterfaceTool::findSuitablePosition(
				       const Vector3D<indexType> position,
				       const Range3D<indexType> interfaceRange,
				       const int direction,
				       const Material *material)
{

    Lattice &lattice = simbox_->getLattice();
    Vector3D<indexType> workVector=position;
    
    uint8_t outOfPlaneDimension = simbox_->getOutOfPlaneDimension();
    uint8_t inPlaneDimension1 = (outOfPlaneDimension+1)%3;
    uint8_t inPlaneDimension2 = (outOfPlaneDimension+2)%3;

    CLOG(TRACE, logName_) << "Starting findSuitablePosition()";
    posJournal_.add("findSuitablePosition: trying " + position.str());

    if (unsuitedPositionsCount_ > unsuitedPositionsLimit_)
        return position;

    unsuitedPositionsCount_++;
    
    if((lattice(position))->getMaterial()->getName() == material->getName()){
	workVector = moveAtom(position, interfaceRange, direction, material);
	return findSuitablePosition(workVector, interfaceRange,
					    direction, material);
    }
    else{
	if (position[outOfPlaneDimension] ==
	    interfaceRange.start[outOfPlaneDimension]){
	    // can not descend any further
	    return position;
	}
	
	//check if all cations 2 layers below are from the same material
	workVector[outOfPlaneDimension] = position[outOfPlaneDimension] -2;
	
	std::vector<int> indices = {0,1,2,3};
	shuffle (indices.begin(), indices.end(), random_);
	
	bool checkedIndexFound = false;
	for(auto index : indices){
    
	    workVector[inPlaneDimension1] = position[inPlaneDimension1] +
		std::get<0>(cationPositionsOutOfPlane[index]);
	    workVector[inPlaneDimension2] = position[inPlaneDimension2] +
		std::get<1>(cationPositionsOutOfPlane[index]);

	    workVector.fitInBox(lattice.getSize());
	    
	    if ( (lattice(workVector))->getMaterial()->getName() !=
		 material->getName() ){
		CLOG(DEBUG, logName_) << "going down to position " <<
		    workVector.str();

		if (checkedIndices_.count(workVector) > 0){
		    checkedIndexFound = true;
		    posJournal_.add("found already checked index at " +
				    workVector.str());
            CLOG(DEBUG, logName_) << "position was already visited. CONTINOUING";
		    continue;
		}

		return findSuitablePosition(workVector, interfaceRange,
					    direction, material);
	    }
	}

	// when all free locations have already been checked mark also current
	// index as checked and proceed to next one
	if (checkedIndexFound){
	    posJournal_.add("marked " + position.str() + " as already checked "+
			    "because all empty spaces below have already been "+
			    "checked!");
        CLOG(DEBUG, logName_) << "marked " + position.str() + " as already checked "+
			    "because all empty spaces below have already been "+
			    "checked!";
	    checkedIndices_.insert(position);
	    workVector =moveAtom(position, interfaceRange, direction, material);
	    return findSuitablePosition(workVector, interfaceRange,
					direction, material);
	}
	    
	//If this line is reached all cation spots 2 layers below are filled so
	//the actual position  is the next valid spot to try.
	return position;
    }
    
}

//------------------------------------------------------------------------------

//! Moves the atom at least once! Checks if the next two cations in the same
//! plane belong to the same material as the specified material. If this is not
//! the case a new location has been found which is returned. In the other case
//! a position two atomic layers above is tried.
//! \param position Starting location in the lattice.
//! \param interfaceRange Range containing upper and lower limit in the
//! out-of-plane direction. the other entries are ignored.
//! \param direction One of four directions the atom is allowed to move over the
//! surface.
//! \param material Material the interface roughness is constructed of.
//! \return A suitable position.
Vector3D<indexType> InterfaceTool::moveAtom(const Vector3D<indexType> position,
				       const Range3D<indexType> interfaceRange,
				       const int direction,
				       const Material *material)
{

    Vector3D<indexType> workVector = position;
    int offset;
    Lattice &lattice = simbox_->getLattice();
    
    uint8_t outOfPlaneDimension = simbox_->getOutOfPlaneDimension();
    uint8_t inPlaneDimension1 = (outOfPlaneDimension+1)%3;
    uint8_t inPlaneDimension2 = (outOfPlaneDimension+2)%3;

    std::uniform_real_distribution<double> directionDistribution(0.0,1.0);

    //! \attention Debug output temporarily commented out!
    posJournal_.add("moveAtom: starting at " + position.str());

    unsuitedPositionsCount_++;
    
    //for each direction two atoms can be jumped to next, since one is the so
    //called "easy direction" it is chosen in relation 3:1
    if (directionDistribution(random_) < 0.75)
	offset=0;
    else
	offset=1;

    //check in plane cations in given direction
    for (int i=0; i<2; i++){
	
	workVector[inPlaneDimension1] = position[inPlaneDimension1] +
	    std::get<0>(cationPositionsInPlane[direction*2 +offset]);
	workVector[inPlaneDimension2] = position[inPlaneDimension2] +
	    std::get<1>(cationPositionsInPlane[direction*2 +offset]);

	workVector.fitInBox(lattice.getSize());

	if((lattice(workVector))->getMaterial()->getName() !=
	   material->getName()){

	    if (checkedIndices_.count(workVector) > 0){
		posJournal_.add("encountered already checked location in plane "
				"at " + workVector.str());
	    }
	    else{
		return workVector;
	    }
	}

	//if first choice failed try other one
	offset = (offset+1)%2;
    }

    //both in plane cations in given direction are already occupied by correct
    //material, therefore go up two atomic layers and try there

    if (workVector[outOfPlaneDimension] !=
	interfaceRange.stop[outOfPlaneDimension]){
	
	std::uniform_int_distribution<int> offsetDistribution(0,1);
	offset = offsetDistribution(random_);
	if (direction>1) offset += 2;
	
	workVector[outOfPlaneDimension] = position[outOfPlaneDimension] + 2;
	workVector[inPlaneDimension1] = position[inPlaneDimension1] +
	    std::get<0>(cationPositionsOutOfPlane[offset]);
	workVector[inPlaneDimension2] = position[inPlaneDimension2] +
	    std::get<1>(cationPositionsOutOfPlane[offset]);

	workVector.fitInBox(lattice.getSize());

	if((lattice(workVector))->getMaterial()->getName() !=
	   material->getName()){

	    if (checkedIndices_.count(workVector) > 0){
		posJournal_.add("encountered already checked location out of "
				"plane at " + workVector.str());
	    }
	    else{
		return workVector;
	    }
	}
    }
    else{
	
	workVector[inPlaneDimension1] = position[inPlaneDimension1] +
	    std::get<0>(cationPositionsInPlane[direction*2 +offset]);
	workVector[inPlaneDimension2] = position[inPlaneDimension2] +
	    std::get<1>(cationPositionsInPlane[direction*2 +offset]);

	workVector.fitInBox(lattice.getSize());

    }

    return moveAtom(workVector, interfaceRange, direction, material);
}
    
//------------------------------------------------------------------------------

//! Checks if the neighboring atoms in the easy directions are of the same
//! material. If that is the case an anion is introduced above them.
//! \param position Starting location in the lattice.
//! \param material Material the interface roughness is constructed of.
//! \param passivationProbability Determines, if an anion is added to matching cations.
void InterfaceTool::addAnions(const Vector3D<indexType> position,
			      const Material *material, double passivationProbability, bool requireMatchingCation)
{
    Vector3D<indexType> workVector;
    Lattice &lattice = simbox_->getLattice();
    uint8_t outOfPlaneDimension = simbox_->getOutOfPlaneDimension();
    uint8_t inPlaneDimension1 = (outOfPlaneDimension+1)%3;
    uint8_t inPlaneDimension2 = (outOfPlaneDimension+2)%3;
    int counter=0;

    CLOG(TRACE, logName_) << "adding anions";

    for (int i=-1; i<2; i+=2){
        workVector[inPlaneDimension1] = position[inPlaneDimension1] + 2*i;
        workVector[inPlaneDimension2] = position[inPlaneDimension2] + 2*i;
        workVector[outOfPlaneDimension] = position[outOfPlaneDimension];

        workVector.fitInBox(lattice.getSize());
	
        if (requireMatchingCation == true) {
            if ( (lattice(workVector))->getMaterial()->getName() !=
                 material->getName())
                continue;
        }

        std::uniform_real_distribution<double> distribution(0.0,1.0);

        if (distribution(random_) < passivationProbability) {

            Atom *atomBelow;
            workVector[inPlaneDimension1] = position[inPlaneDimension1] + i ;
            workVector[outOfPlaneDimension] = position[outOfPlaneDimension] -1;

            // alternate such that interface anions stay roughly at same position
            if ( ( (position[outOfPlaneDimension]/2) % 2) == 0)
                workVector[inPlaneDimension2] = position[inPlaneDimension2] + i - 2 ;
            else
                workVector[inPlaneDimension2] = position[inPlaneDimension2] + i + 2 ;
            
            workVector.fitInBox(lattice.getSize());
            atomBelow = lattice(workVector);

            // If anion in layer below is not of interface material then change
            // only that atom and not the one in the layer above. This way atoms
            // not belonging to the interface material shall be always on top of
            // the interface.
            //if ( atomBelow->getMaterial()->getName() != material->getName()){
            //    atomBelow->setElement(material->getRandomAnion(random_), material, AtomState::ModifiedInterface);
            //    atomBelow->setModificationOrder(modificationIndex_);
            //    continue;
            //}            

            Atom *atomAbove;
            
            workVector[inPlaneDimension1] = position[inPlaneDimension1] + i;
            workVector[inPlaneDimension2] = position[inPlaneDimension2] + i;
            workVector[outOfPlaneDimension] = position[outOfPlaneDimension] +1;

            workVector.fitInBox(lattice.getSize());
            atomAbove = lattice(workVector);
            
            // If an atom that was modified through an exchange reaction was found
            // below, it is copied to the position above - including its
            // modification state.
            // The atom below is set to a random anion element of the initial
            // material. The modification state is first set to Unknown and cleared
            // in the next line.
            if (atomBelow->getState().test(AtomState::ModifiedExchangeReaction)) {
                atomAbove->setElement(atomBelow->getElementId(),
                                                  material, AtomState::ModifiedInterface);

                atomAbove->addModificationState(AtomState::ModifiedExchangeReaction);
                atomBelow->setElement(atomBelow->getMaterial()->getRandomAnion(random_), atomBelow->getMaterial(), 
                        AtomState::ModifiedUnknown);
                atomBelow->clearState();
            } else {
                // Otherwise, a random anion from the given material is chosen.
                atomAbove->setElement(material->getRandomAnion(random_),
                                                  material, AtomState::ModifiedInterface);
            }

            // \todo The modification index is not incremented for anions
            // at the moment. We need to check, which way is better.
            // In its current state, the index is not unique (a cation + potential
            // anions, which are added with the cation) have the same index.
            // Furthermore, if InterfaceTool is run repeatedly (as for multiple
            // interfaces in one structure), the index again is not unique.
            atomAbove->setModificationOrder(modificationIndex_);

            counter++;
        }
    }

    CLOG(DEBUG, logName_) << "added successfully " << counter <<
	" anions for position " << position.str();
}

//------------------------------------------------------------------------------

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
