/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#include "StrainTool.h"

//------------------------------------------------------------------------------

//! The constructor registers initializes the logger.
//! \param logName Optional parameter defining the logger's name.
//! \param simbox Simulation box the tool operates on.
//! is supposed to grow.
StrainTool::StrainTool(SimulationBox &simbox,
		       const char *logName):
    simbox_(&simbox), logName_(logName), strainFieldReady_(false)
{
    // nothing to be done
}

//------------------------------------------------------------------------------

//! The constructor registers initializes the logger.
//! \param logName Optional parameter defining the logger's name.
//! \param simbox Simulation box the tool operates on.
StrainTool::StrainTool(std::shared_ptr<SimulationBox> simbox,
   const char *logName):
    logName_(logName), strainFieldReady_(false)
{
    simbox_ = simbox.get();
}

//------------------------------------------------------------------------------

StrainTool::~StrainTool() {}

//------------------------------------------------------------------------------

//! For all atoms in the simulation box different methods of strain are
//! calculated and written to the output file. The single methods are:
//! - Calculate lattice constant from bond length between two atoms and compare
//! this to the relaxed lattice constant of the material consisting of the two
//! elements.
//! - Compare the distance in growth direction only.
//! \param materials Contains different materials which are required to
//! determine the relaxed lattice constants.
//! \param outputPreamble Preamble of output file names.
//! \param refDistance Single distance that is used on all bonds to determine
//! \param writeToFile Disables file output, if set to false.
//! strain.
LayerStrainInfo StrainTool::calculateStrain( const MaterialCollection &materials,
				    const std::string &outputPreamble,
				    double refDistance, bool writeToFile ) const
{

    LayerStrainInfo result;
    double bondStrain, growthStrain, distanceStrain;
    Lattice const &lattice = simbox_->getLattice();
    uint8_t outOfPlaneDimension = simbox_->getOutOfPlaneDimension();

    CLOG(TRACE, logName_) << "Calculating Strain of simulation box.";
    CLOG(DEBUG, logName_) << "Writing strain results to files " << outputPreamble
			  << "*";

    XmlHandler xmlHandler("strainWriter");
    xmlHandler.set(*simbox_);
    xmlHandler.set(materials);
    xmlHandler.save(outputPreamble + "_information.xml");
    
    std::ofstream bondFile (outputPreamble+"_per_layer_bond.dat",std::ios::out);
    std::ofstream growthFile (outputPreamble+"_per_layer_growth.dat",
			      std::ios::out);
    std::ofstream distanceFile (outputPreamble+"_per_layer_dist.dat",
			      std::ios::out);
    
    bondFile.precision(strainOutputPrecision);
    growthFile.precision(strainOutputPrecision);
    distanceFile.precision(strainOutputPrecision);

    // Take relaxed lattice constant of the material the first atom in the
    // lowermost layer belongs to if no value has been specified.
    if (refDistance < 0){
	refDistance = lattice.getFirstAtomInLayer(0,outOfPlaneDimension)->
	    getMaterial()->getLatticeConstant();
    }
	    
    for (int i=0; i<lattice.getSize()[outOfPlaneDimension]; i++){

        calculateStrainLayer(materials, outputPreamble, i, refDistance,
                     bondStrain,growthStrain, distanceStrain, writeToFile);
        
        bondFile << i << " " << bondStrain << std::endl;
        growthFile << i << " " << growthStrain << std::endl;	
        distanceFile << i << " " << distanceStrain << std::endl;	

        LayerStrainInfoData layerInfo = std::make_tuple((indexType)i, bondStrain, growthStrain, distanceStrain);
        result << layerInfo;
    }

    bondFile.close();
    growthFile.close();
    distanceFile.close();

    CLOG(TRACE, logName_) << "Finished strain calculation successfully.";

    return result;
}


//------------------------------------------------------------------------------

//! For all atoms inside the specified range different methods of strain are
//! calculated and written to the output file. The single methods are:
//! - Calculate lattice constant from bond length between two atoms and compare
//! this to the relaxed lattice constant of the material consisting of the two
//! elements.
//! - Compare the distance in growth direction only.
//! \param materials Contains different materials which are required to
//! determine the relaxed lattice constants.
//! \param outputPreamble Preamble of output file names.
//! \param layer Atomic layer index in growth direction which is processed.
//! \param refDistance Fixed reference distance used to calculate distance
//! strain. 
//! \param bondStrain The average bond strain for the whole layer has to be
//! saved here.
//! \param growthStrain The average growth strain for the whole layer has to be
//! saved here.
//! \param writeToFile Disables output to files, if set to false.
//! \throws MaterialException if a material could not be found.
void StrainTool::calculateStrainLayer(const MaterialCollection &materials,
					 const std::string &outputPreamble,
					 const int layer,
					 const double refDistance,
					 double &bondStrain,
					 double &growthStrain,
					 double &distanceStrain,
                     bool writeToFile) const
{
    uint8_t outOfPlaneDimension = simbox_->getOutOfPlaneDimension();
    uint8_t inPlaneDimension1 = (outOfPlaneDimension+1)%3;
    uint8_t inPlaneDimension2 = (outOfPlaneDimension+2)%3;
    Vector3D<indexType> vector;
    const Atom *atom;
    int atomCount2D=0;
    double bondStrain2D=0, growthStrain2D=0, distanceStrain2D=0;
    Lattice const &lattice = simbox_->getLattice();
    Vector3D<indexType> latticeSize=lattice.getSize();
    Vector3D<spaceType> simboxSize=simbox_->getSize();

    uint32_t reqDigits =  lattice.getSize()[outOfPlaneDimension] > 0 ?
        (int) log10 ((double) lattice.getSize()[outOfPlaneDimension]) + 1 : 1;
    std::stringstream paddedLayerIndex{};
    paddedLayerIndex << std::right << std::setw(reqDigits) << std::setfill('0') << layer;

    std::ofstream bondFile;
    std::ofstream bondFile2D;
    std::ofstream growthFile2D;
    std::ofstream distanceFile2D;

    if (writeToFile == true) {
        bondFile = std::ofstream(outputPreamble+"_per_bond.dat",
                    std::ios::out | std::ios::app);
        bondFile2D = std::ofstream(outputPreamble+"_strain_2D_layer_" +
                  paddedLayerIndex.str() + "_bond.dat", std::ios::out);
        growthFile2D = std::ofstream(outputPreamble+"_strain_2D_layer_" +
                  paddedLayerIndex.str() + "_growth.dat", std::ios::out);
        distanceFile2D = std::ofstream(outputPreamble+"_strain_2D_layer_" +
                  paddedLayerIndex.str() + "_dist.dat", std::ios::out);
        
        bondFile.precision(strainOutputPrecision);
        bondFile2D.precision(strainOutputPrecision);
        growthFile2D.precision(strainOutputPrecision);
        distanceFile2D.precision(strainOutputPrecision);
    }
    
    vector[outOfPlaneDimension] = layer;
    
    for (int i=0; i<latticeSize[inPlaneDimension1]; i++){
	vector[inPlaneDimension1] = i;
	for (int j=0; j<latticeSize[inPlaneDimension2]; j++){
	    vector[inPlaneDimension2] = j;
	    int bondCount = 0;
	    double atomBondStrain=0, atomGrowthStrain=0, atomDistanceStrain=0;

	    CLOG(TRACE, logName_) << "trying atom " << vector.str();
	    
	    atom = lattice(vector);
	    if (atom == nullptr){
            if (writeToFile == true) {
                bondFile2D << "NaN ";
                growthFile2D << "NaN ";
                distanceFile2D << "NaN ";
            }
            continue;
	    }

        if (writeToFile == true)
            bondFile << "processing atom " << atom->getIndex() << std::endl;

	    for (auto neighbor: atom->getNeighbors()){

		bondCount++;
		double relaxedLatticeConstant, bondLatticeConstant;
		double tmpBondStrain, tmpGrowthStrain, tmpDistanceStrain;

		try{
		    relaxedLatticeConstant = materials.getByCationsAnions(
			 {atom->getElementId()},
			 {lattice(neighbor)->getElementId()}
				                     ).getLatticeConstant();
		}catch(MaterialException &e){
		    // If no material found switch cations and anions.
            try{
		    relaxedLatticeConstant = materials.getByCationsAnions(
			 {lattice(neighbor)->getElementId()},
			 {atom->getElementId()}
				                     ).getLatticeConstant();
            }catch(MaterialException &e){
                CLOG(ERROR, logName_) << "Did not find material (" <<
                    atom->getElementId() << "," <<
                    lattice(neighbor)->getElementId() <<
                    ") by cation/anion or anion/cation combination!";
            }
		}
		
		Vector3D<spaceType> distance =
		    atom->getPosition().getMinimalDistance(
				   lattice(neighbor)->getPosition(),simboxSize);
		
		bondLatticeConstant =
		    Zincblende::squaredBondLengthToLatticeConstant(
					 distance.squaredLength());
		    
        if (writeToFile == true) {
            bondFile << "to " << lattice(neighbor)->getIndex() << ": ";
            bondFile << "elementIDs " << atom->getElementId();
            bondFile << " and " << lattice(neighbor)->getElementId()<< ";";
        }

		tmpBondStrain = (bondLatticeConstant-
			 relaxedLatticeConstant)/relaxedLatticeConstant * 100;
		tmpGrowthStrain = (std::abs(4*distance[outOfPlaneDimension])-
			 relaxedLatticeConstant)/relaxedLatticeConstant * 100;
		tmpDistanceStrain = (std::abs(4*distance[outOfPlaneDimension])-
			 refDistance)/refDistance * 100;
		
        if (writeToFile == true) {
            bondFile << " strain B: " << tmpBondStrain;
            bondFile << " G: " << tmpGrowthStrain;
            bondFile << " D: " << tmpDistanceStrain;
            
            bondFile << " lcStrained B: " << bondLatticeConstant;
            bondFile << " G: " <<std::abs(4*distance[outOfPlaneDimension]);
            bondFile << " D: " << refDistance;
            bondFile << " lcRelaxed " << relaxedLatticeConstant;
            bondFile << std::endl;
        }

		atomBondStrain += tmpBondStrain;
		atomGrowthStrain += tmpGrowthStrain;
		atomDistanceStrain += tmpDistanceStrain;
	    }
	    atomCount2D++;
	    
        if (writeToFile == true)
            bondFile2D << atomBondStrain / bondCount << " ";
	    bondStrain2D += atomBondStrain / bondCount;
	    
        if (writeToFile == true)
            growthFile2D << atomGrowthStrain / bondCount << " ";
	    growthStrain2D += atomGrowthStrain / bondCount;

        if (writeToFile == true)
            distanceFile2D << atomDistanceStrain / bondCount << " ";
	    distanceStrain2D += atomDistanceStrain / bondCount;	    
	}

    if (writeToFile == true) {
        bondFile2D << std::endl;
        growthFile2D << std::endl;
        distanceFile2D << std::endl;
    }
	
    }

    if (writeToFile == true) {
        bondFile.close();
        bondFile2D.close();
        growthFile2D.close();
        distanceFile2D.close();
    }
    
    bondStrain = bondStrain2D / atomCount2D;
    growthStrain = growthStrain2D / atomCount2D;
    distanceStrain = distanceStrain2D / atomCount2D;
}

//------------------------------------------------------------------------------

double StrainTool::calculateStrainAtom (const MaterialCollection &materialsInSimbox, 
            const MaterialCollection &otherMaterials,
            const Atom *atom, StrainCalculationModes mode,
            double refDistance) const {
    Lattice const &lattice = simbox_->getLattice();
    Vector3D<spaceType> simboxSize=simbox_->getSize();
    uint8_t outOfPlaneDimension = simbox_->getOutOfPlaneDimension();
    double tmpAtomStrain{};

    uint32_t bondCount{};

    if (atom != nullptr){
        for (auto neighbor: atom->getNeighbors()){
            bondCount++;
            double relaxedLatticeConstant, bondLatticeConstant;
            double tmpBondStrain{};

            try{
                relaxedLatticeConstant = materialsInSimbox.getByCationsAnions(
                 {atom->getElementId()},
                 {lattice(neighbor)->getElementId()}
                                         ).getLatticeConstant();
            }catch(MaterialException &e){
                // If no material found switch cations and anions.
                try{
                relaxedLatticeConstant = materialsInSimbox.getByCationsAnions(
                 {lattice(neighbor)->getElementId()},
                 {atom->getElementId()}
                                         ).getLatticeConstant();
                }catch(MaterialException &e){
                    CLOG(INFO, logName_) << "Did not find material (" <<
                        atom->getElementId() << "," << lattice(neighbor)->getElementId() <<
                        ") in simulation box by cation/anion or anion/cation combination!";

                    try {
                        const Material &tmpMaterial = otherMaterials.getByCationsAnions(
                                                     {atom->getElementId()},
                                                     {lattice(neighbor)->getElementId()}
                                                     );

                        simbox_->addMaterial(tmpMaterial);
                        relaxedLatticeConstant = tmpMaterial.getLatticeConstant();

                        CLOG(INFO, logName_) << "Found material in additional collection";

                    } catch(MaterialException &e){
                        CLOG(ERROR, logName_) << "Did not find material (" << atom->getElementId() <<
                            "," << lattice(neighbor)->getElementId() <<
                            ") in additional materials by cation/anion or anion/cation combination!";
                    }
                }
            }
            
            Vector3D<spaceType> distance =
                atom->getPosition().getMinimalDistance(
                       lattice(neighbor)->getPosition(),simboxSize);
            
            bondLatticeConstant =
                Zincblende::squaredBondLengthToLatticeConstant(
                         distance.squaredLength());
                
            switch (mode) {
            case SCM_BondStrain:
                tmpBondStrain = (bondLatticeConstant -
                     relaxedLatticeConstant) / relaxedLatticeConstant * 100.0;
                break;
            case SCM_GrowthStrain:
                tmpBondStrain = (std::abs(4.0 * distance[outOfPlaneDimension]) -
                     relaxedLatticeConstant) / relaxedLatticeConstant * 100.0;
                break;
            case SCM_DistanceStrain:
                tmpBondStrain = (std::abs(4.0 * distance[outOfPlaneDimension]) -
                     refDistance)/refDistance * 100.0;
                break;
            }
            
            tmpAtomStrain += tmpBondStrain;
        }

        tmpAtomStrain /= (double)bondCount;
    }

    return tmpAtomStrain;
}

//------------------------------------------------------------------------------

//! Calculate strain for each bond
StrainInfo StrainTool::getStrainInfo(const MaterialCollection &materials,
				     double refDistance) const {

    StrainInfo result; //= StrainInfo(); 
    Lattice const &lattice = simbox_->getLattice();
    Vector3D<spaceType> simboxSize=simbox_->getSize();   

    //simbox_->generateNeighbors(1, 0.0, false);

    double strain = 0.0;

    for (auto atom : lattice.getAtomList()) {

        for (auto neighbor: atom->getNeighbors()) {
            const Atom *neighborAtom = lattice(neighbor);

            const Material &mat = materials.getByElementID(atom->getElementId(),
					      neighborAtom->getElementId());

            if (refDistance < 0.0) {
                refDistance = mat.getLatticeConstant();
                refDistance = Zincblende::latticeConstantToBondLength(
								refDistance);
            }

	    CLOG(TRACE, logName_) << "Using ref distance of " << refDistance
				  << " for strain calculation " << mat.str();

            Vector3D<spaceType> distance =
                atom->getPosition().getMinimalDistance(
                       neighborAtom->getPosition(), simboxSize);

            strain = (std::sqrt(distance.squaredLength()) - refDistance)/
		refDistance;

            result.push_back(std::make_tuple(atom->getIndex(),neighbor,strain));

	    CLOG(TRACE, logName_) << "Calculated strain for tuple ( " <<
		atom->getIndex().str() << "," << neighbor.str() << "): " <<
		strain;
        }
    }

    return result;
}

//------------------------------------------------------------------------------

//! Calculate strain state of each atom and return as scalar field.
std::shared_ptr<Field3D<double>> StrainTool::getStrainField (const MaterialCollection &materials,
            double refDistance) {
    Vector3D<indexType> simboxSize = simbox_->getLattice().getSize();
    Vector3D<spaceType> simboxSizePhys = simbox_->getSize();   
    strainField_ = std::make_shared<Field3D<double>>(simboxSize);
    Lattice const &lattice = simbox_->getLattice();

    for (auto atom: simbox_->getLattice().getAtomList()){

        for (auto neighbor: atom->getNeighbors()) {
            const Atom *neighborAtom = lattice(neighbor);

            if(neighborAtom != nullptr) {
                const Material &mat = materials.getByElementID(atom->getElementId(),
                                                               neighborAtom->getElementId());

                if (refDistance < 0.0) {
                    refDistance = mat.getLatticeConstant();
                    refDistance = Zincblende::latticeConstantToBondLength(
                                                                          refDistance);
                }

                CLOG(TRACE, logName_) << "Using ref distance of " << refDistance
                                      << " for strain calculation " << mat.str();

                Vector3D<spaceType> distance =
                    atom->getPosition().getMinimalDistance(
                                                           neighborAtom->getPosition(), simboxSizePhys);

                double strain = (std::sqrt(distance.squaredLength()) - refDistance)/
                    refDistance;

                // calculating scalar strain field --> simply adding up values
                strainField_->addToPoint(atom->getIndex(),strain);
            }
        }

    }

    strainFieldReady_ = true;

    return strainField_;
}

//##############################################################################

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
