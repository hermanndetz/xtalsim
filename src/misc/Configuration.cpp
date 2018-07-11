/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#include "Configuration.h"

//------------------------------------------------------------------------------

Configuration::Configuration(): 
    materialName(""), materialFile(""), materialSearchName(""),
    periodicTableFileName(""), c11(0), c12(0), c44(0),
    latticeType(""),
    latticeConstant(0), size(0,0,0), growthDimension(4),
    neighborRadius(0), neighborLayers(0), latticeTemperature(0),
    startIndex(-1), stopIndex(-1), mmcProbability(-1),
    mmcRunCount(1), minDisplacement(0), maxDisplacement(-1),
    scalingProbability(-1), minScaling(0), maxScaling(-1), runCount(-1),
    checkCount(1), energyDropFactor(2), reductionFactor(2), minEnergy(0),
    anionPassivationProbability(1.0), requireMatchingCation(false), maxThreadCount(1),
    interfaceAtoms(-1), interfacePositions(-1), inputFileName(""), outputFileName(""),
    logFileName(""), tersoffFileName(""), outputPreamble(""), extend(false),
    calculateEnergy(false), print(false), verbose(false), quiet(false),
    renderBonds(false), visualizeBondStrain(false), refLatticeConstant(0.0),
    modifiedOnly(false), filterModificationState(), modifiedNegative(false), 
    modificationIndexMax(0), cameraFocalPoint(0,0,0), cameraDirection(0,0,0), 
    cameraZoom(1.0), modificationAnimation(false), animationStep(1), 
    animationFilePrefix(""), screenshotFileName(""), rotationCount(0),
    dynamicOptimization(false), staticOptimization(false)
{
    //nothing to be done
}

//------------------------------------------------------------------------------

Configuration::Configuration(const Configuration &config) {
    this->operator=(config);
}

//------------------------------------------------------------------------------

Configuration::~Configuration() {}

//------------------------------------------------------------------------------

void Configuration::operator=(const Configuration &config) {
    materialName = config.materialName;
    materialFile = config.materialFile;
    materialSearchName = config.materialSearchName;
    periodicTableFileName = config.periodicTableFileName;
    c11 = config.c11;
    c12 = config.c12;
    c44 = config.c44;
    latticeType = config.latticeType;
    latticeConstant = config.latticeConstant;
    size = config.size;
    growthDimension = config.growthDimension;
    neighborRadius = config.neighborRadius;
    neighborLayers = config.neighborLayers;
    startIndex = config.startIndex;
    stopIndex = config.stopIndex;
    mmcProbability = config.mmcProbability;
    mmcRunCount = config.mmcRunCount;
    minDisplacement = config.minDisplacement;
    maxDisplacement = config.maxDisplacement;
    scalingProbability = config.scalingProbability;
    minScaling = config.minScaling;
    maxScaling = config.maxScaling;
    runCount = config.runCount;
    checkCount = config.checkCount;
    energyDropFactor = config.energyDropFactor;
    reductionFactor = config.reductionFactor;
    minEnergy = config.minEnergy;
    anionPassivationProbability = config.anionPassivationProbability;
    requireMatchingCation = config.requireMatchingCation;
    maxThreadCount = config.maxThreadCount;
    interfaceAtoms = config.interfaceAtoms;
    interfacePositions = config.interfacePositions;
    inputFileName = config.inputFileName;
    outputFileName = config.outputFileName;
    logFileName = config.logFileName;
    tersoffFileName = config.tersoffFileName;
    outputPreamble = config.outputPreamble;
    extend = config.extend;
    calculateEnergy = config.calculateEnergy;
    print = config.print;
    verbose = config.verbose;
    quiet = config.quiet;
    staticOptimization = config.staticOptimization;
    dynamicOptimization = config.dynamicOptimization;
    renderBonds = config.renderBonds;
    visualizeBondStrain = config.visualizeBondStrain;
    refLatticeConstant = config.refLatticeConstant;
    modifiedOnly = config.modifiedOnly;
    filterModificationState = config.filterModificationState;
    modifiedNegative = config.modifiedNegative;
    modificationIndexMax = config.modificationIndexMax;
    cameraFocalPoint = config.cameraFocalPoint;
    cameraDirection = config.cameraDirection;
    cameraZoom = config.cameraZoom;
    modificationAnimation = config.modificationAnimation;
    animationStep = config.animationStep;
    animationFilePrefix = config.animationFilePrefix;
    screenshotFileName = config.screenshotFileName;
    rotationCount = config.rotationCount;
    latticeTemperature = config.latticeTemperature;
}

//------------------------------------------------------------------------------

//! Returns filter definition for atom states
std::bitset<AtomState::StateCount> Configuration::getModificationState(void) const {
    std::bitset<AtomState::StateCount> result{};

    for (auto filter: filterModificationState) {
        if (filter == "filter-modified-interface")
            result.set(AtomState::ModifiedInterface);
        else if (filter == "filter-modified-exchange")
            result.set(AtomState::ModifiedExchangeReaction);
    }

    return result;
}

//------------------------------------------------------------------------------

//! \return Member variables with description as string.
const std::string Configuration::str(void) const
{
    std::ostringstream msg;
    std::string name;
    double share;
   
    msg << "Material Name: '" << materialName << "'" << std::endl;
    msg << "Elastic constants: c11=" << c11 << "; c12=" << c12 <<
        "; c44=" << c44 << std::endl;
    msg << "Lattice type: '" << latticeType << "'" << std::endl;
    msg << "Lattice constant: " << latticeConstant << std::endl;
    
    for (auto el: cations){
        std::tie (name, share) = el;
        msg << "Cation " << name << " with share " << share << std::endl;
    }

    for (auto el: anions){
        std::tie (name, share) = el;
        msg << "Anion " << name << " with share " << share << std::endl;
    }
    
    msg << "Material File: '" << materialFile << "'" << std::endl;
    msg << "Material search name: '" << materialSearchName << "'" << std::endl;
    msg << "---------------------------------------" << std::endl;
    msg << "Size: " << size.str() << std::endl;
    msg << "Growth dimension: " << (int)growthDimension << std::endl;
    msg << "---------------------------------------" << std::endl;
    msg << "Input file name: '" << inputFileName << "'" << std::endl;
    msg << "Output file name: '" << outputFileName << "'" << std::endl;
    msg << "Periodic Table file name: '" << periodicTableFileName << "'" <<
	std::endl;
    msg << "Logging file name: '" << logFileName << "'" << std::endl;
    msg << "Tersoff file name: '" << tersoffFileName << "'" << std::endl;
    msg << "XYZ file name: '" << xyzFileName << "'" << std::endl;
    msg << "Output preamble: '" << outputPreamble << "'" << std::endl;
    msg << "Journal preamble: '" << journalPreamble << "'" << std::endl;
    msg << "---------------------------------------" << std::endl;
    msg << "Neighbor radius: " << neighborRadius << std::endl;
    msg << "Neighbor layers: " << neighborLayers << std::endl;
    msg << "Start index: " << startIndex << std::endl;
    msg << "Stop index: " << stopIndex << std::endl;
    msg << "MMC probability: " << mmcProbability << std::endl;
    msg << "MMC run count: " << mmcRunCount << std::endl;
    msg << "min. displacement: " << minDisplacement << std::endl;
    msg << "max. displacement: " << maxDisplacement << std::endl;
    msg << "Scaling probability: " << scalingProbability << std::endl;
    msg << "min. scaling: " << minScaling << std::endl;
    msg << "max. scaling: " << maxScaling << std::endl;
    msg << "run count: " << runCount << std::endl;
    msg << "check count: " << checkCount << std::endl;
    msg << "energy drop factor: " << energyDropFactor << std::endl;
    msg << "reduction factor: " << reductionFactor << std::endl;
    msg << "minimal energy: " << minEnergy << std::endl;
    msg << "anion passivation probability: " << anionPassivationProbability << std::endl;
    msg << "require matching cation: " << requireMatchingCation << std::endl;
    msg << "maximal parallel threads: " << maxThreadCount << std::endl;
    msg << "---------------------------------------" << std::endl;
    msg << "interface atoms: " << interfaceAtoms << std::endl;
    msg << "interface positions: " << interfacePositions << std::endl;
    msg << "---------------------------------------" << std::endl;
    msg << "Render bonds: " << renderBonds << std::endl;
    msg << "Visualize strain: " << visualizeBondStrain << std::endl;
    msg << "Modified only: " << modifiedOnly << std::endl;
    msg << "Filter modification state: ";

    for (auto state: filterModificationState) {
        msg << state;

        if (state != filterModificationState.back())
            msg << "|";
    }
        
    msg << std::endl;
    msg << "Modified negative: " << modifiedNegative << std::endl;
    msg << "Modification index max: " << modificationIndexMax << std::endl;
    msg << "initial camera focal Point: " << cameraFocalPoint << std::endl;
    msg << "initial camera direction: " << cameraDirection << std::endl;
    msg << "initial camera zoom: " << cameraZoom << std::endl;
    msg << "modification animation: " << modificationAnimation << std::endl;
    msg << "animation step: " << animationStep << std::endl;
    msg << "animation file prefix: " << animationFilePrefix << std::endl;
    msg << "screenshot file name: " << screenshotFileName << std::endl;
    msg << "rotation count: " << rotationCount << std::endl;
    msg << "---------------------------------------" << std::endl;
    msg << "Reference lattice constant: " << refLatticeConstant << std::endl;
    msg << "Lattice temperature[K]: " << latticeTemperature << std::endl;
    msg << "---------------------------------------" << std::endl;
    msg << "Extend: " << extend << std::endl;    
    msg << "Calculate Energy: " << calculateEnergy << std::endl;
    msg << "Quiet: " << quiet << std::endl;
    msg << "Verbose: " << verbose << std::endl;
    msg << "Print statistics: " << print << std::endl;
    msg << "static optimization: " << staticOptimization << std::endl;
    msg << "dynamic optimization: " << dynamicOptimization << std::endl
        << std::endl;
    return msg.str();

}

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
