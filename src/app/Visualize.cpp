/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


/*******************************************************************************

This is simple tool to visualize lattice structures stored in .xml or .xyz files.


*******************************************************************************/

#include "projectConfigure.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdio.h>
#include <math.h>

#include <tclap/CmdLine.h>
#include <easyloggingcpp/easylogging++.h>

#include <xtalsim.h>
#include <physics/Field3D.h>
#include <misc/StrainTool.h>
// #include <misc/Color.h>
// #include <misc/XmlHandler.h>
// #include <misc/Configuration.h>
// #include <misc/ErrorCode.h>
// #include <physics/SimulationBox.h>

// #ifdef __VTK__
// #include <visualization/SimulationBoxRenderer.h>
// #endif

INITIALIZE_EASYLOGGINGPP

std::string programName;

//******************************************************************************

void bailout(const std::string &str, const int errorCode)
{
    CLOG(ERROR, programName) << str;
    exit(errorCode);
}

//******************************************************************************

//! Returns the filename including index for an animation.
//! The index is padded with heading zeros to ensure correct order.
//! The filename is generated using the following pattern: base_index.extension
//! \param basename basename of the output filename
//! \param extension Extension of the output file
//! \param frameIndex Numerical index added to the filename
//! \param lastFrame Index of last frame
//! \param steps The number of frames for the animation, determines zero padding length.
//! \return The generated filename for an animation frame
std::string generateFilenameForAnimationFrame (std::string basename,
        std::string extension, uint32_t frameIndex, uint32_t lastFrame, 
        uint32_t steps) {
    std::stringstream result;

    uint32_t reqDigits =  lastFrame > 0 ? (int) log10 ((double) lastFrame) + 1 : 1;

    result.str(std::string());
    result << basename;
    result << "_" << std::right << std::setw(reqDigits) << std::setfill('0') << frameIndex << '.' << extension;

    return result.str();
}

//******************************************************************************

int main (int argc, char *argv[])
{
    //-------------------------------------------------------------------------
    // configure command line reader tool TCLAP
    //-------------------------------------------------------------------------
    
    TCLAP::ValueArg<std::string> clpInFile("i", "input", 
            "File containing lattice structure in"
            ".xml or .xyz format", false, "", "FILE");

    TCLAP::ValueArg<std::string> clpPeriodicTableFile("p", "periodic-table",
		       "The information of the elements are stored in FILE.",
						      false,"", "FILE");

    TCLAP::ValueArg<std::string> clpConfigFile("", "config",
		       "Use parameters from the configuration file in "
		       "location FILE. If parameters are given both on the "
		       "command line and in the config file the specification "
		       "on the command line dominates.", false,
					"", "FILE");

    TCLAP::ValueArg<std::string> clpMaterialFile("", "material-file",
		       "File containing definition of required materials.",
		       false, "", "FILE");

    TCLAP::ValueArg<std::string> clpLogFile("", "log",
		       "If specified log messages are also written into "
		       "this file.", false, "", "FILE");

    TCLAP::ValueArg<indexType> clpNeighborLayers("", "neighbor-layers",
		       "All atoms at most NUM atomic layers away are considered "
		       "neighbors of a single atom, "
		       "which are afterwards used to calculate the energy. "
		       "If switch -E was activated "
		       "either this or neighbor-radius has to be specified "
		       "(both is also possible). In the latter case only atoms "
		       "which are also within the specified radius are picked.",
		       false, 0, "NUM");

    TCLAP::MultiArg<indexType> clpInterfaceIndex("", "interface-index",
                "Defines an interface layer, which can then be analyzed "
                "within a cut plane.", false, "LAYER");

    TCLAP::ValueArg<indexType> clpInterfaceHeight("", "interface-height",
                "Defines the height of interfaces, as they are generated "
                "using the InterfaceGenerator. The paramter needs to be "
                "speciefied in single atomic layers, i.e. a zincblende "
                "unit cell has a height of 4.", false, 4, "HEIGHT");

    TCLAP::SwitchArg clpRenderBonds("", "bonds",
		       "Render bonds with nearest neighbor atoms.",
		       false);

    TCLAP::SwitchArg clpVisualizeStrain("", "strain",
                "Visualize strain of individual bonds.",
                false);

    TCLAP::ValueArg<double> clpRefLatticeConstant("", "ref-latt-const",
                "Reference lattice constant for strain calculation.",
                false, -1.0, "VAL");

    TCLAP::SwitchArg clpModifiedOnly("", "modified",
                "Visualize only atoms and bonds, which were modified"
                " during optimization.",
                false);

    TCLAP::MultiArg<std::string> clpFilterModificationState("",
                "filter-modified", "Defines filters for possible modification "
                "states of atoms. Possible values: interface, exchange",
                false, "STATE ...");

    TCLAP::SwitchArg clpModifiedNegative("", "modified-negative",
                "Visualize layers, in which atoms were modified, but "
                "filter out modified atoms.", false);

    TCLAP::ValueArg<uint32_t> clpModificationIndexMax("", "max-modification-index",
                "Determines the maximum modification index, which is rendered. "
                "Only effective, if modified parameter is set.", false, 0.0, "VAL");

    TCLAP::SwitchArg clpModificationAnimation("", "animate-modification",
                "Generates individual screen captures with increasing "
                "modification index.", false);

    TCLAP::ValueArg<std::string> clpAnimationFilePrefix("", "animation-prefix",
                "Filename prefix for animation frames. File index + extension "
                "is added automatically.", false, "", "FILE");

    TCLAP::ValueArg<uint32_t> clpAnimationStep("", "animation-step",
                "Number of modifications that are processed into one frame.",
                false, 0, "VAL");

    TCLAP::MultiArg<double> clpCameraFocalPoint("", "camera-focal-point",
		       "initial focal point of camera"
				 , false, "NUM1 NUM2 NUM3");
    
    TCLAP::ValueArg<double> clpCameraDirectionX("", "camera-direction-x",
		"x component of vector from focal point to initial position"
		"of camera", false, 0, "Double");

    TCLAP::ValueArg<double> clpCameraDirectionY("", "camera-direction-y",
		"y component of vector from focal point to initial position"
		"of camera", false, 0, "Double");

    TCLAP::ValueArg<double> clpCameraDirectionZ("", "camera-direction-z",
		"z component of vector from focal point to initial position"
		"of camera", false, 0, "Double");
    
    TCLAP::ValueArg<double> clpCameraZoom("", "camera-zoom",
                "Defines initial zoom value of camera. Higher values means higher "
		"zoom. Depends on values specified by --camera-direction.",
					  false, 1.0, "DOUBLE");
    
    TCLAP::ValueArg<std::string> clpScreenshotFile("", "screenshot",
                "Filename, where screenshot of rendered window is stored.",
                false, "", "FILE");

    TCLAP::ValueArg<uint16_t> clpRotationCount("", "rotation-count",
                "Number of rotational pictures.",
                false, 0, "VAL");
    
    TCLAP::MultiSwitchArg clpQuiet("q", "quiet", "First occurence disables"
                " screen output, second occurence also disables file "
                "logging.", false);

    TCLAP::SwitchArg clpQQuiet("", "qq", "Disables screen output and file "
                "logging.", false);

    TCLAP::SwitchArg clpVerbose("v", "verbose", "Enable debug messages.",
				   false);
    

    std::vector<std::string> usageExamples;

    usageExamples.push_back("EXECNAME -i foo.xml\n"
            "");

    try {
        TCLAP::CmdLine cmd("Visualize lattice structure stored in a .xml or .xyz file."
                ""
                ""
                "", usageExamples, ' ', "0.1");

        cmd.add(clpInFile);
        cmd.add(clpConfigFile);
        cmd.add(clpPeriodicTableFile);
        cmd.add(clpMaterialFile);
        cmd.add(clpLogFile);
        cmd.add(clpVerbose);
        cmd.add(clpNeighborLayers);
        cmd.add(clpInterfaceIndex);
        cmd.add(clpInterfaceHeight);
        cmd.add(clpRenderBonds);
        cmd.add(clpVisualizeStrain);
        cmd.add(clpRefLatticeConstant);
        cmd.add(clpModifiedOnly);
        cmd.add(clpFilterModificationState);
        cmd.add(clpModifiedNegative);
        cmd.add(clpModificationIndexMax);
        cmd.add(clpModificationAnimation);
	cmd.add(clpCameraFocalPoint);
	cmd.add(clpCameraDirectionX);
	cmd.add(clpCameraDirectionY);
	cmd.add(clpCameraDirectionZ);
	cmd.add(clpCameraZoom);
        cmd.add(clpAnimationFilePrefix);
        cmd.add(clpAnimationStep);
        cmd.add(clpScreenshotFile);
        cmd.add(clpRotationCount);
        cmd.add(clpQuiet);
        cmd.add(clpQQuiet);

        cmd.parse(argc, argv);

    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() <<
            std::endl;
        exit(ErrorCode::TCLAP);
    }

    //-------------------------------------------------------------------------
    // initialize logger
    //-------------------------------------------------------------------------
    
    // Load default configuration from file, adapt it according to given command
    // line parameters and set it as default for all loggers.
    el::Configurations conf("../input/easylogging++.conf");
    
    el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
    // To reduce the maintenance effort new loggers are automatically added.
    el::Loggers::addFlag(el::LoggingFlag::CreateLoggerAutomatically);
    
    if (clpQuiet.getValue() > 0){
        conf.set(el::Level::Global, el::ConfigurationType::ToStandardOutput,
             "false");

        if (clpQuiet.getValue() > 1){
            conf.set(el::Level::Global, el::ConfigurationType::ToFile,
                 "false");
        }
    }

    if (clpQQuiet.getValue() == true){
        conf.set(el::Level::Global, el::ConfigurationType::ToStandardOutput,
             "false");
        conf.set(el::Level::Global, el::ConfigurationType::ToFile,
             "false");
    }

    el::Loggers::setDefaultConfigurations(conf, true);

    //-------------------------------------------------------------------------
    // read config file
    //-------------------------------------------------------------------------

    XmlHandler xmlHandler;
    // Stores the user defined parameters necessary to create the lattice.
    Configuration config;
    
    if (clpConfigFile.getValue() != ""){
	try{
        CLOG(DEBUG, programName) << "reading config file " <<
			 clpConfigFile.getValue().c_str();
	    xmlHandler.load(clpConfigFile.getValue());
	    xmlHandler.get(config);
	}
	catch(XmlException &e){
	    std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
	    bailout("Error while reading config file", ErrorCode::XML);
	}	
    }    

    //-------------------------------------------------------------------------
    // Attention! Copied from LatticeGenerator.cpp and modified
    // to a minimum. BEGIN
    if (clpLogFile.getValue() != "")
        config.logFileName = clpLogFile.getValue();

    if (config.logFileName != ""){
	conf.set(el::Level::Global, el::ConfigurationType::ToFile, "true");
	conf.set(el::Level::Global, el::ConfigurationType::Filename,
		 config.logFileName);
    }
    
    if ( clpVerbose.getValue() )
        config.verbose = true;

    if ( config.verbose ){
	conf.set(el::Level::Debug, el::ConfigurationType::Enabled, "true");
	conf.set(el::Level::Trace, el::ConfigurationType::Enabled, "true");
    }

    if (config.quiet){
	conf.set(el::Level::Global, el::ConfigurationType::ToStandardOutput,
		 "false");
    }

    el::Loggers::setDefaultConfigurations(conf, true);

    // Load logger specific settings.
    el::Loggers::configureFromGlobal("../input/loggers.conf");
    
    CLOG(TRACE, programName) << "starting " << programName;

    // Attention! Copied from LatticeGenerator.cpp and modified
    // to a minimum. END
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // fill Configuration object correctly, command line parameter overrule
    // those specified in the configuration file
    //-------------------------------------------------------------------------

    if (clpNeighborLayers.getValue() > 0)
        config.neighborLayers = clpNeighborLayers.getValue();    

    // Other Files
    if (clpPeriodicTableFile.getValue() != "")
	config.periodicTableFileName = clpPeriodicTableFile.getValue();

    if (clpMaterialFile.getValue() != "")
        config.materialFile = clpMaterialFile.getValue();
    
    if (clpInFile.getValue() != "")
	config.inputFileName = clpInFile.getValue();

    if (clpRenderBonds.getValue())
        config.renderBonds = true;

    if (clpVisualizeStrain.getValue())
        config.visualizeBondStrain = true;
    else
        config.visualizeBondStrain = false;

    if (clpRefLatticeConstant.isSet()) {
        config.refLatticeConstant = clpRefLatticeConstant.getValue();
    }

    if (clpModifiedOnly.getValue())
        config.modifiedOnly = true;

    if (clpFilterModificationState.getValue().size() > 0) {
        for (auto state: clpFilterModificationState.getValue())
            config.filterModificationState.push_back("filter-modified-" + state);
    }

    if (clpModifiedNegative.getValue()) {
        // implicitely also enable modified only
        config.modifiedOnly = true;
        config.modifiedNegative = true;
    }

    if (clpModificationIndexMax.getValue() > 0)
        config.modificationIndexMax = clpModificationIndexMax.getValue();

    if (clpCameraFocalPoint.getValue().size() > 0){
	if (clpCameraFocalPoint.getValue().size() < 3)
	    bailout("Error at least three coordinates required for focal point",
		    ErrorCode::Camera);

	for (auto i : {0,1,2})
	    config.cameraFocalPoint[i] = clpCameraFocalPoint.getValue()[i];
    }
    
    if ( std::abs(clpCameraDirectionX.getValue()) > 0)
	config.cameraDirection[0] = clpCameraDirectionX.getValue();

    if ( std::abs(clpCameraDirectionY.getValue()) > 0)
	config.cameraDirection[1] = clpCameraDirectionY.getValue();

    if ( std::abs(clpCameraDirectionZ.getValue()) > 0)
	config.cameraDirection[2] = clpCameraDirectionZ.getValue();
    
    if (config.cameraDirection == 0)
	config.cameraDirection = 1;

    if ((int)clpCameraZoom.getValue() != 1)
	config.cameraZoom = clpCameraZoom.getValue();

    if (clpModificationAnimation.getValue() == true)
	config.modificationAnimation = true;

    if (clpScreenshotFile.getValue() != "")
	config.screenshotFileName = clpScreenshotFile.getValue();

    if (clpAnimationStep.getValue() != 0)
	config.animationStep = clpAnimationStep.getValue();

    if (clpAnimationFilePrefix.getValue() != "")
	config.animationFilePrefix = clpAnimationFilePrefix.getValue();

    if (clpRotationCount.getValue() > 0)
	config.rotationCount = clpRotationCount.getValue();
    
    CLOG(DEBUG, programName) << "Using the following configuration: "
	"\n\n " << config.str();
    
    //-------------------------------------------------------------------------
    // prepare for visualization
    //-------------------------------------------------------------------------
    
    try {
	CLOG(DEBUG, programName) << "loading periodic table file '" <<
	    config.periodicTableFileName << "'";
        xmlHandler.load(config.periodicTableFileName);
        xmlHandler.get(PeriodicTable::getInstance());
    } catch(XmlException &e) {
        std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
        bailout("Error while reading periodic table file", ErrorCode::XML);
    } catch(PeriodicTableException &e) {
        std::cout << e.what() << ": "  << e.getMessage() << std::endl;
        bailout("Error while filling periodic table", ErrorCode::PeriodicTable);
    }

    MaterialCollection materialCollection;

    try {
	CLOG(DEBUG, programName) << "loading material collection from "
	    "file '" <<	config.materialFile << "'";
        xmlHandler.load(config.materialFile);
        xmlHandler.get(materialCollection);
    } catch(XmlException &e) {
        std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
        bailout("Error while reading periodic table file", ErrorCode::XML);
    }

    CLOG(TRACE, programName) << "creating new simulation box";
    std::shared_ptr<SimulationBox> simbox;

    XmlHandler inOutHandler;
	inOutHandler.load(clpInFile.getValue());

    CLOG(DEBUG, programName) << "loading input file " << config.inputFileName;
    
    try {
        if (config.outputFileName == "")
            config.outputFileName = config.inputFileName;

        xmlHandler.load(config.inputFileName);
        xmlHandler.get(simbox);
    } catch(XmlException &e){
        std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
        bailout("Error while reading input file", ErrorCode::XML);
    } catch(LatticeException &e){
        std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
        bailout("Error while reading input file", ErrorCode::Lattice);
    } catch(std::exception &e){
        std::cerr << "caught exception " << e.what() << std::endl;
        bailout("Unknown Error while reading input file", ErrorCode::Unknown);
    }

    std::vector<Range3D<indexType>> ifDefs;

    std::vector<indexType> interfaceLayers = clpInterfaceIndex.getValue();

    for (auto interface: interfaceLayers) {
        Range3D<indexType> ifDef;

        ifDef = Range3D<indexType>(Vector3D<indexType>(0, 0, interface),
                                Vector3D<indexType>(0, 0, interface+clpInterfaceHeight.getValue()),
                                Vector3D<bool>(false, false, true));

        ifDefs.push_back(ifDef);
    }

    // need to generate neighbor list before evaluating strain
    // \todo maybe call this function from strain tool
	simbox->generateNeighbors(1, 0.0, false);

    std::shared_ptr<Field3D<double>> strainField{};

    StrainTool st(simbox);

    if (clpVisualizeStrain.getValue() == true)
        strainField = st.getStrainField(materialCollection, -1.0);

#ifdef __VTK__
    std::shared_ptr<SimulationBoxRenderer> sbr;

    simbox->generateNeighbors(config.neighborLayers, 0.0, false);

    if (config.modificationAnimation == false) {
        sbr = std::make_shared<SimulationBoxRenderer>(config, "SimulationBoxRenderer");

        if (clpVisualizeStrain.getValue() == true) {
            sbr->renderPerspectiveWithInterfaces(*simbox, ifDefs,
						 materialCollection, config, strainField);
        } else
            sbr->renderPerspectiveWithInterfaces(*simbox, ifDefs,
						 materialCollection, config);

        if (config.screenshotFileName != "")
            sbr->saveImage(config.screenshotFileName);
    } else
    {
        uint32_t maxModificationIndex = simbox->getLattice().getMaxModificationIndex();

        std::string screenshotFilewithIndex{};
        uint32_t step = config.animationStep;

        for (uint32_t i = 1; i <= maxModificationIndex; i+=step) {
            sbr = std::make_shared<SimulationBoxRenderer>(config, "SimulationBoxRenderer");

            screenshotFilewithIndex = generateFilenameForAnimationFrame(config.animationFilePrefix,
									"png", i, maxModificationIndex, maxModificationIndex/step);

	    config.modifiedOnly = true;
            config.modificationIndexMax = i;

            if (clpVisualizeStrain.getValue() == true) {
                sbr->renderPerspectiveWithInterfaces(*simbox, ifDefs,
						     materialCollection, config, strainField, false);
            } else
                sbr->renderPerspectiveWithInterfaces(*simbox, ifDefs,
						     materialCollection, config, nullptr, false);

            sbr->saveImage(screenshotFilewithIndex);
            sbr->clear();

        }

	Vector3D<double> savDirection = config.cameraDirection;
	double radius = std::sqrt(config.cameraDirection[0]*config.cameraDirection[0]+
				  config.cameraDirection[1]*config.cameraDirection[1]);
	double angle = std::asin (config.cameraDirection[0] / radius);

	if (config.cameraDirection[1] < 0)
	    angle = 2*M_PI - angle;
	
	config.modifiedOnly = true;
	config.modificationIndexMax = maxModificationIndex;

	for (uint32_t i = 1; i <= config.rotationCount; i++) {
            sbr = std::make_shared<SimulationBoxRenderer>(config, "SimulationBoxRenderer");

	    std::stringstream fileName;

	    fileName << config.animationFilePrefix << "_rotate_" << std::right <<
		std::setw(2) << std::setfill('0') << i << ".png";

	    config.cameraDirection[0] = radius*cos(angle);
	    config.cameraDirection[1] = radius*sin(angle);
	    
            if (clpVisualizeStrain.getValue() == true) {
                sbr->renderPerspectiveWithInterfaces(*simbox, ifDefs,
						     materialCollection, config, strainField, false);
            } else
                sbr->renderPerspectiveWithInterfaces(*simbox, ifDefs,
						     materialCollection, config, nullptr, false);

            sbr->saveImage(fileName.str());
            sbr->clear();
	    angle += 2*M_PI / config.rotationCount;

	    if (angle > 2*M_PI)
		angle -= 2*M_PI;

        }
    }
#endif

    return 0;
}


