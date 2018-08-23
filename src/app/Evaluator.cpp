/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "projectConfigure.h"

#include <string>
#include <algorithm>
#include <vector>
#include <tuple>
#include <cmath>
#include <map>
#include <exception>
#include <signal.h>

#include <tclap/CmdLine.h>
#include <easyloggingcpp/easylogging++.h>

#include <xtalsim.h>
// #include <misc/JsonHandler.h>
// #include <misc/XmlHandler.h>
// #include <misc/Configuration.h>
// #include <misc/ErrorCode.h>
// #include <physics/SimulationBox.h>
// #include <physics/Range3D.h>
// #include <physics/PeriodicTable.h>

#ifdef __VTK__
#include <vtkSmartPointer.h>
#include <vtkTable.h>

// #include <visualization/HistogramPlot.h>
// #include <visualization/XYPlot.h>
#endif

INITIALIZE_EASYLOGGINGPP

std::string programName;

//******************************************************************************

void bailout(const std::string &str, const int errorCode)
{
    CLOG(ERROR, programName.c_str()) << str;
    CLOG(ERROR, programName.c_str()) << "use \"" << programName << " -h\" to "
	"display the help text";
    exit(errorCode);
}

//******************************************************************************

FileHandler* getHandler(std::string fileName, XmlHandler &xmlHandler,
			JsonHandler &jsonHandler)
{   
    std::size_t found = fileName.find_last_of(".");
    std::string ending = fileName.substr(found+1);

    if (ending == "xml"){
	return &xmlHandler;
    }
    else if (ending == "json"){
	return &jsonHandler;
    }
    else{
	CLOG(ERROR, programName.c_str()) << "Do not know how to read file '" <<
	    fileName << "'";
	bailout("Only '.xml' and '.json' supported.", ErrorCode::Config);
    }

    // should be never reached
    return ( nullptr );
}

//******************************************************************************

static void errSignalHandler (int sig) {
    switch (sig) {
    case SIGINT:
        bailout("SIGINT received!", ErrorCode::SigInt);
        break;
    case SIGSEGV:
        bailout("SIGSEGV received!", ErrorCode::SigSegV);
        break;
    default:
        bailout("Unknown signal received!", ErrorCode::SigUnknown);
        break;
    }

    exit(ErrorCode::Unknown);
}

//******************************************************************************

void signalHandleInit (void) {
    signal (SIGINT, errSignalHandler);
    signal (SIGSEGV, errSignalHandler);
}

//******************************************************************************

int main (int argc, char *argv[])
{

    //-------------------------------------------------------------------------
    // configure command line reader tool TCLAP
    //-------------------------------------------------------------------------
       
    TCLAP::ValueArg<std::string> clpInFile("i", "input",
		       "FILE contains the lattice that shall be evaluated."
		       , false, "", "FILE");

    TCLAP::ValueArg<std::string> clpConfigFile("", "config",
		       "Use parameters from the configuration file in "
		       "location FILE. If parameters are given both on the "
		       "command line and in the config file the specification "
		       "on the command line dominates.", false,
					"", "FILE");

    TCLAP::ValueArg<std::string> clpPeriodicTableFile("p", "periodic-table",
		       "The information of the elements are stored in FILE.",
						      false,"", "FILE");
    
    TCLAP::ValueArg<std::string> clpMaterialFile("m", "material-file",
		       "File containing definition of required materials.",
		       false, "", "FILE");

    TCLAP::ValueArg<std::string> clpOutputPreamble("s", "output-preamble",
		       "Preamble of output files.", false,"", "PATH");
    
    TCLAP::ValueArg<std::string> clpOutputFile("o", "output-file",
		       "Specifies output file name. Only effective for gsf "
               "output format.", false,"", "FILE");

    TCLAP::ValueArg<std::string> clpOutputFormat("", "output-format",
		       "Specifies format of output files.", false,"", "FORMAT");
    
    TCLAP::ValueArg<std::string> clpLogFile("", "log",
		       "If specified log messages are also written into "
		       "this file.", false, "", "FILE");

#ifdef __VTK__
    TCLAP::MultiArg<std::string> clpEvaluationMode("", "mode",
                "Defines, which calculations should be performed. "
                "Options are state, strain, composition and render.", false, "MODE");
#else
    TCLAP::MultiArg<std::string> clpEvaluationMode("", "mode",
                "Defines, which calculations should be performed. "
                "Options are state, strain and composition.", false, "MODE");
#endif

    TCLAP::ValueArg<std::string> clpAtomState("",
                "atom-state", "Defines metric for state evaluation. "
                "Possible values are interface, exchange, unknown.",
                false, "", "STATE ...");

    TCLAP::MultiArg<std::string> clpFilterElements("f", "filter",
                "Defines a list of elements that are ignored when "
                "rendering composition data.", false, "ELEMENT");

    TCLAP::SwitchArg clpInterpolate("", "interpolate",
                "Interpolates composition traces of elements "
                "for cation and anion sublattice.", false);

    TCLAP::ValueArg<double> clpNeighborRadius("", "neighbor-radius",
		       "All atoms at most NUM Angstroem away are considered "
		       "neighbors of a single atom, "
		       "which are afterwards used to calculate the energy. "
		       "Either this or neighbor-layers has to be specified, "
		       "both is also possible.",
		       false, 0, "NUM");

    TCLAP::ValueArg<indexType> clpNeighborLayers("", "neighbor-layers",
		       "All atoms at mos NUM atomic layers away are considered "
		       "neighbors of a single atom, "
		       "which are afterwards used to calculate the energy. "
		       "Either this or neighbor-radius has to be specified, "
		       "both is also possible. In the latter case only atoms "
		       "which are also within the specified radius are picked.",
		       false, 0, "NUM");

    TCLAP::ValueArg<indexType> clpLayerIndex("", "layer",
		       "Index of layer to be evaluated.",
		       false, 0, "NUM");
    
    TCLAP::SwitchArg clpVerbose("v", "verbose", "Enable debug messages.",
				   false);

    TCLAP::MultiSwitchArg clpQuiet("q", "quiet", "First occurence disables"
                " screen output, second occurence also disables file "
                "logging.", false);

    TCLAP::SwitchArg clpQQuiet("", "qq", "Disables screen output and file "
                "logging.", false);
    
    std::vector<std::string> examples;

    examples.push_back("EXECNAME -s strain -p periodic_table.xml -i input.xml "
		       "-m material.xml\n"
		       "Calculate strain with preamble 'strain'.");

    examples.push_back("EXECNAME --config foo.xml\n"
		       "Read parameters from configuration file foo.xml.");
    
    try{

	TCLAP::CmdLine cmd("Calculate strain of structure and export data.",
		       examples, ' ', "1.0");

	cmd.add(clpInFile);
	cmd.add(clpLogFile);
	cmd.add(clpConfigFile);
	cmd.add(clpMaterialFile);
	cmd.add(clpPeriodicTableFile);
    cmd.add(clpEvaluationMode);
    cmd.add(clpAtomState);
    cmd.add(clpFilterElements);
    cmd.add(clpInterpolate);
    cmd.add(clpNeighborRadius);
    cmd.add(clpNeighborLayers);
    cmd.add(clpLayerIndex);
	cmd.add(clpVerbose);
	cmd.add(clpOutputPreamble);
    cmd.add(clpOutputFile);
    cmd.add(clpOutputFormat);
	cmd.add(clpQuiet);
    cmd.add(clpQQuiet);
	
	cmd.parse( argc, argv );

	programName = cmd.getProgramName();
	std::size_t position = programName.find_last_of("/\\");
	programName = programName.substr(position+1);	
    }
    catch (TCLAP::ArgException &e){
	std::cerr << "error: " << e.error() << " for arg " << e.argId() <<
		  std::endl;
	exit(ErrorCode::TCLAP);
    }

    //-------------------------------------------------------------------------
    // initialize logger
    //-------------------------------------------------------------------------
    
    // Load default configuration from file, adapt it according to given command
    // line parameters and set it as default for all loggers.
#if defined(EASYLOGGING_CONF_DIR)
    el::Configurations conf(EASYLOGGING_CONF_DIR "/easylogging++.conf");
#else
    el::Configurations conf("./input/easylogging++.conf");
#endif
    
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
    // initialize error signal handler
    //-------------------------------------------------------------------------
    signalHandleInit();
    
    //-------------------------------------------------------------------------
    // read config file
    //-------------------------------------------------------------------------

    // Stores the user defined parameters necessary to create the lattice.
    Configuration config;
    XmlHandler xmlHandler;
    JsonHandler jsonHandler;
    FileHandler *fileHandler;
    
    if (clpConfigFile.getValue() != ""){
	try{
	    fileHandler=getHandler(clpConfigFile.getValue(),
				   xmlHandler, jsonHandler);
	    fileHandler->load(clpConfigFile.getValue());
	    fileHandler->get(config);
	}
	catch(XmlException &e){
	    std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
	    bailout("Error while reading config file", ErrorCode::XML);
	}	
    }    

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
#if defined(EASYLOGGING_CONF_DIR)
    el::Loggers::configureFromGlobal(EASYLOGGING_CONF_DIR "/loggers.conf");
#else
    el::Loggers::configureFromGlobal("./input/loggers.conf");
#endif
    
    CLOG(TRACE, programName.c_str()) << "starting " << programName;
    
    //-------------------------------------------------------------------------
    // fill Configuration object correctly, command line parameter overrule
    // those specified in the configuration file
    //-------------------------------------------------------------------------

    // Files  
    if (clpInFile.getValue() != "")
	config.inputFileName = clpInFile.getValue();

    if (clpMaterialFile.getValue() != "")
	config.materialFile = clpMaterialFile.getValue();

    if (clpPeriodicTableFile.getValue() != "")
	config.periodicTableFileName = clpPeriodicTableFile.getValue();

    if (clpOutputPreamble.getValue() != "")
	config.outputPreamble = clpOutputPreamble.getValue();

    if (clpNeighborRadius.getValue() > 0)
	config.neighborRadius = clpNeighborRadius.getValue();

    if (clpNeighborLayers.getValue() > 0)
	config.neighborLayers = clpNeighborLayers.getValue();    
    
    CLOG(DEBUG, programName.c_str()) << "Using the following configuration: "
	"\n\n " << config.str();

    //-------------------------------------------------------------------------
    // check if all necessary parameters have been specified
    //-------------------------------------------------------------------------

    std::shared_ptr<SimulationBox> simbox;

    //if (config.outputPreamble == ""){
	//bailout("Specify preamble for output files", ErrorCode::Strain);
    //}
    
    if (config.periodicTableFileName == ""){
	bailout("No periodic table file name specified!",
		ErrorCode::PeriodicTable);
    }

    if (config.materialFile == ""){
	bailout("No material file name specified!", ErrorCode::Material);
    } 

    try{
	CLOG(DEBUG, programName.c_str()) << "loading periodic table file '" <<
	    config.periodicTableFileName << "'";
	fileHandler=getHandler(config.periodicTableFileName,
			       xmlHandler, jsonHandler);
	fileHandler->load(config.periodicTableFileName);
	fileHandler->get(PeriodicTable::getInstance());
    }catch(XmlException &e){
	std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
	bailout("Error while reading periodic table file",
		ErrorCode::XML);
    }
    catch(PeriodicTableException &e){
	std::cout << e.what() << ": "  << e.getMessage() << std::endl;
	bailout("Error while filling periodic table",
		ErrorCode::PeriodicTable);
    }
    
    if (config.inputFileName == ""){
	bailout("No input file name specified!", ErrorCode::InOutFiles);
    }
    try {

	CLOG(DEBUG, programName.c_str()) << "loading input file " <<
	    config.inputFileName;
	fileHandler=getHandler(config.inputFileName,
			       xmlHandler, jsonHandler);
	fileHandler->load(config.inputFileName);
	fileHandler->get(simbox);

	simbox->generateNeighbors(config.neighborLayers,
				  config.neighborRadius);

    }catch(XmlException &e){
	std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
	bailout("Error while reading input file", ErrorCode::XML);
    }catch(std::exception &e){
	std::cerr << "caught exception " << e.what() << std::endl;
	bailout("Unknown Error while reading input file",
		ErrorCode::Unknown);
    }
       
    //-------------------------------------------------------------------------
    // Process some input parameters
    //-------------------------------------------------------------------------
    std::transform(clpOutputFormat.getValue().begin(), clpOutputFormat.getValue().end(), clpOutputFormat.getValue().begin(), ::tolower);

    //-------------------------------------------------------------------------
    // Evaluation
    //-------------------------------------------------------------------------

    CLOG(TRACE, programName.c_str()) << "starting evaluation";
    
    try{

	MaterialCollection collection("MaterialCollection");
	
	try{
	    fileHandler=getHandler(config.materialFile,
				   xmlHandler, jsonHandler);
	    fileHandler->load(config.materialFile);
	    fileHandler->get(collection);
	} catch(XmlException &e) {
	    std::cerr << e.what() << std::endl;
	    bailout("Error while reading material file" + std::string(e.what()),
		    ErrorCode::XML);
	}

    std::vector<std::string> modes = clpEvaluationMode.getValue();

    // need to generate neighbor list before evaluating strain
    // \todo maybe call this function from strain tool
    simbox->generateNeighbors(1, 0.0, false);

    if (std::find(modes.begin(), modes.end(), "state") != modes.end()) {
        std::vector<uint32_t> stateCount;
        AtomState tmpState;

        if (clpAtomState.getValue() == "interface")
            tmpState = AtomState::ModifiedInterface;
        else if (clpAtomState.getValue() == "exchange")
            tmpState = AtomState::ModifiedExchangeReaction;
        else if (clpAtomState.getValue() == "unknown")
            tmpState = AtomState::ModifiedUnknown;

        // out of plane dimension hardcoded here
        for (uint32_t i = 0; i < simbox->getLattice().getSize()[2]; i++) {
            uint32_t tmpAtomCount = simbox->getLattice().countAtomsInLayerByState(i, 2, tmpState);
            uint32_t atomsInLayer = simbox->getLattice().getAtomsInLayer(i, 2).size();

            if (tmpAtomCount > 0)
                CLOG(INFO, programName.c_str()) << "Found " << tmpAtomCount << "/" << atomsInLayer << "atoms in layer " << i;
        }
    }

    if (std::find(modes.begin(), modes.end(), "strain") != modes.end()) {
        if (clpOutputFormat.getValue() == "gsf") {
            std::shared_ptr<Field3D<double>> strainField{};

            StrainTool st(simbox);
            strainField = st.getStrainField(collection, -1.0);

            //Field3D<double> compressedField = strainField->downsampling(2, 2, 4, DownsamplingFieldOperation::Avg);
            Field3D<double> compressedField = strainField->flatten().interpolate();

            Vector3D<indexType> res = compressedField.getSize();

            compressedField.saveToGSFFile(clpOutputFile.getValue(), clpLayerIndex.getValue()/4, 2, res[0], res[1]);
        } else {
            LayerStrainInfo layerStrain;

            bool writeToFile = true;

            // disable file output if no file preamble was given
            if (clpOutputPreamble.getValue() == "")
                writeToFile = false;

            layerStrain = simbox->calculateStrain(collection, config.outputPreamble, -1.0, writeToFile);

#ifdef __VTK__

            if (std::find(modes.begin(), modes.end(), "render") != modes.end()) {
                vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

                table = layerStrain.getTable();
                std::vector<DataSeries> dataSeries = layerStrain.get();

                XYPlot xyp = XYPlot(table, "XYPlot");

                axesLimits xAxisLimits = std::make_tuple(0, layerStrain.getLayerCount());
                axesLimits yAxisLimits = std::make_tuple(std::trunc(layerStrain.getMinimumStrain())-1, std::trunc(layerStrain.getMaximumStrain())+1);
                xyp.setAxesLimits(xAxisLimits, yAxisLimits, true, true);

                xyp.setAxesLabels("Layer", "Strain (%)");
                xyp.plot(dataSeries);

            }
#endif
        }
    }
    
        if (std::find(modes.begin(), modes.end(), "composition") != modes.end()) {
#ifdef __VTK__
            if (std::find(modes.begin(), modes.end(), "render") == modes.end()) {
#endif
                simbox->analyzeComposition(config.outputPreamble+"composition.dat");

#ifdef __VTK__
            } else {
                CompositionInfo compInfo;
                compInfo = simbox->analyzeComposition(config.outputPreamble+"composition.dat");

		XmlHandler h;
		try{
		    h.load(config.outputPreamble+"composition.xml");
		} catch (std::exception &e) {}
		
		h.set(compInfo);
		h.save(config.outputPreamble+"composition.xml");

		JsonHandler h2;
		try{
		    h2.load(config.outputPreamble+"composition.json");
		} catch (std::exception &e) {}
		
		h2.set(compInfo);
		h2.save(config.outputPreamble+"composition.json");
		
                std::vector<elementType> filterElements;

                for (auto element: clpFilterElements.getValue()) {
                    filterElements.push_back(PeriodicTable::getInstance().getBySymbol(element).id);
                }

                vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
                table = compInfo.getSparse(filterElements, clpInterpolate.getValue());
                std::vector<DataSeries> dataSeries = compInfo.get();

                XYPlot xyp = XYPlot(table, "XYPlot");
                axesLimits xAxisLimits = std::make_tuple(0, compInfo.getLayerCount());
                axesLimits yAxisLimits = std::make_tuple(0.0, compInfo.getMaxAtomCount()*1.1);
                xyp.setAxesLimits(xAxisLimits, yAxisLimits, true, true);
                xyp.setAxesLabels("Layer", "Composition (atoms/layer)");
                xyp.plot(dataSeries);
                //HistogramPlot hp = HistogramPlot(table, "HistogramPlot");
                //hp.plot();
            }
#endif
        }
	
    } catch(std::exception &e){
	std::cerr << "caught exception " << e.what() << std::endl;
	bailout("Unknown Error while creating analysing lattice: " +
		std::string(e.what()), ErrorCode::Unknown);
    }		


    //-------------------------------------------------------------------------
    
    CLOG(TRACE, programName.c_str()) << "successfully finished " << programName;
    
    return(ErrorCode::NoError);
    
}


/// Local variables:
//  mode: c++
//  indent-tabs-mode: nil
//  c-basic-offset: 4
//  tab-width: 4
//  End:
//  vim:noexpandtab:sw=4:ts=4:
