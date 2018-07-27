/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


/*******************************************************************************

This is a reference implementation which optimizes a given simulation box by
Metropolis Monte Carlo relaxations and lattice scaling.

*******************************************************************************/

#include "projectConfigure.h"

#include <string>
#include <vector>
#include <tuple>
#include <cmath>
#include <map>
#include <exception>
#include <signal.h>

#include <tclap/CmdLine.h>
#include <easyloggingcpp/easylogging++.h>

#include <xtalsim.h>

INITIALIZE_EASYLOGGINGPP

std::string programName;

//******************************************************************************

void bailout(const std::string &str, const int errorCode)
{
    CLOG(ERROR, programName) << str;
    CLOG(ERROR, programName) << "use \"" << programName << " -h\" to "
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
	CLOG(ERROR, programName) << "Do not know how to read file '" <<
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
		       "File containing already "
		       "lattice information in the XML format. If not provided "
		       "or if file does not contain a lattice a new one is "
		       "created", false, "", "FILE");

    TCLAP::ValueArg<std::string> clpOutFile("o", "output",
		       "The generated lattice "
		       "is written to this file. If not specified the "
		       "information are exported to the input file if "
		       "specified. If not the program exits.", false,
					"", "FILE");

    TCLAP::ValueArg<std::string> clpConfigFile("", "config",
		       "Use parameters from the configuration file in "
		       "location FILE. If parameters are given both on the "
		       "command line and in the config file the specification "
		       "on the command line dominates.", false,
					"", "FILE");

    TCLAP::ValueArg<std::string> clpTersoffFile("t", "tersoff",
		       "FILE that contains the tersoff parameters, which are "
		       "required to calculate the energy.", false,"", "FILE");

    TCLAP::ValueArg<std::string> clpPeriodicTableFile("p", "periodic-table",
		       "The information of the elements are stored in FILE.",
						      false,"", "FILE");

    TCLAP::ValueArg<std::string> clpJournalPreamble("j", "journal-preamble",
		       "Preamble of journal files.", false,"", "PATH");
    
    TCLAP::ValueArg<std::string> clpLogFile("", "log",
		       "If specified log messages are also written into "
		       "this file.", false, "", "FILE");

    TCLAP::ValueArg<std::string> clpXYZFile("", "xyz",
	               "FILE xyz data shall be written to. ",
		       false, "", "FILE");
    
    TCLAP::ValueArg<double> clpTemperature("k", "temperature",
		       "Internal lattice temperature in Kelvin.", false,
					   0, "NUM");

    TCLAP::ValueArg<int> clpStartIndex("a", "start-index",
		       "Start index of range that shall be optimized. If not "
		       "specified the first atomic layer is chosen. Valid "
		       "parameter range is [0,#Atomic layers -1].", false,
					   -1, "NUM");

    TCLAP::ValueArg<int> clpStopIndex("z", "stop-index",
		       "Stop index of range that shall be optimized. If not "
		       "specified the top atomic layer is chosen. Valid "
		       "parameter range is [0,#Atomic layers -1].", false,
					   -1, "NUM");

    TCLAP::ValueArg<int> clpRunCount("r", "runs",
		       "Run NUM times a random optimization step and run it. "
		       "This parameter must be defined.", false, -1, "NUM");

    TCLAP::ValueArg<int> clpCheckCount("", "check-runs",
		       "After NUM optimization steps energy in dynamic "
		       "optimization is checked.", false, 0, "NUM");

    TCLAP::ValueArg<double> clpEnergyDrop("", "energy-drop",
		       "If an energy gain in dynamic optimization lower than "
		       "NUM is achieved the parameters are adapted.",
					false, 0, "NUM");

    TCLAP::ValueArg<double> clpReduction("", "reduction",
		       "Parameters in dynamic optimization are divided by NUM "
		       "when they have to be adapted.", false, 0, "NUM");

    TCLAP::ValueArg<double> clpMinEnergy("", "min-energy",
		       "Minimal energy gain that has to be achieved in dynamic "
		       "optimization. Is used as stop condition.", false, 0,
					 "NUM");
    
    TCLAP::ValueArg<double> clpMmcProbability("m", "mmc-probability",
		       "Probability to trigger Metropolis Monte Carlo "
		       "(MMC) Relaxation. This value is set in relation to "
		       "the probabilites of other optimization methods. If "
		       "not specified or set to 0 the action is not carried out"
		       " at all.", false, -1, "NUM");

    TCLAP::ValueArg<int> clpMmcRunCount("", "mmc-run-count",
					"How many Metropolis Monte Carlo "
					"(MMC) Relaxation steps shall be carried out at once.",
					false, 1, "NUM");
    
    TCLAP::ValueArg<double> clpMinDisplacement("", "min-displacement",
		       "Minimal displacement in the Metropolis Monte Carlo "
		       "Relaxation. Is used as stop condition in the dynamic "
		       "optimization.", false, 0, "NUM");
    
    TCLAP::ValueArg<double> clpMaxDisplacement("d", "max-displacement",
		       "Maximal displacement in the Metropolis Monte Carlo "
		       "Relaxation. If MMC relaxation is desired this parameter"
		       " is mandatory.", false, -1, "NUM");

    TCLAP::ValueArg<double> clpScalingProbability("s", "scaling-probability",
		       "Probability to trigger scaling between atomic layers "
		       "during optimization. This value is set in relation to "
		       "the probabilites of other optimization methods. If "
		       "not specified or set to 0 the action is not carried out"
		       " at all.", false, -1, "NUM");

    TCLAP::ValueArg<double> clpMinScaling("", "min-scaling",
		       "Minimal range of the scaling factor. Is used as stop "
		       "condition in the dynamic optimization.",
		       false, 0, "NUM");
    
    TCLAP::ValueArg<double> clpMaxScaling("f", "max-scaling",
		       "Range of the scaling factor around 1 for scaling the "
		       "distance between atomic layers. If scaling is desired "
		       "this parameter is mandatory.", false, -1, "NUM");

    TCLAP::ValueArg<int> clpMaxThreadCount("", "max-thread-count",
				       "Maximal number of parallel threads. ",
				       false, 1, "NUM");
    
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

    TCLAP::MultiSwitchArg clpQuiet("q", "quiet", "First occurence disables"
                " screen output, second occurence also disables file "
                "logging.", false);

    TCLAP::SwitchArg clpQQuiet("", "qq", "Disables screen output and file "
                "logging.", false);
    
    TCLAP::SwitchArg clpVerbose("v", "verbose", "Enable debug messages.",
				   false);

    TCLAP::SwitchArg clpPrint("", "print-statistics", "If set a complete"
		       "statistic of the optimization steps is printed.",
				   false);    

    TCLAP::SwitchArg clpStaticOptimization("", "staticOptimization",
		       "If set use static optimization method.", false);

    TCLAP::SwitchArg clpDynamicOptimization("", "dynamicOptimization",
		       "If set use dynamic optimization method.", false);
    
    std::vector<std::string> examples;

    examples.push_back("EXECNAME -m 100 -d 3e-1 -s 1 -f 0.2 -r 200 -i "
		       "input.xml -o output.xml -t tersoff.xml -p "
		       "periodic_table.xml\n"
		       "Trigger MMC relaxation with a maximal displacement of "
		       " 0.3 100 times more often than scaling between atomic "
		       "layers with a scaling factor in the range [0.9,1.1]. "
		       "200 runs are carried out.");

    examples.push_back("EXECNAME --config foo.xml -d 5e-1\n"
		       "Read parameters from configuration file foo.xml "
		       "and use maximal MMC displacement of 0.5.");
    
    try{

	TCLAP::CmdLine cmd("Optimize existing lattice using Metropolis Monte "
			   "Carlo simulations and/or lattice scaling.",
			   examples, ' ', "1.0");

	cmd.add(clpInFile);
	cmd.add(clpOutFile);
	cmd.add(clpConfigFile);
	cmd.add(clpStartIndex);
	cmd.add(clpStopIndex);
	cmd.add(clpLogFile);
	cmd.add(clpPeriodicTableFile);
	cmd.add(clpTersoffFile);
	cmd.add(clpMmcProbability);
	cmd.add(clpMmcRunCount);
	cmd.add(clpMinDisplacement);
	cmd.add(clpMaxDisplacement);
	cmd.add(clpScalingProbability);
	cmd.add(clpMinScaling);
	cmd.add(clpMaxScaling);
	cmd.add(clpMaxThreadCount);
	cmd.add(clpVerbose);
	cmd.add(clpRunCount);
	cmd.add(clpCheckCount);
	cmd.add(clpEnergyDrop);
	cmd.add(clpReduction);
	cmd.add(clpMinEnergy);
	cmd.add(clpPrint);
	cmd.add(clpNeighborRadius);
	cmd.add(clpNeighborLayers);
	cmd.add(clpTemperature);
	cmd.add(clpXYZFile);
	cmd.add(clpJournalPreamble);
	cmd.add(clpStaticOptimization);
	cmd.add(clpDynamicOptimization);
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
    el::Configurations conf("../input/easylogging++.conf");
    
    //el::Loggers::configureFromGlobal("../easyloggingcpp/easylogging++Global.conf");
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
    el::Loggers::configureFromGlobal("../input/loggers.conf");
    
    CLOG(TRACE, programName) << "starting " << programName;
    
    //-------------------------------------------------------------------------
    // fill Configuration object correctly, command line parameter overrule
    // those specified in the configuration file
    //-------------------------------------------------------------------------

    // Files  
    if (clpInFile.getValue() != "")
	config.inputFileName = clpInFile.getValue();

    if (clpOutFile.getValue() != "")
	config.outputFileName = clpOutFile.getValue();

    if (clpPeriodicTableFile.getValue() != "")
	config.periodicTableFileName = clpPeriodicTableFile.getValue();

    if (clpXYZFile.getValue() != "")
        config.xyzFileName = clpXYZFile.getValue();

    if (clpJournalPreamble.getValue() != "")
	config.journalPreamble = clpJournalPreamble.getValue();

    if (clpTersoffFile.getValue() != "")
        config.tersoffFileName = clpTersoffFile.getValue();
    
    // Optimization

    if (clpStartIndex.getValue() >= 0)
	config.startIndex = clpStartIndex.getValue();

    if (clpStopIndex.getValue() >= 0)
	config.stopIndex = clpStopIndex.getValue();

    if (clpMmcProbability.getValue() >= 0)
	config.mmcProbability = clpMmcProbability.getValue();

    if (clpMinDisplacement.getValue() > 0)
	config.minDisplacement = clpMinDisplacement.getValue();
    
    if (clpMaxDisplacement.getValue() >= 0)
	config.maxDisplacement = clpMaxDisplacement.getValue();
    
    if (clpScalingProbability.getValue() >= 0)
	config.scalingProbability = clpScalingProbability.getValue();

    if (clpMinScaling.getValue() > 0)
	config.minScaling = clpMinScaling.getValue();
    
    if (clpMaxScaling.getValue() >= 0)
	config.maxScaling = clpMaxScaling.getValue();

    if (clpRunCount.getValue() >= 0)
	config.runCount = clpRunCount.getValue();

    if (clpCheckCount.getValue() > 0)
	config.checkCount = clpCheckCount.getValue();    

    if (clpEnergyDrop.getValue() > 0)
	config.energyDropFactor = clpEnergyDrop.getValue();    

    if (clpReduction.getValue() > 0)
	config.reductionFactor = clpReduction.getValue();    

    if (clpMinEnergy.getValue() > 0)
	config.minEnergy = clpMinEnergy.getValue();    
    
    if (clpNeighborRadius.getValue() > 0)
	config.neighborRadius = clpNeighborRadius.getValue();

    if (clpNeighborLayers.getValue() > 0)
	config.neighborLayers = clpNeighborLayers.getValue();    

    if (clpTemperature.getValue() > 0)
	config.latticeTemperature = clpTemperature.getValue();

    if (clpMaxThreadCount.getValue() > 1)
	config.maxThreadCount = clpMaxThreadCount.getValue();

    if (clpMmcRunCount.getValue() > 1)
	config.mmcRunCount = clpMmcRunCount.getValue();
    
    // Switches    
    if ( clpPrint.getValue() == true )
	config.print = true;

    if ( clpStaticOptimization.getValue() == true )
	config.staticOptimization = true;

    if ( clpDynamicOptimization.getValue() == true )
	config.dynamicOptimization = true;
    
    CLOG(DEBUG, programName) << "Using the following configuration: "
	"\n\n " << config.str();

    //-------------------------------------------------------------------------
    // check if all necessary parameters have been specified
    //-------------------------------------------------------------------------

    std::shared_ptr<SimulationBox> simbox;

    if (config.runCount < 0){
	bailout("Invalid number of optimization runs specified.",
		ErrorCode::Runs);
    }
    
    if (config.tersoffFileName == ""){
	bailout("No Tersoff parameter file specified.",
		ErrorCode::Tersoff);
    }

    if (! ((config.mmcProbability > 0 ) || (config.scalingProbability > 0)) ){
	bailout("Specify at least one optimization action!",
		ErrorCode::OptimizationAction);
    }

    if ( (config.mmcProbability > 0 ) && (!(config.maxDisplacement > 0)) ){
	bailout("Specify maximal displacement for Metropolis Monte Carlo "
		"relaxation.", ErrorCode::OptimizationAction);
    }

    if ( (config.scalingProbability > 0 ) && (!(config.maxScaling > 0)) ){
	bailout("Specify maximal swing for scaling optimization.",
		ErrorCode::OptimizationAction);
    }

    // if neither radius nor atomic layer count specified.
    if ( ! ( (config.neighborRadius > 0) || (config.neighborLayers > 0) )  ){
	config.neighborLayers=1;
	CLOG(WARNING, programName) << "No neighbor determination method"
	    " specified! Falling back to default of layer count = 1";
    }

    try{
	CLOG(DEBUG, programName) << "loading periodic table file '" <<
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

	if (config.outputFileName == "")
	    config.outputFileName = config.inputFileName;

	CLOG(DEBUG, programName) << "loading input file '" <<
	    config.inputFileName << "'";
	fileHandler=getHandler(config.inputFileName,
			       xmlHandler, jsonHandler);
	fileHandler->load(config.inputFileName);
	fileHandler->get(simbox);

    }catch(XmlException &e){
	std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
	bailout("Error while reading input file", ErrorCode::XML);
    }
    catch(LatticeException &e){
	std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
	bailout("Error while reading input file", ErrorCode::Lattice);
    }catch(std::exception &e){
	std::cerr << "caught exception " << e.what() << std::endl;
	bailout("Unknown Error while reading input file", 
		ErrorCode::Unknown);
    }

    if (config.latticeTemperature > 0)
	simbox->setLatticeTemperature(config.latticeTemperature);
   
    //-------------------------------------------------------------------------
    // Optimization
    //-------------------------------------------------------------------------
   
    try{

	double energyBefore, energyAfter;
	Range3D<indexType> range;
	uint8_t outOfPlaneDimension = simbox->getOutOfPlaneDimension();

	if (config.startIndex >= 0){
	    range.start[outOfPlaneDimension] = config.startIndex;
	    range.apply[outOfPlaneDimension] = true;
	}
	    
	if (config.stopIndex >= 0){
	    range.stop[outOfPlaneDimension] = config.stopIndex;
	    range.apply[outOfPlaneDimension] = true;
	}
	
	OptimizationParameter optParam;
	Optimization opt;
	
	TersoffPotential potential;

	optParam.set(config);

	if (config.mmcProbability > 0)
	    if (config.maxThreadCount > 1)
		opt.registerAction("MMC relaxation",
				   &SimulationBox::mmcRelaxParallel, config.mmcProbability);
	    else
		opt.registerAction("MMC relaxation",
				   &SimulationBox::mmcRelax, config.mmcProbability);
	
	if (config.scalingProbability > 0)
	    opt.registerAction("lattice scaling",
			       &SimulationBox::scale,config.scalingProbability);

	fileHandler=getHandler(config.tersoffFileName,
			       xmlHandler, jsonHandler);
	fileHandler->load(config.tersoffFileName);
	fileHandler->get(potential);

	simbox->generateNeighbors(config.neighborLayers,
				  config.neighborRadius);
	
	energyBefore = simbox->getEnergy(potential, range);

	if (config.staticOptimization)
	    opt.runStatic(*simbox, range, potential, optParam,
			 config.journalPreamble);
	
	if (config.dynamicOptimization)
	    opt.runDynamic(*simbox, range, potential, optParam,
			 config.journalPreamble);
	
	energyAfter = simbox->getEnergy(potential, range);
	CLOG(INFO, programName) << "Energy before optimization: " <<
	    energyBefore;
	CLOG(INFO, programName) << "\t " <<  energyBefore/simbox->getAtomCount() <<
	    " eV/atom";
	CLOG(INFO, programName) << "Energy after optimization: " <<
	    energyAfter;
	CLOG(INFO, programName) << "\t " <<  energyAfter/simbox->getAtomCount() <<
	    " eV/atom";

	if (config.print)
	    opt.printStatistic();

	fileHandler = getHandler(config.outputFileName, xmlHandler,jsonHandler);
	fileHandler->clear();
        fileHandler->set(*simbox);
        fileHandler->save(config.outputFileName);

	if (config.xyzFileName != "")
	    simbox->writeToXYZ(config.xyzFileName);
	
    }catch(LatticeException &e){
	std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
	bailout("Error during optimization.", ErrorCode::Lattice);
    }catch(std::exception &e){
	std::cerr << "caught exception " << e.what() << std::endl;
	bailout("Unknown Error during optimization.", ErrorCode::Unknown);
    }		

    
    //-------------------------------------------------------------------------

    CLOG(TRACE, programName) << "successfully finished " << programName;
    
    return(ErrorCode::NoError);
    
}


/// Local variables:
//  mode: c++
//  indent-tabs-mode: nil
//  c-basic-offset: 4
//  tab-width: 4
//  End:
//  vim:noexpandtab:sw=4:ts=4:
