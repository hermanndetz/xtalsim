/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


/*******************************************************************************

This is a reference implementation which creates energetically optimal
interfaces.

*******************************************************************************/

#include <string>
#include <vector>
#include <tuple>
#include <cmath>
#include <map>
#include <exception>
#include <stdio.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

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
    
    TCLAP::ValueArg<std::string> clpLogFile("", "log",
		       "If specified log messages are also written into "
		       "this file.", false, "", "FILE");

    TCLAP::ValueArg<std::string> clpXYZFile("", "xyz",
	               "FILE xyz data shall be written to. ",
		       false, "", "FILE");

    TCLAP::ValueArg<std::string> clpJournalPreamble("j", "journal-preamble",
		       "Preamble of journal files.", false,"", "PATH");

    TCLAP::MultiArg<std::string> clpCation("c", "cations", "Defines the cations"
		       " of the introduced material. Multiple cations are "
		       "allowed. Each cation has to be identified by its "
		       "symbol in the periodic table, e.g. Ga, In.",
					   false, "SYMBOL ...");

    TCLAP::MultiArg<double> clpCationShare("", "share-cations",
		       "Determines the share of the single cations in "
		       "the material. If not specified an equal share is "
		       "taken. At least one share for each cation has to be "
		       "defined, otherwise the program exits with an error.",
					   false, "SHARE ...");
    
    TCLAP::MultiArg<std::string> clpAnion("a", "anions", "Defines the anions"
		       " of the introduced material. Multiple anions are "
		       "allowed. Each anion has to be identified by its "
		       "symbol in the periodic table, e.g. As, Sb.",
					   false, "SYMBOL ...");

    TCLAP::MultiArg<double> clpAnionShare("", "share-anions",
		       "Determines the share of the single anions in "
		       "the material. If not specified an equal share is "
		       "taken. At least one share for each anion has to be "
		       "defined, otherwise the program exits with an error.",
					   false, "SHARE ...");
    
    TCLAP::ValueArg<std::string> clpMaterialFile("m", "material-file",
		       "File containing definition of required materials.",
		       false, "", "FILE");

    TCLAP::ValueArg<std::string> clpMaterialName("n", "material-name",
		       "Name of material that shall be used. Only valid in "
		       "combination with --material-file. A specification by "
		       "naming cations and anions overrules this method.",
						 false, "", "NAME");

    TCLAP::ValueArg<double> clpLatticeConstant("l", "lattice-constant",
		       "Lattice constant of the introduced "
		       "material.", false, 0, "VAL");
    
    TCLAP::ValueArg<int> clpStartIndex("x", "start-index",
		       "Start index of range that shall be optimized. If not "
		       "specified the first atomic layer is chosen. Valid "
		       "parameter range is [0,#Atomic layers -1].", false,
					   -1, "NUM");

    TCLAP::ValueArg<int> clpStopIndex("y", "stop-index",
		       "Stop index of range that shall be optimized. If not "
		       "specified the top atomic layer is chosen. Valid "
		       "parameter range is [0,#Atomic layers -1].", false,
					   -1, "NUM");

    TCLAP::ValueArg<double> clpTemperature("k", "temperature",
		       "Internal lattice temperature in Kelvin.", false,
					   0, "NUM");
    
    TCLAP::ValueArg<int> clpRunCount("r", "runs",
		       "Run NUM Metropolis Monte Carlo simulations in a certain"
		       " space around every newly introduced atom. "
		       "This parameter is mandatory.", false, -1, "NUM");
    
    TCLAP::ValueArg<double> clpMaxDisplacement("d", "max-displacement",
		       "Maximal displacement in the Metropolis Monte Carlo "
		       "Relaxation. If MMC relaxation is desired this parameter"
		       " is mandatory.", false, -1, "NUM");

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

    TCLAP::ValueArg<std::string> clpMetric("", "metric", "Metric to use "
                "to determine optimum position of atoms. Possible values "
                "are energy, bondstrain, growthstrain and distancestrain. "
                "If no value is given, energy is taken as default.",
                false, "energy", "STRING");

    TCLAP::ValueArg<int> clpInterfaceAtoms("", "interface-atoms",
		       "Interface is built out of NUM cations. The amount of "
		       "anions may vary.", false, 0, "NUM");

    TCLAP::ValueArg<int> clpInterfacePositions("", "interface-positions",
		       "For each cation NUM different positions are tested "
		       "on the interface.", false, 0, "NUM");

    TCLAP::MultiArg<std::string> clpExchangeReaction("", "exchange-reaction",
                "Defines an exchange reaction between two materials "
                "at an interface. The general format is A_B_x, with A and B "
		"being material names and x the replacement probability. To replace "
                "e.g. InAs with GaSb at a probability of 10% use InAs_GaSb_0.10 .",
                false, "STRING");
    
    TCLAP::SwitchArg clpVerbose("v", "verbose", "Enable debug messages.",
				   false);

    TCLAP::MultiSwitchArg clpQuiet("q", "quiet", "First occurence disables"
                " screen output, second occurence also disables file "
                "logging.", false);

    TCLAP::SwitchArg clpQQuiet("", "qq", "Disables screen output and file "
                "logging.", false);
    
    std::vector<std::string> examples;

    examples.push_back("EXECNAME -x 6 -y 10 -d 0.1 -c Ga -a As -n GaAs -l "
		       "5.6535 --interface-atoms 3 --interface-positions 3 "
		       "-t tersoff.xml -p periodicTable.xml -i input.xml\n"
		       "Creates a GaAs interface consisting of 3 atoms, where "
		       "3 different positions for each atoms are checked, "
		       "between layers 6 and 10");

    examples.push_back("EXECNAME --config foo.xml -d 0.25\n"
		       "Read parameters from configuration file foo.xml and "
		       "use a maximal displacement of 0.25.");
    
    try{

	TCLAP::CmdLine cmd("Introduces energetically optimal interfaces into "
			   "an already existing atomic lattice. Only one "
			   "interface can be created at a time; for multiple "
			   "ones the program has to be rerun.",
			   examples, ' ', "1.0");

	cmd.add(clpInFile);
	cmd.add(clpOutFile);
	cmd.add(clpConfigFile);
	cmd.add(clpStartIndex);
	cmd.add(clpStopIndex);
	cmd.add(clpLogFile);
	cmd.add(clpPeriodicTableFile);
	cmd.add(clpMaxDisplacement);
	cmd.add(clpVerbose);
	cmd.add(clpRunCount);
	cmd.add(clpNeighborRadius);
	cmd.add(clpNeighborLayers);
	cmd.add(clpCation);
	cmd.add(clpCationShare);
	cmd.add(clpAnion);
	cmd.add(clpAnionShare);
	cmd.add(clpMaterialFile);
	cmd.add(clpMaterialName);
	cmd.add(clpLatticeConstant);
	cmd.add(clpTemperature);
	cmd.add(clpXYZFile);
	cmd.add(clpJournalPreamble);
	cmd.add(clpTersoffFile);
    cmd.add(clpMetric);
	cmd.add(clpInterfaceAtoms);
	cmd.add(clpInterfacePositions);
    cmd.add(clpExchangeReaction);
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

#if defined(EASYLOGGING_CONF_DIR)
    el::Loggers::configureFromGlobal(EASYLOGGING_CONF_DIR "/loggers.conf");
#else
    el::Loggers::configureFromGlobal("./input/loggers.conf");
#endif
    
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

    if (clpMaterialFile.getValue() != "")
	config.materialFile = clpMaterialFile.getValue();
    
    // Element and Material
    
    if (clpLatticeConstant.getValue() > 0)
        config.latticeConstant = clpLatticeConstant.getValue();

    if ((clpCation.getValue().size() > 0) &&
        (clpAnion.getValue().size() > 0)) {

        config.cations.clear();
        config.anions.clear();
	
        if (clpCation.getValue().size() > clpCationShare.getValue().size()){
	    double share = 1/clpCation.getValue().size();
	    for (unsigned int i=0; i<clpCation.getValue().size(); i++) {
		config.cations.push_back(std::make_tuple(
						clpCation.getValue()[i],
						share) );
	    }
	}
	else{
	    for (unsigned int i=0; i<clpCation.getValue().size(); i++) {
		config.cations.push_back(std::make_tuple(
						clpCation.getValue()[i],
						clpCationShare.getValue()[i]) );
	    }
	}

        if (clpAnion.getValue().size() > clpAnionShare.getValue().size()){
	    double share = 1/clpAnion.getValue().size();
	    for (unsigned int i=0; i<clpAnion.getValue().size(); i++) {
		config.anions.push_back(std::make_tuple(
						clpAnion.getValue()[i],
						share) );
	    }
	}
	else{
	    for (unsigned int i=0; i<clpAnion.getValue().size(); i++) {
		config.anions.push_back(std::make_tuple(
						clpAnion.getValue()[i],
						clpAnionShare.getValue()[i]) );
	    }
	}
    }

    
    // Optimization

    if (clpStartIndex.getValue() >= 0)
	config.startIndex = clpStartIndex.getValue();

    if (clpStopIndex.getValue() >= 0)
	config.stopIndex = clpStopIndex.getValue();

    if (clpRunCount.getValue() >= 0)
	config.runCount = clpRunCount.getValue();

    if (clpMaxDisplacement.getValue() >= 0)
	config.maxDisplacement = clpMaxDisplacement.getValue();
    
    if (clpNeighborRadius.getValue() > 0)
	config.neighborRadius = clpNeighborRadius.getValue();

    if (clpNeighborLayers.getValue() > 0)
	config.neighborLayers = clpNeighborLayers.getValue();    
    
    if (clpMaterialName.getValue() != "")
	config.materialSearchName = clpMaterialName.getValue();

    if (clpTemperature.getValue() > 0)
	config.latticeTemperature = clpTemperature.getValue();

    if (clpInterfaceAtoms.getValue() > 0)
	config.interfaceAtoms = clpInterfaceAtoms.getValue();    

    if (clpInterfacePositions.getValue() > 0)
	config.interfacePositions = clpInterfacePositions.getValue();    
    
    InterfaceToolMetric itMetric = ITM_Energy;

    if (clpMetric.getValue().length() > 0) {
        if (clpMetric.getValue() == "bondstrain")
            itMetric = ITM_BondStrain;
        else if (clpMetric.getValue() == "growthstrain")
            itMetric = ITM_GrowthStrain;
        else if (clpMetric.getValue() == "distanceStrain")
            itMetric = ITM_DistanceStrain;
    }

    // Switches
   
    CLOG(DEBUG, programName) << "Using the following configuration: "
	"\n\n " << config.str();
    
    //-------------------------------------------------------------------------
    // check if all necessary parameters have been specified
    //-------------------------------------------------------------------------

    if (config.interfaceAtoms < 1){
	bailout("Invalid amount of atoms in the interface specified!",
		ErrorCode::InterfaceAtoms);
    }

    if (config.interfacePositions < 1){
	bailout("Invalid number of positions to check for each aotm specified!",
		ErrorCode::InterfacePositions);
    }
    
    std::shared_ptr<SimulationBox> simbox;
   
    if (config.tersoffFileName == ""){
	bailout("No Tersoff Parameter file specified.", ErrorCode::Tersoff);
    }

    if ( (config.startIndex < 0 ) || (config.stopIndex < 0) ){
	bailout("Specify both start and stop Index.",
		ErrorCode::StartStopIndex);
    }

    // if neither radius nor atomic layer count specified.
    if ( ! ( (config.neighborRadius > 0) || (config.neighborLayers > 0) )  ){
	config.neighborLayers=1;
	CLOG(WARNING, programName) << "No neighbor determination method"
	    " specified! Falling back to default of layer count = 1";
    }

    if ( config.maxDisplacement < 0 ){
	bailout("Specify maximal displacement for Metropolis Monte Carlo "
		"relaxation.", ErrorCode::OptimizationAction);
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
	bailout("Error while reading periodic table file", ErrorCode::XML);
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

    //if (! S_ISDIR(config.journalPreamble))
	//	mkdir(config.journalPreamble);
    
    //-------------------------------------------------------------------------
    // Interface Creation
    //-------------------------------------------------------------------------
   
    try{

	Range3D<indexType> range;
	uint8_t outOfPlaneDimension = simbox->getOutOfPlaneDimension();

	range.start[outOfPlaneDimension] = config.startIndex;
	range.stop[outOfPlaneDimension] = config.stopIndex;
	range.apply[outOfPlaneDimension] = true;
	
	OptimizationParameter optParam;
	optParam.set(config);
	// Carry out all MMC steps in one call of the relaxation function.
	optParam.mmcRunCount = optParam.runCount;
	optParam.runCount = 1;

	TersoffPotential potential;

	fileHandler=getHandler(config.tersoffFileName,
			       xmlHandler, jsonHandler);
	fileHandler->load(config.tersoffFileName);
	fileHandler->get(potential);
	
	simbox->generateNeighbors(config.neighborLayers,
				  config.neighborRadius);
	
	ExchangeReaction exchangeReaction;

	if (clpExchangeReaction.getValue().size() == 1) {
	    exchangeReaction.addReaction(clpExchangeReaction.getValue()[0]);
	}
	else if (clpExchangeReaction.getValue().size() == 3) {
	    exchangeReaction.addReaction(clpExchangeReaction.getValue()[0],
					 clpExchangeReaction.getValue()[1],
					 std::stod(clpExchangeReaction.getValue()[2]) );
	}
	else
	    bailout("wrong number of arguments for exchange-reaction", ErrorCode::Unknown);

    // Material file specified
    MaterialCollection collection;
    
    try{
        fileHandler=getHandler(config.materialFile,
                       xmlHandler, jsonHandler);
        fileHandler->load(config.materialFile);
        fileHandler->get(collection);

    } catch(XmlException &e){
        std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
        bailout("Material file '" + config.materialFile + "' could not "
            "be opened", ErrorCode::XML);
    }
    
	if ( (config.cations.size()==0) && (config.anions.size()==0)){
        const Material &material =
		collection.getByName(config.materialSearchName);

	    // This part is also present in the else branch. Was necessary
	    // because getByName delivers const reference and in the other one
	    // we need a non-const object, making pointers infeasible.
	    CLOG(TRACE, programName) << "creating optimized interface "
		"roughness";
	    
	    InterfaceTool ifTool(*simbox, potential, exchangeReaction);
            ifTool.createOptimizedRoughness(range,
					    config.interfaceAtoms,
					    config.interfacePositions,
					    material,
					    optParam,
                        collection,
                        itMetric,
					    config.journalPreamble);   
	    
        } else {
            // Material directly specified
            MaterialComponents cations, anions;
            std::string symbol;
            double share;

            if (FP_ZERO == std::fpclassify(config.latticeConstant)){
                bailout("Invalid Lattice Constant!",
                    ErrorCode::LatticeConstant);
            }

	    if (config.materialSearchName == ""){
		bailout("No material name specified!",
			ErrorCode::Material);
	    }
	    
            for (auto el: config.cations) {
                std::tie (symbol, share) = el;
                cations.push_back(std::make_tuple
                           (PeriodicTable::getInstance().getBySymbol(symbol).id, share) );
            }

            for (auto el: config.anions) {
                std::tie (symbol, share) = el;
                anions.push_back(std::make_tuple
                          (PeriodicTable::getInstance().getBySymbol(symbol).id, share) );
            }
        
            Material material(config.materialSearchName,
			      cations, anions, config.latticeConstant,
			      config.c11, config.c12, config.c44);

	    CLOG(TRACE, programName) << "creating optimized interface "
		"roughness";
	    InterfaceTool ifTool(*simbox, potential, exchangeReaction);
            ifTool.createOptimizedRoughness(range,
					 config.interfaceAtoms,
					 config.interfacePositions,
					 material,
					 optParam,
                     collection,
                     itMetric,
					 config.journalPreamble);    

        }
	
	fileHandler = getHandler(config.outputFileName, xmlHandler,jsonHandler);
	fileHandler->clear();
        fileHandler->set(*simbox);
        fileHandler->save(config.outputFileName);
	
	if (config.xyzFileName != ""){
	    simbox->writeToXYZ(config.xyzFileName + ".xyz");
	    simbox->writeToXYZ(config.xyzFileName + "Modified.xyz",
			       Range3D<indexType>(),true);
	}
	    
    }catch(LatticeException &e){
	std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
	bailout("Error during Interface creation.", ErrorCode::Lattice);
    }catch(std::exception &e){
	std::cerr << "caught exception " << e.what() << std::endl;
	bailout("Unknown Error during Interface creation",
		ErrorCode::Unknown);
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
