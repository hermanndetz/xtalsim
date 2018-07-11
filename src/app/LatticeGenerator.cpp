
/*******************************************************************************

This is a reference implementation to create a zincblende lattice in the
simulation box.

*******************************************************************************/

#include <string>
#include <vector>
#include <tuple>
#include <cmath>
#include <map>
#include <exception>
#include <signal.h>
#include <memory>

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
    
    
    TCLAP::MultiArg<int> clpSize("s", "size",
		       "NUM1 atomic layers are created in the first, NUM2 "
		       "in the second and NUM3 in the third dimension. If "
		       "an already "
		       "existing lattice is read from an input file it is only "
		       "extended in the growth dimension."
				 , false, "NUM1 NUM2 NUM3");

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

    TCLAP::ValueArg<std::string> clpMaterialFile("", "material-file",
		       "File containing definition of required materials.",
		       false, "", "FILE");

    TCLAP::ValueArg<std::string> clpMaterialName("n", "material-name",
		       "Name of material that shall be used. Also mandatory in "
		       "combination with --material-file. "
		       "Naming cations and anions overrules specifying "
		       "material from file.", false, "", "NAME");
    
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

    TCLAP::ValueArg<std::string> clpPeriodicTableFile("p", "periodic-table",
		       "Read information regarding elements from FILE. This "
		       "parameter is mandatory.", false,"", "FILE");
    
    TCLAP::ValueArg<std::string> clpTersoffFile("t", "tersoff",
		       "FILE that contains the tersoff parameters, which are "
		       "required to calculate the energy. Only required in "
		       "combination with -E", false,"", "FILE");

    TCLAP::ValueArg<std::string> clpXYZFile("", "xyz",
	               "FILE xyz data shall be written to. ",
		       false, "", "FILE");

    TCLAP::ValueArg<std::string> clpLatticeType("", "lattice-type",
		       "Type of lattice that shall be created. Currently "
		       "'zincblende', 'half-heusler' and 'full-heusler' "
		       "are supported.", false, "",
		       "[zincblende|half-heusler|full-heusler]");
    
    TCLAP::ValueArg<double> clpLatticeConstant("l", "lattice-constant",
		       "Lattice constant of the introduced "
		       "material.", false, 0, "VAL");

    TCLAP::ValueArg<double> clpNeighborRadius("", "neighbor-radius",
		       "All atoms at most NUM Angstroem away are considered "
		       "neighbors of a single atom, "
		       "which are afterwards used to calculate the energy. "
		       "Either this or neighbor-layers has to be specified "
		       "(both is also possible) if switch -E is used.",
		       false, 0, "NUM");

    TCLAP::ValueArg<indexType> clpNeighborLayers("", "neighbor-layers",
		       "All atoms at most NUM atomic layers away are considered "
		       "neighbors of a single atom, "
		       "which are afterwards used to calculate the energy. "
		       "If switch -E was activated "
		       "either this or neighbor-radius has to be specified "
		       "(both is also possible). In the latter case only atoms "
		       "which are also within the specified radius are picked.",
		       false, 0, "NUM");
    
    TCLAP::ValueArg<int> clpGrowthDimension("d", "growth-dimension",
		       "Defines dimension DIM as the dimension the lattice "
		       "grows in. Only values 1,2 and 3 are allowed. "
		       "If not specified default value of 3 is chosen.",
		       false, 0, "DIM");
    
    TCLAP::ValueArg<std::string> clpLogFile("", "log",
		       "If specified log messages are also written into "
		       "this file.", false, "", "FILE");
    
    TCLAP::MultiSwitchArg clpVerbose("v", "verbose", "Enable debug messages. "
                "First occurence enables debug and trace, second occurence "
                "also enables info.",
				   false);

    TCLAP::MultiSwitchArg clpQuiet("q", "quiet", "First occurence disables"
                " screen output, second occurence also disables file "
                "logging.", false);

    TCLAP::SwitchArg clpQQuiet("", "qq", "Disables screen output and file "
                "logging.", false);
    
    TCLAP::SwitchArg clpEnergy("E", "calc-energy", "Calculate energy after "
		       "lattice was created.", false);
    
    TCLAP::ValueArg<double> clpTemperature("k", "temperature",
		       "Internal lattice temperature in Kelvin.", false,
					   0, "NUM");
    
    std::vector<std::string> examples;

    examples.push_back("EXECNAME -c Ga -a As -n GaAs -s 4 4 4 "
		       "-l 5.6535 -p periodicTable.xml -o foo.xml\n"
		       "Create GaAs lattice with a size of 4x4x4 atomic"
		       "layers and a lattice constant of 1.4 and store the "
		       "results in foo.xml .");

    examples.push_back("EXECNAME --config foo.xml -o bar.xml\n"
		       "Read parameters from configuration file foo.xml but "
		       "override output file name with bar.xml .");
    
    try{

	TCLAP::CmdLine cmd("Generate a lattice of defined size and introduce"
		       " a zincblende atomic positioning on top of it. Only "
		       "one material can be "
		       "introduced at a time. For additional layers the program"
		       " has to be rerun.",
		       examples, ' ', "1.0");

	cmd.add(clpSize);
	cmd.add(clpCation);
	cmd.add(clpCationShare);
	cmd.add(clpAnion);
	cmd.add(clpAnionShare);
	cmd.add(clpMaterialFile);
	cmd.add(clpMaterialName);
	cmd.add(clpInFile);
	cmd.add(clpOutFile);
	cmd.add(clpConfigFile);
	cmd.add(clpPeriodicTableFile);
	cmd.add(clpTersoffFile);
	cmd.add(clpXYZFile);
	cmd.add(clpLatticeType);
	cmd.add(clpLatticeConstant);
	cmd.add(clpGrowthDimension);
	cmd.add(clpLogFile);
	cmd.add(clpVerbose);
	cmd.add(clpEnergy);
	cmd.add(clpNeighborRadius);
	cmd.add(clpNeighborLayers);
	cmd.add(clpTemperature);
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

    // std::shared_ptr<FileHandler> pHandler;
    // std::shared_ptr<XmlHandler> pXmlHandler =
    // 	std::make_shared<XmlHandler>(distSink);
    // pHandler = pXmlHandler; 
    
    if (clpConfigFile.getValue() != "") {
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

    if ( clpVerbose.getValue() > 0) 
        config.verbose = true;

    if ( config.verbose ){
        conf.set(el::Level::Debug, el::ConfigurationType::Enabled, "true");
        conf.set(el::Level::Trace, el::ConfigurationType::Enabled, "true");

        if (clpVerbose.getValue() > 1)
            conf.set(el::Level::Info, el::ConfigurationType::Enabled, "true");
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

    // Other Files
    if (clpPeriodicTableFile.getValue() != "")
        config.periodicTableFileName = clpPeriodicTableFile.getValue();
    
    if (clpInFile.getValue() != "")
        config.inputFileName = clpInFile.getValue();

    if (clpOutFile.getValue() != "")
        config.outputFileName = clpOutFile.getValue();

    if (clpTersoffFile.getValue() != "")
        config.tersoffFileName = clpTersoffFile.getValue();

    if (clpXYZFile.getValue() != "")
        config.xyzFileName = clpXYZFile.getValue();
    
    // Atoms and Size
    if (clpSize.getValue().size()>0) {
        for (unsigned int i=0; (i<3) && (i<clpSize.getValue().size()); i++)
	    config.size[i] = clpSize.getValue()[i];
    }

    if (clpLatticeConstant.getValue() > 0)
        config.latticeConstant = clpLatticeConstant.getValue();
       
    if (clpGrowthDimension.getValue() > 0)
        config.growthDimension = clpGrowthDimension.getValue();

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

    if (clpMaterialFile.getValue() != "")
        config.materialFile = clpMaterialFile.getValue();

    if (clpMaterialName.getValue() != "")
        config.materialSearchName = clpMaterialName.getValue();  

    if (clpNeighborRadius.getValue() > 0)
        config.neighborRadius = clpNeighborRadius.getValue();

    if (clpNeighborLayers.getValue() > 0)
        config.neighborLayers = clpNeighborLayers.getValue();    

    if (clpTemperature.getValue() > 0)
	config.latticeTemperature = clpTemperature.getValue();

    if (clpLatticeType.getValue() != "")
        config.latticeType = clpLatticeType.getValue();
    
    // Switches
    if ( clpEnergy.getValue() == true )
        config.calculateEnergy = true;
    
    CLOG(DEBUG, programName) << "Using the following configuration: "
	"\n\n " << config.str();

    //-------------------------------------------------------------------------
    // check if all necessary parameters have been specified
    //-------------------------------------------------------------------------

    for (int i=0; i<3; i++){

	if ( ((config.size[i])%2) == 1){
	    CLOG(WARNING, programName) << "Odd size value " <<
		config.size[i] <<" not suited for periodic boundary conditions."
					       << "Increasing it by 1.",
	    config.size[i] += 1;
	}
    }
    
    if ((config.growthDimension < 1) || (config.growthDimension > 3)){
	CLOG(WARNING, programName) <<
	   "Wrong growth dimension specified. Falling back to default value 3.";
        config.growthDimension = 3;
    }
    
    std::shared_ptr<SimulationBox> simbox;

    try {
	CLOG(DEBUG, programName) << "loading periodic table file '" <<
	    config.periodicTableFileName << "'";
	fileHandler=getHandler(config.periodicTableFileName,
			       xmlHandler, jsonHandler);
	fileHandler->load(config.periodicTableFileName);
	fileHandler->get(PeriodicTable::getInstance());
    } catch(XmlException &e) {
        std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
        bailout("Error while reading periodic table file",
            ErrorCode::XML);
    } catch(PeriodicTableException &e) {
        std::cout << e.what() << ": "  << e.getMessage() << std::endl;
        bailout("Error while filling periodic table",
            ErrorCode::PeriodicTable);
    }
    
    // if neither radius nor atomic layer count specified.
    bool radiusOrLayers =  (config.neighborRadius > 0) ||
	(config.neighborLayers > 0);
    if ( config.calculateEnergy && (! radiusOrLayers)  ){
	config.neighborLayers=1;
	CLOG(WARNING, programName) << "No neighbor determination method"
	    "specified! Falling back to default of layer count = 1";
    }
    
    if (config.inputFileName == "") {
        if (config.outputFileName == "") {
            bailout("Neither input nor output file name specified!",
                ErrorCode::InOutFiles);
        }

        if ((config.growthDimension > 3) && (config.growthDimension < 1)) {
            bailout("Invalid growth dimension (-d, --growth-dimension). Is "
                "mandatory if no input file was specified.",
                ErrorCode::GrowthDimension);
        }	

        if ( ! (config.size > 0) ) {
            bailout("Size=0 not allowed when creating new simulation box!",
                ErrorCode::Size);
        }
        
	CLOG(TRACE, programName) << "creating new simulation box";
        simbox = std::make_shared<SimulationBox>(config.growthDimension-1);

    } else {
        if (config.outputFileName == "")
            config.outputFileName = config.inputFileName;
        
	CLOG(DEBUG, programName) << "loading input file " <<
	    config.inputFileName;
	fileHandler=getHandler(config.inputFileName, xmlHandler, jsonHandler);
	fileHandler->load(config.inputFileName);
	fileHandler->get(simbox);

        if (config.size[simbox->getOutOfPlaneDimension()] == 0) {
            bailout("Size=0 in growth dimension not possible!",
                ErrorCode::Size);
        }
    }

    if (config.latticeTemperature > 0)
	simbox->setLatticeTemperature(config.latticeTemperature);
    
    //-------------------------------------------------------------------------
    // start creation of lattice
    //-------------------------------------------------------------------------

    CLOG(TRACE, programName) << "creating Material";
    
    try {

	MaterialCollection collection;
	const Material *material;
	
        if ( (config.cations.size()==0) && (config.anions.size()==0)){
            // Material file specified
                
	    try{
		fileHandler=getHandler(config.materialFile,
				       xmlHandler, jsonHandler);
		fileHandler->load(config.materialFile);
		
	    } catch(XmlException &e){
		std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
		bailout("Material file '" + config.materialFile + "' could not "
			"be opened", ErrorCode::XML);
	    }
	    
            fileHandler->get(collection);
	    
	    material = &(collection.getByName(config.materialSearchName));
            
        } else {
            // Material directly specified
            MaterialComponents cations, anions;
            std::string symbol;
            double share;

            if (FP_ZERO == std::fpclassify(config.latticeConstant)){
                bailout("Invalid Lattice Constant!",ErrorCode::LatticeConstant);
            }

	    if (config.materialSearchName == ""){
		bailout("No material name specified!", ErrorCode::Material);
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

            material = new Material(config.materialSearchName,
			      cations, anions, config.latticeConstant,
			      config.c11, config.c12, config.c44);
	}

	// Create lattice based on defined type.
	if (config.latticeType == "half-heusler"){
	    CLOG(TRACE, programName) << "creating half-Heusler lattice";
	    simbox->createHalfHeusler(config.size, *material);
	}
	else if (config.latticeType == "full-heusler"){
	    CLOG(TRACE, programName) << "creating full-Heusler lattice";
	    simbox->createFullHeusler(config.size, *material);
	}
	else{
	    CLOG(TRACE, programName) << "creating Zincblende lattice";
	    simbox->createZincblende(config.size, *material);
	}
	
	CLOG(TRACE, programName) << "exporting created lattice";

	fileHandler = getHandler(config.outputFileName, xmlHandler,jsonHandler);
	fileHandler->clear();
        fileHandler->set(*simbox);
        fileHandler->save(config.outputFileName);
	
	if (config.xyzFileName != "")
	    simbox->writeToXYZ(config.xyzFileName);
	
        if (config.calculateEnergy) {
            if (config.tersoffFileName == "") {
                bailout("No tersoff parameters specified.",
                    ErrorCode::Tersoff);
            }
    
            TersoffPotential potential;
	    
	    fileHandler=getHandler(config.tersoffFileName,
				   xmlHandler, jsonHandler);
	    fileHandler->load(config.tersoffFileName);
            fileHandler->get(potential);

	    simbox->generateNeighbors(config.neighborLayers,
				      config.neighborRadius);
	    
            double energy = simbox->getEnergy(potential);

	    CLOG(INFO, programName) << "calculated energy of " <<
		energy;
	    CLOG(INFO, programName) << energy/simbox->getAtomCount() <<
		" eV/atom";
        }
    } catch(LatticeException &e) {
        std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
        bailout("Error while initializing lattice",
            ErrorCode::Lattice);
    } catch(XmlException &e) {
        std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
        bailout("Error while exporting newly created lattice",
            ErrorCode::XML);
    } catch(PeriodicTableException &e) {
        std::cerr << e.what() << ": "  << e.getMessage() << std::endl;
        bailout("Error while handling periodic table",
            ErrorCode::PeriodicTable);
    } catch(std::exception &e) {
        std::cerr << "caught exception " << e.what() << std::endl;
        bailout("Unknown Error while creating lattice",
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
