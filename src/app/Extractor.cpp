/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


/*******************************************************************************

This is tool to extract data from .xml files (e.g. journals) and save them
into other file formats.

*******************************************************************************/

#include <iostream>
#include <memory>

#include <tclap/CmdLine.h>
#include <easyloggingcpp/easylogging++.h>

#include <xtalsim.h>
// #include <misc/Color.h>
// #include <misc/XmlHandler.h>
// #include <misc/Configuration.h>
// #include <misc/ErrorCode.h>
// #include <misc/CsvHandler.h>
// #include <physics/SimulationBox.h>

INITIALIZE_EASYLOGGINGPP

std::string programName;

//******************************************************************************

void bailout(const std::string &str, const int errorCode)
{
    CLOG(ERROR, programName) << str;
    exit(errorCode);
}

//******************************************************************************

int main (int argc, char *argv[])
{

    //-------------------------------------------------------------------------
    // configure command line reader tool TCLAP
    //-------------------------------------------------------------------------
    
    TCLAP::ValueArg<std::string> clpInFile("i", "input", 
            "File containing input data in"
            ".xml format", false, "", "FILE");

    TCLAP::ValueArg<std::string> clpOutFile("o", "output",
		       "Specifies the file, where the extracted data are"
               "saved.", true, "", "FILE");

    TCLAP::ValueArg<std::string> clpConfigFile("", "config",
		       "Use parameters from the configuration file in "
		       "location FILE. If parameters are given both on the "
		       "command line and in the config file the specification "
		       "on the command line dominates.", false,
					"", "FILE");

    TCLAP::ValueArg<std::string> clpLogFile("", "log",
		       "If specified log messages are also written into "
		       "this file.", false, "", "FILE");

    TCLAP::MultiArg<std::string> clpJournalName("j", "journal",
            "Specifies the name of the journal to be extracted.", false,
            "NAME ...");

    TCLAP::ValueArg<std::string> clpOutputFormat("f", "format",
            "Specifies the output format.", true, "", "FORMAT");
    
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
        cmd.add(clpOutFile);
        cmd.add(clpConfigFile);
        cmd.add(clpLogFile);
        //cmd.add(clpOutputFormat);
        cmd.add(clpJournalName);
        cmd.add(clpVerbose);
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

    // Load logger specific settings.
    el::Loggers::configureFromGlobal("../input/loggers.conf");
    CLOG(TRACE, programName) << "starting " << programName;
    
    //-------------------------------------------------------------------------
    // Attention! Copied from LatticeGenerator.cpp and modified
    // to a minimum. BEGIN
    if (clpOutFile.getValue() != "")
        config.outputFileName = clpOutFile.getValue();
    else
        return -1;
    
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
    
    CLOG(TRACE, programName) << "starting " << programName;

    // Attention! Copied from LatticeGenerator.cpp and modified
    // to a minimum. END
    //-------------------------------------------------------------------------

    //-------------------------------------------------------------------------
    // fill Configuration object correctly, command line parameter overrule
    // those specified in the configuration file
    //-------------------------------------------------------------------------

    // nothing to be done for now

    // do something useful with the input data
    Journal<UDTuple> journal;

    XmlHandler inputHandler;
    inputHandler.load(clpInFile.getValue());

    std::vector<std::string> journalNames = clpJournalName.getValue();

    CsvHandler outputHandler = CsvHandler(config.outputFileName, "CSV File");

    for (auto journalName: journalNames) {
        inputHandler.get(journal,journalName);
        outputHandler += journal;
    }

    outputHandler.save(config.outputFileName);

    return 0;
}


