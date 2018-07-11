
/*******************************************************************************

This is a reference implementation of how to compute the energy of a crystal.


*******************************************************************************/

#include <string>

#include <tclap/CmdLine.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/dist_sink.h>

#include <misc/XmlHandler.h>

int main (int argc, char *argv[])
{

    auto distSink = std::make_shared<spdlog::sinks::dist_sink_st>();
    
    TCLAP::CmdLine cmd("Evaluate atomic bonds", ' ', "0.1");
    TCLAP::ValueArg<std::string> input("i", "input", "Input File", true, "",
				       "FILE", cmd);

    TCLAP::ValueArg<std::string> output("o", "output", "Output File", false,
					"", "FILE", cmd);
   
    TCLAP::ValueArg<std::string> logFile("l", "log-file",
		       "Path to log file", false, "", "PATH", cmd);
    
    TCLAP::ValueArg<double> cutoff("d", "cutoff", "Cutoff distance",
				 true, -1,"Decimal", cmd);

    TCLAP::ValueArg<int> minIndex("l", "min-index", "Energy calculated only for"
				  " a specific range starting at this layer",
				 false, -1,"Integer", cmd);

    TCLAP::ValueArg<int> maxIndex("u", "max-index", "Energy calculated only for"
				  " a specific range stopping at this layer",
				 false, -1,"Integer", cmd);

    TCLAP::SwitchArg modeABL("a", "mode-avg-bond-length",
			     "Average Bond Length", cmd, false);

    TCLAP::SwitchArg modeBC("b", "mode-bond-cound",
			     "Bond Count", cmd, false);

    TCLAP::SwitchArg modeSM("s", "mode-strain-map",
			    "Strain Map", cmd, false);
    
    TCLAP::SwitchArg range("r", "range", "bond length range", cmd, false);

    TCLAP::SwitchArg csv("c", "csv", "CSV mode", cmd, false);    
    
    TCLAP::SwitchArg printNeighbor("n", "neighbor", "include neighbors",
				   cmd, false);

    try{

	cmd.xorAdd ( range, csv );
	cmd.parse( argc, argv );
	
    }
    catch (TCLAP::ArgException &e){
	std::cerr << "error: " << e.error() << " for arg " << e.argId() <<
		  std::endl;
	exit(1);
    }

#ifdef _WIN32
    auto sinkStdout =std::make_shared<spdlog::sinks::wincolor_stdout_sink_st>();
#else
    auto sinkStdout = std::make_shared<spdlog::sinks::ansicolor_sink>
	(spdlog::sinks::stdout_sink_st::instance());
#endif
    
    distSink->add_sink(sinkStdOut);
    if (logFile.getValue() != ""){
	distSink->add_sink(std::make_shared<spdlog::sinks::simple_file_sink_st>
			   ( logFile.getValue() ) );
    }
    
    auto logger = std::make_shared<spdlog::logger>("ecalc", distSink);
    logger->set_level(spdlog::level::trace);
    logger->set_pattern("%Y/%m/%d %T [%n:%l] %v");    

    logger->info("starting {}", argv[0]);

    //----------------------------------------------------------------------

    XmlHandler inputXml(distSink,"XmlHandler");
    logger->info("reading input file");
    inputXml.load(input.getValue());
    
    if (csv.getValue()){

	logger->info("starting CSV mode");
	
    }

    //----------------------------------------------------------------------

    if (range.getValue()){

	logger->info("starting range mode");
	
    }

    //----------------------------------------------------------------------

    if (modeABL.getValue()){

	logger->info("starting average bond length mode");
	
    }

    //----------------------------------------------------------------------

    if (modeBC.getValue()){

	logger->info("starting bond count mode");
	
    }

    //----------------------------------------------------------------------

    if (modeSM.getValue()){

	logger->info("starting strain map mode");
	
    }

    //----------------------------------------------------------------------

    logger->info("finished {}", argv[0]);

    return 0;
}
