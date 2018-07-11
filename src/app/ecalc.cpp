
/*******************************************************************************

This is a reference implementation of how to compute the energy of a crystal.


*******************************************************************************/

#include <string>

#include <tclap/CmdLine.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/dist_sink.h>

#include <misc/XmlHandler.h>
#include <physics/SimulationBox.h>
#include <simulation/Optimization.h>
#include <simulation/TersoffPotential.h>

int main (int argc, char *argv[])
{

    std::vector<std::string> usageExamples;

    //! \todo define useful examples for ecalc
    usageExamples.push_back("EXECNAME ...");
    
    //----------------------------------------------------------------------
    
    TCLAP::CmdLine cmd("Calculate energy of crystal",
            usageExamples, ' ', "0.1");
    TCLAP::ValueArg<std::string> input("i", "input", "Input File", true, "",
				       "FILE", cmd);

    TCLAP::ValueArg<std::string> output("o", "output", "Output File. If not"
		   " specified input file is taken", false,"", "FILE", cmd);


    
    TCLAP::ValueArg<std::string> logFile("g", "log-file",
		       "Path to log file", false, "", "PATH", cmd);

    TCLAP::ValueArg<std::string> tersoff("t", "tersoff",
	  "Path to file holding tersoff parameters", true, "", "PATH", cmd);
    
    TCLAP::ValueArg<int> maxMove("m", "max-moves", "Maximal moves",
				 false, 5e4,"Integer", cmd);

    TCLAP::ValueArg<int> minIndex("l", "min-index", "Energy calculated only for"
				  " a specific range starting at this layer",
				 false, -1,"Integer", cmd);

    TCLAP::ValueArg<int> maxIndex("u", "max-index", "Energy calculated only for"
				  " a specific range stopping at this layer",
				 false, -1,"Integer", cmd);

    TCLAP::ValueArg<double> temperature("k", "temperature",
		       "Temperature in Kelvin", false, 0, "Decimal", cmd);
    
    TCLAP::SwitchArg printNeighbor("n", "neighbor", "include neighbors",
				   cmd, false);

    try{
	
	cmd.parse( argc, argv );
	
    }
    catch (TCLAP::ArgException &e){
	std::cerr << "error: " << e.error() << " for arg " << e.argId() <<
		  std::endl;
	exit(1);
    }

    //----------------------------------------------------------------------
    
   
#ifdef _WIN32
    auto sinkStdOut =std::make_shared<spdlog::sinks::wincolor_stdout_sink_st>();
#else
    auto sinkStdOut = std::make_shared<spdlog::sinks::ansicolor_sink>
	(spdlog::sinks::stdout_sink_st::instance());
#endif

    auto distSink = std::make_shared<spdlog::sinks::dist_sink_st>();
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

    PeriodicTable periodicTable(distSink);

    // \todo{level not defined}
    //periodicTable.setLoggerLevel(level);

    std::shared_ptr<SimulationBox> simbox = inputXml.get(distSink, periodicTable);
    
    
    XmlHandler tersoffXml(distSink,"XmlHandler");
    logger->info("reading tersoff potential file");
    tersoffXml.load(tersoff.getValue());

    TersoffPotential potential(distSink);
    tersoffXml.get(potential);

    Optimization opt(distSink);
    opt.registerAction(std::string("MMC"), &SimulationBox::mmcRelax, 100);
    opt.registerAction(std::string("Scale"), &SimulationBox::scale, 1);



    if ( ( minIndex.getValue() > -1) && ( maxIndex.getValue() > -1) )
	logger->info("using range mode");
    else
	logger->info("using calculation mode");
    
}
