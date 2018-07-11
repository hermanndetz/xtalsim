
/*******************************************************************************

This is a reference implementation of how to compute the energy of a crystal.


*******************************************************************************/

#include <string>

#include <tclap/CmdLine.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/dist_sink.h>

//------------------------------------------------------------------------------

int main (int argc, char *argv[])
{

    auto distSink = std::make_shared<spdlog::sinks::dist_sink_st>();
    
    TCLAP::CmdLine cmd("Calculate energy of crystal", ' ', "0.1");
    TCLAP::ValueArg<std::string> input("i", "input", "Input File", true, "",
				       "FILE", cmd);

    TCLAP::ValueArg<std::string> output("o", "output", "Output File", true,
					"", "FILE", cmd);

    TCLAP::ValueArg<int> maxMove("m", "max-moves", "Maximal moves",
				 false, 5e4,"Integer", cmd);

    TCLAP::ValueArg<int> temperature("k", "temperature",
		       "Temperature in Kelvin", false, 0, "Integer", cmd);

    TCLAP::ValueArg<std::string> logFile("l", "log-file",
		       "Path to log file", false, "", "PATH", cmd);
    
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

}
