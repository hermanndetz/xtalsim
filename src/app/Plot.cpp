/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


/*******************************************************************************

This is simple tool to plot data files.

*******************************************************************************/

#include "projectConfigure.h"

#include <iostream>
#include <memory>

#include <tclap/CmdLine.h>
#include <easyloggingcpp/easylogging++.h>

#include <xtalsim.h>

#ifdef __VTK__
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#endif

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
    std::string programName;

    //-------------------------------------------------------------------------
    // configure command line reader tool TCLAP
    //-------------------------------------------------------------------------
    
    TCLAP::ValueArg<std::string> clpInFile("i", "input", 
            "File containing data in CSV or XML format",
            false, "", "FILE");

    TCLAP::ValueArg<std::string> clpConfigFile("", "config",
		       "Use parameters from the configuration file in "
		       "location FILE. If parameters are given both on the "
		       "command line and in the config file the specification "
		       "on the command line dominates.", false,
					"", "FILE");

    TCLAP::ValueArg<std::string> clpInputFormat("f", "format",
            "Specifies the format of the input file. "
            "Possible options are csv or xml.", true,
            "csv", "FORMAT");

    TCLAP::MultiArg<std::string> clpJournalName("j", "journal",
            "Specifies the name of the journal to be read.", false,
            "NAME");

    TCLAP::MultiArg<std::string> clpColors("c", "color",
            "Specifies the color of a data series in web format RRGGBB."
            "If specified, one color per journal is required.", 
            false, "COLOR");

    TCLAP::MultiArg<std::string> clpDataSeriesTitles("s", "series-title",
            "Speciefies the title of a data series, which is displayed "
            "in the legend of the plot.", false, "TITLE");

    TCLAP::ValueArg<std::string> clpPlotTitle("t", "title",
            "Speciefies the plot title.", false, "", "TITLE");

    // \todo define data series via parameter list (filename.dat:4@2,color)

    // \todo x y limits

    TCLAP::SwitchArg clpVerbose("v", "verbose", "Enable debug messages.",
				   false);
    

    std::vector<std::string> usageExamples;

    usageExamples.push_back("EXECNAME -i foo.xml\n"
            "");

    try {
        TCLAP::CmdLine cmd("Create plot from data in CSV file."
                ""
                ""
                "", usageExamples, ' ', "0.1");

        cmd.add(clpInFile);
        cmd.add(clpConfigFile);
        cmd.add(clpInputFormat);
        cmd.add(clpJournalName);
        cmd.add(clpColors);
        cmd.add(clpDataSeriesTitles);
        cmd.add(clpPlotTitle);
        cmd.add(clpVerbose);

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
#if defined(EASYLOGGING_CONF_DIR)
    el::Configurations conf(EASYLOGGING_CONF_DIR "/easylogging++.conf");
#else
    el::Configurations conf("./input/easylogging++.conf");
#endif

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

    if (clpVerbose.getValue())
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
    
    CLOG(DEBUG, programName) << "starting " << programName;
    
    //-------------------------------------------------------------------------
    // fill Configuration object correctly, command line parameter overrule
    // those specified in the configuration file
    //-------------------------------------------------------------------------

    
    if (clpInFile.getValue() != "")
	config.inputFileName = clpInFile.getValue();

    //-------------------------------------------------------------------------
    // prepare for visualization
    //-------------------------------------------------------------------------
    
#ifdef __VTK__
    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
    std::vector<DataSeries> dataSeries;
    //! \todo Find a way to determine the correct tuple!
    //! Presently, this is only useful for energy vs MMC steps traces 
    //! and similar data sets.
    Journal<UDTuple> journal;

    char colName[3]{};
    uint32_t curCols = 0;

    std::vector<vtkSmartPointer<vtkTable>> tmpTables = std::vector<vtkSmartPointer<vtkTable>>{};

    if (clpInputFormat.getValue() == "csv") {
        std::shared_ptr<CsvHandler> csvfile;
        csvfile = std::make_shared<CsvHandler>(config.inputFileName, "CSVFile");

        tmpTables.push_back(csvfile->get());

        uint32_t i = 0;

        //! \todo implement iteration over different colors
        std::vector<std::string> colorCodes = clpColors.getValue();

        if (colorCodes.size() < tmpTables.front()->GetNumberOfColumns()) {
            for (uint32_t j = colorCodes.size(); j < tmpTables.front()->GetNumberOfColumns(); j++) {
                colorCodes.push_back("");
            }
        }

        for (i = 1; i < tmpTables.front()->GetNumberOfColumns(); i++) {
            Color col(colorCodes[i]);
            dataSeries.push_back(DataSeries(0, 0, i, col));
        }


    } else if (clpInputFormat.getValue() == "xml") {
        XmlHandler inputHandler;
	    inputHandler.load(clpInFile.getValue());

        std::vector<std::string> journalNames = clpJournalName.getValue();
        std::vector<std::string> colorCodes = clpColors.getValue();

        if (colorCodes.size() < journalNames.size()) {
            for (uint32_t j = colorCodes.size(); j < journalNames.size(); j++) {
                colorCodes.push_back("");
            }
        }

        uint32_t i = 0;

        for (auto journalName: journalNames) {
            try {
                inputHandler.get(journal,journalName);
            } catch (XmlException) {
                continue;
            }

            tmpTables.push_back(vtkSmartPointer<vtkTable>::New());

            std::vector<std::string> colTitles;
            colTitles.push_back("Steps");

            if (i < clpDataSeriesTitles.getValue().size()) {
                colTitles.push_back(clpDataSeriesTitles.getValue()[i]);
            } else {
                colTitles.push_back(journalName);
            }

            journal.get(tmpTables.back(), colTitles);

            //! \todo iterate over different colors here
            Color col(colorCodes[i]);
            //DataSeries ds = DataSeries(i*2, i*2+1, col);
            dataSeries.push_back(DataSeries(i, 0, 1, col));

            i++;
        }

    } else {
        return -1;
    }

    std::shared_ptr<XYPlot> plot;
    plot = std::make_shared<XYPlot>(tmpTables, clpPlotTitle.getValue(), "XYPlot");
    plot->plot(dataSeries);
#endif

    return 0;
}

