
#include <physics/Atom.h>
#include <physics/PeriodicTable.h>
#include <physics/MaterialCollection.h>
#include <physics/Material.h>
#include <physics/Lattice.h>
#include <physics/SimulationBox.h>
#include <misc/XmlHandler.h>
#include <misc/JsonHandler.h>
#include <misc/Journal.h>
#include <physics/Vector3D.h>
#include <physics/Range3D.h>
//#include <physics/Vector3D.cpp>
#include <simulation/Optimization.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/dist_sink.h>
#include <iostream>
#include <memory>
#include <chrono>
#include <ctime>
#include <sstream>
#include <typeinfo>

#include <easyloggingcpp/easylogging++.h>

INITIALIZE_EASYLOGGINGPP

int main (int argc, char *argv[])
{

    //    el::Configurations defaultConf("../easyloggingcpp/easylogging++.config");
    //el::Loggers::setDefaultConfigurations(defaultConf, true);

    std::string a = "JsonHandler";
    const char *name(&a);
    
    el::Loggers::configureFromGlobal("../easyloggingcpp/easylogging++.config");
    
    LOG(INFO) << "First Log";
    CLOG(WARNING, name) << "First Log json";
    CLOG(ERROR, "XmlHandler") << "First Log xml";

    el::Loggers::reconfigureAllLoggers(el::ConfigurationType::ToFile, "true");
    el::Loggers::reconfigureAllLoggers(el::ConfigurationType::Filename, "a.log");
    
    CLOG(INFO, "aaa") << "First Log aaa";
    CLOG(FATAL, "XmlHandler") << "Second Log xml";
    exit(1);
    
    //    std::shared_ptr<Journal> journal;
    std::string s;
    std::time_t now; 
   
    auto distSink = std::make_shared<spdlog::sinks::dist_sink_st>();

    #ifdef _WIN32
    auto sinkStdout =std::make_shared<spdlog::sinks::wincolor_stdout_sink_st>();
    #else
    auto sinkStdout = std::make_shared<spdlog::sinks::ansicolor_sink>
	(spdlog::sinks::stdout_sink_st::instance());
    #endif
    auto sinkFile = std::make_shared<spdlog::sinks::simple_file_sink_st>
	("../tmp/log2");
    distSink->add_sink(sinkStdout);
    distSink->add_sink(sinkFile);

    auto logger = std::make_shared<spdlog::logger>("main", distSink);
    logger->set_level(spdlog::level::trace);

    logger->set_pattern("%Y/%m/%d %T [%n:%l] %v");
    /*
    try{
	journal = std::make_shared<Journal>("../tmp/log1", true);
    }
    catch(JournalException &e){
	std::cout << e.what() << ": "  << e.getMessage() << std::endl;
	exit(1);	
    }
    catch(std::exception &e){
	std::cout << "Journal file tmp/log1 could not be opened" << std::endl;
	exit(1);
    }
    
    journal->setLevel("main",1);
    journal->setLevel("Material",2);
    */

    logger->debug("[START OF LOGGING]");
    logger->debug("creating XmlHandler", "main", 1);

    PeriodicTable table(distSink);
    table.setLoggerLevel(spdlog::level::info);
    XmlHandler handler(distSink,"NoXmlHandler");
    logger->debug("load periodic table");
    handler.load("../xml/elements.xml");
    handler.get(table);
    
    JsonHandler jsonHandler(distSink);
    
    MaterialCollection collection(distSink);
    collection.setLoggerLevel(spdlog::level::info);

    XmlHandler xhandler(distSink,"NoXmlHandler");
    logger->debug("load XML file");
    xhandler.load("../xml/materials.xml");
    xhandler.get(collection, table, distSink);

    xhandler.save("../xml/materials2.xml");
    
    XmlHandler xhandler2(distSink,"NoXmlHandler");
    xhandler2.set(collection, table);
    xhandler2.save("../xml/materials3.xml");

    jsonHandler.set(collection, table);
    jsonHandler.save("../json/materials.json");

    MaterialCollection collection2(distSink);
    collection2.setLoggerLevel(spdlog::level::debug);
    jsonHandler.get(collection2, table, distSink);
    
    JsonHandler jsonHandler2(distSink);
    jsonHandler2.set(collection2, table);
    jsonHandler2.save("../json/materials2.json");
    exit(1);

    // PeriodicTable table(distSink);
    // //distSink->remove_sink(sinkFile);
    
    // // try{
    // // 	XmlHandler handler(distSink,"NoXmlHandler");
    // // 	Journal<int> intJournal("name", "descrption");
    // // 	Journal<std::string> strJournal("name2", "description");

    // // 	intJournal.addEntry(3);
    // // 	intJournal.addEntry(4);
    // // 	intJournal.addEntry(5);

    // // 	strJournal.addEntry("adf");
    // // 	strJournal.addEntry("gsdfg");
    // // 	strJournal.addEntry("ggg");
	
    // // 	handler.set(intJournal);
    // // 	handler.save("../xml/journal.xml");

    // // 	XmlHandler handler2(distSink,"NoXmlHandler");
    // // 	Journal<int> intJournal2;
    // // 	handler2.load("../xml/journal.xml");
    // // 	handler2.get(intJournal2, "name");
    // // 	handler2.get(strJournal, "name");
	
    // // 	XmlHandler handler3(distSink,"NoXmlHandler");
    // // 	handler3.set(intJournal2);
    // // 	handler3.set(intJournal2);
    // // 	handler3.set(strJournal);
    // // 	handler3.save("../xml/journal2.xml");
    // // }
    // // catch(std::exception &e){
    // // 	std::cout << "caught exception " << e.what() << std::endl;
    // // }


    // try{
    // 	XmlHandler handler(distSink,"NoXmlHandler");
    // 	logger->debug("load XML file");
    // 	handler.load("../xml/elements.xml");
    // 	logger->debug("save XML file");
    // 	handler.save("../xml/test_out.xml");

    // 	handler.get(table);
    // 	handler.set(table);
    // 	handler.save("../xml/elements1.xml");

    // 	Element el=table.getBySymbol("Ga");
    // 	logger->info(el.name);

    // 	uint16_t id = el.id;
    // 	logger->info(el.id);
    // 	Element el4=table.getById(id);
    // 	logger->info(el4.name);
	
    // }
    // catch(PeriodicTableException &e){
    // 	std::cout << e.what() << ": "  << e.getMessage() << std::endl;
    // }	
    // catch(XmlException &e){
    // 	std::cout << e.what() << ": "  << e.getMessage() << std::endl;
    // }	
    // catch(std::exception &e){
    // 	std::cout << "caught exception " << e.what() << std::endl;
    // }

    // //sinkStdout->set_color(spdlog::level::info,
    // //			  sinkStdout->magenta + sinkStdout->underline);
    
    // //distSink->add_sink(sinkFile);
    
    // try{
    // 	Element el = table.getByProtonNeutron(14,14);
    // 	logger->info("Got element silicon");
    // 	logger->info("Got element {0} ({1})", el.name, el.symbol);
	
    // 	Element el2 = table.getBySymbol("As");
    // 	logger->info("Got element Arsenide");	
    // }
    // catch(PeriodicTableException &e){
    // 	std::cout << e.what() << ": "  << e.getMessage() << std::endl;
    // }	
    // catch(std::exception &e){
    // 	std::cout << "caught exception " << e.what() << std::endl;
    // }

    // // std::map<uint16_t, double> cations = {
    // // 	{table.getBySymbol("Ga").id, 0.3},
    // // 	{table.getBySymbol("As").id, 0.7}};
    // // std::map<uint16_t, double> anions = {
    // // 	{table.getBySymbol("Ga").id, 0.2},
    // // 	{table.getBySymbol("As").id, 0.8}};

    // std::map<uint16_t, double> cations = {
    // 	{table.getBySymbol("Ga").id, 1}};
    // std::map<uint16_t, double> anions = {
    // 	{table.getBySymbol("As").id, 1}};

    // XmlHandler handler(distSink);
    // handler.load("../xml/materials.xml");
    // MaterialCollection collection22(distSink);
    // handler.get(collection22, table, distSink);

    // Material matm1 = collection22.getByName("InAs");
    // Material matm2 = collection22.getByName("GaSb");
    // SimulationBox simbox22(distSink, 2);
    // simbox22.setLoggerLevel(spdlog::level::trace);
    // simbox22.setLatticeLoggerLevel(spdlog::level::trace);
    // simbox22.createZincblende(Vector3D<indexType>(28,28,8), matm1);
    // simbox22.createZincblende(Vector3D<indexType>(28,28,8), matm2);
    // simbox22.createZincblende(Vector3D<indexType>(28,28,8), matm1);

    // simbox22.writeToXYZ("../xml/scaleStart.xyz", table);
    
    // TersoffPotential tersoff(distSink);
    // OptimizationParameter param;
    // OptimizationStatistic stat;

    // param.maxScaling = 1;
    // Range3D<indexType> range44;
    // range44.start={0,0,10};
    // range44.stop={0,0,5};
    // range44.apply={false,false,true};
    
    // simbox22.scale(tersoff,param,stat, range44);
    // simbox22.writeToXYZ("../xml/scale.xyz", table);

    // return 0;
    
    // Material mat1(distSink,"InGaAs", cations, anions, 1.4, 1.22, 3.45);
    
    // // for (int i=0; i<10; i++){
    // // 	uint16_t id = mat1.getRandomCation();

    // // 	if ( ! mat1.hasCation(id)) continue;
    // // 	Element el3 = table.getById(id);
    // // 	logger->info("chose cation {}", el3.name);
    // // }

    // // for (int i=0; i<10; i++){
    // // 	uint16_t id = mat1.getRandomAnion();

    // // 	if ( ! mat1.hasAnion(id)) continue;
    // // 	Element el3 = table.getById(id);
    // // 	logger->info("chose cation {}", el3.name);
    // // }

    // Vector3D<indexType> vec1(1,2,3);
    // Vector3D<indexType> vec2(2,4,6);

    // Vector3D<indexType> vec3(vec2);

    // vec3 = vec2-vec1;

    // if (vec3==vec1)
    // 	logger->info("== in Vector3D works");

    // indexType arr[2];
    // arr[0]=6;arr[1]=7;arr[2]=8;
    
    // vec3 = arr;

    // logger->info("values are {0}{1}{2}", vec3[0], vec3[1],
    // 		 vec3[2]);

    // logger->info("Simulation Box");

    // SimulationBox simbox(distSink, 2);
    // simbox.setLoggerLevel(spdlog::level::trace);
    // simbox.createZincblende(Vector3D<indexType>(12,12,6), mat1);
    // simbox.createZincblende(Vector3D<indexType>(12,12,6), mat1);

    // XmlHandler handler6(distSink);
    // handler6.load("../xml/materials.xml");
    // MaterialCollection collection(distSink);
    // handler6.get(collection, table, distSink);

    // simbox.generateNeighbors(1);
    // simbox.calculateStrain(collection, "../tmp/strainOut", 3);
    
    // try{
    // 	simbox.generateNeighbors(1,0.672);
    // }
    // catch(LatticeException &e){
    // 	std::cout << e.what() << ": "  << e.getMessage() << std::endl;
    // 	exit(1);
    // }	
    // catch(std::exception &e){
    // 	std::cout << "caught exception " << e.what() << std::endl;
    // }
    

    // try{
    // 	XmlHandler handler(distSink,"NoXmlHandler");
    // 	logger->debug("load XML file");
    // 	handler.load("../xml/elements.xml");
    // 	logger->debug("save XML file");
    // 	handler.set(simbox, table);
    // 	handler.save("../xml/test_out2.xml");
	
    // }
    // catch(XmlException &e){
    // 	std::cout << e.what() << ": "  << e.getMessage() << std::endl;
    // }	
    // catch(std::exception &e){
    // 	std::cout << "caught exception " << e.what() << std::endl;
    // }

    // try{
    // 	XmlHandler handler2(distSink,"NoXmlHandler");
    // 	handler2.load("../xml/test_out2.xml");
    // 	logger->debug("importing Simbox");
    // 	std::shared_ptr<SimulationBox> simbox2 =handler2.get(distSink, table);
    // 	simbox2->createZincblende(Vector3D<indexType>(1,1,4), mat1);
    // 	handler2.set(*simbox2, table);
    // 	handler2.save("../xml/test_out3.xml");

    // 	logger->debug("export Periodic table to get IDs");
    // 	XmlHandler handler3(distSink,"NoXmlHandler");
    // 	handler3.set(table);
    // 	handler3.save("../xml/test_out4.xml");

    // 	logger->debug("loading exported table to check numbers");
    // 	XmlHandler handler4(distSink,"NoXmlHandler");
    // 	PeriodicTable table2(distSink);
    // 	handler4.load("../xml/test_out4.xml");
    // 	handler4.get(table2);

    // 	logger->debug("loading tersoff potentials");
    // 	XmlHandler handler5(distSink,"NoXmlHandler");
    // 	handler5.load("../xml/tersoff_potentials.xml");
    // 	TersoffPotential potential(distSink);
    // 	handler5.get(potential);

    // 	logger->debug("reexporting tersoff potentials");
    // 	XmlHandler handler6(distSink,"NoXmlHandler");
    // 	handler6.set(potential);
    // 	handler6.save("../xml/tersoff_potentials2.xml");

    // 	logger->debug("trying to calculate energy");
    // 	Range3D<indexType> range;
    // 	double energy = simbox.getEnergy(potential, range);
    // 	logger->info("energey = {}",energy);
	
    // }
    // catch(XmlException &e){
    // 	std::cout << e.what() << ": "  << e.getMessage() << std::endl;
    // }	
    // catch(std::exception &e){
    // 	std::cout << "caught exception " << e.what() << std::endl;
    // }

    // Range3D<indexType> range;
    // logger->debug("range start value {}", range.start.str());

    // range.stop = {3,4,5};
    // range.apply = {false,true,true};
    // logger->debug("range stop value {}", range.stop.str());
    // logger->debug("range apply value {}", range.apply.str());
    
    // indexType start[3]= {3,5,7};
    // indexType stop[3]= {445,33,2340};
    // bool apply[3]= {false,false,false};

    // Range3D<indexType> range2(start, stop, apply);
    // logger->debug("range2 stop value {}", range2.stop.str());

    // Range3D<indexType> range3(Vector3D<indexType>(1,2,3),
    // 			     Vector3D<indexType>(666,3,566),
    // 			     Vector3D<bool>(true,false,false));
    // logger->debug("range3 stop value {}", range3.stop.str());

    
    // logger->debug("trying Optimization object");
    // simbox.generateNeighbors(0,0.672);
    // Optimization opt(distSink);
    // opt.setLoggerLevel(spdlog::level::trace);    
    // opt.registerAction(std::string("MMC"),&SimulationBox::mmcRelax, 100);
    // opt.registerAction(std::string("lattice scaling"),
    // 		       &SimulationBox::scale, 100);

    // // OptimizationFunction f = std::bind(&SimulationBox::mmcRelax,&simbox,
    // // 				       std::placeholders::_1,
    // // 				       std::placeholders::_2,
    // // 				       std::placeholders::_3);
    // // OptimizationFunction f2 = std::bind(&SimulationBox::scale,&simbox,
    // // 				       std::placeholders::_1,
    // // 				       std::placeholders::_2,
    // // 				       std::placeholders::_3);
    
    // // opt.registerAction(f, 100);
    // // opt.registerAction(f2, 100);
    
    // XmlHandler handler5(distSink,"NoXmlHandler");
    // handler5.load("../xml/tersoff_potentials.xml");
    // TersoffPotential potential(distSink);
    // handler5.get(potential);
        
    // Range3D<indexType> rangeAll;
    // OptimizationParameter optParam;
    // optParam.mmcMaxDisplacement=1e-5;
    // optParam.maxScaling=0.2;
    // optParam.runCount = 3;
    // simbox.setLoggerLevel(spdlog::level::info);    
    // opt.runStatic(simbox, rangeAll, potential, optParam);
    // opt.printStatistic();
    // opt.clearStatistic();
    // opt.printStatistic();
    
    // try{
    // 	XmlHandler handler6(distSink);
    // 	handler6.load("../xml/materials.xml");
    // 	MaterialCollection collection(distSink);
    // 	handler6.get(collection, table, distSink);

    // 	collection.setLoggerLevel(spdlog::level::debug);
	
    // 	try{
    // 	    collection.getByCationsAnions({8481}, {12593});
    // 	}catch(MaterialException &e){
    // 	    std::cout << e.what() << ": "  << e.getMessage() << std::endl;
    // 	    collection.getByCationsAnions({12593}, {8481});
    // 	}

    // 	XmlHandler handler7(distSink);
    // 	handler7.set(collection, table);
    // 	handler7.save("../xml/materials2.xml");
    // }catch(MaterialException &e){
    // 	std::cout << e.what() << ": "  << e.getMessage() << std::endl;
    // }	
    // catch(std::exception &e){
    // 	std::cout << "caught exception " << e.what() << std::endl;
    // }
    
    // logger->debug("finished test app");
    // logger->trace("[END OF LOGGING]\n\n");

    // //spdlog::set_async_mode(8192);
    // /*    
    // std::vector<spdlog::sink_ptr> sinks;
    // sinks.push_back(std::make_shared<spdlog::sinks::stdout_sink_st>());
    // sinks.push_back(std::make_shared<spdlog::sinks::simple_file_sink_st>
    // 		    ("../tmp/log2"));
    
    // auto logger = std::make_shared<spdlog::logger>("main", begin(sinks),
    // 						      end(sinks));
    // auto xmlLogger = std::make_shared<spdlog::logger>("xml", begin(sinks),
    // 						      end(sinks));
    // auto elLogger = spdlog::stdout_color_mt("Element");
    // logger->info("Welcome");
    // xmlLogger->info("Welcome");
    // elLogger->error("Welcome");    

    // //auto dist_sink = std::make_shared<spdlog::sinks::dist_sink_st>();
    // auto logger1 = std::make_shared<spdlog::logger>("1st", dist_sink);
    // auto logger2 = std::make_shared<spdlog::logger>("2nd", dist_sink);
    // auto logger3 = std::make_shared<spdlog::logger>("2nd", dist_sink);

    // std::cout << typeid(dist_sink).name() << std::endl;
    
    // logger1->set_level(spdlog::level::debug);
    // logger2->set_level(spdlog::level::critical);
    
    // dist_sink->add_sink(sinks[0]);
    // dist_sink->add_sink(sinks[1]);
    
    // logger1->debug("aaa");
    // logger2->debug("bbb");
    // logger2->log(spdlog::level::off,"bbb");
    // logger3->info("bbb");    
    // */
}
