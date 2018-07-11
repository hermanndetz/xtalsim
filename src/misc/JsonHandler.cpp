/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#include "JsonHandler.h"

//! \param logName Optional parameter defining the logger's name.
JsonHandler::JsonHandler(const char *logName): logName_(logName)
{
    // nothing to be done
}

//------------------------------------------------------------------------------

JsonHandler::~JsonHandler() { }

//------------------------------------------------------------------------------

void JsonHandler::clear(void)
{
    jsonDoc_.clear();
}

//------------------------------------------------------------------------------

//! Drops internal DOM structure and creates new one from specified file.
//! \param fileName Path to JSON file that shall be read.
void JsonHandler::load(const std::string &fileName)
{

    std::ifstream jsonFile (fileName);

    jsonDoc_.clear();
    jsonFile >> jsonDoc_;
    jsonFile.close();
}

//------------------------------------------------------------------------------

//! \param fileName Path to JSON file which shall be written.
void JsonHandler::save(const std::string &fileName) const
{

    std::ofstream jsonFile (fileName, std::ios::out);
    
    jsonFile << jsonDoc_;
    jsonFile.close();
    
}

//------------------------------------------------------------------------------

//! \param periodicTable Object storing elements that shall be written.
void JsonHandler::set(const PeriodicTable &periodicTable)
{

    CLOG(TRACE, logName_) << "writing Elements to JSON";

    std::vector<Element> elements = periodicTable.get();
    Json::Value tableNode(Json::arrayValue);
    unsigned int counter=0;
    
    // remove already existing node with the same name
    jsonDoc_.removeMember("periodic_table");
        
    for (auto el: elements){
	Json::Value elementNode(Json::objectValue);
		
	elementNode["name"] = el.name;
	elementNode["symbol"] = el.symbol;
	elementNode["proton"] = el.protonCount;
	elementNode["neutron"] = el.neutronCount;
	elementNode["period"] = el.period;
	elementNode["group"] = el.group;	
	elementNode["mass"] = el.mass;
	elementNode["weight"] = el.weight;
	elementNode["color"] = el.color.str();

	counter++;
	tableNode.append(elementNode);
    }

    jsonDoc_["periodic_table"]=tableNode;

    CLOG(DEBUG, logName_) << counter << " Elements written to JSON";
    CLOG(TRACE, logName_) << "Elements successfully written to JSON";
}

//------------------------------------------------------------------------------

//! For each single element an Element object is created an inserted into the
//! specified object. If a value was not specified in
//! the XML value a default value is stored.
//! \param periodicTable Object accepting Element objects.
void JsonHandler::get(PeriodicTable &periodicTable) const
{

    CLOG(TRACE, logName_) << "generating Periodic Table from JSON";
    
    unsigned int counter=0;
    Json::Value elements = jsonDoc_["periodic_table"];

    for (unsigned int i=0; i < elements.size(); i++, counter++){
	
	Element element(elements[i].get("name","undef").asString(),
			elements[i].get("symbol","undef").asString(),
			(uint8_t)elements[i].get("proton",0).asInt(),
			(uint8_t)elements[i].get("neutron",0).asInt(),
			(uint8_t)elements[i].get("period",0).asInt(),
			(uint8_t)elements[i].get("group",0).asInt(),
			elements[i].get("mass",0).asDouble(),
			elements[i].get("weight",0).asDouble(),
			elements[i].get("color","undef").asString()
			);
	
	periodicTable.add(element);
    }

    CLOG(DEBUG, logName_) << counter <<
	" Elements written to Periodic Table from JSON";
    CLOG(TRACE, logName_) << "Periodic Table successfully generated";
    
}

//------------------------------------------------------------------------------

//! Transfers all Tersoff Parameters stored in the potential to the DOM
//! structure. 
//! \param potential Tersoff potential storing the
//! single parameter objects.
void JsonHandler::set(const TersoffPotential &potential)
{

    const PeriodicTable &pt = PeriodicTable::getInstance();
    
    CLOG(TRACE, logName_) << "writing Tersoff Parameters to JSON";

    Json::Value tersoffNode (Json::arrayValue);
    int counter = 0;

    // remove already existing nodes with the same name
    jsonDoc_.removeMember("tersoff");

    for(auto it = potential.get().begin(); it != potential.get().end();
	++it, counter++){
	Json::Value potentialNode(Json::objectValue);

	potentialNode["name"] = it->second.name;

	const Element *element = &pt.getById(it->second.elementIdLow);
	potentialNode["element1"] = element->symbol;	
	potentialNode["proton1"] = element->protonCount;
	potentialNode["neutron1"] = element->neutronCount;
	
	element = &pt.getById(it->second.elementIdHigh);
	potentialNode["element2"] = element->symbol;	
	potentialNode["proton2"] = element->protonCount;
	potentialNode["neutron2"] = element->neutronCount;
	
	potentialNode["mode"] = it->second.mode;
	potentialNode["beta"] = it->second.beta;
	potentialNode["delta"] = it->second.delta;
	potentialNode["gamma"] = it->second.gamma;
	potentialNode["lambda"] = it->second.lambda;
	potentialNode["lambda3"] = it->second.lambda3;
	potentialNode["mu"] = it->second.mu;
	potentialNode["A"] = it->second.a;
	potentialNode["B"] = it->second.b;
	potentialNode["C"] = it->second.c;
	potentialNode["D"] = it->second.d;
	potentialNode["D0"] = it->second.d0;
	potentialNode["DD"] = it->second.dd;
	potentialNode["H"] = it->second.h;
	potentialNode["N"] = it->second.n;
	potentialNode["Rcut"] = it->second.rCut;
	potentialNode["R0"] = it->second.r0;
	potentialNode["S"] = it->second.s;

	tersoffNode.append(potentialNode);
    }

    jsonDoc_["tersoff"] = tersoffNode;

    CLOG(DEBUG, logName_) << counter << " Tersoff Parameters written to JSON";
    CLOG(TRACE, logName_) << "Tersoff Potential succesfully written";
}

//------------------------------------------------------------------------------
    
//! Read Tersoff Parameters from DOM structure and add them to the specified
//! Tersoff Potential.
//! \param potential Tersoff potential accepting the
//! stored parameter objects.
void JsonHandler::get(TersoffPotential &potential) const
{

    const PeriodicTable &pt = PeriodicTable::getInstance();
    unsigned int counter=0;
    uint8_t proton1, proton2;
    uint8_t neutron1, neutron2;
    elementType elementId1, elementId2;

    CLOG(TRACE, logName_) << "generating Tersoff Parameters from JSON";
    
    Json::Value tersoff = jsonDoc_["tersoff"];
    
    for (unsigned int i=0; i < tersoff.size(); i++, counter++){

	proton1 = (uint8_t)tersoff[i].get("proton1",0).asInt();
	neutron1 = (uint8_t)tersoff[i].get("neutron1",0).asInt();
	elementId1 = pt.getByProtonNeutron(proton1, neutron1).id;

	proton2 = (uint8_t)tersoff[i].get("proton2",0).asInt();
	neutron2 = (uint8_t)tersoff[i].get("neutron2",0).asInt();
	elementId2 = pt.getByProtonNeutron(proton2, neutron2).id;
	
	TersoffParameter parameter(elementId1, elementId2,
				   tersoff[i].get("name","undef").asString(),
				   tersoff[i].get("mode","undef").asString(),
				   tersoff[i].get("A",0).asDouble(),
				   tersoff[i].get("B",0).asDouble(),
				   tersoff[i].get("beta",0).asDouble(),
				   tersoff[i].get("C",0).asDouble(),
				   tersoff[i].get("D",0).asDouble(),
				   tersoff[i].get("D0",0).asDouble(),
				   tersoff[i].get("DD",0).asDouble(),
				   tersoff[i].get("delta",0).asDouble(),
				   tersoff[i].get("gamma",0).asDouble(),
				   tersoff[i].get("H",0).asDouble(),
				   tersoff[i].get("lambda",0).asDouble(),
				   tersoff[i].get("lambda3",0).asDouble(),
				   tersoff[i].get("mu",0).asDouble(),
				   tersoff[i].get("N",0).asDouble(),
				   tersoff[i].get("R0",0).asDouble(),
				   tersoff[i].get("Rcut",0).asDouble(),
				   tersoff[i].get("S",0).asDouble()
				   );
	
	potential.add(parameter);
    }

    CLOG(DEBUG, logName_) << counter <<" Tersoff Parameters generated from JSON";
    CLOG(TRACE, logName_) << "Tersoff Parameters successfully generated";
 
}

//------------------------------------------------------------------------------

//! The periodic table is required here to display the
//! elements by their symbol (e.g. Ga, In, As) and not by the internally used
//! ID.
//! \param collection Material collection storing the single materials.
void JsonHandler::set(const MaterialCollection &collection)
{

    CLOG(TRACE, logName_) << "writing materials to JSON";

    const PeriodicTable &pt = PeriodicTable::getInstance();
    int counter=0;
    Json::Value collectionNode (Json::arrayValue);
    MaterialComponents component;
    
    // remove already existing nodes with the same name
    jsonDoc_.removeMember("material-collection");

    for(auto it = collection.get().begin(); it != collection.get().end();
	++it, counter++){

	Json::Value materialNode(Json::objectValue);
	std::vector<double> elasticConstants = it->getElasticConstants();

	
	materialNode["name"] = it->getName();
	materialNode["lattice-constant"] = it->getLatticeConstant();
	
	materialNode["c11"] = elasticConstants[0];
	materialNode["c12"] = elasticConstants[1];
	materialNode["c44"] = elasticConstants[2];

	component = it->getCations();

	Json::Value cationsNode(Json::arrayValue);

    double lastShare=0;
	for (auto it2 = component.begin(); it2 != component.end(); ++it2){
	    Json::Value cation (Json::objectValue);

        auto element = pt.getById(std::get<0>(*it2));
	    cation["element"] = element.symbol;
	    cation["proton"] = element.protonCount;
	    cation["neutron"] =	element.neutronCount;
	    cation["share"] = std::get<1>(*it2)-lastShare;
        lastShare = std::get<1>(*it2);

	    cationsNode.append(cation);
	}

	materialNode["cations"] = cationsNode;

	component = it->getAnions();

	Json::Value anionsNode(Json::arrayValue);

    lastShare=0;
	for (auto it2 = component.begin(); it2 != component.end(); ++it2){
	    Json::Value anion (Json::objectValue);

        auto element = pt.getById(std::get<0>(*it2));
	    anion["element"] = element.symbol;
	    anion["proton"] = element.protonCount;
	    anion["neutron"] =	element.neutronCount;
	    anion["share"] = std::get<1>(*it2)-lastShare;
        lastShare = std::get<1>(*it2);
        
	    anionsNode.append(anion);
	}

	materialNode["anions"] = anionsNode;

	collectionNode.append(materialNode);
    }

    jsonDoc_["material-collection"] = collectionNode;

    CLOG(DEBUG, logName_) << counter <<" Materials written to JSON";
    CLOG(TRACE, logName_) << "Material succesfully written";
}

//------------------------------------------------------------------------------

//! Creates a Material object and fills it with data from the DOM structure.
//! The periodic table is required here to allow the user to specify the
//! elements by their symbol (e.g. Ga, In, As) and not by the internally used
//! ID.
//! \param collection Material collection that shall be filled from DOM.
//! \throws PeriodicTableException
void JsonHandler::get(MaterialCollection &collection)
    const
{

    int counter=0;
    const PeriodicTable &pt = PeriodicTable::getInstance();
    
    CLOG(TRACE, logName_) << "generating Materials from JSON";

    Json::Value materials = jsonDoc_["material-collection"];
    
    for (unsigned int i=0; i < materials.size(); i++, counter++){
    
	MaterialComponents cations, anions;
	Json::Value jsonCations = materials[i]["cations"];

	for(unsigned int j=0; j< jsonCations.size(); j++){
	    
	    cations.push_back (std::make_tuple (
		      pt.getByProtonNeutron(
			    (uint8_t)jsonCations[j].get("proton",0).asInt(),
			    (uint8_t)jsonCations[j].get("neutron",0).asInt()
					       ).id,
		      jsonCations[j].get("share",1).asDouble()
							  )
			    );
	}
	
	Json::Value jsonAnions = materials[i]["anions"];
	
	for(unsigned int j=0; j< jsonAnions.size(); j++){
	    
	    anions.push_back (std::make_tuple (
		      pt.getByProtonNeutron(
			    (uint8_t)jsonAnions[j].get("proton",0).asInt(),
			    (uint8_t)jsonAnions[j].get("neutron",0).asInt()
					       ).id,
		      jsonAnions[j].get("share",1).asDouble()
							  )
			    );
	}

	Material material(
			  materials[i].get("name","undef").asString(),
			  cations, anions,
			  materials[i].get("lattice-constant",0).asDouble(),
			  materials[i].get("c11",0).asDouble(),
			  materials[i].get("c12",0).asDouble(),
			  materials[i].get("c44",0).asDouble()
			  );
	
	collection.add(material);
	
    }

    CLOG(DEBUG, logName_) << counter <<" Materials generated from JSON";
    CLOG(TRACE, logName_) << "Materials successfully generated";
   
}

//------------------------------------------------------------------------------

//! Transfers content of the lattice into the DOM structure. The ID, (i,j,k) and
//! (x,y,z) values are exported using the PCDATA format.
//! \param lattice Atomic lattice.
void JsonHandler::set(const Lattice &lattice)
{

    CLOG(TRACE, logName_) << "writing Lattice to JSON";

    const PeriodicTable &pt = PeriodicTable::getInstance();
    Json::Value latticeNode(Json::objectValue), sizeNode(Json::arrayValue);
    Json::Value gridNode(Json::arrayValue);
    std::vector<Atom *> atoms;
    Vector3D<indexType> size = lattice.getSize();
    int counter=0;
    
    // remove already existing nodes with the same name
    jsonDoc_.removeMember("lattice");
    
    latticeNode["temperature"] = lattice.getTemperature();

    sizeNode.resize(3);
    for (int i=0; i<3; i++) sizeNode[i] = size[i];
    latticeNode["size"] = sizeNode;

    for (auto atom: lattice.getAtomList()){

	Json::Value atomNode(Json::objectValue);
	Vector3D<spaceType> position = atom->getPosition();
	Vector3D<indexType> index = atom->getIndex();

	const Element &element = pt.getById(atom->getElementId());
	// export as PCDATA nodes
        atomNode["element"] = element.symbol;
	atomNode["proton"] = element.protonCount;
	atomNode["neutron"] = element.neutronCount;
	atomNode["material"] = atom->getMaterial()->getName();
	atomNode["modified"] = atom->wasModified();
    atomNode["modification-index"] = atom->getModificationOrder();

	Json::Value indexNode(Json::arrayValue);
	indexNode.resize(3);
	for (int i=0; i<3; i++) indexNode[i] = index[i];
	atomNode["index"] = indexNode;

	Json::Value positionNode(Json::arrayValue);
	positionNode.resize(3);
	for (int i=0; i<3; i++) positionNode[i] = position[i];
	atomNode["position"] = positionNode;

	counter++;
	gridNode.append(atomNode);
	
    }

    latticeNode["atoms"] = gridNode;
    jsonDoc_["lattice"] = latticeNode;

    CLOG(DEBUG, logName_) << counter <<" Atoms written to JSON";
    CLOG(TRACE, logName_) << "Lattice successfully written";
}

//------------------------------------------------------------------------------

//! Creates a lattice of size specified in the DOM structure and fills the
//! lattice with the stored atoms defined also in the DOM structure.
//! \param lattice Atomic lattice.
//! \param materials Materials used in the simulation box.
void JsonHandler::get(Lattice &lattice, const MaterialCollection &materials) const
{

    unsigned int counter=0;
    Vector3D<spaceType> position;
    Vector3D<indexType> index;
    Json::Value workNode;
    const PeriodicTable &pt = PeriodicTable::getInstance();
    
    CLOG(TRACE, logName_) << "generating Lattice from JSON";

    workNode = jsonDoc_["lattice"];
    lattice.setTemperature(workNode.get("temperature",0).asDouble());

    // 'u' after number required here to tell compiler to use integer parameter
    // function, not char* one.
    lattice.generate((indexType)workNode["size"].get(0u,0).asInt(),
		     (indexType)workNode["size"].get(1u,0).asInt(),
		     (indexType)workNode["size"].get(2u,0).asInt()
		     );
    
    workNode = workNode["atoms"];
    
    for (unsigned int i=0; i < workNode.size(); i++, counter++){

	Json::Value atomNode = workNode[i];

	for (int j=0; j<3; j++){
	    position[j] = atomNode["position"].get(j,0).asDouble();
	    index[j] = (indexType)atomNode["index"].get(j,0).asInt();
	}

	uint8_t proton = (uint8_t)atomNode.get("proton",0).asInt();
	uint8_t neutron =(uint8_t)atomNode.get("neutron",0).asInt();
	std::string name = atomNode.get("material","").asString();
	bool modified = atomNode.get("modified",false).asBool();
    uint32_t modificationIndex = atomNode.get("modification-index",0).asInt();

	lattice(index) =
	    new Atom(pt.getByProtonNeutron(proton,neutron).id, position,
		     index, &(materials.getByName(name)), modified, false, modificationIndex );

	CLOG(DEBUG, logName_) << name << " Material Name " <<
	    materials.getByName(name).getName();

    }

    CLOG(DEBUG, logName_) << counter <<" Atoms generated from JSON";
    CLOG(TRACE, logName_) << "Lattice successfully generated";

}

//------------------------------------------------------------------------------

//! Transfers content of the complete simulation box into the DOM
//! structure. First global properties of the simulation box are imported. Later
//! the function to write the Lattice is called.
//! \param simbox SimulationBox that shall be written.
void JsonHandler::set(const SimulationBox &simbox)
{

    CLOG(TRACE, logName_) << "writing Simulation Box to JSON";

    Json::Value simboxNode(Json::objectValue);

    // remove already existing nodes with the same name
    jsonDoc_.removeMember("simbox");
    
    simboxNode["outOfPlaneDimension"] =	simbox.getOutOfPlaneDimension();
    simboxNode["inPlaneLatticeConstant"] =	simbox.getInPlaneLatticeConstant();
    simboxNode["description"] = simbox.getDescription();
    
    Json::Value sizeNode(Json::arrayValue);
    sizeNode.resize(3);
    for (int i=0; i<3; i++) sizeNode[i] = simbox.getSize()[i];
    simboxNode["size"] = sizeNode;
    
    jsonDoc_["simbox"]=simboxNode;
    
    this->set(simbox.getMaterials());
    this->set(simbox.getLattice());
    
    CLOG(TRACE, logName_) << "Simulation Box successfully written";
    
}

//------------------------------------------------------------------------------

//! Creates a simbox object from the stored data and returns a pointer. This
//! method was chosen as some parameters of the SimulationBox have to be set
//! in the constructor since they are not supposed to change. Therefore it was
//! neglected to create a member function to set these values.
//! \param simbox SimulationBox that shall be filled.
void JsonHandler::get(std::shared_ptr<SimulationBox> &simbox) const
{
    Vector3D<spaceType> size;

    CLOG(TRACE, logName_) << "generating Simulation Box from JSON";
    
    Json::Value simboxNode=jsonDoc_["simbox"];
    
    simbox = std::make_shared<SimulationBox>(
	      (uint8_t)simboxNode.get("outOfPlaneDimension",0).asInt(),
	      Vector3D<spaceType>(
				  simboxNode["size"].get(0u,0).asDouble(),
				  simboxNode["size"].get(1u,0).asDouble(),
				  simboxNode["size"].get(2u,0).asDouble()),
	      simboxNode.get("inPlaneLatticeConstant",-1).asDouble()
                                             );
    simbox->setDescription(simboxNode.get("description","").asString());
    
    this->get(simbox->getMaterials());
    this->get(simbox->getLattice(), simbox->getMaterials());
    // Force creation of neighbor list.
    simbox->generateNeighbors(true);

    CLOG(TRACE, logName_) << "Simulation Box successfully generated";
   
}

//------------------------------------------------------------------------------

//! \param composition CompositionInfo that shall be transferred.
void JsonHandler::set(const CompositionInfo &composition)
{

    CLOG(TRACE, logName_) << "writing Composition Infos to JSON";

    const PeriodicTable &pt = PeriodicTable::getInstance();    
    Json::Value compNode(Json::arrayValue);
    int counter = 0;
  
    for (auto i=0; i<composition.getLayerCount(); i++){

        Json::Value layerNode(Json::objectValue);
        LayerCompositionInfo layer = composition.getLayer(i);

        layerNode["index"] = i;
        
        Json::Value elementNode(Json::arrayValue);
        for (auto entry : layer) {

            Json::Value entryNode(Json::objectValue);
            
            auto key = std::get<0>(entry);
        
            entryNode["material"] = std::get<0>(key).c_str();
            entryNode["id"] = std::get<1>(key);
            entryNode["name"] = pt.getById(std::get<1>(key)).symbol.c_str();
            entryNode["count"] = std::get<1>(entry);

            elementNode.append(entryNode);
        }

        layerNode["elements"] = elementNode;

        counter++;
        compNode.append(layerNode);
    }

    auto workNode = jsonDoc_["compositions"];

    if (workNode.type() == Json::nullValue)
        workNode = Json::Value(Json::arrayValue);
    
    workNode.append(compNode);
    jsonDoc_["compositions"] = workNode;
    
    CLOG(DEBUG, logName_) << counter <<" atomic layers written to JSON";
    CLOG(TRACE, logName_) << "CompositionInfo successfully written";    
}

//------------------------------------------------------------------------------

//! \param compositions Vector where read compositions shall be stored.
void JsonHandler::get(std::vector<CompositionInfo> &compositions) const
{
    unsigned int counter=0;

    CLOG(TRACE, logName_) << "generating composition Infos from JSON";

    Json::Value compNode = jsonDoc_["compositions"];

    for (auto i=0; i<compNode.size(); i++){

        Json::Value layersNode = compNode[i];
        CompositionInfo comp;

        for (auto j=0; j<layersNode.size(); j++){

            LayerCompositionInfo layer;
            Json::Value elements = layersNode[j]["elements"];

            for (auto k=0; k<elements.size(); k++){

                auto key = std::make_tuple(elements[k].get("material","undef").asString(),
                                           (elementType)elements[k].get("id",0).asInt() );

                layer[key] = elements[k].get("count",0).asInt();
            }

            comp << layer;
        }
        
        compositions.push_back(comp);
    }

    CLOG(TRACE, logName_) << "CompositionInfos successfully generated";
    
}

//------------------------------------------------------------------------------

//! \param journal Journal the entries shall be written to.
//! \param entry First entry node of XML <journal> environment
//! \return Number of read entries.
int JsonHandler::readEntries(Journal<int> &journal, Json::Value entries) const
{
    int counter=0;
    
    for(unsigned int i=0; i<entries.size(); i++, ++counter){
	journal.add(entries.get(i,0).asInt());
    }

    return counter;
}

//------------------------------------------------------------------------------

//! \param journal Journal the entries shall be written to.
//! \param entry First entry node of XML <journal> environment
//! \return Number of read entries.
int JsonHandler::readEntries(Journal<double> &journal, Json::Value entries)
    const
{
    int counter=0;
    
    for(unsigned int i=0; i<entries.size(); i++, ++counter){
	journal.add(entries.get(i,0).asDouble());
    }

    return counter;
}

//------------------------------------------------------------------------------

//! \param journal Journal the entries shall be written to.
//! \param entry First entry node of XML <journal> environment
//! \return Number of read entries.
int JsonHandler::readEntries(Journal<std::string> &journal,
			     Json::Value entries) const
{
    int counter=0;
    
    for(unsigned int i=0; i<entries.size(); i++, ++counter){
	journal.add(entries.get(i,"").asString());
    }

    return counter;
}

//------------------------------------------------------------------------------
    
//! Read Configuration file from DOM structure.
//! \param config Configuration object.
void JsonHandler::get(Configuration &config) const
{

    CLOG(TRACE, logName_) << "generating Configuration from JSON";

    Json::Value configNode = jsonDoc_["config"];
    Json::Value workNode = configNode["material"];
    
    for(unsigned int i=0; i< workNode["cations"].size(); i++){
	Json::Value cation = workNode["cations"][i];

	config.cations.push_back( std::make_tuple(
		  cation.get("element","undef").asString(),
		  cation.get("share",0).asDouble() )
				  );
    }
    
    for(unsigned int i=0; i< workNode["anions"].size(); i++){
	Json::Value anion = workNode["anions"][i];

	config.anions.push_back( std::make_tuple(
		  anion.get("element","undef").asString(),
		  anion.get("share",0).asDouble() )
				  );
    }

    config.materialName = workNode.get("name","undef").asString();
    config.c11= workNode.get("c11",0).asDouble();
    config.c12= workNode.get("c12",0).asDouble();
    config.c44= workNode.get("c44",0).asDouble();
    config.latticeConstant= workNode.get("lattice-constant",0).asDouble();
    
    config.materialFile =configNode.get("material-file","").asString();
    config.materialSearchName =	configNode.get("material-name","").asString();
    
    workNode = configNode["size"];
    config.size = {(indexType)workNode.get(0u,0).asInt(),
		   (indexType)workNode.get(1u,0).asInt(),
		   (indexType)workNode.get(2u,0).asInt()};
    
    config.growthDimension =
	(uint8_t)configNode.get("growth-dimension",0).asInt();
    config.inputFileName = configNode.get("input","").asString();
    config.outputFileName = configNode.get("output","").asString();
    config.periodicTableFileName =
	configNode.get("periodic-table","").asString();
    config.logFileName = configNode.get("log","").asString();
    config.tersoffFileName = configNode.get("tersoff","").asString();
    config.xyzFileName = configNode.get("xyz","").asString();
    config.outputPreamble = configNode.get("output-preamble","").asString();
    config.journalPreamble = configNode.get("journal-preamble","").asString();
    config.neighborRadius = configNode.get("neighbor-radius",0).asDouble();
    config.neighborLayers =
	(indexType)configNode.get("neighbor-layers",0).asInt();
    config.latticeTemperature = 
	(indexType)configNode.get("temperature",0).asDouble();
	
    config.startIndex = configNode.get("start-index",-1).asInt();
    config.stopIndex = configNode.get("stop-index",-1).asInt();
    config.mmcProbability = configNode.get("mmc-probability",-1).asDouble();
    config.mmcRunCount = configNode.get("mmc-run-count",1).asInt();
    config.minDisplacement = configNode.get("min-displacement",0).asDouble();
    config.maxDisplacement = configNode.get("max-displacement",-1).asDouble();
    config.scalingProbability =
	configNode.get("scaling-probability",-1).asDouble();
    config.minScaling =	configNode.get("min-scaling",0).asDouble();
    config.maxScaling =	configNode.get("max-scaling",-1).asDouble();
    config.runCount = configNode.get("runs",-1).asInt();
    config.checkCount = configNode.get("check-runs",1).asInt();
    config.energyDropFactor=configNode.get("energy-drop",2).asDouble();
    config.reductionFactor = configNode.get("reduction",2).asDouble();
    config.minEnergy = configNode.get("min-energy",0).asDouble();
    config.anionPassivationProbability = configNode.get("anion-passivation-probability",1).asDouble();
    config.maxThreadCount=(uint16_t)configNode.get("max-thread-count",1).asInt();

    config.modificationIndexMax = configNode.get("max-modification-index",0).asDouble();
    workNode = configNode["camera-focal-point"];
    config.cameraFocalPoint = {workNode.get(0u,0).asDouble(),
		 workNode.get(1u,0).asDouble(),
		 workNode.get(2u,0).asDouble()};
    workNode = configNode["camera-direction"];
    config.cameraDirection = {workNode.get(0u,0).asDouble(),
		 workNode.get(1u,0).asDouble(),
		 workNode.get(2u,0).asDouble()};
    config.cameraZoom = configNode.get("camera-zoom",1).asDouble();
    config.animationStep = (uint32_t)configNode.get("animation-step",1).asInt();
    config.animationFilePrefix = configNode.get("animation-prefix","").asString();
    config.screenshotFileName = configNode.get("screenshot","").asString();
    config.rotationCount = (uint16_t)configNode.get("animation-step",0).asInt();
    
    config.interfaceAtoms = configNode.get("interface-atoms",0).asInt();
    config.interfacePositions =	configNode.get("interface-positions",0).asInt();

    config.refLatticeConstant =configNode.get("ref-lattice-const",0).asDouble();
    
    workNode = configNode["switches"];    

    for(unsigned int i=0; i< workNode.size(); i++){
        std::string text(workNode[i].asString());
        
        if (text == "extend"){
            config.extend=true;
            continue;
        }

        if (text == "calc-energy"){
            config.calculateEnergy=true;
            continue;
        }

        if (text == "quiet"){
            config.quiet=true;
            continue;
        }
        
        if (text == "verbose"){
            config.verbose=true;
            continue;
        }

        if (text == "print-statistics"){
            config.print=true;
            continue;
        }

        if (text == "render-bonds") {
            config.renderBonds=true;
            continue;
        }

        if (text == "visualize-bond-strain") {
            config.visualizeBondStrain=true;
            continue;
        }

        if (text == "static-optimization"){
            config.staticOptimization=true;
            continue;
        }

        if (text == "dynamic-optimization"){
            config.dynamicOptimization=true;
            continue;
        }

        if (text == "require-matching-cation"){
            config.requireMatchingCation=true;
            continue;
        }

        if (text == "modfied-only"){
            config.modifiedOnly=true;
            continue;
        }

        if (text == "modified-negative"){
            config.modifiedNegative=true;
            continue;
        }

        if (text == "animate-modification"){
            config.modificationAnimation=true;
            continue;
        }

        if ((text == "filter-modified-interface") || 
                (text == "filter-modified-exchange")){
            config.filterModificationState.push_back(text);
        }
    }
    
    CLOG(TRACE, logName_) << "Configuration successfully generated";
    
}

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
