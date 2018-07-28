/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#include "XmlHandler.h"

//! \param logName Optional parameter defining the logger's name.
XmlHandler::XmlHandler(const char *logName): logName_(logName)
{
    // nothing to be done
}

//------------------------------------------------------------------------------

XmlHandler::~XmlHandler() { }

//------------------------------------------------------------------------------

void XmlHandler::clear(void)
{
    xmlDoc_.reset();
}

//------------------------------------------------------------------------------

//! Drops internal DOM structure and creates new one from specified file. The
//! result code is catched and translated into Exception-based notification.
//! \param fileName Path to file that shall be read.
//! \throws XmlException
void XmlHandler::load(const std::string &fileName)
{

    CLOG(TRACE, logName_) << "reading XML file '" << fileName << "'";
    
    // use flag pugi::parse_trim_pcdata which trims white spaces??
    pugi::xml_parse_result result = xmlDoc_.load_file(fileName.c_str());

    switch(result.status){
    case (pugi::status_ok): // no error
	break;
    case (pugi::status_file_not_found): // file could not be opened
	throw XmlException(result.description(),fileName,
			   XmlException::Id::FileNotFound);
	break;
    case (pugi::status_out_of_memory): // not enough memory
	throw XmlException(result.description(), fileName,
			   XmlException::Id::OutOfMemory);
	break;
    case (pugi::status_unrecognized_tag): // wrong tag
	throw XmlException(result.description(), fileName,
			   XmlException::Id::Tag);
	break;
    case (pugi::status_bad_pi): // incorrect document declaration/process instr.
	throw XmlException(result.description(),fileName,
			   XmlException::Id::Declaration);
	break;
    case (pugi::status_bad_pcdata): // invalid plain char
	throw XmlException(result.description(), fileName,
			   XmlException::Id::PcData);
	break;
    case (pugi::status_bad_start_element): // invalid start tag
	throw XmlException(result.description(), fileName,
			   XmlException::Id::StartElement);
	break;
    case (pugi::status_bad_attribute): // invalid attribute
	throw XmlException(result.description(), fileName,
			   XmlException::Id::Attribute);
	break;
    case (pugi::status_bad_end_element): // invalid end tag
	throw XmlException(result.description(), fileName,
			   XmlException::Id::EndElement);
	break;
    case (pugi::status_end_element_mismatch): // start and end tag not matching
	throw XmlException(result.description(), fileName,
			   XmlException::Id::Mismatch);
	break;
    case (pugi::status_no_document_element): // no element node
	throw XmlException(result.description(), fileName,
			   XmlException::Id::NoElement);
	break;
    default:
	CLOG(ERROR, logName_) << "import failed";
	throw XmlException(result.description(), fileName,
			   XmlException::Id::Unknown);
	break;
    }

    CLOG(TRACE, logName_) << "import finished successfully";
    
}
//------------------------------------------------------------------------------

//! \param fileName Path to XML file which shall be written.
//! \throws XmlException
void XmlHandler::save(const std::string &fileName) const
{

    CLOG(TRACE, logName_) << "saving XML data to file" << fileName;
    
    bool result = xmlDoc_.save_file(fileName.c_str()); 

    if (!result){
	CLOG(ERROR, logName_) << "writing to file failed";
	throw XmlException("could not write XML file", fileName,
			   XmlException::Id::Save);
    }

    CLOG(DEBUG, logName_) << "XML data successfully written to " << fileName;
    CLOG(TRACE, logName_) << "export finished successfully";
    
}

//------------------------------------------------------------------------------

//! \param periodicTable Object storing elements that shall be written.
void XmlHandler::set(const PeriodicTable &periodicTable)
{

    CLOG(TRACE, logName_) << "writing Elements to XML";

    std::vector<Element> elements = periodicTable.get();
    pugi::xml_node tableNode;
    pugi::xml_node elementNode;
    unsigned int counter=0;
    
    // remove already existing nodes with the same name
    while(xmlDoc_.remove_child("periodic_table")) ;
    
    tableNode = xmlDoc_.append_child("periodic_table");
    
    for (auto el: elements){
	elementNode=tableNode.append_child("element");
		
	// export as PCDATA nodes
	//elementNode.append_child("id").text() = el.id;
	elementNode.append_child("name").text() = el.name.c_str();
	elementNode.append_child("symbol").text() = el.symbol.c_str();
	elementNode.append_child("proton").text() = el.protonCount;
	elementNode.append_child("neutron").text() = el.neutronCount;
	elementNode.append_child("period").text() = el.period;
	elementNode.append_child("group").text() = el.group;	
	elementNode.append_child("mass").text() = el.mass;
	elementNode.append_child("weight").text() = el.weight;
    elementNode.append_child("color").text() = el.color.str().c_str();

	counter++;
    }
    
    CLOG(DEBUG, logName_) << counter << " Elements written to XML";
    CLOG(TRACE, logName_) << "Elements successfully written to XML";
}

//------------------------------------------------------------------------------

//! For each single element an Element object is created an inserted into the
//! specified object. If a value was not specified in
//! the XML value a default value is stored.
//! \param periodicTable Object accepting Element objects.
void XmlHandler::get(PeriodicTable &periodicTable) const
{

    CLOG(TRACE, logName_) << "generating Periodic Table from XML";    
    
    unsigned int counter=0;
       
    pugi::xpath_node_set elNodes =
	xmlDoc_.select_nodes("/periodic_table/element");
    pugi::xpath_node_set::const_iterator it;

    for (it = elNodes.begin(); it != elNodes.end(); ++it, counter++) {
        pugi::xpath_node node=*it;

        Element element(node.node().child("name").text().as_string("undef"),
                node.node().child("symbol").text().as_string("undef"),
                (uint8_t)node.node().child("proton").text().as_int(0),
                (uint8_t)node.node().child("neutron").text().as_int(0),
                (uint8_t)node.node().child("period").text().as_int(0),
                (uint8_t)node.node().child("group").text().as_int(0),
                node.node().child("mass").text().as_double(0),	
                node.node().child("weight").text().as_double(0),
                node.node().child("color").text().as_string("undef"));
        
        periodicTable.add(element);
    }

    CLOG(DEBUG, logName_) << counter <<
	" Elements written to Periodic Table from XML";
    CLOG(TRACE, logName_) << "Periodic Table successfully generated";    
    
}

//------------------------------------------------------------------------------

//! Transfers all Tersoff Parameters stored in the potential to the DOM
//! structure. 
//! \param potential Tersoff potential storing the
//! single parameter objects.
void XmlHandler::set(const TersoffPotential &potential)
{

    CLOG(TRACE, logName_) << "writing Tersoff Parameters to XML";
    const PeriodicTable &pt = PeriodicTable::getInstance();
    
    //    TParam_map parameters = potential.get();
    pugi::xml_node tersoffNode;
    pugi::xml_node potentialNode;
    int counter = 0;

    // remove already existing nodes with the same name
    while(xmlDoc_.remove_child("tersoff")) ;

    tersoffNode=xmlDoc_.append_child("tersoff");

    for (auto it = potential.get().begin(); it != potential.get().end();
        ++it, counter++) {
        potentialNode=tersoffNode.append_child("potential");

        potentialNode.append_child("name").text() = it->second.name.c_str();
        potentialNode.append_child("comment").text() = it->second.comment.c_str();

        const Element *element = &pt.getById(it->second.elementIdLow);
        potentialNode.append_child("element1").text() = element->symbol.c_str();
        potentialNode.append_child("proton1").text() = element->protonCount;
        potentialNode.append_child("neutron1").text() = element->neutronCount;

        element = &pt.getById(it->second.elementIdHigh);
        potentialNode.append_child("element2").text() = element->symbol.c_str();
        potentialNode.append_child("proton2").text() = element->protonCount;
        potentialNode.append_child("neutron2").text() = element->neutronCount;
        
        potentialNode.append_child("mode").text() =
            it->second.mode.c_str();
        potentialNode.append_child("beta").text() = it->second.beta;
        potentialNode.append_child("delta").text() = it->second.delta;
        potentialNode.append_child("gamma").text() = it->second.gamma;
        potentialNode.append_child("lambda").text() = it->second.lambda;
        potentialNode.append_child("lambda3").text() = it->second.lambda3;
        potentialNode.append_child("mu").text() = it->second.mu;
        potentialNode.append_child("A").text() = it->second.a;
        potentialNode.append_child("B").text() = it->second.b;
        potentialNode.append_child("C").text() = it->second.c;
        potentialNode.append_child("D").text() = it->second.d;
        potentialNode.append_child("D0").text() = it->second.d0;
        potentialNode.append_child("DD").text() = it->second.dd;
        potentialNode.append_child("H").text() = it->second.h;
        potentialNode.append_child("N").text() = it->second.n;
        potentialNode.append_child("Rcut").text() = it->second.rCut;
        potentialNode.append_child("R0").text() = it->second.r0;
        potentialNode.append_child("S").text() = it->second.s;
    }

    CLOG(DEBUG, logName_) << counter << " Tersoff Parameters written to XML";
    CLOG(TRACE, logName_) << "Tersoff Potential succesfully written";    
}

//------------------------------------------------------------------------------
    
//! Read Tersoff Parameters from DOM structure and add them to the specified
//! Tersoff Potential.
//! \param potential Tersoff potential accepting the
//! stored parameter objects.
void XmlHandler::get(TersoffPotential &potential)
    const
{

    unsigned int counter=0;
    uint8_t proton1, proton2;
    uint8_t neutron1, neutron2;
    elementType elementId1, elementId2;
    const PeriodicTable &pt = PeriodicTable::getInstance();
    
    CLOG(TRACE, logName_) << "generating Tersoff Parameters from XML";    
    
    pugi::xpath_node_set tersoffNodes =
	xmlDoc_.select_nodes("/tersoff/potential");
    pugi::xpath_node_set::const_iterator it;

    for (it = tersoffNodes.begin(); it != tersoffNodes.end(); ++it, counter++) {
        pugi::xpath_node node=*it;

        proton1 = (uint8_t)node.node().child("proton1").text().as_int(0);
        neutron1 = (uint8_t)node.node().child("neutron1").text().as_int(0);
        elementId1 = pt.getByProtonNeutron(proton1, neutron1).id;

        proton2 = (uint8_t)node.node().child("proton2").text().as_int(0);
        neutron2 = (uint8_t)node.node().child("neutron2").text().as_int(0);
        elementId2 = pt.getByProtonNeutron(proton2, neutron2).id;
        
        TersoffParameter parameter(elementId1, elementId2,
              node.node().child("name").text().as_string("undef"),
              node.node().child("comment").text().as_string("undef"),
              node.node().child("mode").text().as_string("undef"),
              node.node().child("A").text().as_double(0),
              node.node().child("B").text().as_double(0),
              node.node().child("beta").text().as_double(0),
              node.node().child("C").text().as_double(0),
              node.node().child("D").text().as_double(0),
              node.node().child("D0").text().as_double(0),
              node.node().child("DD").text().as_double(0),
              node.node().child("delta").text().as_double(0),
              node.node().child("gamma").text().as_double(0),
              node.node().child("H").text().as_double(0),
              node.node().child("lambda").text().as_double(0),
              node.node().child("lambda3").text().as_double(0),
              node.node().child("mu").text().as_double(0),
              node.node().child("N").text().as_double(0),
              node.node().child("R0").text().as_double(0),
              node.node().child("Rcut").text().as_double(0),
              node.node().child("S").text().as_double(0));
        
        potential.add(parameter);
    }

    CLOG(DEBUG, logName_) << counter <<" Tersoff Parameters generated from XML";
    CLOG(TRACE, logName_) << "Tersoff Parameters successfully generated";    
    
}

//------------------------------------------------------------------------------

//! Creates a Material object and fills it with data from the DOM structure.
//! The periodic table is required here to allow the user to specify the
//! elements by their symbol (e.g. Ga, In, As) and not by the internally used
//! ID.
//! \param collection Material collection that shall be filled from DOM.
//! \throws PeriodicTableException
void XmlHandler::get(MaterialCollection &collection)
    const
{

    int counter=0;
    const PeriodicTable &pt = PeriodicTable::getInstance();

    CLOG(TRACE, logName_) << "generating Materials from XML";    

    pugi::xpath_node_set::const_iterator it;
    pugi::xpath_node_set materialNodes =
	xmlDoc_.select_nodes("/material-collection/material");    

    for (it = materialNodes.begin(); it != materialNodes.end(); ++it, counter++) {
        MaterialComponents cations, anions;
        pugi::xml_node materialNode=it->node();

        for (pugi::xml_node cation = materialNode.child("cation");
            cation; cation=cation.next_sibling("cation")) {
            
            cations.push_back(std::make_tuple (
             pt.getByProtonNeutron(
                  (uint8_t)cation.child("proton").text().as_int(0),
                  (uint8_t)cation.child("neutron").text().as_int(0)
                          ).id,
             cation.child("share").text().as_double(1)
                              ));
        }

        for (pugi::xml_node anion = materialNode.child("anion");
            anion; anion=anion.next_sibling("anion")) {
            
            anions.push_back (std::make_tuple (
             pt.getByProtonNeutron(
                  (uint8_t)anion.child("proton").text().as_int(0),
                  (uint8_t)anion.child("neutron").text().as_int(0)
                          ).id,
              anion.child("share").text().as_double(1)
                              ) );
        }

        Material material(materialNode.child("name").text().as_string("undef"),
                       cations, anions,
              materialNode.child("lattice-constant").text().as_double(0),
              materialNode.child("c11").text().as_double(0),
              materialNode.child("c12").text().as_double(0),
              materialNode.child("c44").text().as_double(0)
                        );
        
        collection.add(material);
    }

    CLOG(DEBUG, logName_) << counter <<" Materials generated from XML";
    CLOG(TRACE, logName_) << "Materials successfully generated";
   
}

//------------------------------------------------------------------------------

//! The periodic table is required here to display the
//! elements by their symbol (e.g. Ga, In, As) and not by the internally used
//! ID.
//! \param collection Material collection storing the single materials.
void XmlHandler::set(const MaterialCollection &collection)
{

    CLOG(TRACE, logName_) << "writing materials to XML";

    const PeriodicTable &pt = PeriodicTable::getInstance();    
    int counter=0;
    pugi::xml_node collectionNode, materialNode, workNode;
    MaterialComponents component;
    
    // remove already existing nodes with the same name
    while(xmlDoc_.remove_child("material-collection")) ;

    collectionNode=xmlDoc_.append_child("material-collection");

    for (auto it = collection.get().begin(); it != collection.get().end();
        ++it, counter++) {

        std::vector<double> elasticConstants = it->getElasticConstants();
        
        materialNode=collectionNode.append_child("material");

        materialNode.append_child("name").text() = it->getName().c_str();
        materialNode.append_child("lattice-constant").text() =
            it->getLatticeConstant();
        
        materialNode.append_child("c11").text() = elasticConstants[0];
        materialNode.append_child("c12").text() = elasticConstants[1];
        materialNode.append_child("c44").text() = elasticConstants[2];

        component = it->getCations();
        
        double lastShare=0;

        for (auto it2 = component.begin();
             it2 != component.end(); ++it2) {
            auto element = pt.getById(std::get<0>(*it2));
            workNode=materialNode.append_child("cation");
            workNode.append_child("element").text() = element.symbol.c_str();
            workNode.append_child("proton").text() = element.protonCount;
            workNode.append_child("neutron").text() = element.neutronCount;
            workNode.append_child("share").text() = std::get<1>(*it2)-lastShare;
            lastShare = std::get<1>(*it2);
        }

        component = it->getAnions();

        lastShare=0;

        for (auto it2 = component.begin(); it2 != component.end(); ++it2) {
            auto element = pt.getById(std::get<0>(*it2));
            workNode=materialNode.append_child("anion");
            workNode.append_child("element").text() = element.symbol.c_str();
            workNode.append_child("proton").text() = element.protonCount;
            workNode.append_child("neutron").text() = element.neutronCount;
            workNode.append_child("share").text() = std::get<1>(*it2)-lastShare;
            lastShare = std::get<1>(*it2);
        }
        
    }

    CLOG(DEBUG, logName_) << counter <<" Materials written to XML";
    CLOG(TRACE, logName_) << "Material succesfully written";    
}

//------------------------------------------------------------------------------

//! Transfers content of the lattice into the DOM structure. The ID, (i,j,k) and
//! (x,y,z) values are exported using the PCDATA format.
//! \param lattice Atomic lattice.
void XmlHandler::set(const Lattice &lattice)
{

    CLOG(TRACE, logName_) << "writing Lattice to XML";

    const PeriodicTable &pt = PeriodicTable::getInstance();    
    pugi::xml_node latticeNode, workNode;
    std::vector<Atom *> atoms;
    Vector3D<indexType> size = lattice.getSize();
    int counter=0;
    
    // remove already existing nodes with the same name
    while(xmlDoc_.remove_child("lattice"));
    
    latticeNode = xmlDoc_.append_child("lattice");

    latticeNode.append_child("temperature").text() = lattice.getTemperature();
    workNode=latticeNode.append_child("size");
    workNode.append_child("x").text() = size[0];
    workNode.append_child("y").text() = size[1];
    workNode.append_child("z").text() = size[2];

    for (auto atom: lattice.getAtomList()) {
        Vector3D<spaceType> position = atom->getPosition();
        Vector3D<indexType> index = atom->getIndex();
        workNode=latticeNode.append_child("atom");
            
        /*
        // export as node attributes
        atomNode.append_attribute("ID") = atom->getElementId();
        atomNode.append_attribute("i") = i;
        atomNode.append_attribute("j") = j;
        atomNode.append_attribute("k") = k;
        atomNode.append_attribute("x") = pos[0];
        atomNode.append_attribute("y") = pos[1];
        atomNode.append_attribute("z") = pos[2];
        */

        const Element &element = pt.getById(atom->getElementId());
        // export as PCDATA nodes
        workNode.append_child("element").text() = element.symbol.c_str();
        workNode.append_child("proton").text() = element.protonCount;
        workNode.append_child("neutron").text() = element.neutronCount;
        workNode.append_child("material").text() =
            atom->getMaterial()->getName().c_str();
        workNode.append_child("modified").text() = atom->wasModified();
        workNode.append_child("modification-index").text() = atom->getModificationOrder();
        workNode.append_child("state").text() = atom->getState().to_string().c_str();

        workNode.append_child("i").text() = index[0];
        workNode.append_child("j").text() = index[1];
        workNode.append_child("k").text() = index[2];
        workNode.append_child("x").text() = position[0];
        workNode.append_child("y").text() = position[1];
        workNode.append_child("z").text() = position[2];

        counter++;
    }

    CLOG(DEBUG, logName_) << counter <<" Atoms written to XML";
    CLOG(TRACE, logName_) << "Lattice successfully written";    
}

//------------------------------------------------------------------------------

//! Creates a lattice of size specified in the DOM structure and fills the
//! lattice with the stored atoms defined also in the DOM structure.
//! \param lattice Atomic lattice.
//! \param materials Materials used in the simulation box.
void XmlHandler::get(Lattice &lattice, const MaterialCollection &materials) const
{

    unsigned int counter=0;
    Vector3D<spaceType> position;
    Vector3D<indexType> index;
    const PeriodicTable &pt = PeriodicTable::getInstance();
    
    CLOG(TRACE, logName_) << "generating Lattice from XML";    

    lattice.setTemperature(
	    xmlDoc_.child("lattice").child("temperature").text().as_double(0));
    pugi::xml_node size=xmlDoc_.child("lattice").child("size");
    lattice.generate((indexType)size.child("x").text().as_int(0),
		     (indexType)size.child("y").text().as_int(0),
		     (indexType)size.child("z").text().as_int(0));
    
    pugi::xpath_node_set atomNodes =
	xmlDoc_.select_nodes("/lattice/atom");
    pugi::xpath_node_set::const_iterator it;

    for (it = atomNodes.begin(); it != atomNodes.end(); ++it, counter++) {
        pugi::xpath_node atom=*it;

        position[0] = atom.node().child("x").text().as_double(0);
        position[1] = atom.node().child("y").text().as_double(0);
        position[2] = atom.node().child("z").text().as_double(0);

        index[0] = (indexType)atom.node().child("i").text().as_int(0);
        index[1] = (indexType)atom.node().child("j").text().as_int(0);
        index[2] = (indexType)atom.node().child("k").text().as_int(0);

        uint8_t proton = (uint8_t)atom.node().child("proton").text().as_int(0);
        uint8_t neutron =(uint8_t)atom.node().child("neutron").text().as_int(0);
        std::string name = atom.node().child("material").text().as_string("");
        bool modified = atom.node().child("modified").text().as_bool(false);
        uint32_t modificationIndex = atom.node().child("modification-index").text().as_int(0);
        std::string state = atom.node().child("state").text().as_string("");

        lattice(index) =
            new Atom(pt.getByProtonNeutron(proton,neutron).id, position,
                 index, &(materials.getByName(name)), modified, false, modificationIndex, std::bitset<AtomState::StateCount>(state));

        CLOG(TRACE, logName_) << name << " Material Name " <<
            materials.getByName(name).getName();
    }

    CLOG(DEBUG, logName_) << counter <<" Atoms generated from XML";
    CLOG(TRACE, logName_) << "Lattice successfully generated";

}

//------------------------------------------------------------------------------

//! Transfers content of the complete simulation box into the DOM
//! structure. First global properties of the simulation box are imported. Later
//! the function to write the Lattice is called.
//! \param simbox SimulationBox that shall be written.
void XmlHandler::set(const SimulationBox &simbox)
{

    CLOG(TRACE, logName_) << "writing Simulation Box to XML";    

    pugi::xml_node simboxNode;

    // remove already existing nodes with the same name
    while(xmlDoc_.remove_child("simbox"));
    
    simboxNode = xmlDoc_.append_child("simbox");
    simboxNode.append_child("description").text() = simbox.getDescription().c_str();
    simboxNode.append_child("outOfPlaneDimension").text() =
        simbox.getOutOfPlaneDimension();
    simboxNode.append_child("inPlaneLatticeConstant").text() =
        simbox.getInPlaneLatticeConstant();
    simboxNode.append_child("x").text() = simbox.getSize()[0];
    simboxNode.append_child("y").text() = simbox.getSize()[1];
    simboxNode.append_child("z").text() = simbox.getSize()[2];

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
void XmlHandler::get(std::shared_ptr<SimulationBox> &simbox) const
{
    Vector3D<spaceType> size;

    CLOG(TRACE, logName_) << "generating Simulation Box from XML";
    
    pugi::xml_node node=xmlDoc_.child("simbox");
    
    simbox = std::make_shared<SimulationBox>(
	      (uint8_t)node.child("outOfPlaneDimension").text().as_int(0),
	      Vector3D<spaceType>(node.child("x").text().as_double(0),
                              node.child("y").text().as_double(0),
                              node.child("z").text().as_double(0)),
	      node.child("inPlaneLatticeConstant").text().as_double(-1)
                                             );

    simbox->setDescription(node.child("description").text().as_string(""));

    this->get(simbox->getMaterials());
    this->get(simbox->getLattice(), simbox->getMaterials());
    // Force creation of neighbor list.
    simbox->generateNeighbors(true);
    
    CLOG(TRACE, logName_) << "Simulation Box successfully generated";
   
}

//------------------------------------------------------------------------------

//! Transfers content of CompositionInfo to XML
//! \param composition Composition information that shall be transferred
void XmlHandler::set(const CompositionInfo &composition)
{

    CLOG(TRACE, logName_) << "writing CompositionInfo to XML";
    
    pugi::xml_node compositionNode, layerNode, entryNode;
    int counter=0;
    const PeriodicTable &pt = PeriodicTable::getInstance();
    
    compositionNode = xmlDoc_.append_child("composition");

    for (auto i=0; i<composition.getLayerCount(); i++) {
        LayerCompositionInfo layer = composition.getLayer(i);
        layerNode=compositionNode.append_child("layer");
        layerNode.append_child("index").text() = i;
        
        for (auto entry : layer) {

            entryNode=layerNode.append_child("element");
            
            auto key = std::get<0>(entry);
        
            entryNode.append_child("material").text() = std::get<0>(key).c_str();
            entryNode.append_child("id").text() = std::get<1>(key);
            entryNode.append_child("name").text() = pt.getById(std::get<1>(key)).symbol.c_str();
            entryNode.append_child("count").text() = std::get<1>(entry);
        }

        counter++;
    }

    CLOG(DEBUG, logName_) << counter <<" atomic layers written to XML";
    CLOG(TRACE, logName_) << "CompositionInfo successfully written";    
}

//------------------------------------------------------------------------------

//! Creates a simbox object from the stored data and returns a pointer. This
//! method was chosen as some parameters of the SimulationBox have to be set
//! in the constructor since they are not supposed to change. Therefore it was
//! neglected to create a member function to set these values.
//! \param simbox SimulationBox that shall be filled.
void XmlHandler::get(std::vector<CompositionInfo> &compositions) const
{

    CLOG(TRACE, logName_) << "generating CompositionInfos from XML";

	for (pugi::xml_node compNode = xmlDoc_.child("composition");
	    compNode; compNode=compNode.next_sibling("composition")) {

        CompositionInfo comp;
        
        for (pugi::xml_node layerNode = compNode.child("layer");
            layerNode; layerNode=layerNode.next_sibling("layer")) {

            LayerCompositionInfo layer;

            for (pugi::xml_node elementNode = layerNode.child("element");
                elementNode; elementNode=elementNode.next_sibling("element")) {

                auto key = std::make_tuple(elementNode.child("material").text().as_string("undef"),
                                           (elementType)elementNode.child("id").text().as_int(0) );

                layer[key] = elementNode.child("count").text().as_int(0);

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
int XmlHandler::readEntries(Journal<int> &journal, pugi::xml_node entry) const
{
    int counter=0;
    
    for (; entry; entry = entry.next_sibling("entry"), ++counter) {
        journal.add(entry.text().as_int(0));
    }

    return counter;
}

//------------------------------------------------------------------------------

//! \param journal Journal the entries shall be written to.
//! \param entry First entry node of XML <journal> environment
//! \return Number of read entries.
int XmlHandler::readEntries(Journal<double> &journal, pugi::xml_node entry)const
{
    int counter=0;
    
    for (; entry; entry = entry.next_sibling("entry"), ++counter) {
        journal.add(entry.text().as_double(0));
    }

    return counter;
}

//------------------------------------------------------------------------------

//! \param journal Journal the entries shall be written to.
//! \param entry First entry node of XML <journal> environment
//! \return Number of read entries.
int XmlHandler::readEntries(Journal<std::string> &journal,
			    pugi::xml_node entry) const
{
    int counter=0;
    
    for (; entry; entry = entry.next_sibling("entry"), ++counter) {
        journal.add(entry.text().as_string(""));
    }

    return counter;
}

//------------------------------------------------------------------------------

//! \param journal Journal the entries shall be written to.
//! \param entry First entry node of XML <journal> environment
//! \return Number of read entries.
int XmlHandler::readEntries(Journal<UDTuple> &journal,
                pugi::xml_node entry) const
{
    int counter = 0;

    for (; entry; entry = entry.next_sibling("entry"), ++counter) {
        journal.add(std::make_tuple(entry.child("U").text().as_ullong(), entry.child("D").text().as_double(0)));
    }

    return counter;
}

//------------------------------------------------------------------------------

//! Specific write function for std::string. Is required because Pugi XML can
//! not handle std::string but only char*.
//! \param journal Journal containing entries that shall be written.
//! \param journalNode XML Node which stores Journal data.
//! \return Number of written entries.
int XmlHandler::writeEntries(const Journal<std::string> &journal,
			      pugi::xml_node &journalNode)
{
    int counter = 0;

    for (auto entry: journal.getEntries()) {
        journalNode.append_child("entry").text() = entry.c_str();
        counter ++;
    }

    return counter;
}

//------------------------------------------------------------------------------

int XmlHandler::writeEntries(const Journal<UDTuple> &journal,
                pugi::xml_node &journalNode)
{
    int counter = 0;
    pugi::xml_node tmpNode;

    for (auto entry: journal.getEntries()) {
        tmpNode = journalNode.append_child("entry");
        tmpNode.append_child("U").text() = std::get<0>(entry);
        tmpNode.append_child("D").text() = std::get<1>(entry);
        counter++;
    }

    return counter;
}

//------------------------------------------------------------------------------
    
//! Read Configuration file from DOM structure.
//! \param config Configuration object.
//! \todo read cations and anions also by proton/neutron count
void XmlHandler::get(Configuration &config) const
{
    
    CLOG(TRACE, logName_) << "generating Configuration from XML";

    pugi::xml_node workNode, configNode=xmlDoc_.child("config");
    pugi::xpath_node_set nodeSet;
    pugi::xpath_node_set::const_iterator it;

    nodeSet = xmlDoc_.select_nodes("/config/material/cation");    

    for (it = nodeSet.begin(); it != nodeSet.end(); ++it) {
        pugi::xml_node cationNode=(*it).node();

        config.cations.push_back( std::make_tuple(
            cationNode.child("element").text().as_string("undef"),
            cationNode.child("share").text().as_double(0) )
                   );
    }
    
    nodeSet = xmlDoc_.select_nodes("/config/material/anion");    

    for (it = nodeSet.begin(); it != nodeSet.end(); ++it) {
        pugi::xml_node anionNode=(*it).node();

        config.anions.push_back( std::make_tuple(
            anionNode.child("element").text().as_string("undef"),
            anionNode.child("share").text().as_double(0) )
                   );
    }

    workNode = configNode.child("material");
    config.materialName = workNode.child("name").text().as_string("undef");
    config.c11= workNode.child("c11").text().as_double(0);
    config.c12= workNode.child("c12").text().as_double(0);
    config.c44= workNode.child("c44").text().as_double(0);
    config.latticeConstant=
	workNode.child("lattice-constant").text().as_double(0);
	config.latticeType = workNode.child("lattice-type").text().as_string("");
    config.materialFile =configNode.child("material-file").text().as_string("");
    config.materialSearchName =
	configNode.child("material-name").text().as_string("");
    
    workNode = configNode.child("size");
    config.size = {(indexType)workNode.child("dimension1").text().as_int(0),
		 (indexType)workNode.child("dimension2").text().as_int(0),
		 (indexType)workNode.child("dimension3").text().as_int(0)};
    config.growthDimension =
	(uint8_t)configNode.child("growth-dimension").text().as_int(0);
    config.inputFileName = configNode.child("input").text().as_string("");
    config.outputFileName = configNode.child("output").text().as_string("");
    config.periodicTableFileName =
	configNode.child("periodic-table").text().as_string("");
    config.logFileName = configNode.child("log").text().as_string("");
    config.tersoffFileName = configNode.child("tersoff").text().as_string("");
    config.xyzFileName = configNode.child("xyz").text().as_string("");
    config.outputPreamble =
	configNode.child("output-preamble").text().as_string("");
    config.journalPreamble =
	configNode.child("journal-preamble").text().as_string("");
    config.neighborRadius =
	configNode.child("neighbor-radius").text().as_double(0);
    config.neighborLayers =
	(indexType)configNode.child("neighbor-layers").text().as_int(0);
    config.latticeTemperature = 
	(indexType)configNode.child("temperature").text().as_double(0);
	
    config.startIndex = configNode.child("start-index").text().as_int(-1);
    config.stopIndex = configNode.child("stop-index").text().as_int(-1);
    config.mmcProbability =
	configNode.child("mmc-probability").text().as_double(-1);
	config.mmcRunCount = configNode.child("mmc-run-count").text().as_int(1);
    config.minDisplacement =
	configNode.child("min-displacement").text().as_double(0);
    config.maxDisplacement =
	configNode.child("max-displacement").text().as_double(-1);
    config.scalingProbability =
	configNode.child("scaling-probability").text().as_double(-1);
    config.minScaling =
	configNode.child("min-scaling").text().as_double(0);
    config.maxScaling =
	configNode.child("max-scaling").text().as_double(-1);
    config.runCount = configNode.child("runs").text().as_int(-1);
    config.checkCount = configNode.child("check-runs").text().as_int(1);
    config.energyDropFactor=configNode.child("energy-drop").text().as_double(2);
    config.reductionFactor = configNode.child("reduction").text().as_double(2);
    config.minEnergy = configNode.child("min-energy").text().as_double(0);
    config.anionPassivationProbability = configNode.child("anion-passivation-probability").text().as_double(1);
    config.maxThreadCount = (uint16_t)configNode.child("max-thread-count").text().as_int(1);

    config.modificationIndexMax = configNode.child("max-modification-index").text().as_double(0);
    workNode = configNode.child("camera-focal-point");
    config.cameraFocalPoint = {workNode.child("dimension1").text().as_double(0),
		 workNode.child("dimension2").text().as_double(0),
		 workNode.child("dimension3").text().as_double(0)};
    workNode = configNode.child("camera-direction");
    config.cameraDirection = {workNode.child("dimension1").text().as_double(0),
		 workNode.child("dimension2").text().as_double(0),
		 workNode.child("dimension3").text().as_double(0)};
    config.cameraZoom = configNode.child("camera-zoom").text().as_double(1);
    config.animationStep = (uint32_t)configNode.child("animation-step").text().as_int(1);
    config.animationFilePrefix = configNode.child("animation-prefix").text().as_string("");
    config.screenshotFileName = configNode.child("screenshot").text().as_string("");
    config.rotationCount = (uint16_t)configNode.child("rotation_count").text().as_int(0);
    
    config.interfaceAtoms =
	configNode.child("interface-atoms").text().as_int(0);
    config.interfacePositions =
	configNode.child("interface-positions").text().as_int(0);

    config.refLatticeConstant =
	configNode.child("ref-lattice-const").text().as_double(0);

    nodeSet = xmlDoc_.select_nodes("/config/switch");    

    for (it = nodeSet.begin(); it != nodeSet.end(); ++it) {
        std::string text((*it).node().text().as_string(""));
        
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

        if (text == "modified-only") {
            config.modifiedOnly=true;
            continue;
        }

        if (text == "modified-negative") {
            config.modifiedNegative=true;
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

//##############################################################################

XmlException::XmlException(const std::string &message,
			   const std::string &fileName,
			   const XmlException::Id id):
    message_(message), fileName_(fileName), id_(id) {}

//------------------------------------------------------------------------------

//! Returns a short description of the exception. This is the implementation of
//! the virtual function inherited from std::exception
//! \remark This function does not throw exceptions!
const char * XmlException::what () const throw ()
{
    static std::string text;
    text = "XML Exception -- File '" + fileName_ + "': " +message_;
    return text.c_str();
}

//------------------------------------------------------------------------------

//! Returns a more detailled description of the exception, making it possible to
//! use a single exception class for multiple errors.
//! \remark This function does not throw exceptions!
const std::string XmlException::getMessage () const noexcept
{
    return message_;
}

//------------------------------------------------------------------------------

//! Returns the ID of the exception source, making possible to easily
//! distinguish different errors.
//! \remark This function does not throw exceptions!
const std::string XmlException::getFileName () const noexcept
{
    return fileName_;
}

//------------------------------------------------------------------------------

//! Returns the ID of the exception source, making possible to easily
//! distinguish different errors.
//! \remark This function does not throw exceptions!
const XmlException::Id XmlException::getId () const noexcept
{
    return id_;
}

//##############################################################################

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
