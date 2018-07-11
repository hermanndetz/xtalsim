/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __XML_HANDLER_INL__
#define __XML_HANDLER_INL__

//------------------------------------------------------------------------------

//! Arbitrary types can be stored in the Journal class. Therefore at first the
//! values present in every object are written. The single entries are then
//! processed in function writeEntries which is implemented for each data type
//! separately.
//! \param journal Journal that shall be transferred to XML.
//! \sa writeEntries
template <class T>
void XmlHandler::set(const Journal<T> &journal)
{

    int counter=0;
    pugi::xml_node journalNode=xmlDoc_.append_child("journal");

    CLOG(TRACE, logName_) << "writing Journal to XML";
    journalNode.append_child("name").text() = journal.getName().c_str();
    journalNode.append_child("description").text() =
	journal.getDescription().c_str();

    counter = writeEntries(journal, journalNode);

    CLOG(DEBUG, logName_) << counter << " Entries written to XML";
    CLOG(TRACE, logName_) << "Material succesfully written";
}

//------------------------------------------------------------------------------

//! The name of the Journal entry to read has to be specified.
//! Arbitrary types can be stored in the Journal class. Therefore at first the
//! values present in every object are read. The single entries are then
//! processed in function readEntries which is implemented for each data type
//! separately.
//! \warning The first occurance of "name" is read.
//! \param journal Journal that shall be filled.
//! \param name Name of journal that shall be read.
//! \sa readEntries
template <class T>
void XmlHandler::get(Journal<T> &journal, const std::string name) const
{
    int counter=0;
    bool journalFound=false;

    CLOG(TRACE, logName_) << "generating Journal with name " << name <<
	" from XML";

    journal.clear();
    
    pugi::xpath_node_set::const_iterator it;
    pugi::xpath_node_set journalNodes =
	xmlDoc_.select_nodes("/journal");

    for(it = journalNodes.begin(); it != journalNodes.end(); ++it){
	pugi::xml_node journalNode=it->node();

	if (journalNode.child("name").text().as_string("") != name)
	    continue;

	journalFound = true;

	journal.setName(journalNode.child("name").text().as_string(""));
	journal.setDescription(
		      journalNode.child("description").text().as_string(""));

	counter = readEntries(journal, journalNode.child("entry"));

	break;
	
    }

    if (journalFound == false){
	CLOG(WARNING, logName_) << "Journal with name " << name <<" not found!";
	throw XmlException("Journal with name " + name + " not found!", "",
			   XmlException::Id::Journal);
    }

    CLOG(DEBUG, logName_) << counter << " Journal entries generated from XML";
    CLOG(TRACE, logName_) << "Journal successfully generated";
}

//------------------------------------------------------------------------------

//! \param journal Journal containing entries that shall be written.
//! \param journalNode XML Node which stores Journal data.
//! \return Number of written entries.
template <class T>
int XmlHandler::writeEntries(const Journal<T> &journal,
			      pugi::xml_node &journalNode)
{
    int counter = 0;
    
    for (auto entry: journal.getEntries()){
	journalNode.append_child("entry").text() = entry;
	counter ++;
    }

    return counter;
}

//------------------------------------------------------------------------------

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
