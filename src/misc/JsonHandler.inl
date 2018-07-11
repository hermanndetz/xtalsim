/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __JSON_HANDLER_INL__
#define __JSON_HANDLER_INL__

//------------------------------------------------------------------------------

//! Arbitrary types can be stored in the Journal class. Therefore at first the
//! values present in every object are written. The single entries are then
//! processed in function writeEntries which is implemented for each data type
//! separately.
//! \param journal Journal that shall be transferred to XML.
//! \sa writeEntries
template <class T>
void JsonHandler::set(const Journal<T> &journal)
{

    int counter=0;

    Json::Value journalNode(Json::objectValue);
    Json::Value entryNode(Json::arrayValue);

    CLOG(TRACE, logName_) << "writing Journal to JSON";
    journalNode["name"] = journal.getName();
    journalNode["description"] = journal.getDescription();

    counter = writeEntries(journal, entryNode);
    journalNode["entries"] = entryNode;

    if (jsonDoc_["journals"].type() == Json::nullValue)
	jsonDoc_["journals"] = Json::Value(Json::arrayValue);
    
    jsonDoc_["journals"].append(journalNode);

    CLOG(DEBUG, logName_) << counter << " Entries written to JSON";
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
void JsonHandler::get(Journal<T> &journal, const std::string name) const
{
    int counter=0;
    bool journalFound=false;

    CLOG(TRACE, logName_) << "generating Journal with name " << name <<
	" from JSON";
    
    Json::Value startNode = jsonDoc_["journals"];

    for (int i=0; i< startNode.size(); i++){
	Json::Value journalNode=startNode[i];

	if (journalNode.get("name","").asString() != name)
	    continue;

	journalFound = true;

	journal.setName(journalNode.get("name","").asString());
	journal.setDescription(journalNode.get("description","").asString());

	counter = readEntries(journal, journalNode["entries"]);

	break;
	
    }

    if (journalFound == false){
	CLOG(WARNING, logName_) << "Journal with name " << name <<" not found!";
    }

    CLOG(DEBUG, logName_) << counter << " Journal entries generated from JSON";
    CLOG(TRACE, logName_) << "Journal successfully generated";
}

//------------------------------------------------------------------------------

//! \param journal Journal containing entries that shall be written.
//! \param journalNode XML Node which stores Journal data.
//! \return Number of written entries.
template <class T>
int JsonHandler::writeEntries(const Journal<T> &journal,
			      Json::Value &entryNode)
{
    int counter = 0;
    
    for (auto entry: journal.getEntries()){
	entryNode[entryNode.size()] = entry;
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
