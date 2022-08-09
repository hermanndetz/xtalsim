/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "CsvHandler.h"

//! Constructor
//! \param logName Optional parameter defining the logger's name.
CsvHandler::CsvHandler (const char *logName): logName_(logName), colCount_(0),
					      rowCount_(0)
{
    data_.resize(0);
}

//------------------------------------------------------------------------------

//! Constructor
//! \param fileName Path to the file, which should be loaded.
//! \param logName Optional parameter defining the logger's name.
CsvHandler::CsvHandler (const std::string &fileName, const char *logName):
    logName_(logName), colCount_(0), rowCount_(0)
{
    loadFile(fileName);
}

//------------------------------------------------------------------------------

//! Destructor
CsvHandler::~CsvHandler () {

}

//------------------------------------------------------------------------------

//! \param str The string that has to be analyzed.
//! \return True, if string consists of whitespace only or if it is empty.
bool CsvHandler::isWhitespace (const std::string str) const {

    // quick check first
    if (str.empty())
        return true;

    return std::all_of(str.begin(), str.end(), isspace);

}

//------------------------------------------------------------------------------

//! Loads data from a CSV or TSV file into memory.
//! \warning This function is not fully typesafe yet. The data types of the
//! first valid row determine the data type for all remaining ones. Conflicting
//! data are disarded and set to 0.
//! \param fileName Path to the file, which should be loaded.
//!
//! Example:
//!
//! # col1  col2
//! 0       1
//! 1       1.2
//!
//! The value in the second row of col2 is discarded since the
//! column is assumed to be an integer column due to the value of
//! 1 in the first row.
//!
//! \todo Find a way to be more flexible with respect to data types.
//! \todo {Implement a more sophisticated way to read the column names 
//!       to allow strings with blanks.}
void CsvHandler::loadFile (const std::string &fileName) {

    data_.resize(0);

    std::ifstream ifHandle(fileName, std::ios::binary);

    if (ifHandle.is_open()) {
        std::string line;
        std::string prevLine;

        // ignore header lines starting with #
        do {
            prevLine = line;
            std::getline(ifHandle, line);
        } while (line[0] == '#');

        // process last header line
        // treat whatever is there as column names
        char delimiter;
        std::stringstream lineStream(prevLine);
        std::string field;

        if (prevLine != "") {
            if (prevLine.find('\t') != std::string::npos) {
                delimiter = '\t';
            } else if (prevLine.find(',') != std::string::npos) {
                delimiter = ',';
            } else {
                delimiter = ' ';
            }

            while (std::getline(lineStream, field, delimiter)) {
                if (field == "#")
                    continue;

                if (isWhitespace(field) == false) {
                    colNames_.push_back(field);
                }
            }
        }

        // process first data line
        // determine delimiting character
        if (line.find(',') != std::string::npos)
            delimiter = ',';

        if (line.find(' ') != std::string::npos)
            delimiter = ' ';

        if (line.find('\t') != std::string::npos)
            delimiter = '\t';

        lineStream.str(line);
        dataField tmpData;

        unsigned long rowCount_ = 0;
        unsigned long curCol = 0;

        while (std::getline(lineStream, field, delimiter)) {
            if (isWhitespace(field) == false) {
                data_.push_back(std::vector<dataField>());
                colCount_++;

                if ((field.find('.') != std::string::npos) ||
                    (field.find('e') != std::string::npos)) {
                    colTypes_.push_back(ColumnDataType::Float);

                    tmpData.d = std::stold(field);
                } else {
                    colTypes_.push_back(ColumnDataType::Integer);
                    tmpData.n = std::stoll(field);
                }

                data_[curCol].push_back(tmpData);
                curCol++;
            }

            rowCount_++;
        }

        // now process all other lines
        while (std::getline(ifHandle, line)) {
            lineStream.str(line);
            curCol = 0;

            while (std::getline(lineStream, field, delimiter)) {
                if (isWhitespace(field) == false) {
                    switch (colTypes_[curCol]) {
                    case ColumnDataType::Integer:
                        if ((field.find('.') != std::string::npos) ||
                            (field.find('e') != std::string::npos)) {
                            tmpData.n = 0;
                        } else {
                            tmpData.n = std::stoll(field);
                        }
                        break;
                    case ColumnDataType::Float:
                        tmpData.d = std::stold(field);
                        break;
                    }

                    data_[curCol].push_back(tmpData);
                    curCol++;
                }

                rowCount_++;
            }
        }

        // fill up empty column names
        // the default value contains the column index as
        // e.g. vtkTable requires unique column names.
        for (unsigned long i = colNames_.size(); i < colCount_; i++) {
            colNames_.push_back("col" + std::to_string(i));
        }
    }
}

//------------------------------------------------------------------------------

//! Returns the number of columns in the CSV file loaded.
//! \return Number of columns.
unsigned long CsvHandler::getColCount (void) const {
    return colCount_;
}

//------------------------------------------------------------------------------

//! Returns the number or rows in the CSV file loaded.
//! \return Number of rows.
unsigned long CsvHandler::getRowCount (void) const {
    return rowCount_;
}

//------------------------------------------------------------------------------

//! Clear data representation in memory.
//! Does not affect data in file.
void CsvHandler::clear () {
    for (auto col : data_) {
        col.clear();
    }

    data_.clear();
}

//------------------------------------------------------------------------------

//! Add column of a given type.
uint32_t CsvHandler::addCol (const ColumnDataType dt, const std::string & colName, const dataField defaultValue) {
    data_.push_back(std::vector<dataField>());
    colTypes_.push_back(dt);
    colNames_.push_back(colName);

    if (rowCount_ > 0) {
        // \todo JM: use "auto item&" here, because this way nothing is stored I predict!
        for (auto item: data_[colCount_]) {
            switch (dt) {
                case ColumnDataType::Integer:
                    item.n = defaultValue.n;
                    break;
                case ColumnDataType::Float:
                    item.d = defaultValue.d;
                    break;
                default:
                    throw CsvException("Unknown column data type!", "", CsvException::Id::UnknownColumnType);
                    break;
            }
        }
    }

    return colCount_ ++;
}

//------------------------------------------------------------------------------

//! Add row.
void CsvHandler::addRow (std::vector<dataField> df) {
    uint32_t i = 0;

    for (auto col: data_) {
        try {
            col.push_back(df[i]);
            i++;
        } catch (std::exception &e) {
            throw CsvException("Wrong number of colunns!", "", CsvException::Id::ColumnCountMismatch);
        }
    }

    rowCount_++;
}

//------------------------------------------------------------------------------

//! Add row.
void CsvHandler::addRow (UDTuple df) {
    dataField tmp{};

    try {
        tmp.n = std::get<0>(df);
        data_[0].push_back(tmp);
        tmp.d = std::get<1>(df);
        data_[1].push_back(tmp);
    } catch (std::exception &e) {
        throw CsvException("Wrong number or type of columns!", "", CsvException::Id::ColumnCountMismatch);
    }
}

//------------------------------------------------------------------------------

//! Save data to file.
void CsvHandler::save (const std::string &fileName, const std::string separator, const std::string &header) {
    save(fileName, -1, separator, header);
}

//------------------------------------------------------------------------------

//! Save data fo file, replace element id with element name.
//! \remark This function will be deprecated in future.
void CsvHandler::save (const std::string &fileName, int32_t substCol, const std::string separator, const std::string &header) {
    save(fileName, std::vector<int32_t>{substCol}, separator, header);
    /* std::ofstream out(fileName, std::ofstream::out); */

    /* const PeriodicTable &pt = PeriodicTable::getInstance(); */
    
    /* if (out.is_open()) { */
    /*     if (header != "") { */
    /*         out << "# " << header << std::endl; */
    /*     } */

    /*     for (uint32_t i = 0; i < data_[0].size(); i++) { */
    /*         for (uint32_t j = 0; j < data_.size(); j++) { */
    /*             // j == substCol is sufficient here */
    /*             // if substitution is disabled, substCol is -1, which will not */
    /*             // be reached here (within reasonable limits) */
    /*             // worst case scenario: saving 4,294,967,295 columns will */
    /*             // result in a substitution, if a periodic table reference */ 
    /*             // is provided. */
    /*             if ((j == (unsigned) substCol) ) { */
    /*                 out << pt.getById(data_[j][i].n).symbol; */
    /*             } else { */
    /*                 switch (colTypes_[j]) { */
    /*                 case ColumnDataType::Integer: */
    /*                     out << data_[j][i].n; */
    /*                     break; */
    /*                 case ColumnDataType::Float: */
    /*                     out << data_[j][i].d; */
    /*                     break; */
    /*                 default: */
    /*                     throw CsvException("Unknown column data type!", "", CsvException::Id::UnknownColumnType); */
    /*                     break; */
    /*                 } */
    /*             } */

    /*             out << separator; */
    /*         } */
            
    /*         out << std::endl; */
    /*     } */

    /*     out.close(); */
    /* } */
}

//------------------------------------------------------------------------------
//! Save data fo file, replace element id with element name.
void CsvHandler::save (const std::string &fileName, std::vector<int32_t> substCols,
           const std::string separator,
           const std::string &header) {
    std::ofstream out(fileName, std::ofstream::out);

    const PeriodicTable &pt = PeriodicTable::getInstance();
    
    if (out.is_open()) {
        if (header != "") {
            out << "# " << header << std::endl;
        }

        for (uint32_t i = 0; i < data_[0].size(); i++) {
            for (uint32_t j = 0; j < data_.size(); j++) {
                if (std::find(substCols.begin(), substCols.end(), j) != substCols.end()) {
                    out << pt.getById(data_[j][i].n).symbol;
                } else {
                    switch (colTypes_[j]) {
                    case ColumnDataType::Integer:
                        out << data_[j][i].n;
                        break;
                    case ColumnDataType::Float:
                        out << data_[j][i].d;
                        break;
                    default:
                        throw CsvException("Unknown column data type!", "", CsvException::Id::UnknownColumnType);
                        break;
                    }
                }

                out << separator;
            }
            
            out << std::endl;
        }

        out.close();
    }
}

//------------------------------------------------------------------------------

//! Print data to std out
void CsvHandler::dump () {
    if (data_.size() > 0) {
        for (uint32_t i = 0; i < data_[0].size(); i++) {
            for (uint32_t j = 0; j < data_.size(); j++) {
                switch (colTypes_[j]) {
                case ColumnDataType::Integer:
                    std::cout << data_[j][i].n;
                    break;
                case ColumnDataType::Float:
                    std::cout << data_[j][i].d;
                    break;
                default:
                    throw CsvException("Unknown column data type!", "", CsvException::Id::UnknownColumnType);
                    break;
                }

                std::cout << "\t";
            }
            
            std::cout << std::endl;
        }
    }
}

//------------------------------------------------------------------------------

//! Loads a journal of UDTuples (unsigned long long, double) into object.
//! \param src Journal of UDTuples to be loaded.
//! \todo return type
CsvHandler & CsvHandler::operator = (const Journal<UDTuple> &src) {
    std::vector<UDTuple> entries = src.getEntries();

    this->clear();
    this->addCol(ColumnDataType::Integer, "x");
    this->addCol(ColumnDataType::Float, "y");

    for (auto entry: entries) {
        this->addRow(entry);
    }

    return *this;
}

//------------------------------------------------------------------------------

//! Adds a journal of UDTuples (unsigned long long, double) to object.
//! \param src Journal of UDTuples to be added.
//! \todo return type
CsvHandler & CsvHandler::operator += (const Journal<UDTuple> &src) {
    std::vector<UDTuple> entries = src.getEntries();

    uint32_t ix{};
    uint32_t iy{};

    ix = this->addCol(ColumnDataType::Integer, "x");
    iy = this->addCol(ColumnDataType::Float, "y");


    if (rowCount_ == 0) {
        dataField tmp{};

        for (auto entry: entries) {
            tmp.n = std::get<0>(entry);
            data_[ix].push_back(tmp);
            tmp.d = std::get<1>(entry);
            data_[iy].push_back(tmp);
        }
    } else {
        uint32_t i{};
        for (auto entry: entries) {
            if (i >= rowCount_)
                break;

            data_[ix][i].n = std::get<0>(entry);
            data_[iy][i].d = std::get<1>(entry);
            i++;
        }
    }

    return *this;
}

//------------------------------------------------------------------------------

//! Loads CompositionInfo data into object.
//! \param src CompositionInfo structure to be loaded.
//! \todo does actually not return anything
CsvHandler & CsvHandler::operator = (const CompositionInfo &src) {
    this->clear();
    this->addCol(ColumnDataType::Integer, "layer");
    this->addCol(ColumnDataType::Integer, "element");
    this->addCol(ColumnDataType::Integer, "count");

    //! \todo implement iterator in CompositionInfo to make this loop simpler
    for (uint32_t i = 0; i < src.getLayerCount(); i++) {
        LayerCompositionInfo layer = src.getLayer(i);
        dataField tmp{};

        for (auto entry: layer) {
            tmp.n = i;
            data_[0].push_back(tmp);
            tmp.n = std::get<1>(std::get<0>(entry));
            data_[1].push_back(tmp);
            tmp.n = std::get<1>(entry);
            data_[2].push_back(tmp);
        }
    }

    return *this;
}

//------------------------------------------------------------------------------

//! Loads BondInfo data into object.
//! \param src BondInfo structure to be loaded.
//! \todo does actually not return anything
CsvHandler & CsvHandler::operator = (const BondInfo &src) {
    this->clear();
    this->addCol(ColumnDataType::Integer, "layer");
    this->addCol(ColumnDataType::Integer, "element1");
    this->addCol(ColumnDataType::Integer, "element2");
    this->addCol(ColumnDataType::Integer, "count");
    //
    //! \todo implement iterator in CompositionInfo to make this loop simpler
    for (uint32_t i = 0; i < src.getLayerCount(); i++) {
        LayerBondInfo layer = src.getLayer(i);
        dataField tmp{};

        for (auto entry: layer) {
            tmp.n = i;
            data_[0].push_back(tmp);
            tmp.n = std::get<0>(std::get<0>(entry));
            data_[1].push_back(tmp);
            tmp.n = std::get<1>(std::get<0>(entry));
            data_[2].push_back(tmp);
            tmp.n = std::get<1>(entry);
            data_[3].push_back(tmp);
        }
    }

    return *this;
}

//------------------------------------------------------------------------------

#ifdef __VTK__

vtkSmartPointer<vtkVariantArray> CsvHandler::getCol (unsigned long col) {
    vtkSmartPointer<vtkVariantArray> column = vtkSmartPointer<vtkVariantArray>::New();

    if (col < colCount_) {
        for (auto val : data_[col]) {
            switch (colTypes_[col]) {
            case ColumnDataType::Integer:
                column->InsertNextValue(vtkVariant(val.n));
                break;
            case ColumnDataType::Float:
                column->InsertNextValue(vtkVariant(val.d));
                break;
            }
        }
    }

    return column;
}

//------------------------------------------------------------------------------

vtkSmartPointer<vtkTable> CsvHandler::get() {
    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
    vtkSmartPointer<vtkLongLongArray> coli;
    vtkSmartPointer<vtkDoubleArray> cold;

    for (unsigned long i = 0; i < colCount_; i++) {
        vtkSmartPointer<vtkVariantArray> col = vtkSmartPointer<vtkVariantArray>::New();
        col = getCol(i);

        switch (colTypes_[i]) {
        case ColumnDataType::Integer:
            coli = vtkSmartPointer<vtkLongLongArray>::New();

            for (uint32_t j = 0; j < col->GetNumberOfValues(); j++) {
                coli->InsertNextValue(data_[i][j].n);
            }

            coli->SetName(colNames_[i].c_str());
            table->AddColumn(coli);
            break;
        case ColumnDataType::Float:
            cold = vtkSmartPointer<vtkDoubleArray>::New();

            for (uint32_t j = 0; j < col->GetNumberOfValues(); j++) {
                cold->InsertNextValue(data_[i][j].d);
            }

            cold->SetName(colNames_[i].c_str());
            table->AddColumn(cold);
            break;
        default:
            break;
        }
    }

    return table;
}

#endif

//##############################################################################

CsvException::CsvException(const std::string &message,
			   const std::string &fileName,
			   const CsvException::Id id):
    message_(message), fileName_(fileName), id_(id) {}

//------------------------------------------------------------------------------

//! Returns a short description of the exception. This is the implementation of
//! the virtual function inherited from std::exception
//! \remark This function does not throw exceptions!
const char * CsvException::what () const throw ()
{
    static std::string text;
    text = "XML Exception -- File '" + fileName_ + "': " +message_;
    return text.c_str();
}

//------------------------------------------------------------------------------

//! Returns a more detailled description of the exception, making it possible to
//! use a single exception class for multiple errors.
//! \remark This function does not throw exceptions!
const std::string CsvException::getMessage () const noexcept
{
    return message_;
}

//------------------------------------------------------------------------------

//! Returns the ID of the exception source, making possible to easily
//! distinguish different errors.
//! \remark This function does not throw exceptions!
const std::string CsvException::getFileName () const noexcept
{
    return fileName_;
}

//------------------------------------------------------------------------------

//! Returns the ID of the exception source, making possible to easily
//! distinguish different errors.
//! \remark This function does not throw exceptions!
const CsvException::Id CsvException::getId () const noexcept
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
