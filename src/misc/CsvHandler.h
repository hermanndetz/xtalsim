/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#ifndef __CSV_HANDLER_H__
#define __CSV_HANDLER_H__

#include <memory>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "projectConfigure.h"

#ifdef __VTK__
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkDoubleArray.h>
#include <vtkLongLongArray.h>
#include <vtkVariant.h>
#include <vtkVariantArray.h>
#endif

#include <physics/PeriodicTable.h>
#include <misc/BondInfo.h>
#include <misc/CompositionInfo.h>
#include <misc/Journal.h>

class CsvHandler {
private:
    union dataField {
        long long int n;
        double d;
    };

    enum class ColumnDataType{
        Integer,
        Float
    } ;
    
    const std::string logName_;

    std::vector<std::vector<dataField>> data_;
    std::vector<ColumnDataType> colTypes_;
    std::vector<std::string> colNames_;
    unsigned long colCount_;
    unsigned long rowCount_;

    //! Checks, if string consists of whitespace only.
    bool isWhitespace (const std::string str) const;

    //! Empty strings are also treated as whitespace.
    void loadFile (const std::string &fileName);
    
public:
    //! Constructor
    CsvHandler (const char *logName="CsvHandler");

    //! Constructor
    CsvHandler (const std::string &fileName, const char *logName="CsvHandler");

    //! Destructor
    ~CsvHandler ();

    //! Returns column count.
    unsigned long getColCount (void) const;
    //! Returns row count.
    unsigned long getRowCount (void) const;

    //! Clear data representation in memory.
    void clear ();

    //! Add column of a given type.
    uint32_t addCol (const ColumnDataType dt, const std::string &colName,
                     const dataField defaultValue=dataField{});

    //! Add row.
    void addRow (std::vector<dataField> df);
    //! Add row.
    void addRow (UDTuple df);

    //! Save data to file.
    void save (const std::string &fileName, const std::string separator="\t",
               const std::string &header="");
    //! Save data fo file, replace element id with element name.
    void save (const std::string &fileName, int32_t substCol,
               const std::string separator="\t",
               const std::string &header="");
    //! Save data fo file, replace element id with element name.
    void save (const std::string &fileName, std::vector<int32_t> substCols,
               const std::string separator="\t",
               const std::string &header="");

    //! Print data to std out
    void dump ();

    //! Loads a journal of UDTuples (unsigned long long, double) into object.
    CsvHandler & operator = (const Journal<UDTuple> &src);
    //! Adds a journal of UDTuples (unsigned long long, double) to object.
    CsvHandler & operator += (const Journal<UDTuple> &src);

    //! Loads CompositionInfo data into object.
    CsvHandler & operator = (const CompositionInfo &src);

    //! Loads BondInfo data into object.
    CsvHandler & operator = (const BondInfo &src);

#ifdef __VTK__
    vtkSmartPointer<vtkVariantArray> getCol (unsigned long col);
    vtkSmartPointer<vtkTable> get ();
#endif
};

//------------------------------------------------------------------------------

//! \brief CSV exception

//! Exceptions thrown by the CSV Handler class.
class CsvException : public std::exception {
public:

    //! Identifies source of the Exception.
    enum class Id{
    FileNotFound,
    UnknownColumnType,
    ColumnCountMismatch,
    Unknown
    };

    //! Constructor
    CsvException(const std::string &message="",
            const std::string &fileName="",
            const CsvException::Id id=CsvException::Id::Unknown);

    //! Return description of execution.
    const char * what () const throw ();
    //! Return more detailled error message.
    const std::string getMessage() const noexcept;
    //! Return name of file which caused the exception.
    const std::string getFileName () const noexcept;
    //! Return Exception identification number.
    const CsvException::Id getId() const noexcept;

private:
    std::string message_; //!< error message
    std::string fileName_; //!< name of file which caused the exception
    CsvException::Id id_; //!< Unique ID identifying the exception source

};

#endif


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
