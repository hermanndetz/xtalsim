/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#ifndef __FIELD_3D_H__
#define __FIELD_3D_H__

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <math.h>

#include <easyloggingcpp/easylogging++.h>

#include <physics/Vector3D.h>

typedef enum {
    Max,
    AbsMax,
    Sum,
    Avg} DownsamplingFieldOperation;

//! Stores data in a three-dimensional field.
//! Defined as template to serve for basic data types, but e.g. also Vector3D.
template <class T>
class Field3D{
private:
    std::vector<std::vector<std::vector<T>>> values_{};

    //! Logger name.
    const std::string logName_;

    //! Generates field with given dimensions.
    void generate (const Vector3D<indexType> &size);
    //! Extends field in one specific dimension.
    void extend (indexType amount, uint8_t dimension);
    //! Extends field in three dimensions.
    void extend (const Vector3D<indexType> size);

    //! Returns the number of negative data points in given layer.
    uint32_t negValuesInLayer (indexType layer, indexType dimension);
    //! Returns the number of positive data points in given layer.
    uint32_t posValuesInLayer (indexType layer, indexType dimension);

public:
    //! Constructor
    Field3D (const indexType x, const indexType y, const indexType z, const char *logName="Field");
    //! Constructor
    Field3D (const Vector3D<indexType> v, const char *logName="Field");
    //! Destructor
    ~Field3D ();

    //! Returns single element
    T & operator()(indexType x, indexType y, indexType z);
    //! Returns single element
    T & operator()(Vector3D<indexType> coordinates);

    //! Assignment operator
    T & operator=(const Field3D<T> &field);

    //! Add two fields
    Field3D<T> operator+(const Field3D<T> &field);
    //! Add two fields
    Field3D<T> & operator+=(const Field3D<T> &field);
    //! Add value point by point
    Field3D<T> operator+(const T value);
    //! Add value point by point
    Field3D<T> & operator+=(const T value);
    //! Subtract two fields
    Field3D<T> operator-(const Field3D<T> &field);
    //! Subtract two fields
    Field3D<T> & operator-=(const Field3D<T> &field);
    //! Subtract value point by point
    Field3D<T> operator-(const T value);
    //! Subtract value point by point
    Field3D<T> & operator-=(const T value);
    //! Multiply with value point by point
    Field3D<T> operator*(const T value);
    //! Multiply with value point by point
    Field3D<T> & operator*=(const T value);

    //! Returns the maximum value found in the field
    T max(void) const;
    //! Returns the maximum value in a particular layer
    T max(indexType layer, indexType dimension) const;
    //! Returns the maximum absolute value found in the field
    T maxAbs(void) const;
    //! Returns the minimum value found in the field
    T min(void) const;
    //! Returns the minimum value in a particular layer
    T min(indexType layer, indexType dimension) const;
    //! Returns the minimum absolute value found in the field
    T minAbs(void) const;

    //! Adds a value to one point in the field
    void addToPoint(const indexType x, const indexType y, const indexType z, const T value);
    //! Adds a value to one point in the field
    void addToPoint(const Vector3D<indexType> coordinates, const T value);

    //! Returns size of field in three dimensions.
    Vector3D<indexType> getSize(void) const;

    //! Return values as string
    const std::string str(int maxX=-1, int maxY=-1, int maxZ=-1, int digits=4, char separator='\t') const;

    //! Transforms data into a field with lower resolution
    Field3D<T> downsampling(uint32_t factorX, uint32_t factorY, uint32_t factorZ, DownsamplingFieldOperation op);

    //! Transforms data into a field where a unit cell is compressed into one layer
    Field3D<T> flatten();

    //! Interpolates elements in a field.
    Field3D<T> interpolate();

    //! Returns general information on all layers (slices in index=2 direction)
    std::string getReport();

    //! Returns general information on one given layer
    std::string getReport(indexType layer, indexType dimension);

    //! Save a single layer of the field structue into a Gwyddion Simple Field file.
    void saveToGSFFile(std::string fileName, uint32_t layer, uint32_t dimension, 
            uint32_t xRes=1, uint32_t yRes=1, double xOff=0.0, double yOff=0.0, 
            double xReal=1.0, double yReal=1.0, std::string xyUnits="m", 
            std::string zUnits="m", std::string title="");

    //! Generates a test pattern of straight lines.
    void generateLines(indexType layer, indexType dimension, indexType lineOrientation, T valueA, T valueB, uint32_t widthA, uint32_t widthB);

    //! Generates a checkerboard test pattern.
    void generateCheckerboard(indexType layer, indexType dimension, T valueA, T valueB, uint32_t widthA, uint32_t widthB);
};

// Unfortunately it is necessary to include the code of template member
// functions in the header file.
// Other solutions are documented at
// http://www.codeproject.com/Articles/48575/How-to-define-a-template-class-in-a-h-file-and-imp 
#include "Field3D.inl"

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
