/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#ifndef __FIELD_3D_INL__
#define __FIELD_3D_INL__

//------------------------------------------------------------------------------

//! Constructor
//! \param x, y, z Initial size in three dimensions
//! \param logName Optional parameter defining the logger's name.
template <class T>
Field3D<T>::Field3D(const indexType x, const indexType y, const indexType z, const char *logName) :
    logName_(logName)
{
    generate(Vector3D<indexType>(x, y, z));
}

//------------------------------------------------------------------------------

//! Constructor
//! \param v Initial size in three dimensions
//! \param logName Optional parameter defining the logger's name.
template <class T>
Field3D<T>::Field3D(const Vector3D<indexType> v, const char *logName) :
    logName_(logName)
{
    generate(v);
}

//------------------------------------------------------------------------------

//! Destructor
template <class T>
Field3D<T>::~Field3D() {

}

//------------------------------------------------------------------------------

//! Returns an element of the field, given by the indices in three dimensions.
//! \param x, y, z Index in field.
//! \return Rererence to specified field element.
template <class T>
T & Field3D<T>::operator()(indexType x, indexType y, indexType z) {
    return operator()(Vector3D<indexType>(x, y, z));
}

//------------------------------------------------------------------------------

//! Returns an element of the field, given by the indices in three dimensions.
//! \param coordinates Vector3D<indexType> defining indices of element in field.
//! \return Rererence to specified field element.
template <class T>
T & Field3D<T>::operator()(Vector3D<indexType> coordinates) {
    if ((coordinates[0] > values_.size()) ||
        (coordinates[1] > values_[0].size() ) ||
        (coordinates[2] > values_[0][0].size())) {
        // \todo throw an exception
        // is there any reason not to use exceptions in template classes?
    }

    return values_[coordinates[0]][coordinates[1]][coordinates[2]];
}

//------------------------------------------------------------------------------

//! Assignment operator
//! \param field Field that shall be assigned to current one.
//! \returns Reference to resulting field after assignment.
template <class T>
T & Field3D<T>::operator=(const Field3D<T> &field) {
    this->generate(field.getSize());

    for (uint32_t i = 0; i < field.getSize()[0]; i++) {
        for (uint32_t j = 0; j < field.getSize()[1]; j++) {
            for (uint32_t k = 0; k < field.getSize()[2]; k++) {
                this->values_(i,j,k) = field.values_(i,j,k);
            }
        }
    }
}

//------------------------------------------------------------------------------

//! Adds a field to the current field element by element.
//! \param field The field to be added.
//! \return Resulting field after addition operation.
template <class T>
Field3D<T> Field3D<T>::operator+(const Field3D<T> &field) {
    Field3D<T> result = Field3D<T>(values_.size(), values_[0].size(), values_[0][0].size());

    if (field.getSize() == this->getSize()) {
        for (uint32_t i = 0; i < values_.size(); i++) {
            for (uint32_t j = 0; j < values_[0].size(); j++) {
                for (uint32_t k = 0; k < values_[0][0].size(); k++) {
                    result(i,j,k) = values_[i][j][k] + field.values_[i][j][k];
                }
            }
        }

    } else {
        //! throw exception
    }

    return result;
}

//------------------------------------------------------------------------------

//! Adds a field to the current field element by element.
//! \param field The field to be added.
//! \return Reference to the resulting field after addition operation.
template <class T>
Field3D<T> & Field3D<T>::operator+=(const Field3D<T> &field) {
    if (field.getSize() == this->getSize()) {
        for (uint32_t i = 0; i < values_.size(); i++) {
            for (uint32_t j = 0; j < values_[0].size(); j++) {
                for (uint32_t k = 0; k < values_[0][0].size(); k++) {
                    values_[i][j][k] += field.values_[i][j][k];
                }
            }
        }

    } else {
        //! throw exception
    }

    return *this;
}

//------------------------------------------------------------------------------

//! Adds a value to each element of the current field.
//! \param value The value to be added.
//! \return Resulting field after addition operation.
template <class T>
Field3D<T> Field3D<T>::operator+(const T value) {
    Field3D<T> result = Field3D<T>(values_.size(), values_[0].size(), values_[0][0].size());

    for (uint32_t i = 0; i < values_.size(); i++) {
        for (uint32_t j = 0; j < values_[0].size(); j++) {
            for (uint32_t k = 0; k < values_[0][0].size(); k++) {
                result(i,j,k) = values_[i][j][k] + value;
            }
        }
    }

    return result;
}

//------------------------------------------------------------------------------

//! Adds a value to each element of the current field.
//! \param value The value to be added.
//! \return Referencto to the resulting field after addition operation.
template <class T>
Field3D<T> & Field3D<T>::operator+=(const T value) {
    for (uint32_t i = 0; i < values_.size(); i++) {
        for (uint32_t j = 0; j < values_[0].size(); j++) {
            for (uint32_t k = 0; k < values_[0][0].size(); k++) {
                values_[i][j][k] += value;
            }
        }
    }

    return *this;
}

//------------------------------------------------------------------------------

//! Subtracts a field from the current field element by element.
//! \param field The field to be subtracted
//! \return Resulting field after subtraction operation.
template <class T>
Field3D<T> Field3D<T>::operator-(const Field3D<T> &field) {
    Field3D<T> result = Field3D<T>(values_.size(), values_[0].size(), values_[0][0].size());

    if (field.getSize() == this->getSize()) {
        for (uint32_t i = 0; i < values_.size(); i++) {
            for (uint32_t j = 0; j < values_[0].size(); j++) {
                for (uint32_t k = 0; k < values_[0][0].size(); k++) {
                    result(i,j,k) = values_[i][j][k] - field.values_[i][j][k];
                }
            }
        }

    } else {
        //! throw exception
    }

    return result;
}

//------------------------------------------------------------------------------

//! Subtracts a field from the current field element by element.
//! \param field The field to be subtracted
//! \return Reference to the esulting field after subtraction operation.
template <class T>
Field3D<T> & Field3D<T>::operator-=(const Field3D<T> &field) {
    if (field.getSize() == this->getSize()) {
        for (uint32_t i = 0; i < values_.size(); i++) {
            for (uint32_t j = 0; j < values_[0].size(); j++) {
                for (uint32_t k = 0; k < values_[0][0].size(); k++) {
                    values_[i][j][k] -= field.values_[i][j][k];
                }
            }
        }

    } else {
        //! throw exception
    }

    return *this;
}

//------------------------------------------------------------------------------

//! Subtracts a value from each element of the current field.
//! \param value The value to be subtracted.
//! \return Resulting field after subtraction operation.
template <class T>
Field3D<T> Field3D<T>::operator-(const T value) {
    Field3D<T> result = Field3D<T>(values_.size(), values_[0].size(), values_[0][0].size());

    for (uint32_t i = 0; i < values_.size(); i++) {
        for (uint32_t j = 0; j < values_[0].size(); j++) {
            for (uint32_t k = 0; k < values_[0][0].size(); k++) {
                result(i,j,k) = values_[i][j][k] - value;
            }
        }
    }

    return result;
}

//------------------------------------------------------------------------------

//! Subtracts a value from each element of the current field.
//! \param value The value to be subtracted.
//! \return Referencto to the resulting field after subtraction operation.
template <class T>
Field3D<T> & Field3D<T>::operator-=(const T value) {
    for (uint32_t i = 0; i < values_.size(); i++) {
        for (uint32_t j = 0; j < values_[0].size(); j++) {
            for (uint32_t k = 0; k < values_[0][0].size(); k++) {
                values_[i][j][k] -= value;
            }
        }
    }

    return *this;
}

//------------------------------------------------------------------------------

//! Multiply field with value point by point.
//! \param value The value each field element is multiplied with.
//! \return Resulting field after multiplication.
template <class T>
Field3D<T> Field3D<T>::operator*(const T value) {
    Field3D<T> result = Field3D<T>(values_.size(), values_[0].size(), values_[0][0].size());

    for (uint32_t i = 0; i < values_.size(); i++) {
        for (uint32_t j = 0; j < values_[0].size(); j++) {
            for (uint32_t k = 0; k < values_[0][0].size(); k++) {
                result(i,j,k) = values_[i][j][k] * value;
            }
        }
    }

    return result;
}

//------------------------------------------------------------------------------

//! Multiply field with value point by point.
//! \param value The value each field element is multiplied with.
//! \return Reference to resulting field after multiplication.
template <class T>
Field3D<T> & Field3D<T>::operator*=(const T value) {
    for (uint32_t i = 0; i < values_.size(); i++) {
        for (uint32_t j = 0; j < values_[0].size(); j++) {
            for (uint32_t k = 0; k < values_[0][0].size(); k++) {
                values_[i][j][k] *= value;
            }
        }
    }

    return *this;
}

//------------------------------------------------------------------------------

//! Returns the maximum value found in the field
//! \return The maximum value of the field.
template <class T>
T Field3D<T>::max(void) const {
    T result{};

    for (uint32_t i = 0; i < values_.size(); i++) {
        for (uint32_t j = 0; j < values_[0].size(); j++) {
            for (uint32_t k = 0; k < values_[0][0].size(); k++) {
                if (values_[i][j][k] > result)
                    result = values_[i][j][k];
            }
        }
    }

    return result;
}

//------------------------------------------------------------------------------

//! Returns the maximum value in a particular layer
//! \param layer index of layer, where maximum value is searched.
//! \param dimension Normal to the layer plane. (0 = (100), 1 = (010), 2 = (001))
//! \return The maximum element in the given layer.
template <class T>
T Field3D<T>::max(indexType layer, indexType dimension) const {
    T result{};

    if (dimension < 3) {
        indexType iMin[3] {0, 0, 0};
        indexType iMax[3] {values_.size(), values_[0].size(), values_[0][0].size()};

        iMin[dimension] = layer;
        iMax[dimension] = layer+1;

        for (indexType i = iMin[0]; i < iMax[0]; i++) {
            for (indexType j = iMin[1]; j < iMax[1]; j++) {
                for (indexType k = iMin[2]; k < iMax[2]; k++) {
                    if (values_[i][j][k] > result)
                        result = values_[i][j][k];
                }
            }
        }
    }

    return result;
}

//------------------------------------------------------------------------------

//! Returns the minimum value found in the field
//! \return The mininum value of the field.
template <class T>
T Field3D<T>::min(void) const {
    T result{};

    for (uint32_t i = 0; i < values_.size(); i++) {
        for (uint32_t j = 0; j < values_[0].size(); j++) {
            for (uint32_t k = 0; k < values_[0][0].size(); k++) {
                if (values_[i][j][k] < result)
                    result = values_[i][j][k];
            }
        }
    }

    return result;
}

//------------------------------------------------------------------------------

//! Returns the maximum value in a particular layer
//! \param layer index of layer, where minimum value is searched.
//! \param dimension Normal to the layer plane. (0 = (100), 1 = (010), 2 = (001))
//! \return The minimum element in the given layer.
template <class T>
T Field3D<T>::min(indexType layer, indexType dimension) const {
    T result{};

    if (dimension < 3) {
        indexType iMin[3] {0, 0, 0};
        indexType iMax[3] {values_.size(), values_[0].size(), values_[0][0].size()};

        iMin[dimension] = layer;
        iMax[dimension] = layer+1;

        for (indexType i = iMin[0]; i < iMax[0]; i++) {
            for (indexType j = iMin[1]; j < iMax[1]; j++) {
                for (indexType k = iMin[2]; k < iMax[2]; k++) {
                    if (values_[i][j][k] < result)
                        result = values_[i][j][k];
                }
            }
        }
    }

    return result;
}


//------------------------------------------------------------------------------

//! Generates field with given dimensions.
//! \param size Vector3D<indexType> Size of field ot be generated.
//! \remark Undefined behavior does not delete values, if called on non-empty field.
template <class T>
void Field3D<T>::generate (const Vector3D<indexType> &size) {
    values_.resize(size[0]);

    for (uint32_t i = 0; i < values_.size(); i++) {
        values_[i].resize(size[1]);

        for (uint32_t j = 0; j < values_[0].size(); j++) {
            values_[i][j].resize(size[2]);
        }
    }
}

//------------------------------------------------------------------------------

//! Extends a field along one dimension.
//! \param amount Number of data points by which field is extended.
//! \param dimension Direction in which field is extended.
//!         (0 = (100), 1 = (010), 2 = (001))
template <class T>
void Field3D<T>::extend (indexType amount, uint8_t dimension) {
    uint32_t oldSize = 0;

    switch (dimension) {
    case 0:
        oldSize = values_.size();
        values_.resize(oldSize+amount);

        for (uint32_t i = 0; i < amount; i++) {
            values_[oldSize+i].resize(values_[0].size());

            for (uint32_t j = 0; j < values_[0][0].size(); j++) {
                values_[oldSize+i][j].resize(values_[0][0].size());
            }
        }
        break;
    case 1:
        oldSize = values_[0].size();

        for (uint32_t i = 0; i < values_[0].size(); i++) {
            values_[i].resize(oldSize+amount);

            for (uint32_t j = 0; j < amount; j++) {
                values_[i][oldSize+j].resize(values_[0][0].size());
            }
        }
        break;
    case 2:
        oldSize = values_[0][0].size();

        for (uint32_t i = 0; i < values_.size(); i++) {
            for (uint32_t j = 0; j < values_[0].size(); j++) {
                values_[i][j].resize(oldSize+amount);
            }
        }
        break;
    default:
        CLOG(WARNING, logName_) << "Illegal dimension " << (int)dimension <<
            " specified.";
        // do nothing
        break;
    }
}

//------------------------------------------------------------------------------

//! Extends field in three dimensions.
//! \param size Vector3D<indexType> Size, by which field is extended in 3D.
template <class T>
void Field3D<T>::extend (const Vector3D<indexType> size) {
    extend(size[0], 0);
    extend(size[1], 1);
    extend(size[2], 2);
}

//------------------------------------------------------------------------------

//! Returns the number of negative data points in given layer.
//! \param layer index of layer, where minimum value is searched.
//! \param dimension Normal to the layer plane. (0 = (100), 1 = (010), 2 = (001))
//! \return Number of negative values in given layer.
template <class T>
uint32_t Field3D<T>::negValuesInLayer (indexType layer, indexType dimension) {
    uint32_t result = 0;

    if (dimension < 3) {
        indexType iMin[3] {0, 0, 0};
        indexType iMax[3] {values_.size(), values_[0].size(), values_[0][0].size()};

        iMin[dimension] = layer;
        iMax[dimension] = layer+1;

        for (indexType i = iMin[0]; i < iMax[0]; i++) {
            for (indexType j = iMin[1]; j < iMax[1]; j++) {
                for (indexType k = iMin[2]; k < iMax[2]; k++) {
                    if (values_[i][j][k] < 0)
                        result++;
                }
            }
        }
    }

    return result;
}

//------------------------------------------------------------------------------
//! Returns the number of positive data points in given layer.
//! \param layer index of layer, where maximum value is searched.
//! \param dimension Normal to the layer plane. (0 = (100), 1 = (010), 2 = (001))
//! \return Number of positive values in given layer.
template <class T>
uint32_t Field3D<T>::posValuesInLayer (indexType layer, indexType dimension) {
    uint32_t result = 0;

    if (dimension < 3) {
        indexType iMin[3] {0, 0, 0};
        indexType iMax[3] {values_.size(), values_[0].size(), values_[0][0].size()};

        iMin[dimension] = layer;
        iMax[dimension] = layer+1;

        for (indexType i = iMin[0]; i < iMax[0]; i++) {
            for (indexType j = iMin[1]; j < iMax[1]; j++) {
                for (indexType k = iMin[2]; k < iMax[2]; k++) {
                    if (values_[i][j][k] > 0)
                        result++;
                }
            }
        }
    }

    return result;
}

//------------------------------------------------------------------------------

//! Adds a value to one point in the field
//! \param x, y, z Index in field.
//! \param value The value, which is added to the given element.
//! \remark Useful to construct a field step by step.
template <class T>
void Field3D<T>::addToPoint(const indexType x, const indexType y, 
                        const indexType z, const T value) {
    if ((x < values_.size()) && 
        (y < values_[0].size()) && 
        (z < values_[0][0].size())) {
        values_[x][y][z] += value;
    } else
    {
        // exception
    }
}

//------------------------------------------------------------------------------

//! Adds a value to one point in the field
//! \param coordinates Index in field.
//! \param value The value, which is added to the given element.
//! \remark Useful to construct a field step by step.
template <class T>
void Field3D<T>::addToPoint(const Vector3D<indexType> coordinates, const T value) {
    if ((coordinates[0] < values_.size()) && 
        (coordinates[1] < values_[0].size()) && 
        (coordinates[2] < values_[0][0].size())) {
        values_[coordinates[0]][coordinates[1]][coordinates[2]] += value;
    } else
    {
        // exception
    }
}

//------------------------------------------------------------------------------

//! Returns the field size in three dimensions in form of a Vector3D<indexType<.
//! \ŗeturn Vector3D<indexType> defining the size of the field in 3 dimensions.
template <class T>
Vector3D<indexType> Field3D<T>::getSize(void) const {
    return Vector3D<indexType>(values_.size(), values_[0].size(), values_[0][0].size());
}

//------------------------------------------------------------------------------

//! Prints values of the field to a stream, separated by tabs.
//! \param maxX maximum index in x direction (index 0) to be printed
//! \param maxY maximum index in y direction (index 1) to be printed
//! \param maxZ maximum index in z direction (index 2) to be printed
//! \param digits defines number of digits for floating point numbers
//! \param separator character to separate individual columns
//! \return String containing values up to specified maximum indices.
template <class T>
const std::string Field3D<T>::str(int maxX, int maxY, int maxZ, int digits, char separator) const {
    std::ostringstream msg;

    msg << std::fixed << std::setprecision(digits);

    if (maxX == -1)
        maxX = values_.size();

    if (maxY == -1)
        maxY = values_[0].size();

    if (maxZ == -1)
        maxZ = values_[0][0].size();

    for (uint32_t k = 0; k < values_[0][0].size(), k < maxZ; k++) {
        msg << '{' << std::endl;

        for (uint32_t j = 0; j < values_[0].size(), j < maxY; j++) {
            msg << '[';

            for (uint32_t i = 0; i < values_.size(), i < maxX; i++) {
                msg << (values_[i][j][k] >= 0 ? " ":"") << values_[i][j][k] << separator;
            }
            
            msg << ']' << std::endl;
        }

        msg << '}' << std::endl;
    }

    return msg.str();
}

//------------------------------------------------------------------------------

//! Transforms data into a field with lower resolution
//! \param factorX factor, by which resolution is downsampled in x (index 0) direction
//! \param factorY factor, by which resolution is downsampled in y (index 1) direction
//! \param factorZ factor, by which resolution is downsampled in z (index 2) direction
//! \param op method that is used to transform resolution.
//! \return Returns downsampled field.
//! \note This is a basic version of a downsampling function. The data points are reduced
//!         with left-sided operations. For future versions, centered operations would
//!         be better.
template <class T>
Field3D<T> Field3D<T>::downsampling(uint32_t factorX, uint32_t factorY, uint32_t factorZ, DownsamplingFieldOperation op) {
    uint32_t dim[3];

    dim[0] = values_.size() / factorX;
    dim[1] = values_[0].size() / factorY;
    dim[2] = values_[0][0].size() / factorZ;

    Field3D<T> result = Field3D<T>(dim[0], dim[1], dim[2]);

    for (uint32_t i = 0; i < dim[0]; i++) {
        for (uint32_t j = 0; j < dim[1]; j++) {
            for (uint32_t k = 0; k < dim[2]; k++) {
                for (uint32_t ii = 0; ii < factorX; ii++) {
                    for (uint32_t jj = 0; jj < factorY; jj++) {
                        for (uint32_t kk = 0; kk < factorZ; kk++) {
                            switch (op) {
                            case DownsamplingFieldOperation::Max:
                                if (values_[i*factorX+ii][j*factorY+jj][k*factorZ+kk] > result(i,j,k))
                                    result(i,j,k) = values_[i*factorX+ii][j*factorY+jj][k*factorZ+kk];
                                break;
                            case DownsamplingFieldOperation::Avg:
                                [[fallthrough]]
                            case DownsamplingFieldOperation::Sum:
                                result.addToPoint(i,j,k, values_[i*factorX+ii][j*factorY+jj][k*factorZ+kk]);
                                break;

                            }

                        }
                    }
                }
            }
        }
    }

    if (op == DownsamplingFieldOperation::Avg) {
        for (uint32_t i = 0; i < dim[0]; i++) {
            for (uint32_t j = 0; j < dim[1]; j++) {
                for (uint32_t k = 0; k < dim[2]; k++) {
                    result(i,j,k) = result(i,j,k) / (factorX * factorY * factorZ);
                }
            }
        }
    }

    return result;
}

//------------------------------------------------------------------------------

//! Transforms data into a field where a unit cell is compressed into one layer
//! \return Returns flattened field.
template <class T>
Field3D<T> Field3D<T>::flatten() {
    uint32_t dim[3];

    dim[0] = values_.size();
    dim[1] = values_[0].size();
    dim[2] = values_[0][0].size();

    Field3D<T> result = Field3D<T>(dim[0], dim[1], dim[2]/4);

    for (uint32_t i = 0; i < dim[0]; i++) {
        for (uint32_t j = 0; j < dim[1]; j++) {
            for (uint32_t k = 0; k < dim[2]; k++) {
                if (std::fpclassify(values_[i][j][k]) == FP_ZERO)
                    continue;

                result(i,j,k/4) = values_[i][j][k];
            }
        }
    }

    return result;
}

//------------------------------------------------------------------------------

//! Interpolates elements in a field.
//! \return Returns the interpolated field.
//! \remark This is a very basic function that currently works with zincblende 
//!         lattices.
//! \todo   Interpolates values, which equal 0, no matter, if they were never
//!         set or actually set to 0.
template <class T>
Field3D<T> Field3D<T>::interpolate() {
    int32_t dim[3];

    dim[0] = values_.size();
    dim[1] = values_[0].size();
    dim[2] = values_[0][0].size();

    Field3D<T> result = Field3D<T>(dim[0], dim[1], dim[2]);

    for (int32_t i = 0; i < dim[0]; i++) {
        for (int32_t j = 0; j < dim[1]; j++) {
            for (int32_t k = 0; k < dim[2]; k++) {
                if (std::fpclassify(values_[i][j][k]) == FP_ZERO) {
                    // found unassigned value --> calculate average wrt neighbors
                    T tmp{};
                    uint32_t cnt = 0;
                    
                    if ((i-1) >= 0) {
                        tmp += values_[i-1][j][k];
                        cnt ++;
                    }

                    if ((i+1) < dim[0]) {
                        tmp += values_[i+1][j][k];
                        cnt ++;
                    }

                    if ((j-1) >= 0) {
                        tmp += values_[i][j-1][k];
                        cnt ++;
                    }

                    if ((j+1) < dim[1]) {
                        tmp += values_[i][j+1][k];
                        cnt ++;
                    }

                    if (cnt > 0) {
                        result(i,j,k) = tmp/cnt;
                    }

                } else {
                    result(i,j,k) = values_[i][j][k];
                }
            }
        }
    }

    return result;
}

//------------------------------------------------------------------------------

//! Returns general information on all layers
//! \return string containing information about the field.
template <class T>
std::string Field3D<T>::getReport () {
    std::stringstream result{};

    result << "Field size: " << values_.size() << "x" << values_[0].size() << "x" << values_[0][0].size() << std::endl;

    for (indexType i = 0; i < values_[0][0].size(); i++) {
        result << "Information on layer " << i << std::endl;
        result << getReport(i, 2);
    }

    return result.str();
}

//------------------------------------------------------------------------------

//! Returns general information on one given layer.
//! \return string containing information about a given layer in the field.
template <class T>
std::string Field3D<T>::getReport (indexType layer, indexType dimension) {
    std::stringstream result{};

    result << "Max: " << std::scientific << this->max(layer, dimension) << "\tMin: " << std::scientific << this->min(layer, dimension) << std::endl;
    result << "Values > 0: " << posValuesInLayer(layer, dimension) << "\tValues < 0: " << negValuesInLayer(layer, dimension) << std::endl;

    return result.str();
}

//------------------------------------------------------------------------------

//! Save a single layer of the field structue into a Gwyddion Simple Field file.
//! \param fileName Specifies the output file name.
//! \param layer Index of the layer to be printed.
//! \param dimension normal to the layer plane. (0 = (100), 1 = (010), 2 = (001))
//! \param xRes resolution (no. of data points) of output field in x direction
//! \param yRes resolution (no. of data points) of output field in y direction
//! \param xOff offset of output field in x direction
//! \param yOff offset of output field in y direction
//! \param xReal size of field in x direction in real coordinates
//! \param yReal size of field in y direction in real coordinates
//! \param xyUnits specifies the length unit of x and y scale 
//!         (without prefix, i.e. m but not nm or km)
//! \param zUnits specifies the length unit of x and y scale 
//!         (without prefix, i.e. m but not nm or km)
//! \param title name of data set.
//! \remark Field elements have to be written as float (not double!) to conform
//!         with the file format definition.
//! \remark File format documentation: http://gwyddion.net/documentation/user-guide-en/gsf.html
template <class T>
void Field3D<T>::saveToGSFFile(std::string fileName, uint32_t layer, 
        uint32_t dimension, uint32_t xRes, uint32_t yRes, double xOff, 
        double yOff, double xReal, double yReal, std::string xyUnits, 
        std::string zUnits, std::string title) {
    std::ofstream gsfFile (fileName, std::ios::out);

    CLOG(TRACE, logName_) << "Writing atoms to GSF file " << fileName;

    // file format description
    gsfFile << "Gwyddion Simple Field 1.0" << std::endl;

    // mandatory fields
    gsfFile << "XRes = " << xRes << std::endl;
    gsfFile << "YRes = " << yRes << std::endl;

    // optional fields
    if (! (FP_ZERO == std::fpclassify(xOff)) )
        gsfFile << "XOffset = " << xOff << std::endl;

    if (! (FP_ZERO == std::fpclassify(yOff)) )
        gsfFile << "YOffset = " << yOff << std::endl;

    if (! (FP_ZERO == std::fpclassify(xReal - 1.0)) )
        gsfFile << "XReal = " << xReal << std::endl;

    if (! (FP_ZERO == std::fpclassify(yReal - 1.0)) )        
        gsfFile << "YReal = " << yReal << std::endl;

    if (xyUnits != "")
        gsfFile << "XYUnits = " << xyUnits << std::endl;

    if (zUnits != "")
        gsfFile << "ZUnits = " << zUnits << std::endl;

    if (title != "")
        gsfFile << "Title = " << title << std::endl;

    // padding before data (min one, max four lines containing 0x00
    uint32_t p = 0;
    char padding=0x00;

    for (p = 0; p < ((gsfFile.tellp() % 4) + 1); p++)
        gsfFile.write((char*)&padding,sizeof(padding));

    float f = 0.0;

    indexType iMin[3] {0, 0, 0};
    ulong iMax[3] {values_.size(), values_[0].size(), values_[0][0].size()};

    switch (dimension) {
    case 0:
        for (indexType i = iMin[1]; i < iMax[1]; i++) {
            for (indexType j = iMin[2]; j < iMax[2]; j++) {
                f = values_[layer][i][j];
                gsfFile.write( reinterpret_cast<const char*>( &f ), sizeof( float ));
            }
        }
        break;
    case 1:
        for (indexType i = iMin[0]; i < iMax[0]; i++) {
            for (indexType j = iMin[2]; j < iMax[2]; j++) {
                f = values_[i][layer][j];
                gsfFile.write( reinterpret_cast<const char*>( &f ), sizeof( float ));
            }
        }
        break;
    case 2:
        for (indexType i = iMin[0]; i < iMax[0]; i++) {
            for (indexType j = iMin[1]; j < iMax[1]; j++) {
                f = values_[i][j][layer];
                gsfFile.write( reinterpret_cast<const char*>( &f ), sizeof( float ));
            }
        }
        break;
    default:
        break;
    }

    gsfFile.close();

    CLOG(TRACE, logName_) << "Field successfully written.";
}

//------------------------------------------------------------------------------

//! Generates a test pattern of straight lines.
//! \param layer layer index, where pattern is generated.
//! \param dimension normal to the layer plane. (0 = (100), 1 = (010), 2 = (001))
//! \param lineOrientation defines the orientation of the pattern.
//! \param valueA, valueB numerical values, which are arranged in the pattern.
//! \param widthA, widthB lengths of the half-periods of the pattern.
template <class T>
void Field3D<T>::generateLines(indexType layer, indexType dimension, 
        indexType lineOrientation, T valueA, T valueB, uint32_t widthA, 
        uint32_t widthB) {
    uint32_t tmpWidth = 0;

    if (dimension < 3) {
        indexType iMin[3] {0, 0, 0};
        indexType iMax[3] {values_.size(), values_[0].size(), values_[0][0].size()};

        iMin[dimension] = layer;
        iMax[dimension] = layer+1;

        for (indexType i = iMin[0]; i < iMax[0]; i++) {
            for (indexType j = iMin[1]; j < iMax[1]; j++) {
                for (indexType k = iMin[2]; k < iMax[2]; k++) {
                    switch (dimension) {
                    case 0:
                        switch (lineOrientation) {
                        case 1:
                            if ((j % (widthA + widthB)) < widthA) 
                                values_[i][j][k] = valueA;
                            else
                                values_[i][j][k] = valueB;
                            break;
                        case 2:
                            if ((k % (widthA + widthB)) < widthA) 
                                values_[i][j][k] = valueA;
                            else
                                values_[i][j][k] = valueB;
                            break;
                            break;
                        default:
                            break;
                        }
                        break;
                    case 1:
                        switch (lineOrientation) {
                        case 0:
                            if ((i % (widthA + widthB)) < widthA) 
                                values_[i][j][k] = valueA;
                            else
                                values_[i][j][k] = valueB;
                            break;
                        case 2:
                            if ((k % (widthA + widthB)) < widthA) 
                                values_[i][j][k] = valueA;
                            else
                                values_[i][j][k] = valueB;
                            break;
                            break;
                        default:
                            break;
                        }
                        break;
                    case 2:
                        switch (lineOrientation) {
                        case 0:
                            if ((i % (widthA + widthB)) < widthA) 
                                values_[i][j][k] = valueA;
                            else
                                values_[i][j][k] = valueB;
                            break;
                        case 1:
                            if ((j % (widthA + widthB)) < widthA) 
                                values_[i][j][k] = valueA;
                            else
                                values_[i][j][k] = valueB;
                            break;
                            break;
                        default:
                            break;
                        }
                        break;
                    default:
                        break;
                    }
                }
            }
        }
    }
}

//------------------------------------------------------------------------------

//! Generates a checkerboard test pattern.
//! \param layer layer index, where pattern is generated.
//! \param dimension normal to the layer plane. (0 = (100), 1 = (010), 2 = (001))
//! \param valueA, valueB numerical values, which are arranged in the pattern.
//! \param widthA, widthB lengths of the half-periods of the pattern.
template <class T>
void Field3D<T>::generateCheckerboard(indexType layer, indexType dimension, 
        T valueA, T valueB, uint32_t widthA, uint32_t widthB) {
    uint32_t tmpWidth = 0;

    if (dimension < 3) {
        indexType iMin[3] {0, 0, 0};
        indexType iMax[3] {values_.size(), values_[0].size(), values_[0][0].size()};
        indexType inPlane[2];

        switch (dimension) {
        case 0:
            inPlane[0] = 1;
            inPlane[1] = 2;
            break;
        case 1:
            inPlane[0] = 0;
            inPlane[1] = 2;
            break;
        case 2:
            inPlane[0] = 0;
            inPlane[1] = 1;
            break;
        }

        iMin[dimension] = layer;
        iMax[dimension] = layer+1;

        indexType ijk[3]{};

        for (ijk[0] = iMin[0]; ijk[0] < iMax[0]; ijk[0]++) {
            for (ijk[1] = iMin[1]; ijk[1] < iMax[1]; ijk[1]++) {
                for (ijk[2] = iMin[2]; ijk[2] < iMax[2]; ijk[2]++) {
                    if ((ijk[inPlane[0]] % (widthA + widthB)) < widthA) 
                        if ((ijk[inPlane[1]] % (widthA + widthB)) < widthB)
                            values_[ijk[0]][ijk[1]][ijk[2]] = valueA;
                        else
                            values_[ijk[0]][ijk[1]][ijk[2]] = valueB;
                    else
                        if ((ijk[inPlane[1]] % (widthA + widthB)) < widthB)
                            values_[ijk[0]][ijk[1]][ijk[2]] = valueB;
                        else
                            values_[ijk[0]][ijk[1]][ijk[2]] = valueA;
                }
            }
        }
    }
}

#endif

//##############################################################################

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
