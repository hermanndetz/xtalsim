/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "Field3D.h"

//------------------------------------------------------------------------------

//! Returns the maximum absolute value found in the field
template <>
double Field3D<double>::maxAbs(void) const {
    double result{};

    for (uint32_t i = 0; i < values_.size(); i++) {
        for (uint32_t j = 0; j < values_[0].size(); j++) {
            for (uint32_t k = 0; k < values_[0][0].size(); k++) {
                if (fabs(values_[i][j][k]) > result)
                    result = fabs(values_[i][j][k]);
            }
        }
    }

    return result;
}

//------------------------------------------------------------------------------

//! Returns the minimum absolute value found in the field
template <>
double Field3D<double>::minAbs(void) const {
    double result{};

    for (uint32_t i = 0; i < values_.size(); i++) {
        for (uint32_t j = 0; j < values_[0].size(); j++) {
            for (uint32_t k = 0; k < values_[0][0].size(); k++) {
                if (fabs(values_[i][j][k]) > result)
                    result = fabs(values_[i][j][k]);
            }
        }
    }

    return result;
}

//------------------------------------------------------------------------------

//! Returns the number of negative data points in given layer.
template <>
uint32_t Field3D<double>::negValuesInLayer (indexType layer, indexType dimension) {
    uint32_t result = 0;

    if (dimension < 3) {
        indexType iMin[3] {0, 0, 0};
        indexType iMax[3] {(indexType)values_.size(),
                (indexType)values_[0].size(), (indexType)values_[0][0].size()};

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
//! Returns the number of positive data points in given layer.
template <>
uint32_t Field3D<double>::posValuesInLayer (indexType layer, indexType dimension) {
    uint32_t result = 0;

    if (dimension < 3) {
        indexType iMin[3] {0, 0, 0};
        indexType iMax[3] {(indexType)values_.size(),
                (indexType)values_[0].size(), (indexType)values_[0][0].size()};

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

//! Transforms data into a field with lower resolution
//! \param factorX factor, by which resolution is downsampled in x (index 0) direction
//! \param factorY factor, by which resolution is downsampled in y (index 1) direction
//! \param factorZ factor, by which resolution is downsampled in z (index 2) direction
//! \param op method that is used to transform resolution.
//! \note This is a basic version of a downsampling function. The data points are reduced
//!         with left-sided operations. For future versions, centered operations would
//!         be better.
template <>
Field3D<double> Field3D<double>::downsampling(uint32_t factorX, uint32_t factorY, uint32_t factorZ, DownsamplingFieldOperation op) {
    uint32_t dim[3];

    dim[0] = values_.size() / factorX;
    dim[1] = values_[0].size() / factorY;
    dim[2] = values_[0][0].size() / factorZ;

    Field3D<double> result = Field3D<double>(dim[0], dim[1], dim[2]);

    for (uint32_t i = 0; i < dim[0]; i++) {
        for (uint32_t j = 0; j < dim[1]; j++) {
            for (uint32_t k = 0; k < dim[2]; k++) {
                for (uint32_t ii = 0; ii < factorX; ii++) {
                    for (uint32_t jj = 0; jj < factorX; jj++) {
                        for (uint32_t kk = 0; kk < factorZ; kk++) {
                            switch (op) {
                            case DownsamplingFieldOperation::Max:
                                if (values_[i*factorX+ii][j*factorY+jj][k*factorZ+kk] > result(i,j,k))
                                    result(i,j,k) = values_[i*factorX+ii][j*factorY+jj][k*factorZ+kk];
                                break;
                            case DownsamplingFieldOperation::AbsMax:
                                if (fabs(values_[i*factorX+ii][j*factorY+jj][k*factorZ+kk]) > fabs(result(i,j,k)))
                                    result(i,j,k) = fabs(values_[i*factorX+ii][j*factorY+jj][k*factorZ+kk]);
                                break;
                            case DownsamplingFieldOperation::Sum:
                                result(i,j,k) += values_[i*factorX+ii][j*factorY+jj][k*factorZ+kk];
                                break;

                            }

                        }
                    }
                }
            }
        }
    }

    return result;
}

//##############################################################################


// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
