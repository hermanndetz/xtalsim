/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __TERSOFF_PARAMETER_H__
#define __TERSOFF_PARAMETER_H__

#include <stdint.h>
#include <string>
#include <cmath>

#include <physics/Element.h>
#include <physics/Vector3D.h>

//! \todo add description of constants!
const std::string TERSOFF_POTENTIAL_S="S";
const std::string TERSOFF_POTENTIAL_A="A";

//! \brief Stores Tersoff parameter for bond between elements.

//! Specific parameters for a unique combination of two elements are stored  in
//! this class. These parameters are used to calculate the bond energy.

class TersoffParameter{

    //! Introduces a smooth transition between 0 and 1.
    double cutoffRaisedCosine_(const double distance) const;
    
 public:

    //! Name of parameter set.
    std::string name;
    //! Comment (e.g. to specify origin of parameter set)
    std::string comment;
    //! Mode of parameter set. Valid options are S/??.
    std::string mode;
    //! lower ID of involved elements
    elementType elementIdLow;
    //! higher ID of involved elements
    elementType elementIdHigh;
    double a, b, beta, c, d, d0, dd, delta, gamma;
    double h, lambda, lambda3, mu, n, r0, rCut, s;
    
    //! Constructor
    TersoffParameter(elementType elementId1, elementType elementId2,
		   std::string _name, std::string _comment, std::string _mode,
		   double _a, double _b, double _beta,double _c,
		   double _d, double _d0, double _dd, double _delta,
		   double _gamma, double _h, double _lambda, double _lambda3,
		   double _mu, double _n, double _r0, double _rCut, double _s);

    //! Copy Constructor
    TersoffParameter(const TersoffParameter &parameter);
    //! Destructor
    ~TersoffParameter();

    //! Check if given distance is above the cutoff region.
    bool isAboveCutoff(const double distance) const;

    //! Calculates angular terms for the energy calculations.
    double getAngularCoefficient(const Vector3D<spaceType> vector1,
				const double lengthVector1,
				const Vector3D<spaceType> vector2) const;

    //! Calculate energy for single atom.
    double getEnergy(const double distance,
		    const double angular) const;
    
};

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
