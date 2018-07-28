/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "TersoffParameter.h"

//------------------------------------------------------------------------------

TersoffParameter::TersoffParameter(elementType elementId1,
		   elementType elementId2, std::string _name, std::string _comment,
           std::string _mode,
		   double _a, double _b, double _beta,double _c,
		   double _d, double _d0, double _dd, double _delta,
		   double _gamma, double _h, double _lambda, double _lambda3,
           double _mu, double _n, double _r0, double _rCut, double _s)
    :name(_name), comment(_comment), mode(_mode),
     a(_a), b(_b), beta(_beta), c(_c), d(_d), d0(_d0), dd(_dd), delta(_delta),
     gamma(_gamma), h(_h), lambda(_lambda), lambda3(_lambda3), mu(_mu),
     n(_n), r0(_r0), rCut(_rCut), s(_s)
{
    if (elementId1 < elementId2) {
	elementIdLow = elementId1;
	elementIdHigh = elementId2;
    }
    else{
	elementIdLow = elementId2;
	elementIdHigh = elementId1;
    }
}

//------------------------------------------------------------------------------

TersoffParameter::~TersoffParameter() {}

//------------------------------------------------------------------------------

TersoffParameter::TersoffParameter(const TersoffParameter &parameter)
{
    elementIdLow = parameter.elementIdLow;
    elementIdHigh = parameter.elementIdHigh;
    name = parameter.name;    
    comment = parameter.comment;
    mode = parameter.mode;
    a = parameter.a;
    b = parameter.b;
    beta = parameter.beta;
    c = parameter.c;
    d = parameter.d;
    d0 = parameter.d0;
    dd = parameter.dd;
    delta = parameter.delta;
    gamma = parameter.gamma;
    h = parameter.h;
    lambda = parameter.lambda;
    lambda3 = parameter.lambda3;
    mu = parameter.mu;
    n = parameter.n;
    r0 = parameter.r0;
    rCut = parameter.rCut;
    s = parameter.s;
    
}
    
//------------------------------------------------------------------------------

//! Introduces a smooth transition between 0 and 1 following a cosine shape in
//! the range from [rCut-d,rCut+d].
//! \param distance Value for which a value shall be computed.
//! \return Value of smooth transition curve for the given value.
double TersoffParameter::cutoffRaisedCosine_(const double distance) const
{
    
    if (distance < (rCut - dd))
	return(1.0);
    else if (distance > (rCut + dd))
	return(0.0);
    else{
        double tmp =cos( M_PI/(4*dd) * (distance - (rCut - dd)));
	return(tmp*tmp);
    }
}

//------------------------------------------------------------------------------

//! Checks if the given parameter is above the cutoff region.
//! \sa cutoffRaisedCosine_
//! \param distance Value that shall be checked.
//! \return True if above, false otherwise.
bool TersoffParameter::isAboveCutoff(const double distance) const
{
    if (distance > (rCut+dd) )
	return true;
    else
	return false;
}

//------------------------------------------------------------------------------

//! The angular term is calculated between two vector both starting at the same
//! location.
//! \note The parameter lengthVector1 can of course be calculated also within
//! the function however it uses sqrt() which is costly.
//! \param vector1 Vector to neighbor atom whose energy share is calculated.
//! \param lengthVector1 Absolute length of vector1.
//! \param vector2 Vector to neighbor atom whose angular term is calculated.
//! \return Angular term for energy calculations.
double TersoffParameter::getAngularCoefficient(const Vector3D<spaceType>vector1,
				const double lengthVector1,
				const Vector3D<spaceType> vector2) const
{

    double angle, g, lengthVector2;

    lengthVector2 = sqrt(vector2.squaredLength());

    if (isAboveCutoff(lengthVector2)) return 0.0;
    
    angle = vector1.angle(vector2);

    g = 1.0 + pow(c, 2.0) / pow(d, 2.0) -
	pow(c,2.0) / (pow(d,2.0) + pow(h-cos(angle),2.0));

    if (mode == TERSOFF_POTENTIAL_S){

	return(cutoffRaisedCosine_(lengthVector2) * g *
	       exp( pow(lambda,3.0) * pow(lengthVector1-lengthVector2,3.0)) );

    }
    else{
	if ( FP_ZERO == std::fpclassify(lambda3)){
	    return(cutoffRaisedCosine_(lengthVector2) * g);
	}
	else{
	    return(cutoffRaisedCosine_(lengthVector2) * g *
		   exp( pow(lambda3,3.0) *pow(lengthVector1-lengthVector2,3.0)) );
	}
    }

    return 0.0;
}

//------------------------------------------------------------------------------

//! The energy for an atom with a certain distance and specific angular terms is
//! calculated.
//! \param distance Distance of the atom to calculate.
//! \param angular Previously calculated angular term.
//! \return Calculated energy.
double TersoffParameter::getEnergy(const double distance,
				  const double angular) const
{

    double b, Vr, Va, fc;
    
    b = pow(1.0 + pow(gamma * angular, n), (-1.0) / (2.0 * n));

    if (mode == TERSOFF_POTENTIAL_S){
	    
	Vr = d0 / (s-1.0) * exp((-1.0) * beta *sqrt(2.0 * s) *(distance-r0));
	Va = s * d0 / (s-1.0) * exp((-1.0) * beta *sqrt(2.0 / s) * (distance-r0));

    }
    else{

	Vr = a * exp((-1.0)* lambda * distance);
	Va = b * exp((-1.0)* mu * distance);
	    
    }

    fc = cutoffRaisedCosine_(distance);
    return 0.5 * (fc * (Vr -b *Va));

}
    
//------------------------------------------------------------------------------

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
