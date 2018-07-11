/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include <string>
#include <map>
#include <random>
#include <memory>
#include <exception>
#include <cmath>
#include <algorithm>
#include <tuple>

#include <easyloggingcpp/easylogging++.h>

#include <physics/Element.h>

typedef std::vector<std::tuple<elementType, double>> MaterialComponents;

//! \brief Stores information to one specific III-V material.

//! It is assumed that materials consits of cations and anions. Their share is
//! automatically balanced, so only their relationship count, e.g. 0.47 and 0.13
//! lead to the same result as 47 and 13.
//! The unique identification of a material is its name. Materials with the same
//! name may lead to undesired behaviour.

class Material{
 private:
    
    std::string name_; //!< Full name of Material.
    double latticeConstant_; //!< lattice constant
    double c11_, c12_, c44_; //!< elastic constants

    //! Cation and anion elements with corresponding share.
    MaterialComponents cations_, anions_;

    //! Logger name.
    const std::string logName_;

    //! Check and set shares appropriately.
    void balanceShares(MaterialComponents &vector);

    //! Pick a random entry of the vector
    elementType getRandomEntry(const MaterialComponents &vector,std::mt19937 &rng) const;
    
 public:

    //! Constructor
    Material(const std::string name, const MaterialComponents cations,
	     const MaterialComponents anions, double lc=0, double c11=0,
	     double c12=0, double c44=0, const char *logName="Material");
    //! Destructor
    ~Material();
    
    //! Check if element is part of cations.
    bool hasCation(const elementType elementID) const;
    //! Check if element is part of anions.
    bool hasAnion(const elementType elementID) const;
    //! Check if element is part of material.
    bool hasElement(const elementType elementID) const;
    //! Check if element is the only cation.
    bool isOnlyCation(const elementType elementID) const;
    //! Check if element is the only anion.
    bool isOnlyAnion(const elementType elementID) const;

    //! Determine randomly one of the cation elements.
    elementType getRandomCation(std::mt19937 &rng) const;
    //! Determine randomly one of the anion elements.
    elementType getRandomAnion(std::mt19937 &rng) const;
    
    //! Return lattice constant as stored in the object.
    double getLatticeConstant(void) const;

    //! Calculate out of plane lattice constant based on given in plane one.
    double getLatticeConstant(double inPlaneLatticeConstant) const;

    //! Compare name of material to given name.
    bool operator==(const std::string &name) const;

    //! Return name of material.
    std::string getName(void) const;
    //! Get elastic constants of material.
    std::vector<double> getElasticConstants(void) const;
    //! Return cations of material.
    const MaterialComponents & getCations() const;
    //! Return anions of material.
    const MaterialComponents & getAnions() const;
    
    //! Return information of object.
    std::string str(bool verbose=false) const;
};

//------------------------------------------------------------------------------

//! \brief Material exception

class MaterialException : public std::exception {
   
 public:

    //! Identifies source of the Exception.
    enum class Id{
    	CationAnionMissing,
    	CationAnionCollection,
    	NameCollection,
    	Unknown
    } ;
    
    //! Constructor
    MaterialException(const std::string &message="",
	    const MaterialException::Id id = MaterialException::Id::Unknown ); 

    //! Return description of execution.
    const char * what () const throw ();
    //! Return more detailled error message.
    const std::string getMessage() const noexcept;
    //! Return Exception identification number.
    const MaterialException::Id getId() const noexcept;

 private:

    std::string message_; //!< error message
    MaterialException::Id id_; //!< unique ID identifying source
    
};


#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
