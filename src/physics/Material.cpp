/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#include "Material.h"

#include <sstream>

//! Initializes the random number generator.
//! \throws MaterialException If either anions or cations not specified an
//! exception is thrown.
Material::Material(const std::string name,
		   MaterialComponents cations, MaterialComponents anions, double lc,
		   double c11, double c12, double c44,
		   const char *logName):
    name_(name), latticeConstant_(lc), c11_(c11), c12_(c12), c44_(c44),
    cations_(cations), anions_(anions), logName_(logName)
{

    if (cations_.empty() || anions_.empty()){
	CLOG(ERROR, logName_) << "Either cations or anions not set";
	throw MaterialException("Cations or Anions missing",
				MaterialException::Id::CationAnionMissing);
    }

    CLOG(TRACE, logName_) << "Creating material " << name;
    balanceShares(cations_);
    balanceShares(anions_);
}

//------------------------------------------------------------------------------

Material::~Material() {}

//------------------------------------------------------------------------------

//! The shares are set by the user but for the program to run smoothely they
//! have to sum up to 1. This is assured in this method by summing up all of
//! them and scaling each share appropriately.
//! \param map The std::map to check.
void Material::balanceShares(MaterialComponents &components)
{

    double sum=0, factor;

    CLOG(TRACE, logName_) << "starting to balance element share";
    
    for (auto el: components) sum += std::get<1>(el);
    factor = 1/sum;

    // If sum == 1 then nothing has to be done.
    if ( FP_ZERO != std::fpclassify(sum-1)){

        CLOG(DEBUG, logName_) << "shares are not balanced correctly";
        CLOG(DEBUG, logName_) << "have to apply factor " << factor;

        sum = 0;
        for(auto it= components.begin(); it != components.end(); ++it){
            sum += std::get<1>(*it) * factor;
            std::get<1>(*it)=sum;
        }
    }
    
    CLOG(TRACE, logName_) << "finished balancing element share successfully";
}

//------------------------------------------------------------------------------

//! Checks if given element is either part of the anions or cations of the
//! material.
//! \param elementID Identification number of search element.
//! \return True if element is either in cations or anions, otherwise False.
bool Material::hasElement(const elementType elementID) const
{

    if ( hasAnion(elementID) || hasCation(elementID))
	return true;

    return false;
}

//------------------------------------------------------------------------------

//! Check if given element is part of the cations of the material.
//! \param elementID Element to check.
//! \return True if in cations, False otherwise.
bool Material::hasCation(const elementType elementID) const
{
    for(auto it= cations_.begin(); it != cations_.end(); ++it)
        if (std::get<0>(*it) == elementID) return true;

    return false;
}


//------------------------------------------------------------------------------

//! Check if given element is part of the anions of the material.
//! \param elementID Element to check.
//! \return True if in anions, False otherwise.
bool Material::hasAnion(const elementType elementID) const
{
    for(auto it= anions_.begin(); it != anions_.end(); ++it)
        if (std::get<0>(*it) == elementID) return true;

    return false;
}

//------------------------------------------------------------------------------

//! Check if element is the only cation.
bool Material::isOnlyCation(const elementType elementID) const {

    if ((cations_.size() == 1) && (std::get<0>(cations_[0]) == elementID))
        return true;

    return false;
}

//------------------------------------------------------------------------------

//! Check if element is the only anion.
bool Material::isOnlyAnion(const elementType elementID) const {

    if ((anions_.size() == 1) && (std::get<0>(anions_[0]) == elementID))
        return true;

    return false;
}

//------------------------------------------------------------------------------

//! \return Element ID which was chosen.
elementType Material::getRandomCation(std::mt19937 &rng) const
{   
    return getRandomEntry(cations_, rng);
}

//------------------------------------------------------------------------------

//! \return Element ID which was chosen.
elementType Material::getRandomAnion(std::mt19937 &rng) const
{
    return getRandomEntry(anions_, rng);
}

//------------------------------------------------------------------------------

//! The choice is based on the share of the single elements, meaning
//! that ones with a higher share are more probable picked than those with a
//! lower share.
elementType Material::getRandomEntry(const MaterialComponents &components,
				     std::mt19937 &rng) const
{
    
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double choice;

    CLOG(TRACE, logName_) << "starting to get a random entry";

    choice = distribution(rng);
    for(auto it= components.begin(); it != components.end(); ++it){
        if (choice <= std::get<1>(*it)) return std::get<0>(*it);
    }

    CLOG(WARNING, logName_) << "vector of Material object not filled correctly";
    CLOG(WARNING, logName_) << "falling back to first material";
    
    // If vector not filled correctly return first element.
    return std::get<0>(components[0]);

}

//------------------------------------------------------------------------------

//! \return Lattice constant of material.
double Material::getLatticeConstant(void) const
{
    return latticeConstant_;
}

//------------------------------------------------------------------------------

//! \return Lattice constant of material.
double Material::getLatticeConstant(double inPlaneLatticeConstant) const
{
    double result = 0.0;

    // if elastic parameters are not set, simply interpolate
    if ((c11_ > 0.0) && (c12_ > 0.0))
        result = (inPlaneLatticeConstant +
	    (latticeConstant_ -inPlaneLatticeConstant) *(1+2*c12_/c11_));
    else
        result = (inPlaneLatticeConstant + 
                (latticeConstant_ -inPlaneLatticeConstant) * 0.5);

    return result;
}

//------------------------------------------------------------------------------

//! \return True if material name is equal to provided name, false otherwise.
bool Material::operator==(const std::string &name) const
{
    return (name_ == name);
}

//------------------------------------------------------------------------------

//! \return Name of material.
std::string Material::getName(void) const
{
    return name_;
}

//------------------------------------------------------------------------------

//! \return [c11,c12,c44]
std::vector<double> Material::getElasticConstants(void) const
{
    return {c11_, c12_, c44_};
}

//------------------------------------------------------------------------------

//! \return ID and share of cations in material.
const MaterialComponents & Material::getCations() const
{
    return cations_;
}

//------------------------------------------------------------------------------

//! \return ID and share of anions in material.
const MaterialComponents & Material::getAnions() const
{
    return anions_;
}

//------------------------------------------------------------------------------

//! \param verbose If true besides the name also the cations and anions of the
//! material are included.
//! \return Description of the material.
std::string Material::str(bool verbose) const
{
    std::stringstream result;
    result << name_;

    if (verbose == true) {
        result << " (";

        for (auto cation : cations_) {
            result << " " << std::get<0>(cation);
        }

        result << " |";

        for (auto anion : anions_) {
            result << " " << std::get<0>(anion);
        }

        result << " )";

        result << " c11: " << c11_ << " c12: " << c12_;
    }

    return result.str();
}

//##############################################################################

MaterialException::MaterialException(const std::string &message,
				     const MaterialException::Id id)
    :message_(message), id_(id) {}

//------------------------------------------------------------------------------

//! Returns a short description of the exception. This is the implementation of
//! the virtual function inherited from std::exception
//! \remark This function does not throw exceptions!
const char * MaterialException::what () const throw ()
{
    static std::string text;
    text = "Material Exception -- " + message_;
    return text.c_str();
}

//------------------------------------------------------------------------------

//! Returns a more detailled description of the exception, making it possible to
//! use a single exception class for multiple errors.
//! \remark This function does not throw exceptions!
const std::string MaterialException::getMessage () const noexcept
{
    return message_;
}

//------------------------------------------------------------------------------

//! Returns the ID of the exception source, making possible to easily
//! distinguish different errors.
//! \remark This function does not throw exceptions!
const MaterialException::Id MaterialException::getId () const noexcept
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
