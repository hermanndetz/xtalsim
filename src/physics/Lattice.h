/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __LATTICE_H__
#define __LATTICE_H__

#include <cmath>
#include <vector>
#include <iterator>
#include <iostream>

#include <stdio.h>
#include <stdint.h>

#include <easyloggingcpp/easylogging++.h>

#include <physics/Atom.h>
#include <physics/Material.h>
#include <physics/Range3D.h>
#include <physics/Vector3D.h>

//! Manages a three dimensional lattice structure.

//! The lattice is managed such that in all circumstances it is a cuboid. If
//! other structures are required later this has to be extended.

class Lattice{
private:

    //! 3D array containing atoms.
    std::vector< std::vector < std::vector<Atom *>>> lattice_;
    //! Number of layers in each direction.
    Vector3D<indexType> size_;
    //! Lattice temperature in Kelvin.
    double temperature_;

    //! random number generator
    std::mt19937 random_;
    
    //! Logger name.
    const std::string logName_;

    //! Extends lattice in one specific dimension.
    void extend(indexType amount, uint8_t dimension);

    //! Checks if value is within [0,limit-1] and adapts value if needed.
    inline void checkOverflow(indexType &value, const indexType limit) const;
    
public:

    enum Layer{
        CationRegular=0,
        CationRegularSwitched,
        CationShifted,
        CationShiftedSwitched,
        AnionRegular,
        AnionRegularSwitched,
        AnionShifted,
        AnionShiftedSwitched
    } ;

    enum class Type{
        Zincblende,
        HalfHeusler,
        FullHeusler
    };
    
    //! Constructor
    Lattice(const double temperature =0,
	    const char *logName="Lattice");
    //! Destructor
    ~Lattice();

    //! Create or extend lattice in specified size.
    void generate(const indexType x, const indexType y, const indexType z,
		  const uint8_t dimension=3);
    //! Create or extend lattice in specified size.
    void generate(const Vector3D<indexType> &size, const uint8_t dimension=3);
    //! Introduce Zincblende Lattice structure in a single layer.
    void createZincblendeLayer(const indexType i,
                               Vector3D<spaceType> position,
                               const Material *material,
                               const Lattice::Layer layerType,
                               uint8_t outOfPlaneDimension,
                               double inPlaneLatticeConstant);
    //! Introduce Rocksalt Lattice structure in a single layer.    
    void createRocksaltLayer(const indexType i,
                             Vector3D<spaceType> position,
                             const Material *material,
                             const Lattice::Layer layerType,
                             uint8_t outOfPlaneDimension,
                             double inPlaneLatticeConstant);
    
    //! Return arbitrary atom around (0,0) in specified layer.
    const Atom* getFirstAtomInLayer(const indexType layerNumber,
				    const uint8_t outOfPlaneDimension) const;
    
    //! Returns single element of lattice.
    Atom* & operator()(indexType x, indexType y, indexType z);
    //! Returns single element of lattice.
    Atom* & operator()(Vector3D<indexType> position);

    //! Returns single element of lattice.
    const Atom * operator()(indexType x, indexType y, indexType z) const;
    //! Returns single element of lattice.
    const Atom * operator()(Vector3D<indexType> position) const;

    //! Returns size of lattice.
    Vector3D<indexType> getSize(void) const;
    //! Returns maximum coordinate of atoms in each spatial direction.
    Vector3D<spaceType> getMaxCoordinates(void) const;
    //! Returns lattice temperature.
    double getTemperature(void) const;
    //! Sets lattice temperature.
    void setTemperature(const double temperature);
    //! Get list of existing atoms inside specified range.
    std::vector<Atom *> getAtomList(
			       Range3D<indexType> range = Range3D<indexType>(),
			       const bool checkModification=false,
			       const std::string &materialName="") const;

    //! Get the list of atoms in a given layer.
    std::vector<Atom *> getAtomsInLayer(
                    indexType layerID, indexType dimension) const;
    
    //! Count atoms according to their state
    uint32_t countAtomsByState(AtomState metric);

    //! Count atoms within a layer according to their state
    uint32_t countAtomsInLayerByState(indexType layerID, indexType dimension,
                    AtomState metric);

    //! Backup range of atoms.
    std::vector<Atom> backupAtoms(const Range3D<indexType> range=
				    Range3D<indexType>()) const;

    //! Restore atoms.
    void restoreAtoms(const std::vector<Atom> atoms);

    //! Count atoms (not lattice sites)
    uint16_t getAtomCount(Range3D<indexType> range =Range3D<indexType>()) const;

    //! Get average coordinate
    double getAvgCoordinate (indexType layer, uint8_t dimension, bool checkModification=false) const;

    //! Returns the maximum modification index of all atoms in object
    uint32_t getMaxModificationIndex(void) const;
};

//##############################################################################

//! \brief Lattice exception

class LatticeException : public std::exception {
   
 public:

    //! Identifies source of the Exception.
    enum class Id{
    	Exists,
    	NotExists,
    	SizeZero,
	WrongDimension,
	PositionTooBig,
	RangeTooBig,
    PositionNAN,
    	Unknown
    } ;
    
    //! Constructor
    LatticeException(const std::string &message="",
	    const LatticeException::Id id = LatticeException::Id::Unknown ); 

    //! Return description of execution.
    const char * what () const throw ();
    //! Return more detailled error message.
    const std::string getMessage() const noexcept;
    //! Return Exception identification number.
    const LatticeException::Id getId() const noexcept;
    
 private:

    std::string message_; //!< error message
    LatticeException::Id id_; //!< unique ID identifying source
    
};


//##############################################################################



#endif



// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
