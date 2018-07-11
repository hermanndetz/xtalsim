/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __ERROR_CODE_H__
#define __ERROR_CODE_H__


namespace ErrorCode{
    enum {
	NoError=0,
	Unknown,
	Config,
	GrowthDimension,
	CationShare,
	AnionShare,
	LatticeConstant,
	Size,
	InOutFiles,
	PeriodicTable,
	Lattice,
	Tersoff,
	Material,
	Neighbor,
	OptimizationAction,
	Runs,
	Strain,
	StartStopIndex,
    Camera,
	TCLAP,
	XML,
	InterfaceAtoms,
	InterfacePositions,
	SigSegV,
	SigInt,
	SigUnknown
    };
}  

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
