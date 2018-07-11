/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/

#ifndef __VTKDEFS_H__
#define __VTKDEFS_H__

// this is necessary to allow build with VTK 6.3
// for VTK 7.1, VTK_OVERRIDE is defined in vtkConfigure.h
#ifndef VTK_OVERRIDE
    #define VTK_OVERRIDE    override
#endif

#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
