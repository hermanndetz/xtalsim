/**
Copyright (C) 2018 Hermann Detz and Juergen Maier

This software may be modified and distributed under the terms
of the MIT license.  See the LICENSE file for details.
*/


#ifndef __ATOMSTATE_H__
#define __ATOMSTATE_H__

#include <bitset>

#define AtomState_Default 0

enum AtomState {
    // Default = 0
    ModifiedInterface, // = 0x01
    ModifiedExchangeReaction, // = 0x02
    ModifiedUnknown, // = 0x04

    StateCount // dummy item - must be last entry
};


#endif

// Local variables:
// mode: c++
// indent-tabs-mode: nil
// tab-width: 4
// End:
// vim:noexpandtab:sw=4:ts=4:
