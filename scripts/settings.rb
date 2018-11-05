#!/bin/ruby
###
# Copyright (C) 2018 Hermann Detz and Juergen Maier
# 
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.
###


#-------------------------------------------------------------------------------
# Paths
#-------------------------------------------------------------------------------

$xtalsimBinDir = "./bin"
$xtalsimInputDir = "./input"

#-------------------------------------------------------------------------------
# Simulation settings
#-------------------------------------------------------------------------------

$maxDisplacement = 1e-3
$mmcRuns = 50 # MMC runs during interface generation (runs per atom within +- unit cells)
$interfacePositions = 10
$neighborsCutoff = 1
$temperature = 0

#-------------------------------------------------------------------------------


# Local variables:
# mode: ruby
# indent-tabs-mode: nil
# tab-width: 4
# End:
# vim:noexpandtab:sw=4:ts=4:
