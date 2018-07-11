#!/usr/bin/ruby
###
# Copyright (C) 2018 Hermann Detz and Juergen Maier
# 
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.
###


require 'singleton'
require './xtalsimconfig.rb'

#-------------------------------------------------------------------------------
# This class parses command line options for the XtalSim Shell.
# Settings which only affect the shell but not XtalSim programs directly
# are also stored here.
#-------------------------------------------------------------------------------
class XtalSimCLI
    include Singleton

    @mShowCommands = false

    def initialize
        parse_arguments
    end

    def show_commands?
        return @mShowCommands
    end

    private def parse_arguments
        xtalsim = XtalSim.instance

        # Parse command line options
        OptionParser.new do |opt|
            opt.on('-a', '--all', "Perform all tasks") do
                xtalsim.enable_all
            end

            opt.on('-l', '--lattice', "Generate lattice") do |lat|
                xtalsim.enable_task(:lattice)
            end

            opt.on('--interface', "Optimize interface") do |int|
                xtalsim.enable_task(:interface)
            end

            opt.on('--optimize', "Optimize structure") do
                xtalsim.enable_task(:optimize)
            end

            opt.on('--render', "Render structure") do |int|
                xtalsim.enable_task(:render)
            end

            opt.on('-fFILE', '--file=FILE', "Path to .xml file") do |file|
                xtalsim.set_project_file file
            end

            opt.on('--[no-]timestamp', "Add timestamp to project directories and files") do |enabled|
                xtalsim.add_timestamp! enabled
            end

            opt.on('-bDIR', '--base-dir=DIR', "Base directory for all project directories and files") do |dir|
                xtalsim.set_basedir! dir
            end

            opt.on('-pNAME', '--project-name=NAME', "Set project name") do |name|
                xtalsim.set_project_name! name
            end

            opt.on('-v', '--verbose', "Detailed output. -q or --quiet overrides -v or --verbose!") do
                xtalsim.set_verbose
            end

            opt.on('-q', '--quiet', "No output to terminal. -q or --quiet overrides -v or --verbose!") do
                xtalsim.set_quiet
            end

            opt.on('--qq', "Output to terminal and file logging disabled.") do
                xtalsim.set_quiet
                xtalsim.set_logging false
            end

            opt.on('--disable-logging', "Disable generation of log files.") do
                xtalsim.set_logging false
            end

            opt.on('--show-commands', "Print commands including parameter lists before they are executed") do
                @mShowCommands = true
            end

            opt.on('--structure-file=FILE', "Path to ruby script containing layer definitions") do |file|
                xtalsim.set_file :structure,file
            end

            opt.on('--settings-file=FILE', "Path to ruby script containing simulation settings") do |file|
                xtalsim.set_file :settings,file
            end

            opt.on('--user-config=FILE', "Path to XML file containing a configuration for XtalSim applications") do |file|
                ConfigFiles.instance.override_filename "config",file
            end
        end.parse!
    end
end

# Local variables:
# mode: ruby
# indent-tabs-mode: nil
# tab-width: 4
# End:
# vim:noexpandtab:sw=4:ts=4:
