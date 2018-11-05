#!/usr/bin/ruby
###
# Copyright (C) 2018 Hermann Detz and Juergen Maier
# 
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.
###


require 'fileutils'
require 'optparse'
require './scripts/ioutils.rb'
require './scripts/xtalsim.rb'
require './scripts/xtalsimcli.rb'
require './scripts/xtalsimconfig.rb'

#-------------------------------------------------------------------------------
# Prepare core objects
#-------------------------------------------------------------------------------

xtalsim = XtalSim.instance # stores all settings
XtalSimCLI.instance # parses command line arguments
ConfigFiles.instance.set_filename xtalsim.get_project_file

#-------------------------------------------------------------------------------
# Prepare input files
#-------------------------------------------------------------------------------

fileName = xtalsim.get_project_file # .xml is added later

if xtalsim.task_enabled? :lattice
    if xtalsim.get_file(:structure).to_s.empty? == false
        require xtalsim.get_file :structure
    else
        puts "No structure file given!".red.bold
        abort
    end
end

if xtalsim.get_file(:settings).to_s.empty? == false
    require xtalsim.get_file :settings

    # update path in ConfigFiles
    ConfigFiles.instance.set_base_path($xtalsimInputDir)
end

if xtalsim.user_set_file? == false
    xtalsim.prepare_work_dir
end

#-------------------------------------------------------------------------------
# Step through tasks
#-------------------------------------------------------------------------------
puts "Starting project #{xtalsim.get_project_name} using file #{fileName}"

$inFile = ""

if xtalsim.task_enabled? :lattice
    puts "Generating lattice".blue.bold

    $layers.each_with_index do |l,index|

        if index > 0
            $inFile = "-i #{fileName}.xml "
        end

        tmpSize = [$size[0]*4, $size[1]*4, l[:thickness]*4]
        tmpName = l[:cation].to_s + l[:anion].to_s

        puts "Layer #{l[:name]}".blue

        cmd = "#{$xtalsimBinDir}/LatticeGenerator"
        cmd << " #{$inFile}"
        cmd << " -o #{fileName}.xml"
        cmd << " -n #{tmpName}"
        cmd << " -c #{l[:cation].to_s}"
        cmd << " -a #{l[:anion].to_s}"
        cmd << " -s #{tmpSize[0]} #{tmpSize[1]} #{tmpSize[2]}"
        cmd << " -l #{$latticeConstant}"
        cmd << " -d #{$hsdimension}"
        #cmd << " -E"
        #cmd << " --neighbor-layers #{$neighborsCutoff}"
        cmd << " #{ConfigFiles.instance.get_file_parameters("periodic","materials","tersoff")}"

        if xtalsim.logging?
            cmd << " #{ConfigFiles.instance.get_file_parameters("log")}"
        else
            cmd << " --qq"
        end

        if xtalsim.quiet?
            cmd << " -q"
        else
            if xtalsim.verbose?
                cmd << " -v"
            end
        end

        if XtalSimCLI.instance.show_commands?
            puts "#{cmd}"
        end

        %x(#{cmd})
    end
end

#-------------------------------------------------------------------------------

if xtalsim.task_enabled? :interface
    $interfaces.each_with_index do |interface,index|
        puts "Optimizing interface at layer #{interface[:layer]} with height #{interface[:height]} using #{interface[:cation]} and #{interface[:anion]}".blue

        # start index and stop index
        cmd = "#{$xtalsimBinDir}/InterfaceGenerator"
        cmd << " -i #{fileName}.xml"
        cmd << " -o #{fileName}.xml"
        cmd << " -c #{interface[:cation]}"
        cmd << " -a #{interface[:anion]}"
        cmd << " -n #{interface[:material]}"
        cmd << " --interface-atoms #{interface[:atoms]}"
        cmd << " -l #{$latticeConstant}"
        cmd << " --start-index #{interface[:layer]}"
        cmd << " --stop-index #{interface[:layer]+interface[:height]}"
        cmd << " --neighbor-layers 1"
        cmd << " --exchange-reaction #{interface[:exchange]}"
        cmd << " --metric #{interface[:metric]}"
        cmd << " #{ConfigFiles.instance.get_file_parameters("config","periodic","materials","tersoff","journal")}"

        if xtalsim.logging?
            cmd << " #{ConfigFiles.instance.get_file_parameters("log")}"
        else
            cmd << " --qq"
        end

        if xtalsim.quiet?
            cmd << " -q"
        else
            if xtalsim.verbose?
                cmd << " -v"
            end
        end

        if XtalSimCLI.instance.show_commands?
            puts "#{cmd}"
        end

        result = 0
        result = %x(#{cmd})
        puts "Result: #{result}" # this is a work around to wait for completion

    end

end

#-------------------------------------------------------------------------------

if xtalsim.task_enabled? :optimize
    cmd = "#{$xtalsimBinDir}/Optimizer"
    cmd << " -i #{fileName}.xml"
    cmd << " -o #{fileName}.xml"
    cmd << " --xyz #{xtalsim.get_project_dir}"
    cmd << " #{ConfigFiles.instance.get_file_parameters("config","periodic","tersoff","journal")}"

    if xtalsim.logging?
        cmd << " #{ConfigFiles.instance.get_file_parameters("log")}"
    else
        cmd << " --qq"
    end

    if xtalsim.quiet?
        cmd << " -q"
    else
        if xtalsim.verbose?
            cmd << " -v"
        end
    end

    if XtalSimCLI.instance.show_commands?
        puts "#{cmd}"
    end

    %x(#{cmd})
end

#-------------------------------------------------------------------------------

if xtalsim.task_enabled? :render
    if File.exists?("../bin//Visualize")
        fileParams = ConfigFiles.instance.get_file_parameters("config","periodic","materials")
        puts "File params: #{fileParams}"

        cmd = "../bin/Visualize"
        cmd << " -i #{fileName}.xml"
        cmd << " #{fileParams}"
        cmd << " --neighbor-layers 1"

        if xtalsim.quiet?
            cmd << " -q"
        else
            if xtalsim.verbose?
                cmd << " -v"
            end
        end

        #cmd << " --bonds"
        #cmd << " --strain"
        cmd << " --modified"

        if XtalSimCLI.instance.show_commands?
            puts "#{cmd}"
        end


        %x(#{cmd})
    end
end

puts "Done".blue.bold

#-------------------------------------------------------------------------------


# Local variables:
# mode: ruby
# indent-tabs-mode: nil
# tab-width: 4
# End:
# vim:noexpandtab:sw=4:ts=4:
