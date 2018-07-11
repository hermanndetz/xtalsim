#!/usr/bin/ruby
###
# Copyright (C) 2018 Hermann Detz and Juergen Maier
# 
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.
###


require 'date'
require 'singleton'

#-------------------------------------------------------------------------------
# This class manages projecct related directories and settings.
#-------------------------------------------------------------------------------
class XtalSim
    include Singleton

    # task switches
    @mTasks = {}
    @mFiles = {}

    @mAddTimeStamp = true

    @mBaseDir = ""
    @mProjectDir = nil
    @mProjectName = ""
    @mTimestamp = DateTime.new
    @mUserSetFilename = false
    @mVerbose = false
    @mQuiet = false
    @mLogging = true
    @mStructureFile = ""

    def initialize
        @mTimeStamp = DateTime.new
        @mUserSetFilename = false

        @mTasks = [
            {:id => :lattice, :name => 'Lattice Generator', :enabled => false},
            {:id => :interface, :name => 'Interface Optimizer', :enabled => false},
            {:id => :optimize, :name => 'Optimizer', :enabled => false},
            {:id => :render, :name => 'Visualizer', :enabled => false}
        ]

        @mFiles = [
            {:id => :structure, :path => ''},
            {:id => :settings, :path => ''}
        ]

        @mLogging = true
    end

    def enable_all
        @mTasks.each do |t|
            t[:enabled] = true
        end
    end

    def enable_task (id)
        @mTasks.detect {|t| t[:id] == id}[:enabled] = true
    end

    def task_enabled? (id)
        return @mTasks.detect {|t| t[:id] == id}[:enabled]
    end

    def list_enabled_tasks
        @mTasks.each do |t|
            if t[:enabled] == true
                puts "#{t[:name]} enabled"
            end
        end
    end

    def add_timestamp! (enabled)
        @mAddTimeStamp = enabled
    end

    def add_timestamp?
        return @mAddTimeStamp
    end

    def set_basedir! (dir)
        @mBaseDir = dir
    end

    def get_basedir
        return @mBaseDir
    end

    def set_project_name! (name)
        @mProjectName = name
    end

    def get_project_name
        return @mProjectName
    end

    def get_project_dir
        if @mProjectDir == nil
            timeString = ""

            if (add_timestamp? == true)
                timeString = @mTimeStamp.strftime("%Y%m%d-%H%M")
            end
            
            @mProjectDir = "#{@mBaseDir}#{timeString}#{@mProjectName}"
        end

        return @mProjectDir
    end

    def get_project_file
        dir = get_project_dir
        
        tmpFile = ""
        
        if @mUserSetFilename == false
            tmpFile = "#{dir}/#{@mProjectName}"
        else
            tmpFile = @mProjectFile
        end

        return tmpFile
    end

    def set_project_file (file)
       @mProjectFile = File.join(File.dirname(file),File.basename(file, ".xml"))
       @mUserSetFilename = true
       @mProjectDir = File.dirname(file)
    end

    def prepare_work_dir
        get_project_dir

        if Dir.exists?(@mProjectDir) == false
            FileUtils.mkdir(@mProjectDir)
        end

        FileUtils.rm_f(get_project_file + ".xml")
    end

    def user_set_file?
        return @mUserSetFilename
    end

    def verbose?
        return @mVerbose
    end

    def set_verbose
        @mVerbose = true
    end

    def quiet?
        return @mQuiet
    end

    def set_quiet
        @mQuiet = true
    end

    def set_logging (enabled)
        @mLogging = enabled
    end

    def logging?
        return @mLogging
    end

    def set_file (file,path)
        @mFiles.detect {|f| f[:id] == file}[:path] = path
    end

    def get_file (file)
        return @mFiles.detect {|f| f[:id] == file}[:path]
    end
end

#-------------------------------------------------------------------------------


# Local variables:
# mode: ruby
# indent-tabs-mode: nil
# tab-width: 4
# End:
# vim:noexpandtab:sw=4:ts=4:
