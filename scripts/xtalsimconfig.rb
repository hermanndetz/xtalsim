###
# Copyright (C) 2018 Hermann Detz and Juergen Maier
# 
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.
###

require 'singleton'

class ConfigFiles
    include Singleton

    @mFileName = ""
    @mBasePath = ""

    def initialize
        @mEntries = {
		"config" => {"param" => "--config", "path" => "config.xml", :type => :user},
            "periodic" => {"param" => "--periodic-table", "path" => "elements.xml", :type => :input},
            "tersoff" => {"param" => "--tersoff", "path" => "tersoff.xml", :type => :input},
            "materials" => {"param" => "--material-file", "path" => "materials.xml", :type => :input},
            "log" => {"param" => "--log", "path" => "#{@mFileName}.log", :type => :output},
            "journal" => {"param" => "-j", "path" => "#{@mFileName}", :type => :output}
        }
    end

    def set_filename (fileName)
        @mFileName = fileName
        @mEntries["log"]["path"] = "#{@mFileName}.log"
        @mEntries["journal"]["path"] = "#{@mFileName}"
    end

    def override_filename (id,filename)
        @mEntries[id]["path"] = filename
    end

    def get_file_parameters (*args)
        paramString = ""

        args.each do |item|
			if @mEntries[item][:type] == :input
				paramString << " #{@mEntries[item]["param"]} #{@mBasePath}/#{@mEntries[item]["path"]}"
			else
				paramString << " #{@mEntries[item]["param"]} #{@mEntries[item]["path"]}"
			end
        end

        return paramString
    end

    def set_base_path (path = "")
		@mBasePath = path
    end
end


# Local variables:
# mode: ruby
# indent-tabs-mode: nil
# tab-width: 4
# End:
# vim:noexpandtab:sw=4:ts=4:
