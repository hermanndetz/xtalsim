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

    def initialize
        @mEntries = {
            "config" => {"param" => "--config", "path" => "../input/config.xml"},
            "periodic" => {"param" => "--periodic-table", "path" => "../input/elements.xml"},
            "tersoff" => {"param" => "--tersoff", "path" => "../input/tersoff.xml"},
            "materials" => {"param" => "--material-file", "path" => "../input/materials.xml"},
            "log" => {"param" => "--log", "path" => "#{@mFileName}.log"},
            "journal" => {"param" => "-j", "path" => "#{@mFileName}"}
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
            paramString << " #{@mEntries[item]["param"]} #{@mEntries[item]["path"]}"
        end

        return paramString
    end
end


# Local variables:
# mode: ruby
# indent-tabs-mode: nil
# tab-width: 4
# End:
# vim:noexpandtab:sw=4:ts=4:
