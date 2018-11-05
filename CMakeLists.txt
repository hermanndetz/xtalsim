###
# Copyright (C) 2018 Hermann Detz and Juergen Maier
# 
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.
###

cmake_minimum_required (VERSION 2.6)
project (Xtalsim CXX)

# Version number
set (Xtalsim_MAJOR_VERSION 0)
set (Xtalsim_MINOR_VERSION 1)

# seems to work only for VERSION > 3
SET(CMAKE_CXX_STANDARD 11)

SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin")
#SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" CACHE PATH "")
#mark_as_advanced(EXECUTABLE_OUTPUT_PATH)

SET(LIBRARY_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin")
#SET(LIBRARY_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" CACHE PATH "")
#mark_as_advanced(LIBRARY_OUTPUT_PATH)

SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

IF( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
  EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_SOURCE_DIR}/src/compile_commands.json
  )
ENDIF()

################################################################################
## DOXYGEN  ####################################################################
################################################################################

find_program(DOXYGEN_BIN doxygen)

message("-- doxygen")

if (NOT ${DOXYGEN_BIN} MATCHES DOXYGEN_BIN-NOTFOUND)
  # add doxygen documentation
  message("-- doxygen - Success")

  set (DOXYGEN_INPUT_DIR ${CMAKE_SOURCE_DIR}/src)
  set (DOXYGEN_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/doc)
  set (DOXYGEN_EXAMPLE_PATH
        "${CMAKE_SOURCE_DIR}/src/doc")
  configure_file(
    ${CMAKE_SOURCE_DIR}/src/doc/Doxyfile.in
    ${CMAKE_SOURCE_DIR}/src/doc/Doxyfile)
  add_custom_target ( doc
    ${DOXYGEN_BIN} ${PROJECT_BINARY_DIR}/Doxyfile
    DEPENDS ${PROJECT_BINARY_DIR}/Doxyfile
    )
else()
  message("-- doxygen - Fail")
endif()

################################################################################
##  USER CONFIGURATION  ########################################################
################################################################################
set (USERCONF "${CMAKE_SOURCE_DIR}/CMakeUser.txt")
  
message("-- checking for User configuration file")

if (EXISTS ${USERCONF})
  message("-- checking for User configuration file - Found")
  include (${USERCONF})
endif()

################################################################################
##  VTK  #######################################################################
################################################################################

find_package(VTK QUIET NO_MODULE)

message("-- checking for Visualization Toolkit (VTK)")
if (VTK_FOUND)
  message("-- checking for Visualization Toolkit (VTK) - Success")
  # from VTK manual
  #IF(NOT VTK_USE_RENDERING)
  #  MESSAGE(FATAL_ERROR
  #    "Example ${PROJECT_NAME} requires VTK_USE_RENDERING.")
  #ENDIF(NOT VTK_USE_RENDERING)
  set (__VTK__ True)
  #add_definitions(-D__VTK__)
  include (${VTK_USE_FILE})
else()
  set (__VTK__ False)
  message("-- checking for Visualization Toolkit (VTK) - Fail")
endif()

configure_file(
  ${CMAKE_SOURCE_DIR}/src/projectConfigure.h.in
  ${PROJECT_BINARY_DIR}/projectConfigure.h)


################################################################################
##  BUILD  #####################################################################
################################################################################

if (NOT ${CMAKE_VERSION} VERSION_GREATER 3.0)
   # check C++11 support
   include(CheckCXXCompilerFlag)
   CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
   CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
   if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
   elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
   else()
	message(FATAL_ERROR "Compiler ${CMAKE_CXX_COMPILER} has no C++11 support.")
   endif()
endif()

# compile flags
set (PEDANTIC OFF CACHE BOOL "activate pedantic mode and all warnings")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msseregparm -fstack-protector-all -DELPP_THREAD_SAFE")
if (NOT APPLE)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--no-as-needed")
endif()

if (PEDANTIC)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wfloat-equal")
endif()

set (CMAKE_CXX_FLAGS_DEBUG "-g")

add_subdirectory(src)
