###
# Copyright (C) 2018 Hermann Detz and Juergen Maier
# 
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.
###

include_directories(
  .
  app
  easyloggingcpp
  jsoncpp
  misc
  physics
  pugixml
  simulation
  tclap
  visualization
  ${PROJECT_BINARY_DIR}
  )

set (sourceFiles "")
set (srcFolders
  easyloggingcpp
  jsoncpp
  misc
  physics
  pugixml
  simulation
  visualization
  )
set (execFolders app)

# foreach (folder ${folders})

#   file (GLOB SRC_FILES ${folder}/*.cpp)
#   set (sourceFiles ${sourceFiles} ${SRC_FILES})

# endforeach()

foreach (folder ${srcFolders})
  include (${folder}/CMakeLists.txt)
  foreach (FILE ${FILES})
    set (sourceFiles ${sourceFiles} ${folder}/${FILE})
  endforeach()
endforeach()

add_library(XtalsimLib ${sourceFiles})

#message(${sourceFiles})
foreach (folder ${execFolders})
  include (${folder}/CMakeLists.txt)
  foreach (FILE ${FILES})
    add_executable(${FILE} ${folder}/${FILE}.cpp)

    if (VTK_FOUND)
      target_link_libraries(${FILE} m pthread XtalsimLib ${VTK_LIBRARIES})
    else()
      target_link_libraries(${FILE} m pthread XtalsimLib)
    endif()
    endforeach()
endforeach()


#install (TARGETS LatticeGenerator DESTINATION bin)
