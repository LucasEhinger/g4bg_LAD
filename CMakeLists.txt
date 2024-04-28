#----------------------------------------------------------------------------

# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(G4RC_LAD)

enable_language (Fortran)

#ROOT path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
find_package(ROOT REQUIRED)

# this is to ensure we find the header we need
include_directories(${CMAKE_CURRENT_BINARY_DIR}/../include)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})


#debugging symbols, complains
set(CMAKE_CXX_FLAGS "-ggdb -Wall -O0")
set(CMAKE_Fortran_FLAGS "-C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132")

#----------------------------------------------------------------------------
# Generate dictionaries as well (so we can write objects to ROOT files

include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/include ${PROJECT_SOURCE_DIR}/src/)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/g4rc*.cc ${PROJECT_SOURCE_DIR}/src/model/*.f)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Set CERN libraries
# set(cernlibs -L/share/apps/cernlib/x86_64-slc5/2006b/lib -lmathlib -lpacklib -lc -lm)
set(cernlibs -L/apps/cernlib/x86_64_rhel5/2005/lib/ -lmathlib -lpacklib -lc -lm)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
#add_executable(g4rc g4rc.cc ${sources} ${headers} g4rcDict.cxx)
add_executable(g4rc g4rc.cc ${sources} ${headers})
#add_library(g4rcroot SHARED ${libsources} ${libheaders} g4rcDict.cxx)
target_link_libraries(g4rc ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ${cernlibs})
#target_link_libraries(g4rcroot ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build g4rc. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

#FILE( GLOB G4RC_GEO "geometry/*" )

set(G4RC_SCRIPTS
    macros/vis.mac
    vis/OIX.mac
    vis/Qt.mac
    vis/vis.mac
    vis/X.mac
    macros/runexample.mac
    macros/gui.mac
    macros/vis.mac
    macros/batchrun.mac
    macros/runscript.sh
    macros/handler.sh
	macros/merger.sh
	cdbn.qwave
)

foreach(_script ${G4RC_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(G4RC DEPENDS g4rc)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS g4rc DESTINATION bin)


