#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.16...3.21)
project(B4a)

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
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exampleB4a exampleB4a.cc ${sources} ${headers})
target_link_libraries(exampleB4a ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B4a. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEB4A_SCRIPTS
    exampleB4a.out
    exampleB4.in
    gui.mac
    init_vis.mac
    plotHisto.C
    plotNtuple.C
    run0.mac
    Run_C4_2_0.mac
    Run_C4_2_1.mac
    Run_C4_2_2.mac
    Run_C4_2_3.mac
    Run_C4_2_4.mac
    Run_C4_2_5.mac
    Run_C4_2_6.mac
    Run_C4_2_7.mac
    Run_C4_2_8.mac
    Run_C4_2_9.mac
    Run_C4_2_10.mac
    Run_C4_2_11.mac
    Run_C4_2_12.mac
    Run_C4_2_13.mac
    Run_C4_2_14.mac
    Run_C4_2_15.mac
    Run_C4_2_16.mac
    Run_C4_2_17.mac
    Run_C4_2_18.mac
    Run_C4_2_19.mac
    Run_C4_2_20.mac
    Run_PC3_PIP_0.mac
    Run_PC3_PIP_1.mac
    Run_PC3_PIP_2.mac
    Run_PC3_PIP_3.mac
    Run_PC3_PIP_4.mac
    Run_PC3_PIP_5.mac
    Run_PC3_PIP_6.mac
    Run_PC3_PIP_7.mac
    Run_PC3_PIP_8.mac
    Run_PC3_PIP_9.mac
    Run_PC3_PIP_10.mac
    Run_PC3_PIP_11.mac
    Run_PC3_PIP_12.mac
    Run_PC3_PIP_13.mac
    Run_PC3_PIP_14.mac
    Run_PC3_PIP_15.mac
    Run_PC3_PIP_16.mac
    Run_PC3_PIP_17.mac
    Run_PC3_Flu_0.mac
    Run_PC3_Flu_1.mac
    Run_PC3_Flu_2.mac
    Run_PC3_Flu_3.mac
    Run_PC3_Flu_4.mac
    Run_PC3_Flu_5.mac
    vis.mac
    jobScript_C4_2.sh
    jobScript_PC3_PIP.sh
    jobScript_PC3_Flu.sh
    jobScript0.sh
  )

foreach(_script ${EXAMPLEB4A_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS exampleB4a DESTINATION bin)
