cmake_minimum_required( VERSION 2.8 )

project( SPhyR )

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH})

set( python_mod_src
  src/matrix.cpp
  src/utils.cpp
  src/columngen.cpp
  src/python.cpp
)

set( python_mod_hdr
  src/matrix.h
  src/utils.h
  src/columngen.h
)

set (kDP_src
  src/kdpmain.cpp
  src/matrix.cpp
  src/utils.cpp
  src/columngen.cpp
)

set (kDP_hdr
  src/matrix.h
  src/utils.h
  src/columngen.h
)

set (kDPFC_src
  src/kdpfcmain.cpp
  src/matrix.cpp
  src/utils.cpp
  src/coordinateascent.cpp
  src/columngenflipclustered.cpp
  src/columngenflip.cpp
  src/columngen.cpp
  src/cluster.cpp
)

set (kDPFC_hdr
  src/matrix.h
  src/utils.h
  src/coordinateascent.h
  src/cluster.h
  src/dkm/dkm.hpp
  src/dkm/dkm_utils.hpp
  src/columngenflipclustered.h
  src/columngenflip.h
  src/columngen.h
)

set( analyze_src
  src/analyzemain.cpp
  src/utils.cpp
  src/matrix.cpp
  src/comparison.cpp
  src/phylogenetictree.cpp
  src/dollophylogenetictree.cpp
)

set( analyze_hdr
  src/utils.h
  src/matrix.h
  src/comparison.h
  src/phylogenetictree.h
  src/dollophylogenetictree.h
)

set( perturb_src
  src/perturbmain.cpp
  src/utils.cpp
  src/matrix.cpp
)

set( perturb_hdr
  src/utils.h
  src/matrix.h
)

set( visualize_src
  src/visualizemain.cpp
  src/utils.cpp
  src/matrix.cpp
  src/phylogenetictree.cpp
  src/dollophylogenetictree.cpp
)

set( visualize_hdr
  src/utils.h
  src/matrix.h
  src/phylogenetictree.h
  src/dollophylogenetictree.h
)

set( simulate_hdr
  src/utils.h
  src/matrix.h
  src/phylogenetictree.h
  src/dollophylogenetictree.h
)

set( simulate_src
  src/simulatemain.cpp
  src/utils.cpp
  src/matrix.cpp
  src/phylogenetictree.cpp
  src/dollophylogenetictree.cpp
)

# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

if (CMAKE_BUILD_TYPE MATCHES "^[Rr]elease")
  # check if Doxygen is installed
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
      # set input and output files
      set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
      set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

      # request to configure the file
      configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
      message("Doxygen build started")

      # note the option ALL which allows to build the docs together with the application
      add_custom_target( doc_doxygen ALL
          COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          COMMENT "Generating API documentation with Doxygen"
          VERBATIM )
  else (DOXYGEN_FOUND) 
    message("Doxygen need to be installed to generate the doxygen documentation")
  endif (DOXYGEN_FOUND)
endif()

if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR)
  set( CPLEX_AUTODETECTION_PREFIXES "$ENV{HOME}/ILOG;/opt/ibm/ILOG;$ENV{HOME}/Applications/IBM/ILOG" CACHE STRING "CPLEX automatic lookup directories" )
  foreach( CplexInstallPrefix ${CPLEX_AUTODETECTION_PREFIXES} )
    file( GLOB_RECURSE CplexLibraryPaths ${CplexInstallPrefix}/*/cplex/lib/*/libcplex.a )
    if( CplexLibraryPaths )
      list( LENGTH CplexLibraryPaths NumPath )
      if( NumPath GREATER 1 )
        message( WARNING "Found multiple CPLEX installations in ${CplexInstallPrefix}, automatically selecting the first one, please set {CPLEX|CONCERT}_{INC|LIB}_DIR manually if this is not the right selection" )
      else()
        message( "Found CPLEX's installation in ${CPlexInstallPrefix}" )
      endif()
      list( GET CplexLibraryPaths 0 CplexLibraryPath )

      set( MyRegex "${CplexInstallPrefix}/(.*)/cplex/?/lib/(.*)/libcplex.a" )
      string( REGEX MATCH "${MyRegex}" ValidatedCplexLibPath ${CplexLibraryPath} )
      if( NOT ValidatedCplexLibPath )
        message( WARNING "Unable to recognize ${CplexLibraryPath} as a valid CPLEX path" )
      else()
        string( REGEX REPLACE "${MyRegex}" "\\1" CplexHome ${ValidatedCplexLibPath} )
        string( REGEX REPLACE "${MyRegex}" "\\2" CplexArch ${ValidatedCplexLibPath} )

        set( CPLEX_INC_DIR   "${CplexInstallPrefix}/${CplexHome}/cplex/include"            CACHE PATH "CPLEX's include dir" )
        set( CPLEX_LIB_DIR   "${CplexInstallPrefix}/${CplexHome}/cplex/lib/${CplexArch}"   CACHE PATH "CPLEX's library dir" )
        set( CONCERT_INC_DIR "${CplexInstallPrefix}/${CplexHome}/concert/include"          CACHE PATH "Concert's include dir" )
        set( CONCERT_LIB_DIR "${CplexInstallPrefix}/${CplexHome}/concert/lib/${CplexArch}" CACHE PATH "Concert's library dir" )

        break()
      endif()
    else()
      message( "Unable to find CPLEX and Concert installation in ${CplexInstallPrefix}" )
    endif()
  endforeach()

  if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR )
    message( "Unable to autodeted CPLEX's and Concert's installation directories, please set {CPLEX|CONCERT}_{INC|LIB}_DIR manually" )
  else()
    message( "Autodetected CPLEX's include dir   : ${CPLEX_INC_DIR}" )
    message( "Autodetected CPLEX's library dir   : ${CPLEX_LIB_DIR}" )
    message( "Autodetected Concert's include dir : ${CONCERT_INC_DIR}" )
    message( "Autodetected Concert's library dir : ${CONCERT_LIB_DIR}" )
  endif()
endif()

set( Boost_USE_STATIC_LIBS ON )
find_package( Boost REQUIRED COMPONENTS python2.7 )
MESSAGE( STATUS "found boost:"
"\nINCLUDE: ${Boost_INCLUDE_DIRS}"
"\nLIB: ${Boost_LIBRARIES}"
)

FIND_PACKAGE(PythonInterp REQUIRED)
FIND_PACKAGE(PythonLibsNew REQUIRED)

# custom macro with most of the redundant code for making a python example module
MACRO( addPythonExe _name _srccpp )
ADD_EXECUTABLE(${_name} ${_srccpp} ${python_mod_src} ${python_mod_hdr} )
# make the pyd library link against boost_numpy python and boost
TARGET_INCLUDE_DIRECTORIES(${_name} PUBLIC "${PROJECT_SOURCE_DIR}/src" ${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(${_name} ${CommonLibs} boost_python2.7)
ENDMACRO()

MACRO( addPythonMod _name _srccpp )
PYTHON_ADD_MODULE(${_name} ${_srccpp} ${python_mod_src} ${python_mod_hdr} )
# make the pyd library link against boost_numpy python and boost
TARGET_INCLUDE_DIRECTORIES(${_name} PUBLIC "${PROJECT_SOURCE_DIR}/src" ${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(${_name} ${CommonLibs})
ENDMACRO()

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")
add_definitions( -DIL_STD )

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wextra -Wno-long-long -Wno-unused-parameter -Wno-unknown-pragmas -g -ggdb")

include_directories( ${CPLEX_INC_DIR} ${CONCERT_INC_DIR} "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} )
link_directories( ${CPLEX_LIB_DIR} ${CONCERT_LIB_DIR} "${LIBLEMON_ROOT}/lib" )

set( CommonLibs
  m
  ilocplex
  cplex
  concert
  pthread
  emon
  dl
  ${Boost_LIBRARIES}
)

add_executable( kDP ${kDP_src} ${kDP_hdr} )
target_link_libraries( kDP ${CommonLibs} )

add_executable( kDPFC ${kDPFC_src} ${kDPFC_hdr} )
target_link_libraries( kDPFC ${CommonLibs} )

add_executable( visualize ${visualize_src} ${visualize_hdr} )
target_link_libraries( visualize ${CommonLibs} )

add_executable( simulate ${simulate_src} ${simulate_hdr} )
target_link_libraries( simulate ${CommonLibs} )

add_executable( analyze ${analyze_src} ${analyze_hdr} )
target_link_libraries( analyze ${CommonLibs} )

add_executable( perturb ${perturb_src} ${perturb_hdr} )
target_link_libraries( perturb ${CommonLibs} )

addPythonMod( sphyr src/python.cpp )
