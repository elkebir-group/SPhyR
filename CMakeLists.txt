cmake_minimum_required( VERSION 2.8 )

project( SPhyR )

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH})

set( Dollo_src
  src/dollomain.cpp
  src/matrix.cpp
  src/utils.cpp
  src/ilpsolverdollo.cpp
  src/ilpsolverdolloflip.cpp
  src/ilpsolverdolloflipcluster.cpp
  src/ilpsolverdolloflipclustered.cpp
  src/coordinateascent.cpp
  src/phylogenetictree.cpp
)

set( Dollo_hdr
  src/utils.h
  src/ilpsolverdollo.h
  src/ilpsolverdolloflip.h
  src/ilpsolverdolloflipcluster.h
  src/ilpsolverdolloflipclustered.h
  src/matrix.h
  src/dollocallback.h
  src/coordinateascent.h
  src/dkm/dkm.hpp
  src/dkm/dkm_utils.hpp
  src/phylogenetictree.h
)

set( analyze_src
  src/analyzemain.cpp
  src/utils.cpp
  src/matrix.cpp
)

set( analyze_hdr
  src/utils.h
  src/matrix.h
)

set( visualize_hdr
  src/utils.h
  src/matrix.h
  src/phylogenetictree.h
)

set( visualize_src
  src/visualizemain.cpp
  src/utils.cpp
  src/matrix.cpp
  src/phylogenetictree.cpp
)

set( visualize_hdr
  src/utils.h
  src/matrix.h
  src/phylogenetictree.h
)

set( simulate_hdr
  src/utils.h
  src/matrix.h
  src/phylogenetictree.h
)

set( simulate_src
  src/simulatemain.cpp
  src/utils.cpp
  src/matrix.cpp
  src/phylogenetictree.cpp
)

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
find_package( Boost REQUIRED )

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
)

add_executable( sphyr ${Dollo_src} ${Dollo_hdr} )
target_link_libraries( sphyr ${CommonLibs} )

add_executable( visualize ${visualize_src} ${visualize_hdr} )
target_link_libraries( visualize ${CommonLibs} )

add_executable( simulate ${simulate_src} ${simulate_hdr} )
target_link_libraries( simulate ${CommonLibs} )

add_executable( analyze ${analyze_src} ${analyze_hdr} )
target_link_libraries( analyze ${CommonLibs} )
