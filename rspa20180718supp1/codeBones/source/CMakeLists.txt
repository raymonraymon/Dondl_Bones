cmake_minimum_required(VERSION 3.0)

# Project name
project(CPPLIBSHAPEDESIGNBONES)

# C++11 is required to build this project
set(CMAKE_CXX_STANDARD 11)

# Position Independent Code (need for python)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

###########################################
#Boost
############################################ 
option( USE_BOOST "Boost" ON )
if ( USE_BOOST )
    set (Boost_USE_STATIC_LIBS OFF) #must be OFF if one wants to build shared libraries
    set (Boost_USE_MULTITHREADED ON)  
    set (Boost_USE_STATIC_RUNTIME OFF) 
    find_package (Boost COMPONENTS system filesystem program_options log log_setup ) 
    include_directories (${Boost_INCLUDE_DIRS}) 
#     set_property(TARGET ${Boost_INCLUDE_DIRS} PROPERTY POSITION_INDEPENDENT_CODE ON)
endif ( USE_BOOST )


###########################################
#SuiteSparse
############################################ 

find_path(CHOLMOD_INCLUDE_DIR NAMES cholmod.h PATHS "${PROJECT_SOURCE_DIR}/externalLibs/SuiteSparse/include" NO_DEFAULT_PATH )
if ( CHOLMOD_INCLUDE_DIR )
  list ( APPEND SUITESPARSE_INCLUDE_DIRS ${CHOLMOD_INCLUDE_DIR} )
endif( CHOLMOD_INCLUDE_DIR )

find_path( SUITESPARSE_LIBRARY_DIR NAMES libcholmod.so libcholmod.a PATHS "${PROJECT_SOURCE_DIR}/externalLibs/SuiteSparse/lib" NO_DEFAULT_PATH )
if ( SUITESPARSE_LIBRARY_DIR )
  list ( APPEND SUITESPARSE_LIBRARIES amd)
  list ( APPEND SUITESPARSE_LIBRARIES btf)
  list ( APPEND SUITESPARSE_LIBRARIES camd)
  list ( APPEND SUITESPARSE_LIBRARIES ccolamd)
  list ( APPEND SUITESPARSE_LIBRARIES cholmod)
  list ( APPEND SUITESPARSE_LIBRARIES colamd)
  #list ( APPEND SUITESPARSE_LIBRARIES csparse)
  list ( APPEND SUITESPARSE_LIBRARIES cxsparse)
  list ( APPEND SUITESPARSE_LIBRARIES klu)
  #list ( APPEND SUITESPARSE_LIBRARIES spqr)
  list ( APPEND SUITESPARSE_LIBRARIES umfpack)

  # from version 4 on an additional config library is needed
#   find_library( SUITESPARSE_CONFIG_LIB NAMES suitesparseconfig PATHS ${SUITESPARSE_LIBRARY_DIR})
  find_library( SUITESPARSE_CONFIG_LIB NAMES libsuitesparseconfig.so PATHS "${PROJECT_SOURCE_DIR}/externalLibs/SuiteSparse/lib" NO_DEFAULT_PATH )
  if ( EXISTS ${SUITESPARSE_CONFIG_LIB} )
    list ( APPEND SUITESPARSE_LIBRARIES ${SUITESPARSE_CONFIG_LIB} )
  endif ( EXISTS ${SUITESPARSE_CONFIG_LIB} )

  # Metis and spqr are optional
#   find_library( SUITESPARSE_METIS_LIBRARY NAMES metis PATHS ${SUITESPARSE_LIBRARY_DIR} )
#   if (SUITESPARSE_METIS_LIBRARY)                   
#     list ( APPEND SUITESPARSE_LIBRARIES metis)
#   endif(SUITESPARSE_METIS_LIBRARY)

#   if(EXISTS  "${CHOLMOD_INCLUDE_DIR}/SuiteSparseQR.hpp")
#     set(SUITESPARSE_SPQR_VALID TRUE CACHE BOOL "SuiteSparseSPQR valid")
#   else()
#     set(SUITESPARSE_SPQR_VALID false CACHE BOOL "SuiteSparseSPQR valid")
#   endif()
# 
#   if(SUITESPARSE_SPQR_VALID)
#     find_library( SUITESPARSE_SPQR_LIBRARY NAMES spqr PATHS ${SUITESPARSE_LIBRARY_DIR} )
#     if (SUITESPARSE_SPQR_LIBRARY)                 
#       list ( APPEND SUITESPARSE_LIBRARIES spqr)
#     endif (SUITESPARSE_SPQR_LIBRARY)
#   endif()
else( SUITESPARSE_LIBRARY_DIR )  
  MESSAGE ( WARNING "SUITESPARSE_LIBRARY_DIR not found!" )
endif( SUITESPARSE_LIBRARY_DIR )  

if (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)
    set(SUITESPARSE_FOUND TRUE)
else (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)
    set( SUITESPARSE_FOUND FALSE )
endif (SUITESPARSE_INCLUDE_DIRS AND SUITESPARSE_LIBRARIES)

include_directories(${SUITESPARSE_INCLUDE_DIRS})
link_directories(${SUITESPARSE_LIBRARY_DIR})

###########################################
#Eigen
############################################ 
include_directories ("${PROJECT_SOURCE_DIR}/externalLibs/eigen")
                
 
###########################################
#IPOPT
###########################################
option( USE_IPOPT "Interior Point Solver" ON )
if( USE_IPOPT )

 find_path(IPOPT_INCLUDE_DIR IpTNLP.hpp PATHS "${PROJECT_SOURCE_DIR}/externalLibs/Ipopt/build/include/coin" NO_DEFAULT_PATH )  
 find_library(IPOPT_LIBRARY ipopt PATHS "${PROJECT_SOURCE_DIR}/externalLibs/Ipopt/build/lib" NO_DEFAULT_PATH ) 
#  find_path(IPOPT_INCLUDE_DIR IpTNLP.hpp PATHS "${PROJECT_SOURCE_DIR}/externalLibs/Ipopt_WithHSL/build/include/coin" NO_DEFAULT_PATH )  
#  find_library(IPOPT_LIBRARY ipopt PATHS "${PROJECT_SOURCE_DIR}/externalLibs/Ipopt_WithHSL/build/lib" NO_DEFAULT_PATH ) 
 
 set(IPOPT_INCLUDE_DIRS "${IPOPT_INCLUDE_DIR}" )
 set(IPOPT_LIBRARIES "${IPOPT_LIBRARY}")    
 add_definitions(-DUSE_IPOPT)
 #wrong config under Debian workaround    
 add_definitions( -DHAVE_CSTDDEF )    
 
 include(FindPackageHandleStandardArgs)   
 # handle the QUIETLY and REQUIRED arguments and set LIBIPOPT_FOUND to TRUE if all listed variables are TRUE   
 find_package_handle_standard_args(IPOPT  DEFAULT_MSG  IPOPT_LIBRARY IPOPT_INCLUDE_DIR) 
 mark_as_advanced(IPOPT_INCLUDE_DIR IPOPT_LIBRARY )  
 
 include_directories(SYSTEM ${IPOPT_INCLUDE_DIR})
endif( USE_IPOPT )


###########################################
#OpenMP
############################################  
option( USE_OPENMP "OpenMP for parallel computing" ON )
if ( USE_OPENMP )
 find_package(OpenMP)
endif( USE_OPENMP )


###########################################
#VTK
############################################  
option( USE_VTK "VTK Library" ON )
if ( USE_VTK )
    set(VTK_DIR "${PROJECT_SOURCE_DIR}/externalLibs/VTK/build" CACHE PATH "VTK directory override" FORCE)
    find_package(VTK)
    #find_path(VTK_DIR PATHS "${PROJECT_SOURCE_DIR}/externalLibs/VTK/build" NO_DEFAULT_PATH )  
    include(${VTK_USE_FILE})
    add_definitions( -DUSE_VTK )
endif( USE_VTK )


############################################
#Libraries
############################################  
include_directories(modules/aol)
include_directories(modules/Energies)
include_directories(modules/materialOptimization)
include_directories(modules/Solver)

include_directories(modules/finiteElements/quocElements)

include_directories(external/Eigen)
include_directories(external/IPOPT)
include_directories(external/VTK)


############################################
#Projects
############################################ 
add_subdirectory(projects/elasticBodies/PeriodicHomogenizationBones)




#  To let the script know where the source and thus the error parse scripts are sneak in the path as argument.
set_property ( GLOBAL PROPERTY RULE_LAUNCH_COMPILE "bash ${CMAKE_CURRENT_SOURCE_DIR}/util/cmakeParseError.sh ${CMAKE_CURRENT_SOURCE_DIR}" )
