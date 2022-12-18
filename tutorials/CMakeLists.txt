## -------------------------------------------------------------------------- ##
#
#   Project: Example
#
#   Description:
#       CMake & Fortran project for graduate students
#
#   Author:
#       Young-Myung CHOI
#         Dept. of Naval Arch. and Ocean Eng.,
#         Pusan National University, Republic of Korea
#
## -------------------------------------------------------------------------- ##

## Default Settings --------------------------------------------------------- ##

##... Set the minimum version of cmake
cmake_minimum_required(VERSION 2.9 FATAL_ERROR)

##... CMake Project Name
project(FortranProjectBase)

##... Set Default Fortran Compiler & Eable the Fortran Compilation Environment
enable_language(Fortran)
if (NOT DEFINED CMAKE_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER gfortran)
endif()

## CMake Configuration ------------------------------------------------------ ##

##... Load Default CMake Configuration
include( config/config.cmake )

##... Set Fortran Compiler Flags
include( config/fortranFlag.cmake )

##... Auto Compilation Functions
include( config/Func_AutoCompile.cmake )

## Print CMake Project Information ------------------------------------------ ##
message("")
message("${BoldYellow}Project: ${BoldWhite}${CMAKE_PROJECT_NAME}")
message("")
message("  - ${Green}Projct Path       ${White}: ${PROJECT_COMPILE_PATH} ${Green}")
message("  - ${Green}GIT_BRANCH        ${White}: ${GIT_BRANCH} ${Green}")
message("  - ${Green}GIT_COMMIT_HASH   ${White}: ${GIT_COMMIT_HASH} ${Green}")
message("")
message("  - ${Green}Lib. Install Path ${White}: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${Green}")
message("  - ${Green}Exe. Install Path ${White}: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} ${Green}")
message("")

# set(FORTRAN_COMPILE_MODE "Release")
# set(FORTRAN_COMPILE_MODE "Debug")

if( ${FORTRAN_COMPILE_MODE} STREQUAL "Release" )
message("  - ${Green}Compilation Mode  ${White}: Release ${Green}")
elseif( ${FORTRAN_COMPILE_MODE} STREQUAL "Debug" )
message("  - ${Green}Compilation Mode  ${White}: Debug ${Green}")
endif()

message("${Green}")
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
message("  - ${Green}Fortran Compiler  ${White}: Intel Fortran ${Green}")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
message("  - ${Green}Fortran Compiler  ${White}: GNU Fortran ${Green}")
endif()
message("  - ${Green}Global Flag       ${White}: ${CMAKE_Fortran_FLAGS} ${Green}")
message("${ColourReset}")
## Compile Library ---------------------------------------------------------- ##

set(COMPILE_DYNAMIC_LIB False)
# set(COMPILE_DYNAMIC_LIB True)

message(${CMAKE_CURRENT_SOURCE_DIR})

set(LAPACK_FLAG
    ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/lapack/liblapack.a
    ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/lapack/librefblas.a )

#... Library
set(LIB_NAME_mathPack mathPack)
# set(LIB_PATH_mathPack )
# set(LIB_DEPEND_mathPack )
set(LIB_EXTLIB_mathPack ${LAPACK_FLAG})
# set(LIB_EXTINC_mathPack )
# set(LIB_EXTFLAG_mathPack )
add_subdirectory( src )

## Compile Applications ----------------------------------------------------- ##

#... Application
set(APP_NAME_appTutorial appTutorial)
# set(APP_PATH_appTutorial )
set(APP_DEPEND_appTutorial mathPack)
set(APP_EXTLIB_appTutorial ${LAPACK_FLAG})
# set(APP_EXTINC_appTutorial )
# set(APP_EXTFLAG_appTutorial  )
add_subdirectory(applications/appTutorial)
