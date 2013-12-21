#  Copyright &copy; 2013, TU Muenchen, Chair of Structural Analysis,
#  Stefan Sicklinger, Tianyang Wang, Munich
#
#  All rights reserved.
#
#  This file is part of EMPIRE.
#
#  EMPIRE is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  EMPIRE is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with EMPIRE.  If not, see http://www.gnu.org/licenses/.
#------------------------------------------------------------------------------------#
project (EMPIRE)
#------------------------------------------------------------------------------------#
cmake_minimum_required (VERSION 2.8.10)
#------------------------------------------------------------------------------------#
IF( (${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")) 
    enable_language(Fortran OPTIONAL)
    MESSAGE("EMPIRE is going to be configured for Linux x86_64!")
ELSEIF( (${CMAKE_SYSTEM_NAME} MATCHES "Windows" AND ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "AMD64")) 
    MESSAGE("EMPIRE is going to be configured for Windows x64!")	
ELSE() 
    MESSAGE(FATAL_ERROR "EMPIRE is supported on Linux x86_64 and Windows x64 only")	
ENDIF()
#------------------------------------------------------------------------------------#
# 1. Options
# 2. Macros
# 3. Compiler
# 4. Subdirectories
# 5. Configure
# 6. MPI
#------------------------------------------------------------------------------------#
######################################################################################
# THESE GLOBAL VARIABLES ARE SET DURING THE CONFIGURATION
######################################################################################
#EMPIRE_thirdparty_SOURCES
#EMPIRE_thirdparty_INCLUDES
#Emperor_SOURCES
#Emperor_INCLUDES
#Emperor_API_SOURCES
#Emperor_API_INCLUDES
######################################################################################
#1. Option
######################################################################################
#option(USE_INTEL_COMPILERS_MKL_IMPI      "Use Intel Compilers C/C++,     Intel MKL and Intel MPI"  ON )
#option(USE_MICROSOFT_COMPILERS_MKL_IMPI  "Use Microsoft Compilers C/C++, Intel MKL and Intel MPI"  OFF )
option(BUILD_FORTRAN_CLIENTS           "This builds FORTRAN test clients"     OFF )
######################################################################################
#2. Macros
######################################################################################
# Add EMPIRE cmake modules
set(CMAKE_MODULE_PATH ${EMPIRE_SOURCE_DIR}/cmake/modules)
include(MacroOutOfSourceBuild)
include(MacroAppendGlobalVariable)
include(ExternalProject)
include(GetGitRevisionDescription)
include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
######################################################################################
#3. Compiler
######################################################################################

#------------------------------------------------------------------------------------#
# Avoid in source builds
macro_ensure_out_of_source_build("In source build is not allowed for ${PROJECT_NAME}. Please run CMake from your build directory. You may want to remove CMakeCache.txt.")
#------------------------------------------------------------------------------------#
# Check compilers and their versions
IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" AND "${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C/C++
  SET(USE_INTEL_COMPILERS_MKL_IMPI_LINUX ON)
  MESSAGE("C compiler is ${CMAKE_C_COMPILER_ID} and has version ${CMAKE_C_COMPILER_VERSION} \nC++ compiler is ${CMAKE_CXX_COMPILER_ID} and has version ${CMAKE_CXX_COMPILER_VERSION} " ) 
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC" AND "${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
  # using MSVC
  SET(USE_INTEL_COMPILERS_MKL_MSMPI_WINDOWS ON)
  MESSAGE("C compiler is ${CMAKE_C_COMPILER_ID} and has version ${CMAKE_C_COMPILER_VERSION} \nC++ compiler is ${CMAKE_CXX_COMPILER_ID} and has version ${CMAKE_CXX_COMPILER_VERSION} " ) 
  MESSAGE("As OpenMP 3.0 and VLA (C99) is used in EMPIRE you need to switch to Intel Compiler in Visual Studio" ) 
ELSE()
  MESSAGE(FATAL_ERROR "Your CXX compiler is ${CMAKE_CXX_COMPILER_ID} and your C compiler is ${CMAKE_C_COMPILER_ID} they are not supported by the build system yet")	
ENDIF()
#------------------------------------------------------------------------------------#
# Set Intel Compiler & MKL
IF (${USE_INTEL_COMPILERS_MKL_IMPI_LINUX})
  SET (CMAKE_C_COMPILER "mpiicc")
  SET (CMAKE_CXX_COMPILER "mpiicpc")
  #SET (MKL_INCLUDE_DIR "")
  #SET (MKL_LIB_DIR     "")
  find_package(MKL REQUIRED)
  #Let the source code know
  add_definitions(-DUSE_INTEL_MKL) 
ENDIF()

IF (${USE_INTEL_COMPILERS_MKL_MSMPI_WINDOWS})
  # Require MPI for this project:
  SET(MPI_ROOT_PATH "C:/Program Files/Microsoft HPC Pack 2012")
  SET(MPI_CXX_INCLUDE_PATH "${MPI_ROOT_PATH}/Inc")
  SET(MPI_CXX_LINK_FLAGS "${MPI_ROOT_PATH}/Lib/amd64/msmpi.lib")
  #SET (MKL_INCLUDE_DIR "")
  #SET (MKL_LIB_DIR     "")
  find_package(MKL REQUIRED)
  #Let the source code know
  add_definitions(-DUSE_INTEL_MKL)
ENDIF()


#------------------------------------------------------------------------------------#
# Global valid defines
add_definitions(-DFLANN) # TEMP
add_definitions(-DMPICH_IGNORE_CXX_SEEK)
#------------------------------------------------------------------------------------#
# Set Compiler flags and Linker flags (CMAKE_BUILD_TYPE is not used!) 
## Here INTEL specific flags are taken care for
IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" AND "${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
## This is for all compilers
    check_c_compiler_flag(-xHost   CHECKXHOST)
    check_cxx_compiler_flag(-xHost CHECKXHOST)
	SET (CMAKE_C_FLAGS "-xHost -std=c99")
	SET (CMAKE_CXX_FLAGS "-xHost")
ENDIF()
## Here Intel Windows specific flags are taken care for
IF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC" AND "${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
	SET (CMAKE_C_FLAGS "/QxHost /Qstd=c99")
	SET (CMAKE_CXX_FLAGS "/QxHost /Qstd=c++11")
ENDIF()
######################################################################################
#4. Configure
######################################################################################
IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
    configure_file("${PROJECT_SOURCE_DIR}/configure/bashrc_config" "${PROJECT_SOURCE_DIR}/etc/bashrc.sh")
ENDIF()
IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
    configure_file("${PROJECT_SOURCE_DIR}/configure/batch_config" "${PROJECT_SOURCE_DIR}/etc/batch.bat")
ENDIF()
######################################################################################
#5. Subdirectories THE ORDER IS IMPORTANT
######################################################################################
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../EMPIRE-Libraries" "${CMAKE_CURRENT_BINARY_DIR}/EMPIRE-Libraries")
add_subdirectory(Emperor)
#add_subdirectory(EMPIRE_API)
#add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../EMPIRE-Testing" "${CMAKE_CURRENT_BINARY_DIR}/EMPIRE-Testing")

