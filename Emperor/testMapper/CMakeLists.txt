#-------------------------------------------------------------------------------
# Copyright (c) 2013, TU Muenchen, Chair of Structural Analysis,
# Stefan Sicklinger, Tianyang Wang, Munich. All rights reserved.
# 
# This file is part of EMPIRE.
# 
# EMPIRE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# EMPIRE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with EMPIRE.  If not, see http://www.gnu.org/licenses/.
#-------------------------------------------------------------------------------
SET(Emperor_MapperTest_SOURCES "${SOURCES}")
#------------------------------------------------------------------------------------#
SET(Emperor_MapperTest_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}" )
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
FILE(GLOB Emperor_src_subdirectories "*")
FOREACH(f ${Emperor_src_subdirectories})
  IF(NOT f MATCHES "\\.txt" AND IS_DIRECTORY ${f} )
    add_subdirectory(${f})
  ENDIF()
ENDFOREACH()
#------------------------------------------------------------------------------------#
get_property(Emperor_INCLUDES GLOBAL PROPERTY Emperor_INCLUDES)
get_property(EMPIRE_thirdparty_INCLUDES GLOBAL PROPERTY EMPIRE_thirdparty_INCLUDES) 
#------------------------------------------------------------------------------------#
include_directories(${Emperor_INCLUDES})
include_directories(${EMPIRE_thirdparty_INCLUDES})
#------------------------------------------------------------------------------------#
add_executable(EmperorMapperTest ${Emperor_MapperTest_SOURCES})
target_link_libraries(EmperorMapperTest EmperorLib ${Emperor_LIBS})
#------------------------------------------------------------------------------------#
add_dependencies(EmperorMapperTest EmperorLib)
#------------------------------------------------------------------------------------#
#MESSAGE("C Compiler Flags are:\t${CMAKE_C_FLAGS} \nC++ Compiler Flags are:\t${CMAKE_CXX_FLAGS}")
#------------------------------------------------------------------------------------#
