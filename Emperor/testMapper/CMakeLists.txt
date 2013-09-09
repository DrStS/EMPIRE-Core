#-------------------------------------------------------------------------------
file(GLOB SOURCES *.cpp)
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
