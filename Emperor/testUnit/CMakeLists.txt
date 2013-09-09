#------------------------------------------------------------------------------------#
file(GLOB SOURCES *.cpp)
MACRO_APPEND_GLOBAL_VARIABLE(Emperor_UnitTest_SOURCES "${SOURCES}")
#------------------------------------------------------------------------------------#
MACRO_APPEND_GLOBAL_VARIABLE(Emperor_UnitTest_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}" )
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
get_property(Emperor_UnitTest_INCLUDES GLOBAL PROPERTY Emperor_UnitTest_INCLUDES) 
get_property(Emperor_UnitTest_SOURCES GLOBAL PROPERTY Emperor_UnitTest_SOURCES) 
get_property(EMPIRE_thirdparty_INCLUDES GLOBAL PROPERTY EMPIRE_thirdparty_INCLUDES) 
#------------------------------------------------------------------------------------#
include_directories(${Emperor_INCLUDES})
include_directories(${Emperor_UnitTest_INCLUDES})
include_directories(${EMPIRE_thirdparty_INCLUDES})
#------------------------------------------------------------------------------------#
add_executable(EmperorUnitTest ${Emperor_UnitTest_SOURCES})
target_link_libraries(EmperorUnitTest EmperorLib ${Emperor_LIBS})
#------------------------------------------------------------------------------------#
add_dependencies(EmperorUnitTest EmperorLib)
#------------------------------------------------------------------------------------#
#MESSAGE("C Compiler Flags are:\t${CMAKE_C_FLAGS} \nC++ Compiler Flags are:\t${CMAKE_CXX_FLAGS}")
#MESSAGE("Emperor_UnitTest_LIBS: ${Emperor_UnitTest_LIBS}")
#MESSAGE("Emperor_UnitTest_INCLUDES: ${Emperor_UnitTest_INCLUDES}")
#MESSAGE("Emperor_UnitTest_SOURCES: ${Emperor_UnitTest_SOURCES}")
#------------------------------------------------------------------------------------#