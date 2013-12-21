#------------------------------------------------------------------------------------#
file(GLOB SOURCES *.cpp)
MACRO_APPEND_GLOBAL_VARIABLE(Emperor_SOURCES "${SOURCES}")
#------------------------------------------------------------------------------------#
MACRO_APPEND_GLOBAL_VARIABLE(Emperor_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}")
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
get_property(Emperor_SOURCES GLOBAL PROPERTY Emperor_SOURCES) 
get_property(EMPIRE_thirdparty_INCLUDES GLOBAL PROPERTY EMPIRE_thirdparty_INCLUDES) 
#------------------------------------------------------------------------------------#
# MPI
IF(MPI_CXX_INCLUDE_PATH)
include_directories(${MPI_CXX_INCLUDE_PATH})
ENDIF()
#------------------------------------------------------------------------------------#
include_directories(${Emperor_INCLUDES})
include_directories(${EMPIRE_thirdparty_INCLUDES})
LIST(REMOVE_ITEM Emperor_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_library(EmperorLib STATIC ${Emperor_SOURCES} $<TARGET_OBJECTS:EMPIRE_thirdparty>)
#------------------------------------------------------------------------------------#
add_executable(Emperor ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
target_link_libraries(Emperor EmperorLib ${Emperor_LIBS} ${MPI_CXX_LINK_FLAGS})
#------------------------------------------------------------------------------------#
add_dependencies(Emperor EmperorLib EMPIRE_thirdparty)
#------------------------------------------------------------------------------------#
MESSAGE("C Compiler Flags are:\t${CMAKE_C_FLAGS} \nC++ Compiler Flags are:\t${CMAKE_CXX_FLAGS}")
#MESSAGE("Emperor_LIBS: ${Emperor_LIBS}")
#MESSAGE("Emperor_INCLUDES: ${Emperor_INCLUDES}")
#MESSAGE("Emperor_SOURCES: ${Emperor_SOURCES}")
#------------------------------------------------------------------------------------#