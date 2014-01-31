#------------------------------------------------------------------------------------#
# Openmp parallelism
find_package(OpenMP 3.0 REQUIRED)
IF(OPENMP_FOUND)
  IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()
  IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
    SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Qopenmp")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Qopenmp")
ENDIF()
ENDIF(OPENMP_FOUND)
#------------------------------------------------------------------------------------#
IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -C")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -C")
ENDIF()
IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /C")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /C")
ENDIF()
#------------------------------------------------------------------------------------#
include_directories(${MKL_INCLUDE_DIR})
SET(Emperor_LIBS ${Emperor_LIBS} ${MKL_LIBRARIES})
#------------------------------------------------------------------------------------#
ExternalProject_Get_Property(flann INSTALL_DIR)
LINK_DIRECTORIES("${INSTALL_DIR}/libflann_s.a")
#------------------------------------------------------------------------------------#
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#------------------------------------------------------------------------------------#
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_describe(GIT_TAG --tags)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/tools/AuxiliaryParameters.cpp.in" 
"${CMAKE_CURRENT_SOURCE_DIR}/src/tools/AuxiliaryParameters.cpp" @ONLY)
#------------------------------------------------------------------------------------#
add_subdirectory(src)
add_subdirectory(testUnit)
add_subdirectory(testMapper)