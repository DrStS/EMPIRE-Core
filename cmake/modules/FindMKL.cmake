# a cmake module to find Intel Math Kernel Library
#------------------------------------------------------------------------------------#
# MKL includes and libraries are searched for in MKL_INCLUDE_DIR and MKL_LIB_DIR.
# If MKL_INCLUDE_DIR and MKL_LIB_DIR are not set, the module searches under 
# the environment variable MKLROOT and /opt/intel/mkl subdirectories.
#------------------------------------------------------------------------------------#

set(MKLROOT_DIR $ENV{MKLROOT})
if (NOT MKLROOT_DIR)
   if (EXISTS "/opt/intel/mkl")
      set(MKLROOT_DIR "/opt/intel/mkl")
   endif (EXISTS "/opt/intel/mkl")
endif ()

#------------------------------------------------------------------------------------#
# Stage 1: find the include directory
#------------------------------------------------------------------------------------#
if (NOT MKL_INCLUDE_DIR)
  find_path(MKL_INCLUDE_DIR
    NAMES mkl.h
    HINTS $MKLROOT_DIR
    /opt/intel/mkl
    /usr/intel/mkl
    PATH_SUFFIXES include
    )  
endif ()

#------------------------------------------------------------------------------------#
# Stage 2: find the lib directory
#------------------------------------------------------------------------------------#	
if (NOT MKL_LIB_DIR)
  if (MKLROOT_DIR)
    if (CMAKE_SYSTEM_NAME MATCHES "Linux")
      if (CMAKE_SIZEOF_VOID_P MATCHES 8)
	set(EXPECT_MKL_LIBPATH "${MKLROOT_DIR}/lib/intel64")
      else ()
	set(EXPECT_MKL_LIBPATH "${MKLROOT_DIR}/lib/ia32")
      endif ()
    endif ()	
    if (IS_DIRECTORY ${EXPECT_MKL_LIBPATH})
      set(MKL_LIB_DIR ${EXPECT_MKL_LIBPATH})
    endif ()
  endif ()
endif ()

#------------------------------------------------------------------------------------#
# Stage 3: find the libraries
#------------------------------------------------------------------------------------#	
set (CMAKE_FIND_LIBRARY_SUFFIXES .a)
if (MKL_LIB_DIR)
  find_library(MKL_INTEL_LP64_LIBRARY     mkl_intel_lp64     ${MKL_LIB_DIR})
  find_library(MKL_INTEL_THREAD_LIBRARY   mkl_intel_thread   ${MKL_LIB_DIR})
  find_library(MKL_CORE_LIBRARY           mkl_core           ${MKL_LIB_DIR})
  if (MKL_INTEL_LP64_LIBRARY AND MKL_INTEL_THREAD_LIBRARY AND MKL_CORE_LIBRARY)
    set (MKL_LIBRARIES "-Wl,--start-group ${MKL_INTEL_LP64_LIBRARY} ${MKL_INTEL_THREAD_LIBRARY} ${MKL_CORE_LIBRARY} -Wl,--end-group -lpthread")
  else ()
    set (MKL_LIBRARIES "")
  endif ()
endif ()

# set MKL_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INTEL_LP64_LIBRARY MKL_INTEL_THREAD_LIBRARY MKL_CORE_LIBRARY)
#mark_as_advanced(LIB_PTHREAD MKL_INCLUDE_DIR MKL_LIBRARY_DIR)

