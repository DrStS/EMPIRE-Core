CC  = gcc
CXX = g++
FC  = gfortran
LINKER = $(CXX)

CFLAGS    = -g
CXXFLAGS  = $(CFLAGS)
CPPFLAGS  = -DMATLAB_MEX_FILE
FCFLAGS   = 
LFLAGS    =
DEFINES   =

# EMPIRE API
INCLUDES = -I$(EMPIRE_API_INC_ON_MACHINE)
LIBS      = $(EMPIRE_API_LIBSO_ON_MACHINE)


# LINK MATLIB INTERFACE LIBRARIES
MATLIB_HOME = /opt/MATLAB/R2012b
MATLIB_ARC = glnxa64
INCLUDES += -I$(MATLIB_HOME)/extern/include
#RPATH = -Wl,-rpath-link,$(MATLIB_HOME)/bin/glnxa64



CPPFLAGS += #-ansi -D_GNU_SOURCE
CPPFLAGS += -fPIC #-fno-omit-frame-pointer -pthread
LIBS += $(RPATH) -L$(MATLIB_HOME)/bin/$(MATLIB_ARC) -lmx -lmex -lmat -lm 

LFLAGS += -shared -pthread -Wl,--version-script,mexFunction.map
LFLAGS += -Wl,--no-undefined

