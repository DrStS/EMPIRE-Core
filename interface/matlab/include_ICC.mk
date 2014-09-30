CC  = icc
CXX = icpc
FC  = ifort
LINKER = $(CXX)

CFLAGS    = -g
CXXFLAGS  = $(CFLAGS)
CPPFLAGS  = -DMATLAB_MEX_FILE
FCFLAGS   = 
LFLAGS    =
DEFINES   =

# EMPIRE API
INCLUDES  = -I$(EMPIRE_API_INC_ON_MACHINE)
INCLUDES += -I./GiDFileIO
LIBS      = $(EMPIRE_API_LIBSO_ON_MACHINE)
LIBS     += -Wl,-whole-archive GiDFileIO/GiDFileIO.a -Wl,-no-whole-archive


# LINK MATLIB INTERFACE LIBRARIES
MATLIB_HOME = /opt/MATLAB/R2012a
MATLIB_ARC = glnxa64
INCLUDES += -I$(MATLIB_HOME)/extern/include
#RPATH = -Wl,-rpath,$(MATLIB_HOME)/bin/glnxa64
#RPATHL = -Wl,-rpath-link,$(MATLIB_HOME)/bin/glnxa64



CPPFLAGS += #-ansi -D_GNU_SOURCE
CPPFLAGS += -fPIC #-fno-omit-frame-pointer -pthread
LIBS += $(RPATH) $(RPATHL) -L$(MATLIB_HOME)/bin/$(MATLIB_ARC) -lmx -lmex -lmat -lm

LFLAGS += -shared -pthread -Wl,--version-script,mexFunction.map
LFLAGS += -Wl,--no-undefined

