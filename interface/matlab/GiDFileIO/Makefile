# LINK MATLIB INTERFACE LIBRARIES
MATLIB_HOME = /opt/MATLAB/R2012a
MATLIB_ARC = glnxa64
INCLUDES += -I$(MATLIB_HOME)/extern/include
#RPATH = -Wl,-rpath,$(MATLIB_HOME)/bin/glnxa64

# libGiDFileIO.so : GiDFileIO.cpp
# 	icpc -c -fPIC GiDFileIO.cpp -o GiDFileIO.o 
# 	icpc $(RPATH) -shared -Wl,-soname,libGiDFileIO.so -o libGiDFileIO.so GiDFileIO.o

GiDFileIO.a : GiDFileIO.cpp
	icpc $(RPATH) -c -fPIC GiDFileIO.cpp -o GiDFileIO.o #without -fPIC cannot include .a into .so (.mex)
	ar -cvq GiDFileIO.a GiDFileIO.o
