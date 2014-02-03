libGiDFileIO.so : GiDFileIO.cpp
	icpc -c -fPIC GiDFileIO.cpp -o GiDFileIO.o 
	icpc -shared -Wl,-soname,libGiDFileIO.so -o libGiDFileIO.so GiDFileIO.o

GiDFileIO.a : GiDFileIO.cpp
	icpc -c -fPIC GiDFileIO.cpp -o GiDFileIO.o #without -fPIC cannot include .a into .so (.mex)
	ar -cvq GiDFileIO.a GiDFileIO.o