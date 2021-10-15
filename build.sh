#!/bin/bash

if [ -d build ] ; then rm -rf build ; fi

mkdir build 
cd build


# Need to find your compiler. I know clang works, but not sure about g++
cxx="/usr/bin/g++"
# Need to point to your version of Hadrons
cxxflags="-I/Users/kshitij/Desktop/Software/HadronsR/install/include  "`/Users/kshitij/Desktop/Software/HadronsR/install/bin/hadrons-config --cxxflags`
# Need to point to your version of Hadrons, omp, gmp, mpfr
ldflags="-L/usr/local/Cellar/libomp/lib -L/usr/local/lib -L/usr/local/Cellar/gmp/6.2.1/lib -L/usr/local/Cellar/mpfr/4.1.0/lib -L/Users/kshitij/Desktop/Software/HadronsR/install/lib "`/Users/kshitij/Desktop/Software/HadronsR/install/bin/hadrons-config --ldflags` 
# Need to point to your version of Hadrons
libs=`/Users/kshitij/Desktop/Software/HadronsR/install/bin/hadrons-config --libs`


cmake \
    -DCMAKE_CXX_COMPILER=$cxx \
    -DMPI_CXX_COMPILER=/usr/local/bin/mpicxx \
    -DCMAKE_CXX_FLAGS="$cxxflags $libs $ldflags -I/Users/kshitij/Desktop/Software/GridR/install/include" \
    -DCMAKE_C_FLAGS=-I/usr/local/Cellar/libomp/12.0.1/include \
    ../


make -j 6
make -j 6 install



