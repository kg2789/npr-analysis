# The version here is the version I had when writing this; I have not tested
# any other versions
cmake_minimum_required (VERSION 3.13)
project(analysis CXX C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-O2")

SET (CMAKE_CXX_COMPILER "/usr/bin/g++")


# SET (CMAKE_CXX_COMPILER "/usr/bin/clang++")
# Try compiling with clang++


find_package(HDF5 REQUIRED)
# find_package(MPI REQUIRED)
find_package(OpenSSL REQUIRED)

#Try compiling openmpi if doesn't work without

include_directories(${HDF5_INCLUDE_DIRS})
# include_directories(${MPI_CXX_INCLUDE_DIRS})
include_directories(${OPENSSL_INCLUDE_DIR})

# Sets ${EVALUATE_SOURCES}
add_subdirectory(evaluate)
# sets ${DATA_CONVERSION_SOURCES}, ${DATA_CONVERSION_INCLUDES}
add_subdirectory(data)


## From Bigeng's CMake
include_directories("/usr/local/Cellar/tclap/1.2.4/include")
include_directories("/usr/local/Cellar/hdf5/1.12.1/include")
include_directories("/Users/kshitij/Desktop/Software/Grid/install/include")
include_directories("/Users/kshitij/Desktop/Software/lime/install/include")
# include_directories("/home/diaogen/Software/install/openmpi/include")
include_directories("/usr/local/Cellar/hdf5/1.12.1/lib")
include_directories("/Users/kshitij/Desktop/Software/Grid/install/lib")
#include_directories("/Users/kshitij/Desktop/Software/HadronsR/install/include")
#include_directories("/Users/kshitij/Desktop/Software/HadronsR/install/lib")


set(ANALYSIS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/amputate.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/evaluate.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/projector.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/parity.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/subtraction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tree.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/util.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/data.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/msbar.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Zq.cpp
    ${DATA_CONVERSION_SOURCES} ${EVALUATE_SOURCES})


add_executable(analysis main.cpp ${ANALYSIS_SOURCES})
add_executable(convert data/main.cpp ${DATA_CONVERSION_SOURCES})
add_executable(compare data/compare.cpp ${DATA_CONVERSION_SOURCES})
add_executable(merge data/merge.cpp ${DATA_CONVERSION_SOURCES})

include_directories(/Users/kshitij/Desktop/Software/GridR/install/include) 


# set(GRID_LIBS "-lGrid" hdf5 hdf5_cpp "-lz" "-lcrypto" "-llime" "-lmpfr" "-lgmp" "-lstdc++" "-lm" "-lz")

##From Bigeng's CMake
set(GRID_LIBS "-L/usr/local/Cellar/hdf5/1.12.1/lib -L/Users/kshitij/Desktop/Software/Grid/install/lib -lGrid"  hdf5 hdf5_cpp )

target_link_libraries(analysis ${GRID_LIBS})
target_link_libraries(convert ${GRID_LIBS})
target_link_libraries(compare ${GRID_LIBS})
target_link_libraries(merge ${GRID_LIBS})

target_include_directories(convert PRIVATE ${DATA_CONVERSION_INCLUDES})
target_include_directories(compare PRIVATE ${DATA_CONVERSION_INCLUDES})
target_include_directories(merge PRIVATE ${DATA_CONVERSION_INCLUDES})

enable_testing()
add_subdirectory(tests)
