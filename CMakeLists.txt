cmake_minimum_required(VERSION 2.8.10)
project(DigitalImagesLibrary)

set(SOURCES
    src/fourier_transform/CosFourier.cpp
    src/fourier_transform/Fourier.cpp
    src/fourier_transform/InverseCosFourier.cpp
    src/fourier_transform/InverseFourier.cpp
    src/fourier_transform/InverseSinFourier.cpp
    src/fourier_transform/SinFourier.cpp
    src/registration/FluidCurvatureRegistration.cpp
    src/solver/RKNystroem.cpp
    src/solver/RKV43.cpp
    src/templates/VectorArray2D.cpp
    src/utilities/BracketMethod.cpp
    src/utilities/ImageDifference.cpp
    src/utilities/Utilities.cpp
)

include_directories(
    src
    lib
    /opt/local/include
)

add_custom_target(
    Version
    COMMAND /bin/bash version.sh
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_executable(fluidReg
    src/fluidReg.cpp
    ${SOURCES}
)
add_dependencies(fluidReg Version)
find_library(INTEL_LIBRARIES
    NAMES tbb
)
find_library(FFTW3_LIBRARIES NAMES fftw3)
find_library(X11_LIBRARIES NAMES X11)
target_link_libraries(fluidReg pthread ${X11_LIBRARIES} ${INTEL_LIBRARIES} ${FFTW3_LIBRARIES})

if(UNIX AND NOT APPLE)
    message("\nEntering configuration for UNIX computers:\n")
    set(CMAKE_CXX_FLAGS "-O4 -pipe -Wall")
    set(CMAKE_C_FLAGS "-O4 -pipe -Wall")
    set_target_properties(fluidReg
        PROPERTIES
        LINK_FLAGS ""
     )
endif()

if(APPLE)
    message("\nEntering configuration for APPLE computers:\n")
    set(CMAKE_CXX_FLAGS "-O4 -pipe -Wall -I/opt/local/include")
    set(CMAKE_C_FLAGS "-O4 -pipe -Wall")
    set_target_properties(fluidReg
        PROPERTIES
        LINK_FLAGS "-v"
     )
endif()
