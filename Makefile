SRC_DIR=src
BUILD_DIR=build

CXX=icpc
CXXFLAGS=-parallel -O3 -fopenmp -I/opt/intel/tbb/include -Ilib -Ilib/CImg-1.5.0 -Isrc/filter -Isrc/fourier_transform -Isrc/registration -Isrc/segmentation -Isrc/templates -Isrc/utilities
LFLAGS=-L/opt/intel/composer_xe_2011_sp1/lib/ -lirc -limf -liomp5 -L/opt/intel/tbb/lib -ltbb -L/usr/X11R6/lib -lpthread -lX11 -Llib/MacOSX-x86-64 -lfftw3 -lstdc++

OBJECTS=${BUILD_DIR}/fluidReg.o ${BUILD_DIR}/FluidCurvatureRegistration.o ${BUILD_DIR}/Fourier.o \
	${BUILD_DIR}/InverseFourier.o ${BUILD_DIR}/SinFourier.o ${BUILD_DIR}/InverseSinFourier.o ${BUILD_DIR}/CosFourier.o \
	${BUILD_DIR}/InverseCosFourier.o ${BUILD_DIR}/RKNystroem.o ${BUILD_DIR}/RKV43.o ${BUILD_DIR}/VectorArray2D.o \
	${BUILD_DIR}/Utilities.o ${BUILD_DIR}/ImageDifference.o ${BUILD_DIR}/BracketMethod.o

all : ${OBJECTS}
	${CXX} ${LFLAGS} -o fluidReg ${OBJECTS}

${BUILD_DIR}/fluidReg.o : ${SRC_DIR}/fluidReg.cpp
	./version.sh
	${CXX} ${CXXFLAGS} -c ${SRC_DIR}/fluidReg.cpp -o ${BUILD_DIR}/fluidReg.o

${BUILD_DIR}/FluidCurvatureRegistration.o : ${SRC_DIR}/registration/FluidCurvatureRegistration.cpp ${SRC_DIR}/registration/FluidCurvatureRegistration.hpp
	${CXX} ${CXXFLAGS} -c ${SRC_DIR}/registration/FluidCurvatureRegistration.cpp -o ${BUILD_DIR}/FluidCurvatureRegistration.o

${BUILD_DIR}/Fourier.o : ${SRC_DIR}/fourier_transform/Fourier.cpp ${SRC_DIR}/fourier_transform/Fourier.hpp
	${CXX} ${CXXFLAGS} -c ${SRC_DIR}/fourier_transform/Fourier.cpp -o ${BUILD_DIR}/Fourier.o

${BUILD_DIR}/InverseFourier.o : ${SRC_DIR}/fourier_transform/InverseFourier.cpp ${SRC_DIR}/fourier_transform/InverseFourier.hpp
	${CXX} ${CXXFLAGS} -c ${SRC_DIR}/fourier_transform/InverseFourier.cpp -o ${BUILD_DIR}/InverseFourier.o

${BUILD_DIR}/CosFourier.o : ${SRC_DIR}/fourier_transform/CosFourier.cpp ${SRC_DIR}/fourier_transform/CosFourier.hpp
	${CXX} ${CXXFLAGS} -c ${SRC_DIR}/fourier_transform/CosFourier.cpp -o ${BUILD_DIR}/CosFourier.o

${BUILD_DIR}/InverseCosFourier.o : ${SRC_DIR}/fourier_transform/InverseCosFourier.cpp ${SRC_DIR}/fourier_transform/InverseCosFourier.hpp
	${CXX} ${CXXFLAGS} -c ${SRC_DIR}/fourier_transform/InverseCosFourier.cpp -o ${BUILD_DIR}/InverseCosFourier.o

${BUILD_DIR}/SinFourier.o : ${SRC_DIR}/fourier_transform/SinFourier.cpp ${SRC_DIR}/fourier_transform/SinFourier.hpp
	${CXX} ${CXXFLAGS} -c ${SRC_DIR}/fourier_transform/SinFourier.cpp -o ${BUILD_DIR}/SinFourier.o

${BUILD_DIR}/InverseSinFourier.o : ${SRC_DIR}/fourier_transform/InverseSinFourier.cpp ${SRC_DIR}/fourier_transform/InverseSinFourier.hpp
	${CXX} ${CXXFLAGS} -c ${SRC_DIR}/fourier_transform/InverseSinFourier.cpp -o ${BUILD_DIR}/InverseSinFourier.o

${BUILD_DIR}/RKNystroem.o : ${SRC_DIR}/solver/RKNystroem.cpp ${SRC_DIR}/solver/RKNystroem.hpp
	${CXX} ${CXXFLAGS} -c ${SRC_DIR}/solver/RKNystroem.cpp -o ${BUILD_DIR}/RKNystroem.o

${BUILD_DIR}/RKV43.o : ${SRC_DIR}/solver/RKV43.cpp ${SRC_DIR}/solver/RKV43.hpp
	${CXX} ${CXXFLAGS} -c ${SRC_DIR}/solver/RKV43.cpp -o ${BUILD_DIR}/RKV43.o

${BUILD_DIR}/VectorArray2D.o : ${SRC_DIR}/templates/VectorArray2D.cpp ${SRC_DIR}/templates/VectorArray2D.hpp
	${CXX} ${CXXFLAGS} -c ${SRC_DIR}/templates/VectorArray2D.cpp -o ${BUILD_DIR}/VectorArray2D.o

${BUILD_DIR}/Utilities.o : ${SRC_DIR}/utilities/Utilities.cpp ${SRC_DIR}/utilities/Utilities.hpp
	${CXX} ${CXXFLAGS} -c ${SRC_DIR}/utilities/Utilities.cpp -o ${BUILD_DIR}/Utilities.o

${BUILD_DIR}/ImageDifference.o : ${SRC_DIR}/utilities/ImageDifference.cpp ${SRC_DIR}/utilities/ImageDifference.hpp
	${CXX} ${CXXFLAGS} -c ${SRC_DIR}/utilities/ImageDifference.cpp -o ${BUILD_DIR}/ImageDifference.o

${BUILD_DIR}/BracketMethod.o : ${SRC_DIR}/utilities/BracketMethod.cpp ${SRC_DIR}/utilities/BracketMethod.hpp
	${CXX} ${CXXFLAGS} -c ${SRC_DIR}/utilities/BracketMethod.cpp -o ${BUILD_DIR}/BracketMethod.o

clean :
	rm build/*

