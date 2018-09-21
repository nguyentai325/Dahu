override CXXFLAGS += -std=c++0x -I ../.. -W -Wall
override LDFLAGS += -lfreeimage -ltbb -L ~/lib

all: ../lib/libtos.so ../lib/libtos_debug.so ../lib/libtos.a
	cp ../lib/libtos.so ~/lib
	cp ../lib/libtos_debug.so ~/lib

bin: marginalfilter simplify thicken lexico_grainfilter k12k0 rgbgrainfilter tos grainfilter
	[ -e bin ] || mkdir bin
	cp $^ bin


bin: CXXFLAGS += -O3 -DNDEBUG -march=native


pgm_objdetection_mser: LDFLAGS += -ltos
ppm_objdetection_mser: LDFLAGS += -ltos
ppm_objdetection_mser: ppm_objdetection_mser.o colorToSGrad.o
pgm_objdetection: LDFLAGS += -ltos
ppm_objdetection: LDFLAGS += -ltos
ppm_objdetection: ppm_objdetection.o colorToSGrad.o


objdetection: objdetection.o colorToSGrad.o
distances: distances.o colorToSGrad.o

tos: LDFLAGS += -lboost_program_options


../lib/libtos.so: libtos.hpp routines.hpp
../lib/libtos_debug.so: libtos.hpp routines.hpp

../lib/libtos.so: libtos.cpp routines.cpp
	g++ -shared -fPIC $(CXXFLAGS) $^ -o ../lib/libtos.so

../lib/libtos_debug.so: libtos.cpp routines.cpp
	g++ -shared -fPIC $(CXXFLAGS) $^ -o ../lib/libtos_debug.so

../lib/libtos.a: libtos.o routines.o
	ar -cvq ../lib/libtos.a  $^

../lib/libtos.so: CXXFLAGS += -DNDEBUG -O2
../lib/libtos_debug.so: CXXFLAGS += -g -ggdb -O1

../lib/libtos.a: CXXFLAGS += -DNDEBUG -O2