override CXXFLAGS += -I ../.. -I/work/carlinet/boost_1_55_0  -std=c++11 -W -Wall \
					-DBOOST_RESULT_OF_USE_DECLTYPE
LDFLAGS=-lfreeimage

all: g2 satmaxtree-cli


g2-maxdepth: g2-maxdepth.o routines.o compute_g2.o
g2: g2.o compute_g2.o routines.o
satmaxtree-cli: satmaxtree-cli.o satmaxtree.o

g2.o: reconstruct.hpp remove_parent_relation.hpp