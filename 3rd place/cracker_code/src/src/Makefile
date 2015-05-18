# Build targets
QUERY4_V1_O=./src/icde14_delta_pfs.cpp ./src/time.cpp
QUERY4_V2_O=./src/final_exact_2-hop.cpp ./src/time.cpp
QUERY4_V3_O=./src/random_okamoto_alg.cpp ./src/time.cpp
QUERY4_V4_O=./src/random+delta_pfs.cpp ./src/time.cpp
QUERY4_V5_O=./src/query4-2-hop.cpp ./src/time.cpp


USER_FLAGS=
USER_LIBS=

# Compiler flags
CC  = gcc
CXX = g++
CFLAGS=$(USER_FLAGS) -O3 -fPIC -Wall -g -I. -I./include
CXXFLAGS=$(CFLAGS)
LDFLAGS=-lpthread

# The name of the excution that will be built
QUERY4_V1=icde14_delta_pfs
QUERY4_V2=final_exact_2-hop
QUERY4_V3=random_okamoto_alg
QUERY4_V4=random+delta_pfs
QUERY4_V5=query4-2-hop

all: $(QUERY4_V1) $(QUERY4_V2) $(QUERY4_V3) $(QUERY4_V4) $(QUERY4_V5)


$(QUERY4_V1):
	$(CXX) $(CXXFLAGS)  -pthread $(LDFLAGS) -o $(QUERY4_V1) $(QUERY4_V1_O) $(USER_LIBS)

$(QUERY4_V2):
	$(CXX) $(CXXFLAGS)  -pthread $(LDFLAGS) -o $(QUERY4_V2) $(QUERY4_V2_O) $(USER_LIBS)

$(QUERY4_V3):
	$(CXX) $(CXXFLAGS)  -pthread $(LDFLAGS) -o $(QUERY4_V3) $(QUERY4_V3_O) $(USER_LIBS)
	
$(QUERY4_V4):
	$(CXX) $(CXXFLAGS)  -pthread $(LDFLAGS) -o $(QUERY4_V4) $(QUERY4_V4_O) $(USER_LIBS)

$(QUERY4_V5):
	$(CXX) $(CXXFLAGS)  -pthread $(LDFLAGS) -o $(QUERY4_V5) $(QUERY4_V5_O) $(USER_LIBS)

clean:
	rm -f $(QUERY4_V1) $(QUERY4_V2) $(QUERY4_V3) $(QUERY4_V4) $(QUERY4_V5)
	find . -name '*.o' -print | xargs rm -f
