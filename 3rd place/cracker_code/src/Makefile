# Build targets
QUERY1_O=./src/query1.cpp
QUERY1_SEQ_O=./src/query1-seq.cpp
QUERY2_O=./src/query2.cpp
QUERY3_O=./src/query3.cpp 
QUERY4_O=./src/query4.cpp
PRE_O=./src/splitfile.cpp
POST_O=./src/mergeans.cpp
RUN_O=./src/runquery.cpp

USER_FLAGS=
USER_LIBS=

# Compiler flags
CC  = gcc
CXX = g++
CFLAGS=$(USER_FLAGS) -O3 -fPIC -Wall -g -I. -I./include
CXXFLAGS=$(CFLAGS)
LDFLAGS=-lpthread

# The name of the excution that will be built
QUERY1=query1
QUERY1_SEQ=query1-seq
QUERY2=query2
QUERY3=query3
QUERY4=query4
PREPROCESS=splitfile
POSTPROCESS=mergeans
RUNQUERY=runquery

all: $(QUERY1) $(QUERY1_SEQ) $(QUERY2) $(QUERY3) $(QUERY4) $(PREPROCESS) $(POSTPROCESS) $(RUNQUERY)

$(QUERY1):
	$(CXX) $(CXXFLAGS)  -pthread $(LDFLAGS) -o $(QUERY1) $(QUERY1_O) $(USER_LIBS)

$(QUERY1_SEQ):
	$(CXX) $(CXXFLAGS)  -pthread $(LDFLAGS) -o $(QUERY1_SEQ) $(QUERY1_SEQ_O) $(USER_LIBS)

$(QUERY2):
	$(CXX) $(CXXFLAGS)  -pthread $(LDFLAGS) -o $(QUERY2) $(QUERY2_O) $(USER_LIBS)
	
$(QUERY3):
	$(CXX) $(CXXFLAGS)  -pthread $(LDFLAGS) -o $(QUERY3) $(QUERY3_O) $(USER_LIBS)
	
$(QUERY4):
	$(CXX) $(CXXFLAGS)  -pthread $(LDFLAGS) -o $(QUERY4) $(QUERY4_O) $(USER_LIBS)
	
$(PREPROCESS):
	$(CXX) $(CXXFLAGS)  -pthread $(LDFLAGS) -o $(PREPROCESS) $(PRE_O) $(USER_LIBS)
	
$(POSTPROCESS):
	$(CXX) $(CXXFLAGS)  -pthread $(LDFLAGS) -o $(POSTPROCESS) $(POST_O) $(USER_LIBS)

$(RUNQUERY):
	$(CXX) $(CXXFLAGS)  -pthread $(LDFLAGS) -o $(RUNQUERY) $(RUN_O) $(USER_LIBS)
	
clean:
	rm -f $(QUERY1) $(QUERY1_SEQ) $(QUERY2) $(QUERY3) $(QUERY4) $(PREPROCESS) $(POSTPROCESS) $(RUNQUERY)
	find . -name '*.o' -print | xargs rm -f
