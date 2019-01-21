CXX=$(BM_ROOT)/local/bin/g++
CXXFLAGS=-std=c++17 -I$(BM_ROOT)/local/include
LDFLAGS=-L$(BM_ROOT)/local/lib -L$(BM_ROOT)/local/lib64
LIBS=-lmm -ldl

pagerank: main.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

clean:
	rm -f pagerank
