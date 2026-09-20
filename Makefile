CXX=g++ -std=c++17
CPPFLAGS=-Ivendor
INCLUDEPATH=-Iinclude
LIBPATH=-lz
OPTFLAGS=-O3 -fopenmp

all:bin_dir uniqsketch querysketch comparesketch

uniqsketch_SOURCES = src/uniqsketch.cpp

querysketch_SOURCES = src/querysketch.cpp

comparesketch_SOURCES = src/comparesketch.cpp

bin_dir:
	mkdir -p bin

uniqsketch:$(uniqsketch_SOURCES)
	$(CXX) $(OPTFLAGS) $(CPPFLAGS) $(INCLUDEPATH) -o bin/$@ $^ $(LIBPATH)

querysketch:$(querysketch_SOURCES)
	$(CXX) $(OPTFLAGS) $(CPPFLAGS) $(INCLUDEPATH) -o bin/$@ $^ $(LIBPATH)

comparesketch:$(comparesketch_SOURCES)
	$(CXX) $(OPTFLAGS) $(CPPFLAGS) $(INCLUDEPATH) -o bin/$@ $^ $(LIBPATH)

uniqsketchtest: test/uniqsketchtest.cpp
	$(CXX) $(INCLUDEPATH) $(CPPFLAGS) -o test/$@ $^

querysketchtest: test/querysketchtest.cpp
	$(CXX) $(INCLUDEPATH) $(CPPFLAGS) -o test/$@ $^ $(LIBPATH)

bloomfiltertest: test/bloomfiltertest.cpp
	$(CXX) $(INCLUDEPATH) $(CPPFLAGS) -pthread -o test/$@ $^

comparesketchtest: test/comparesketchtest.cpp
	$(CXX) $(INCLUDEPATH) $(CPPFLAGS) -o test/$@ $^ $(LIBPATH)

treesketchtest: test/treesketchtest.cpp
	$(CXX) $(INCLUDEPATH) $(CPPFLAGS) -o test/$@ $^ $(LIBPATH)

headerstest: test/headerstest.cpp
	$(CXX) $(INCLUDEPATH) $(CPPFLAGS) -o test/$@ $^ $(LIBPATH)

# End-to-end: builds an index and queries it in-process, so it needs OpenMP like
# the real tools (the pipeline functions carry omp pragmas).
e2etest: test/e2etest.cpp
	$(CXX) $(INCLUDEPATH) $(CPPFLAGS) -fopenmp -o test/$@ $^ $(LIBPATH)

# Suites are chained with && so the first failing one fails the target; joined
# with ';' (as before) a failed assertion left `make test` exiting 0.
test: uniqsketchtest querysketchtest bloomfiltertest comparesketchtest treesketchtest headerstest e2etest
	cd test && ./uniqsketchtest > /dev/null && ./querysketchtest > /dev/null && ./bloomfiltertest > /dev/null && ./comparesketchtest > /dev/null && ./treesketchtest > /dev/null && ./headerstest > /dev/null && ./e2etest > /dev/null

clean:
	rm bin/uniqsketch bin/querysketch bin/comparesketch 

test_clean: 
	cd test && \
	rm uniqsketchtest \
	querysketchtest \
	bloomfiltertest \
	ref_test_1.tsv \
	ref_test_2.tsv \
	test_out_uniqsketch.txt \
	out_querysketch.tsv