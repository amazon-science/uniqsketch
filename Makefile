CXX=g++ -std=c++17
# -pthread: the query engine parses reads on std::threads (QuerySketchUtil.hpp).
CPPFLAGS=-Ivendor -pthread
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

# Runs the partitioned group scan against the serial one, so it needs OpenMP.
uniqsketchtest: test/uniqsketchtest.cpp
	$(CXX) $(INCLUDEPATH) $(CPPFLAGS) -fopenmp -o test/$@ $^ $(LIBPATH)

querysketchtest: test/querysketchtest.cpp
	$(CXX) $(INCLUDEPATH) $(CPPFLAGS) -o test/$@ $^ $(LIBPATH)

bloomfiltertest: test/bloomfiltertest.cpp
	$(CXX) $(INCLUDEPATH) $(CPPFLAGS) -pthread -o test/$@ $^

# Checks the row order of identifyDifference on several threads, so it needs OpenMP.
comparesketchtest: test/comparesketchtest.cpp
	$(CXX) $(INCLUDEPATH) $(CPPFLAGS) -fopenmp -o test/$@ $^ $(LIBPATH)

# Builds a within-group sketch in-process (parallel job lists), so it needs OpenMP.
treesketchtest: test/treesketchtest.cpp
	$(CXX) $(INCLUDEPATH) $(CPPFLAGS) -fopenmp -o test/$@ $^ $(LIBPATH)

headerstest: test/headerstest.cpp
	$(CXX) $(INCLUDEPATH) $(CPPFLAGS) -o test/$@ $^ $(LIBPATH)

# End-to-end: builds an index and queries it in-process, so it needs OpenMP like
# the real tools (the pipeline functions carry omp pragmas).
e2etest: test/e2etest.cpp
	$(CXX) $(INCLUDEPATH) $(CPPFLAGS) -fopenmp -o test/$@ $^ $(LIBPATH)

# Rule-stack test of the tree query, through the binaries: runTreeQuery lives in
# src/querysketch.cpp, where no in-process suite can call it. Generates 16
# synthetic genomes and their reads, builds a tree index and checks 69 queries
# against test/tree_rules/expected.tsv (a few seconds; needs python3).
# The two tools are file targets here, so they are rebuilt only when a source is
# newer than the binary (the uniqsketch and querysketch targets always compile).
bin/%: src/%.cpp $(wildcard include/*.hpp) $(wildcard vendor/*)
	$(MAKE) bin_dir $*

treerulestest: bin/uniqsketch bin/querysketch
	python3 test/tree_rules/run_tree_rules.py bin test/tree_rules_work

# Suites are chained with && so the first failing one fails the target; joined
# with ';' (as before) a failed assertion left `make test` exiting 0.
test: uniqsketchtest querysketchtest bloomfiltertest comparesketchtest treesketchtest headerstest e2etest treerulestest
	cd test && ./uniqsketchtest > /dev/null && ./querysketchtest > /dev/null && ./bloomfiltertest > /dev/null && ./comparesketchtest > /dev/null && ./treesketchtest > /dev/null && ./headerstest > /dev/null && ./e2etest > /dev/null

clean:
	rm bin/uniqsketch bin/querysketch bin/comparesketch 

test_clean: 
	cd test && \
	rm uniqsketchtest \
	querysketchtest \
	bloomfiltertest \
	ref_test_1.cand expected_ref_test_1.cand \
	ref_test_2.cand expected_ref_test_2.cand \
	test_out_uniqsketch.txt \
	out_querysketch.tsv