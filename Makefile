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
	$(CXX) $(INCLUDEPATH) $(CPPFLAGS) -o test/$@ $^

test: uniqsketchtest querysketchtest bloomfiltertest
	cd test; ./uniqsketchtest > /dev/null; ./querysketchtest > /dev/null; ./bloomfiltertest > /dev/null

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