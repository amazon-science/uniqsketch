#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <string>
#include <vector>
#include <cassert>

#include "BloomFilter.hpp"


namespace opt {
size_t size;
unsigned hash;
unsigned k;
}

// initialize Bloom filter parameters
void initialize() {
    opt::size = 65536;
    opt::hash = 3;
    opt::k = 32; 
}

// check if insert and contain functions correctly work on Bloom filter
void testInsertContain() {
    std::cerr << "START: Test insert, contain... \t";

    BloomFilter bFilter(opt::size, opt::hash, opt::k);
    uint64_t hVec_1[] = {1, 3000, 50000000};
    uint64_t hVec_2[] = {2000000, 700000000, 4000000000};

    bFilter.insert(hVec_1);

    assert(bFilter.contains(hVec_1)==true);
    assert(bFilter.contains(hVec_2)==false);


    std::string kmer_1("GCCTAGCTAGCTTTTAGCTGGGATTTTTT");
    std::string kmer_2("ATTTAGCTAGCGCTAAGCTGGGACCCTGC");

    bFilter.insert(kmer_1.c_str());
    
    assert(bFilter.contains(kmer_1.c_str())==true);
    assert(bFilter.contains(kmer_2.c_str())==false);

    std::cerr << "PASSED: Test insert, contain" << std::endl;   
}

// check if insert_make_change correctly works on Bloom filter
void testInsertMakeChange() {
    std::cerr << "START: Test insert_make_change... \t";

    BloomFilter bFilter(opt::size, opt::hash, opt::k);
 
    uint64_t hVec_1[] = {34987645632, 89, 37692833};
    uint64_t hVec_2[] = {34359738368, 4398046511238, 678};

    bFilter.insert(hVec_1);

    bool query_1 = bFilter.insert_make_change(hVec_1);
    bool query_2 = bFilter.insert_make_change(hVec_2);

    assert(query_1==false);
    assert(query_2==true);

    std::cerr << "PASSED: Test insert_make_change" << std::endl;   
}

int main() {

    initialize();
    testInsertContain();
    testInsertMakeChange();

    std::cerr << "BloomFilter: All tests PASSED!" << std::endl << std::endl;
    return 0;
}

