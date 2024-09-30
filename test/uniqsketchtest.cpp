#include <fstream>
#include <iterator>
#include <algorithm>
#include <string>
#include <vector>
#include <cassert>

#include "SequenceUtil.hpp"
#include "TestUtil.hpp"
#include "UniqSketchUtil.hpp"

// check if identifyUniqKmers function correctly computes unique kmers from fatsa references ref_test_1 and ref_test_1
void testIdentifyUniqKmers(const std::vector<std::string> &refFiles) {
    std::cerr << "START: Test identifyUniqKmers... \t";

    /*
    ref_test_1 has two sequences {CAGGCTACAAA, GAGACTGAACT} of lenght 11bp with total four kmers of 10bp
    all_kmer_ref_test_1 = {CAGGCTACAA, AGGCTACAAA, GAGACTGAAC, AGACTGAACT}

    ref_test_2 has two sequences {CTTAGCTGAGG, CAGGCTACAAC} of lenght 11bp with total four kmers of 10bp
    all_kmer_ref_test_2 = {CTTAGCTGAG, TTAGCTGAGG, CAGGCTACAA, AGGCTACAAC}

    kmer CAGGCTACAA is available in both referencers and hence is not unique, so the canonical unique sets are:
    uniq_kmer_ref_test_1 = {AGGCTACAAA, GAGACTGAAC, AGACTGAACT}
    uniq_kmer_ref_test_2 = {CTCAGCTAAG, CCTCAGCTAA, AGGCTACAAC}
    */
    identifyUniqKmers(refFiles);

    assert(compareFiles("ref_test_1.tsv", "expected_ref_test_1.tsv")==true);
    assert(compareFiles("ref_test_2.tsv", "expected_ref_test_2.tsv")==true);

    std::cerr << "PASSED: Test identifyUniqKmers" << std::endl;
}

// check if buildSketch function generates the expected unique sketch output
void testBuildSketch() {
    std::cerr << "START: Test buildSketch... \t";
    
    /*
    each reference has three uniqe kmers:
    uniq_kmer_ref_test_1 = {AGGCTACAAA, GAGACTGAAC, AGACTGAACT}
    uniq_kmer_ref_test_2 = {CTCAGCTAAG, CCTCAGCTAA, AGGCTACAAC}
    
    and given opt::sketchnum = 2, each ref will be represented by 2 unique kmers uniformly selected
    uniq_sketch_ref_test_1 = {AGACTGAACT, AGGCTACAAA}
    uniq_sketch_ref_test_2 = {AGGCTACAAC, CTCAGCTAAG} 
    
    so the total number of uniqe kmers in output should be 4:
    uniq_sketch = {AGACTGAACT, AGGCTACAAA, AGGCTACAAC, CTCAGCTAAG} 
    */
    buildSketch("expected_db_uniq_count.tsv");
    assert(compareFiles(opt::outfile, "expected_out_uniqsketch.txt")==true);
    
    std::cerr << "PASSED: Test buildSketch" << std::endl;
}

// check if lowComplexity function correctly identifies low complexity signatures
void testLowComplexity() {
    opt::maxEntropy = 24.56;

    std::cerr << "START: Test lowComplexity... \t";

    std::string signature_1("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCT");
    std::string signature_2("CGTTGCATGGCATGCGGTTCGTCTGCGCGTGCGTGCTGATGCAGTAAAAAAGCTGCATAAAATCGATCGAAAGTACATGCG");

    assert(lowComplexity(signature_1)==true);
    assert(lowComplexity(signature_2)==false);

    std::cerr << "PASSED: Test lowComplexity" << std::endl;
}

// check if get_canonical function correctly computes canonical kmer
void testGetCanonical() {
    std::cerr << "START: Test get_canonical... \t";

    std::string kmer_1("GTTTATTTACC");
    std::string kmer_2("AGGTTAGCTTTC");

    assert(get_canonical(kmer_1)=="GGTAAATAAAC");
    assert(get_canonical(kmer_2)=="AGGTTAGCTTTC");
                                   
    std::cerr << "PASSED: Test get_canonical" << std::endl;
}

// check if getBaseId function correctly extracts basename
void testGetBaseId() {
    std::cerr << "START: Test getBaseId... \t";

    std::string fPath_1("/home/user/workspace/test/ref1.fasta");
    std::string fPath_2("ref2.fa");

    assert(getBaseId(fPath_1)=="ref1");
    assert(getBaseId(fPath_2)=="ref2");

    std::cerr << "PASSED: Test getBaseId" << std::endl;
}

int main() {
    opt::kmerLen = 10;
    opt::bits = 128;
    opt::dbfSize = 7;
    opt::sbfSize = 1;
    opt::m1 = opt::bits*opt::dbfSize;
    opt::m2 = opt::bits*opt::sbfSize;
    opt::numRef = 2;
    opt::sketchnum = 2;
    opt::outfile = "test_out_uniqsketch.txt";
    opt::outdir = ".";

    std::vector<std::string> refFiles;
    refFiles.push_back("ref_test_1.fa");
    refFiles.push_back("ref_test_2.fa");

    testIdentifyUniqKmers(refFiles);
    testBuildSketch();
    testLowComplexity();
    testGetCanonical();
    testGetBaseId();

    std::cerr << "UniqSketch: All tests PASSED!" << std::endl << std::endl;
    return 0;
}
