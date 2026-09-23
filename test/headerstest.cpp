// Every *Util header in a single translation unit.
//
// Until the config refactor this file could not compile: each header defined
// the same file-scope globals (threads, bits, kmerLen, ...) inside `namespace
// opt`, so including two of them produced redefinition errors. That is why each
// tool was a single .cpp and nothing could be tested across tool boundaries.
//
// Nothing here is asserted at runtime. The test is that this file compiles and
// links; if someone reintroduces a file-scope definition in one of these
// headers, the build breaks here first.
#include <iostream>

#include "UniqSketchUtil.hpp"
#include "QuerySketchUtil.hpp"
#include "CompareSketchUtil.hpp"
#include "TreeSketchUtil.hpp"
#include "GroupSketchUtil.hpp"
#include "ClusterUtil.hpp"

int main() {
    std::cerr << "START: Test headers compose in one TU... \t";

    // Instantiate one config per tool side by side, and exercise a derived
    // member from each so the structs are genuinely distinct types in scope.
    UniqSketchConfig uniq;
    QuerySketchConfig query;
    CompareSketchConfig cmp;

    const bool ok = uniq.kmerLen == 81 && uniq.m1() == 0        // dbfSize defaults to 0
                 && query.sketchHit == 10 && query.kmerLen == 0  // set by loadSketch, not a default
                 && cmp.nhash == 3 && cmp.memBudgetGB == 64;

    if (!ok) {
        std::cerr << "FAILED: unexpected config defaults\n";
        return 1;
    }
    std::cerr << "PASSED: Test headers compose in one TU\n";
    std::cerr << "Headers: All tests PASSED!\n\n";
    return 0;
}
