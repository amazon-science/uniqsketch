# Margin-based signature selection (Hamming-distance hardening against read errors)

**Status:** Design note / proposal (candidate for v2.0). Not yet implemented.

## Motivation

When selecting signature candidates for a reference, we want to prefer signatures
that are at least **2 bases away** from every k-mer in other references. This
closes a false-positive channel driven by sequencing noise: a single substitution
error in a read can flip one base of a foreign k-mer so that it collides with one
of our signatures.

## The property we actually want

UniqSketch signatures are already *unique* — Hamming distance ≥ 1 from every other
reference's k-mers, by the solid-Bloom-filter construction. So "≥ 2 bases from
cross-references" reduces to a single condition:

> A signature `S` of reference A is **2-safe** iff **no other reference contains a
> k-mer at Hamming distance exactly 1** from `S`.

Why this removes the noise channel: if reference B contains `T` with
`Hamming(S, T) = 1`, then a read truly from B carrying `T`, plus one substitution
error at that single differing base, becomes `S` → false hit for A. Require
margin ≥ 2 and a *single* error can no longer move any foreign k-mer onto `S`; you
would need ≥ 2 coincident errors at specific positions
(≈ `error_rate² · (1/3)²` per occurrence), which is orders of magnitude rarer.

This is minimum-distance / error-correcting-code design: we pick signatures whose
minimum distance to the *foreign* k-mer set is ≥ `2t + 1`, with `t = 1`
correctable read error.

**Scope:** this models **substitution** errors (the dominant Illumina mode).
Indels would require edit distance, which is much costlier and is out of scope for
now. The margin filter also stacks with the existing solid-read filter, which
already suppresses low-coverage single-error k-mers on the query side.

## One implementation reality first

The Bloom filters (`dbFilter`, `sbFilter`) store only **ntHash values**, so we
cannot enumerate Hamming neighbors on the BF itself. We do not need to: `checkRef`
already writes the actual **canonical k-mer strings** to the per-reference
candidate files, so the sequences are in hand at selection time. And ntHash is
canonical, so hashing a mutated neighbor string and querying the BF is
automatically reverse-complement-correct — no manual canonicalization needed.

Distance primitive, using the existing `b2f` 2-bit encoding (an 81-mer packs into
3×64-bit words):

```cpp
// number of differing bases between two 2-bit-packed k-mers
inline unsigned hamming2bit(const uint64_t* a, const uint64_t* b, int words) {
    unsigned d = 0;
    for (int w = 0; w < words; ++w) {
        uint64_t x = a[w] ^ b[w];                     // bits differ where bases differ
        x = (x | (x >> 1)) & 0x5555555555555555ULL;   // collapse to 1 bit/base
        d += __builtin_popcountll(x);
    }
    return d;
}
```

## Method A — neighbor enumeration vs. the distinct BF (recommended default)

For k = 81 there are only `81 × 3 = 243` Hamming-1 neighbors. Generate them and
query the **distinct BF** (all universe k-mers):

```cpp
// true  => S has NO neighbor present anywhere in the universe  => 2-safe
bool isTwoSafe(const std::string& S, const BloomFilter& dbFilter) {
    std::string m = S;
    uint64_t h[MAXHASH];                                  // MAXHASH >= nhash1
    for (size_t i = 0; i < S.size(); ++i) {
        char orig = m[i];
        for (char b : {'A', 'C', 'G', 'T'}) {
            if (b == orig) continue;
            m[i] = b;
            NTMC64(m.c_str(), opt::kmerLen, opt::nhash1, h);   // canonical hash
            if (dbFilter.contains(h)) { m[i] = orig; return false; }
        }
        m[i] = orig;
    }
    return true;
}
```

Three things make this work well here:

- **BF false positives are a non-issue.** Naively, 243 lookups amplify FPR, but the
  *index* distinct BF uses `bits = 128`/element with `nhash1 = 5` → per-query FPR
  ≈ 8e-8, so P(any of 243 spuriously hits) ≈ 2e-5. (This is very different from
  `comparesketch`'s 16-bit BF, where this would *not* be safe.)
- **Be deliberately conservative: do not subtract A's own k-mers.** If a neighbor
  occurs in *both* A and B, `S` is still genuinely unsafe, and a presence-only BF
  cannot distinguish "A only" from "A and B." Treating *any* present neighbor as
  unsafe never keeps a bad signature; it only occasionally drops a good one (when
  the neighbor lives solely inside A, which requires a near-repeat in A and is rare
  for 81-mers). Erring toward rejection is the safe direction for FPR.
- **Do it lazily.** Only test candidates actually considered in `getUniqSet`'s
  spacing slots, so cost is ~`sketchnum × oversample × 243`, not all unique k-mers.
  This requires plumbing `dbFilter` through to the selection stage (today it is
  local to `identifyUniqKmers`).

Optimization: the vendored `nthash.hpp` does **not** expose a named
single-substitution update — its rolling functions only slide the window (drop the
leftmost base, append a rightmost one). But because ntHash is linear over XOR and
the header exports its shifted-seed tables, an O(1) arbitrary-position substitution
is straightforward to build; see the reference implementation below.

### O(1) substitution update (reference implementation)

ntHash exposes only full recompute (`NTF64`/`NTR64`/`NTC64`/`NTMC64(seq, k, …)`)
and window-slide updates (`…(hVal, k, charOut, charIn)` and the left-slide `…L`
variants). None take an internal position, so the library API alone forces a full
`NTMC64` recompute per neighbor.

It is implementable in O(1) anyway. In `NTF64`, base `i` contributes
`srol^(k-1-i)(seedTab[base])`, and those shifted seeds are exactly the header's
precomputed `msTab31l` / `msTab33r` tables (the slide functions use them to remove
the outgoing base). A substitution is therefore: XOR out the old base's shifted
seed and XOR in the new one, on both strands, then recombine. This reproduces
`NTMC64` bit-for-bit — canonical is `fhVal + rhVal` (this ntHash sums the two
strands; it does not take the min), and `NTE64` regenerates hashes `1..m-1`.

```cpp
#include "nthash.hpp"   // exports seedTab, msTab31l, msTab33r, cpOff, NTC64, NTE64

// forward contribution of base c at position i:  srol^(k-1-i)(seedTab[c])
static inline uint64_t fContrib(unsigned char c, unsigned i, unsigned k) {
    unsigned s = k - 1 - i;
    return msTab31l[c][s % 31] | msTab33r[c][s % 33];
}
// reverse-strand contribution: complemented base (c & cpOff), shift i
static inline uint64_t rContrib(unsigned char c, unsigned i, unsigned k) {
    return msTab31l[c & cpOff][i % 31] | msTab33r[c & cpOff][i % 33];
}

// Fill hVal[m] for the k-mer with position i changed oldc -> newc; no O(k) rehash.
static inline void subMultiHash(uint64_t fhVal, uint64_t rhVal,
                                unsigned i, unsigned char oldc, unsigned char newc,
                                unsigned k, unsigned m, uint64_t* hVal) {
    uint64_t f = fhVal ^ fContrib(oldc, i, k) ^ fContrib(newc, i, k);
    uint64_t r = rhVal ^ rContrib(oldc, i, k) ^ rContrib(newc, i, k);
    uint64_t bVal = f + r;                       // canonical = fhVal + rhVal
    hVal[0] = bVal;
    for (unsigned j = 1; j < m; j++) hVal[j] = NTE64(bVal, k, j);
}
```

The faster `isTwoSafe` computes the base hash once, then walks the 243 neighbors
without rehashing:

```cpp
bool isTwoSafe(const std::string& S, const BloomFilter& dbFilter) {
    uint64_t fhVal, rhVal;
    NTC64(S.c_str(), opt::kmerLen, fhVal, rhVal);          // O(k), once per candidate
    uint64_t h[MAXHASH];                                    // MAXHASH >= nhash1
    for (unsigned i = 0; i < opt::kmerLen; ++i)
        for (char nc : {'A', 'C', 'G', 'T'})
            if ((unsigned char)nc != (unsigned char)S[i]) {
                subMultiHash(fhVal, rhVal, i, S[i], nc, opt::kmerLen, opt::nhash1, h);
                if (dbFilter.contains(h)) return false;
            }
    return true;
}
```

Cost drops from `243 × O(k)` to one `O(k)` base hash plus `243 × O(m)` (~15× fewer
word-ops at k = 81, m = 5); the "remove old base" term can be factored once per
position for a little more. The 243 BF probes are unchanged and usually dominate,
so this is worth the ~15 lines at scale (very-sensitive, thousands of references);
when Method A runs lazily on a handful of slot candidates, the plain `NTMC64`
recompute shown earlier is correct and likely fast enough.

Note: `seedTab`, `msTab31l`, `msTab33r`, and `cpOff` have internal linkage in the
header, so this must live in a translation unit that includes `nthash.hpp`
(e.g. `UniqSketchUtil.hpp`). Correctness of the tables was verified against the
vendored `nthash.hpp` (`NTF64` at line 224, canonical `NTC64` = `fhVal + rhVal` at
line 264, `NTMC64` expansion at line 361, `NTE64` at line 341).

## Method B — pigeonhole blocking (exact, when the true margin matters)

If conservative over-rejection costs too many signatures on close strains, compute
the margin exactly. For distance ≤ 1, split each k-mer into 2 blocks
(⌈k/2⌉, ⌊k/2⌋); any two k-mers within distance 1 must share ≥ 1 identical block.
Index every k-mer by `(block_id, block_value) → refId`, then for each candidate
look up its 2 blocks, gather colliders from *other* references, and confirm with
`hamming2bit`. This yields the exact foreign margin and knows the reference
identity, so there is no over-rejection.

The scaling trick: **only run it within clusters / among references above a
similarity threshold.** Cross-species 81-mer near-collisions are astronomically
unlikely; the danger is entirely among close strains, which `comparesketch
--cluster` already groups. That bounds an otherwise universe-wide join to small
neighborhoods.

## Turning detection into "optimal" selection

Detection yields a per-candidate margin `m(S)` (or a 2-safe boolean). Fold it into
`getUniqSet` as a **tiered greedy** choice rather than a hard filter:

1. Within each evenly-spaced slot, prefer a candidate that is non-low-complexity
   **and** margin ≥ 2.
2. Fall back to a margin-1 candidate only if the slot has no 2-safe option, so
   `sketchnum` and the genomic spacing are preserved.
3. Expose `--min-margin=N` (default 1 = current behavior; 2 = the hardening),
   mirroring `--entropy`, and report per-reference how many signatures met the
   margin as a sensitivity diagnostic.

Because signatures are scored independently, this greedy per-slot choice is
effectively optimal for the objective (maximize the number of high-margin,
well-spaced, complex signatures); a full ILP would be overkill.

## The tradeoff to watch

For near-identical strains, many unique k-mers exist *because* they sit one base
from a sibling strain, so a hard margin-≥2 rule can starve those references and
create false negatives. Therefore the margin should be (a) a soft preference with
fallback, and (b) applied *after* clustering collapses the near-identical strains.

Viewed this way it is a third hardening axis, complementary to the two in the
ablation study:

- **Clustering** → accurate quantification of closely related strains.
- **Entropy / low-complexity filter** → fewer cross-reactions from repeat regions.
- **Margin selection** → fewer cross-reactions from single-base read errors.
