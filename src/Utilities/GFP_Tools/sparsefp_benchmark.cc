// Benchmarking to examine performance characteristics of things supporting Sparse_Fingerprints.

#include <cstdint>
#include <random>
#include <vector>

#include <benchmark/benchmark.h>

#include "Foundational/iwmisc/sparse_fp_creator.h"

#include "Utilities/GFP_Tools/sparsefp.h"

namespace {

void
BM_BinarySearch(benchmark::State& state) {
  std::default_random_engine rng;
  int nbits = state.range(0);
  std::uniform_int_distribution<int> count_distribution(1, 255);

  Sparse_Fingerprint_Creator sfc;
  for (int i = 0; i < nbits; ++i) {
    uint32_t b = rng();
    int count = count_distribution(rng);
    sfc.hit_bit(b, count);
  }
  Sparse_Fingerprint sfp;
  sfp.build_from_sparse_fingerprint_creator(sfc);
  int missing = 0;
  for (auto _ : state) {
    uint32_t b = rng();
    int c = sfp.is_set(b);
    if (c == 0) {
      missing++;
    }
  }
}
BENCHMARK(BM_BinarySearch)
->DenseRange(20, 180, 10);

}  // namespace

BENCHMARK_MAIN();
