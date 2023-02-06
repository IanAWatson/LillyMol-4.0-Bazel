// Question of whether integer or float properties are faster

#include <random>

#include <benchmark/benchmark.h>

#include "Utilities/GFP_Tools/gfp.h"
#include "Utilities/GFP_Tools/tanimoto_float.h"

namespace {

void
BM_TestInt(benchmark::State& state) {
  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_int_distribution<int> property(0, 255);

  constexpr const int kNproperties = 8;
  int prop1[kNproperties];
  int prop2[kNproperties];
  for (int i = 0; i < kNproperties; ++i) {
    prop1[i] = property(rng);
    prop2[i] = property(rng);
  }

  Molecular_Properties_Integer mpr1, mpr2;
  mpr1.build_from_array(prop1, kNproperties);
  mpr2.build_from_array(prop2, kNproperties);

  double total = 0.0;
  for (auto _ : state) {
    total += mpr1.similarity(mpr2);
  }

  // Do something useless to make sure the compiler does not short circuit things.
  if (total == -3.0) {
    exit(0);
  }
}
BENCHMARK(BM_TestInt);

void
BM_TestFloat(benchmark::State& state) {
  std::random_device rd;
  std::default_random_engine rng(rd());
  std::uniform_real_distribution<float> property(0.0, 255.0);

  constexpr const int kNproperties = 8;
  float prop1[kNproperties];
  float prop2[kNproperties];
  for (int i = 0; i < kNproperties; ++i) {
    prop1[i] = property(rng);
    prop2[i] = property(rng);
  }

  Molecular_Properties_Continuous mpr1, mpr2;
  mpr1.build_from_array(prop1, kNproperties);
  mpr2.build_from_array(prop2, kNproperties);
  set_continuous_property_distance_metric(ContinuousPropertyDistanceMetric::kDice);

  double total = 0.0;
  for (auto _ : state) {
    total += mpr1.similarity(mpr2);
  }

  // Do something useless to make sure the compiler does not short circuit things.
  if (total == -3.0) {
    exit(0);
  }
}
BENCHMARK(BM_TestFloat);

}  // namespace

BENCHMARK_MAIN();
