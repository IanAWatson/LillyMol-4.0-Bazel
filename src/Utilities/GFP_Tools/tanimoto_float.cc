// Use AVX instructions to compute a tanimoto between arrays of floats
// Adapted from
// https://stackoverflow.com/questions/52415188/sse-avx-choose-from-two-m256-float-vectors-based-on-per-element-min-and-max-a
// And for the reduction
// https://photolens.tech/fastest-way-to-do-horizontal-sse-vector-sum-or-other-reduction/
// Needs options:
// -mavx512vl -mavx -mssse3

#include <immintrin.h>

namespace tanimoto_float {

float hsum_ps_sse3(__m128 v) {
    __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
    __m128 sums = _mm_add_ps(v, shuf);
    shuf        = _mm_movehl_ps(shuf, sums); // high half -> low half
    sums        = _mm_add_ss(sums, shuf);
    return        _mm_cvtss_f32(sums);
}

float
AvxTanimoto(const float * v1, const float* v2, int n) {
  __m256 m1 = _mm256_load_ps(v1);
  __m256 m2 = _mm256_load_ps(v2);
  __m256 mask = _mm256_cmp_ps(m1, m2, _CMP_GE_OS);
  __m256 mins = _mm256_blendv_ps(m1, m2, mask);
  __m256 maxs = _mm256_blendv_ps(m2, m1, mask);
  __m256 ratio = _mm256_div_ps(mins, maxs);

  // Reduction.
  __m128 vlow = _mm256_castps256_ps128(ratio);
  __m128 vhigh = _mm256_extractf128_ps(ratio, 1);  // high 128
  vlow = _mm_add_ps(vlow, vhigh);
  return hsum_ps_sse3(vlow) * 0.125;
}

}  // namespace tanimoto_float
