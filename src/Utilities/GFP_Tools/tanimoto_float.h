#ifndef UTILITIES_GFP_TOOLS_TANIMOTO_FLOAT_H
#define UTILITIES_GFP_TOOLS_TANIMOTO_FLOAT_H

namespace tanimoto_float {

// compute sum(min/max) over `v1` and `v2`.
// for now `n` is ignored and assumed to be 8.
float AvxTanimoto(const float * v1, const float* v2, int n);

}  // namespace tanimoto_float

#endif  // UTILITIES_GFP_TOOLS_TANIMOTO_FLOAT_H
