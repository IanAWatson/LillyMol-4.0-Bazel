#ifndef FOUNDATIONAL_IWMISC_ABSL_HASH_H
#define FOUNDATIONAL_IWMISC_ABSL_HASH_H

#include <memory>

#include "absl/hash/hash.h"
#include "highwayhash/highwayhash_target.h"
#include "highwayhash/instruction_sets.h"

#include "Foundational/iwstring/iwstring.h"

// namespace iw_absl_hash {

using namespace highwayhash;  // not good style.

template <typename H>
H AbslHashValue(H h, const IWString& s) {
  if (s.empty()) {
    return H::combine(std::move(h), 709902);   // arbitrary number.
  }

  //static const HHKey key HH_ALIGNAS(32) = { 14123, 665242, 2, 902362};
  static const HHKey key HH_ALIGNAS(32) = { 1, 2, 3, 4};

  HHResult64 result;
  InstructionSets::Run<HighwayHash>(key, s.data(), s.length(), &result);

  // Very small, and possibly zero, speed improvement by including s.size.
  return H::combine(std::move(h), result, s.size());
}

// }  // namespace iw_absl_hash

#endif  // FOUNDATIONAL_IWMISC_ABSL_HASH_H
