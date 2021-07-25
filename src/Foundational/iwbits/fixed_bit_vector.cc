#include <algorithm>
#include <nmmintrin.h>

#include "fixed_bit_vector.h"

namespace fixed_bit_vector {

/* Each single bit. Generate via
a = one(UInt64)
for o in 0:63
  b = IOBuffer()
  show(b, a << o)
  println(String(take!(b)))
end
*/

uint64_t one_bit_64[] =  {
0x0000000000000001,
0x0000000000000002,
0x0000000000000004,
0x0000000000000008,
0x0000000000000010,
0x0000000000000020,
0x0000000000000040,
0x0000000000000080,
0x0000000000000100,
0x0000000000000200,
0x0000000000000400,
0x0000000000000800,
0x0000000000001000,
0x0000000000002000,
0x0000000000004000,
0x0000000000008000,
0x0000000000010000,
0x0000000000020000,
0x0000000000040000,
0x0000000000080000,
0x0000000000100000,
0x0000000000200000,
0x0000000000400000,
0x0000000000800000,
0x0000000001000000,
0x0000000002000000,
0x0000000004000000,
0x0000000008000000,
0x0000000010000000,
0x0000000020000000,
0x0000000040000000,
0x0000000080000000,
0x0000000100000000,
0x0000000200000000,
0x0000000400000000,
0x0000000800000000,
0x0000001000000000,
0x0000002000000000,
0x0000004000000000,
0x0000008000000000,
0x0000010000000000,
0x0000020000000000,
0x0000040000000000,
0x0000080000000000,
0x0000100000000000,
0x0000200000000000,
0x0000400000000000,
0x0000800000000000,
0x0001000000000000,
0x0002000000000000,
0x0004000000000000,
0x0008000000000000,
0x0010000000000000,
0x0020000000000000,
0x0040000000000000,
0x0080000000000000,
0x0100000000000000,
0x0200000000000000,
0x0400000000000000,
0x0800000000000000,
0x1000000000000000,
0x2000000000000000,
0x4000000000000000,
0x8000000000000000
};

// Return the number of 64 bit words needed for `nbits` bits.
int
words_for_bits(int nbits) {
  int result = nbits / 64;
  if (nbits % 64 == 0)
    return result;
  return result + 1;
}

// Initialise our data structures to handle `nb` bits.
// If it is not a multiple of 64, it is rounded up.
void
FixedBitVector::_allocate_bits(int nb) {
  _nwords = words_for_bits(nb);
  _bits = new uint64_t[_nwords];
  std::fill_n(_bits, _nwords, 0);
  _nbits = _nwords * 64;
}


void
FixedBitVector::_default_values() {
  _bits = nullptr;
  _nwords = 0;
  _nbits = 0;
}

FixedBitVector::FixedBitVector() {
  _default_values();
}

FixedBitVector::FixedBitVector(int nb) {
  if (nb == 0) {
    _default_values();
    return;
  }

  _nwords = words_for_bits(nb);
  _allocate_bits(nb);
}

FixedBitVector::~FixedBitVector() {
  if (_bits != nullptr)
    delete [] _bits;
}

void
FixedBitVector::resize(int nb) {
  if (_bits != nullptr)
    delete [] _bits;

  if (nb == 0) {
    _default_values();
    return;
  }

  _allocate_bits(nb);
}

bool
FixedBitVector::is_set(int b) {
  return _bits[b / 64] & one_bit_64[b % 64];
}

void
FixedBitVector::set_bit(int b) {
  _bits[b / 64] |= one_bit_64[b % 64];
}

void
FixedBitVector::unset_bit(int b) {
  _bits[b / 64] &= ~one_bit_64[b % 64];
}

int
FixedBitVector::nset() const {
  int rc = 0;
  for (int i = 0; i < _nwords; ++i) {
    rc +=  _mm_popcnt_u64(_bits[i]);
  }
  return rc;
}

// Could possibly be made more efficient with loop unrolling, see gfp_standard.cc
int
FixedBitVector::BitsInCommon(const FixedBitVector& rhs) const {
  int rc = 0;
  for (int i = 0; i < _nwords; ++i) {
    rc +=  _mm_popcnt_u64(_bits[i] & rhs._bits[i]);
  }
  return rc;
}

}  // namespace fixed_bit_vector
