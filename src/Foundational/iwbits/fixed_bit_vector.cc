#include <algorithm>
#include <limits>
#include <iostream>
#include <iomanip>
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

static uint64_t one_bit_64[] = {
  0x1,
  0x2,
  0x4,
  0x8,
  0x10,
  0x20,
  0x40,
  0x80,
  0x100,
  0x200,
  0x400,
  0x800,
  0x1000,
  0x2000,
  0x4000,
  0x8000,
  0x10000,
  0x20000,
  0x40000,
  0x80000,
  0x100000,
  0x200000,
  0x400000,
  0x800000,
  0x1000000,
  0x2000000,
  0x4000000,
  0x8000000,
  0x10000000,
  0x20000000,
  0x40000000,
  0x80000000,
  0x100000000,
  0x200000000,
  0x400000000,
  0x800000000,
  0x1000000000,
  0x2000000000,
  0x4000000000,
  0x8000000000,
  0x10000000000,
  0x20000000000,
  0x40000000000,
  0x80000000000,
  0x100000000000,
  0x200000000000,
  0x400000000000,
  0x800000000000,
  0x1000000000000,
  0x2000000000000,
  0x4000000000000,
  0x8000000000000,
  0x10000000000000,
  0x20000000000000,
  0x40000000000000,
  0x80000000000000,
  0x100000000000000,
  0x200000000000000,
  0x400000000000000,
  0x800000000000000,
  0x1000000000000000,
  0x2000000000000000,
  0x4000000000000000,
  0x8000000000000000,
};

#ifdef NOT_USED_MAYBE_USEFUL_SOMETIME
  static const uint64_t bx[] = {
    0xFFFFFFFF00000000,
    0xFFFF0000,
    0xFF00,
    0xF0,
    0xC,
    0x2,
    };

static const uint8_t BitReverseTable256[256] = 
{
#   define R2(n)     n,     n + 2*64,     n + 1*64,     n + 3*64
#   define R4(n) R2(n), R2(n + 2*16), R2(n + 1*16), R2(n + 3*16)
#   define R6(n) R4(n), R4(n + 2*4 ), R4(n + 1*4 ), R4(n + 3*4 )
    R6(0), R6(2), R6(1), R6(3)
};
#endif

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

void
FixedBitVector::reset() {
  if (_bits == nullptr) {
    return;
  }
  std::fill_n(_bits, _nwords, 0);
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

// Taken from http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogObvious. Amazing!

int
most_significant_bit(uint64_t v)
{
  static const uint64_t b[] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000, 0xFFFFFFFF00000000};
  static const uint64_t S[] = {1, 2, 4, 8, 16, 32};

  uint64_t r = 0; // result of first_bit_set(v) will go here
  for (int i = 5; i >= 0; i--) // unroll for speed...
  {
    if (v & b[i])
    {
      v >>= S[i];
      r |= S[i];
    }
  }

  return r;
}

int
first_bit_set(uint64_t v) {
  static const uint64_t b[] = {
    0x00000000FFFFFFFF,
    0x0000FFFF,
    0x00FF,
    0x0F,
    0x3,
    0x1,
  };
  static const uint64_t S[] = {32, 16, 8, 4, 2, 1};

  if (v == 0) {
    return -1;
  }

  uint64_t result = 0;
  for (int i = 0; i < 6; ++i)
  {
    const uint64_t vnbi = v & b[i];
    if (vnbi == 0)
    {
      v >>= S[i];
      result |= S[i];
    }
    else if (vnbi == b[i]) {
      return result;
    }
  }

  return result;
}

// Note this could be optimized if we knew low order bits were preferentially set.
int
first_unset_bit(uint64_t v)
{

  static const uint64_t b[] = {
    0x00000000FFFFFFFF,
    0x0000FFFF,
    0x00FF,
    0x0F,
    0x3,
    0x1,
  };
  static const uint64_t S[] = {32, 16, 8, 4, 2, 1};

  if (v == 0) {
    return 0;
  }
  if (v == std::numeric_limits<uint64_t>::max()) {
    return -1;
  }

  uint64_t result = 0;
  for (int i = 0; i < 6; ++i)
  {
    const uint64_t vnbi = v & b[i];
    if (vnbi == 0) {
      return result;
    }
    if (vnbi == b[i])  // All bits set, shift
    {
      v >>= S[i];
      result |= S[i];
    }
  }

  return result;
}

int
FixedBitVector::FirstBitSet() const {
  for (int i = 0; i < _nwords; ++i) {
    if (_bits[i] == 0)
      continue;
    return i * 64 + first_bit_set(_bits[i]);
  }
  return -1;
}

std::ostream&
operator<<(std::ostream& output, const FixedBitVector& b) {
  output << "FixedBitVector " << b.nbits();
  output << std::hex;
  for (int i = 0; i < b._nwords; ++i) {
    output << ' ' << b._bits[i];
  }
  output << std::dec;

  return output;
}

}  // namespace fixed_bit_vector
