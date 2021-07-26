#ifndef FOUNDATIONAL_IWBITS_FIXED_BIT_VECTOR_H
#define FOUNDATIONAL_IWBITS_FIXED_BIT_VECTOR_H

#include <cstdint>
namespace fixed_bit_vector {
// A BitVector class that only operates on 64 bit words.
// Deliberately minimal functionality. Designed for speed.
class FixedBitVector {
  protected:
    uint64_t * _bits;
    int _nwords;
    int _nbits;

  // Private functions
  private:
    void _default_values();
    void _allocate_bits(int nb);

  public:
    FixedBitVector();
    FixedBitVector(int nb);
    ~FixedBitVector();

    void resize(int nb);

    int nbits() const { return _nbits;}

    void set_bit(int b);
    void unset_bit(int b);
    bool is_set(int b);

    // Turn off all bits.
    void reset();

    // The number of bits set.
    int nset() const;

    // Bits in common with another vector. Must be same size.
    int BitsInCommon(const FixedBitVector& rhs) const;

    // Returns -1 if no bits are set.
    int FirstBitSet() const;
    // Returns -1 if all bits are set.
    int FirstUnsetBit() const;
};

// Returns the first bit set in a uint64_t. Exposed just for testing.
// Note that this is not the most significant bit, but the least significant bit.
int first_bit_set(uint64_t);
// Returns the first unset bit in a uint64_t. Note that it returns the first
// unset least significant bit.
int first_unset_bit(uint64_t);

}  // namespace fixed_bit_vector
#endif  // FOUNDATIONAL_IWBITS_FIXED_BIT_VECTOR_H
