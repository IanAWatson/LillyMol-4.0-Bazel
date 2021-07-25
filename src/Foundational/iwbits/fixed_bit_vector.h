#ifndef FOUNDATIONAL_IWBITS_FIXED_BIT_VECTOR_H
#define FOUNDATIONAL_IWBITS_FIXED_BIT_VECTOR_H

#include <cstdint>
namespace fixed_bit_vector {
// A BitVector class that only operates on 64 bit words.
// Minimal functionality.
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

    int nset() const;

    int BitsInCommon(const FixedBitVector& rhs) const;
};

}  // namespace fixed_bit_vector
#endif  // FOUNDATIONAL_IWBITS_FIXED_BIT_VECTOR_H
