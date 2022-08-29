#ifndef GFP_STANDARD_H
#define GFP_STANDARD_H

#include <array>
#include <optional>

/*
  We so often deal with MPR MK MK2 IW we want to make a very efficient
  implementation
*/

#include "gfp.h"

inline constexpr int kMkFp = 1;
inline constexpr int kMk2Fp = 2;
inline constexpr int kIwFp = 0;

class GFP_Standard
{
  private:
    int _molecular_properties[8];
    unsigned char _iw[256];
    unsigned char _mk[32];
    unsigned char _mk2[32];

    int _nset_mk;
    int _nset_mk2;
    int _nset_iw;

  public:
    void build_molecular_properties (const Molecular_Properties_Integer &);
    void build_iw  (IWDYFP &);
    void build_mk  (IWDYFP &);   // only non const because it calls nset()
    void build_mk2 (IWDYFP &);

    int * molecular_properties () { return _molecular_properties;}

    void build_iwfp(const unsigned char *, int nset);             // note last param is bits set
    template <typename T> void build_iwfp(const T *, int nset);   // note last param is bits set

    void build_mk(const int *);
    void build_mk2(const int *);

    int natoms () const { return _molecular_properties[0];}
    int nrings () const { return _molecular_properties[1];}
    int aromatic_atoms () const { return _molecular_properties[4];}

    float tanimoto (const GFP_Standard &) const;
    float tanimoto_distance (const GFP_Standard & rhs) const { return 1.0f - tanimoto(rhs);}
    void tanimoto_distance_2 (const GFP_Standard & rhs, const GFP_Standard & rhs2, float *) const;

    // Only return a result if the value will be <= `must_be_closer_than`.
    std::optional<float> tanimoto_distance_if_less(const GFP_Standard& rhs, float must_be_closer_than) const;
    // Only return a result if the value will be >= `must_be_further_than`.
    std::optional<float> tanimoto_distance_if_greater(const GFP_Standard& rhs, float must_be_further_than) const;
};

extern int can_be_compared (const GFP_Standard & fp1, const GFP_Standard & fp2);
extern int standard_fingerprints_present ();

namespace gfp {
// Fingerprints can be present in any order in the input, and
// we need a mapping from fixed fingerprint number to the
// components of a GFP_Standard.
// function GetStandardFingerprintIndices fills an std::array
// and this enum describes which index holds what.
// Deliberately not an enum class because the primary
// purpose of this is to use as an array index.
enum StdFpIndex {
  kIWfp = 0,
  kMK = 1,
  kMK2 = 2
};

extern int GetStandardFingerprintIndices(std::array<int, 3>& xref);
}  // namespace gfp

extern void set_bits_in_mk (int s);

#endif
