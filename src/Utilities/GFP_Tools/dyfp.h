#ifndef DYFP_H
#define DYFP_H

/*
  This is really a re-implementation of iw_dy_fp which had grown too
  large
*/

typedef float tversky_coeff_t;

class Tversky;

#define IWDYFP_NSET_NOT_COMPUTED -9

#include "Foundational/iwbits/iwbits.h"
#include "Foundational/iwbits/fixed_bit_vector.h"

class IWDYFP : public fixed_bit_vector::FixedBitVector
{
  private:
    int _nset;

//  private functions

    void _default_values();

    int _bits_in_common(const IWDYFP & f2) const;

    similarity_type_t _tanimoto(IWDYFP &);
    similarity_type_t _tanimoto_multiplier(IWDYFP &);

  public:
    IWDYFP();
    IWDYFP(int);
    ~IWDYFP();

    int ok() const;
    int debug_print(std::ostream &) const;
 
    int allocate_space_for_bits(int);

    int construct_from_array_of_ints(const int *, int nb);    // we overload this method
    int ConstructFromArrayOfInts(const int * ii, int nb);

    int construct_from_tdt_record(const IWString &);
    int construct_from_tdt_record(const const_IWSubstring &);
    int construct_from_daylight_ascii_representation(const const_IWSubstring &);
    int construct_from_sparse_representation(const const_IWSubstring &);
    int construct_from_ascii_01_representation(const char *, int);

    int construct_from_hex(const const_IWSubstring &);

    int construct_from_descriptor_tdt_record(const IWString &);        // includes the tag and <> chars
    int construct_from_descriptor_record(const const_IWSubstring &);   // no tags
    int ConstructFromDaylightAsciiBitRep(const const_IWSubstring & s);   // no tags

    int write_daylight_ascii_representation(std::ostream & os,
                                const const_IWSubstring & data_item_name);
    int daylight_ascii_representation(IWString & result);
    int daylight_ascii_tdt(IWString & result, const const_IWSubstring &);

    int nset() { if (_nset >= 0) return _nset;
                        else return _nset = FixedBitVector::nset();}

    int compute_nset() { _nset = FixedBitVector::nset(); return _nset;}

    int set_nset(int);

//  Many of the common operators are overloaded because if _nset changes,
//  we need to be able to record that. Should use the base object if you
//  want to use these operators

    void iwand(const IWDYFP &);
    void iwand(const IWDYFP &, int &);
    void iwor(const IWDYFP &);
    void iwor(const IWDYFP &, int &);
    void iwxor(const IWDYFP &);

    // We can fairly quickly see if a fingerprint can be compared with another
    // to achieve a given similarity or larger.
    int SimilarityMightBeGreaterThan(const IWDYFP& rhs, similarity_type_t threshold) const;

//  Sometimes this is useful

    int recompute_nset() { return _nset = FixedBitVector::nset();}

    similarity_type_t tanimoto         (IWDYFP &);
    similarity_type_t fraction_matched (IWDYFP &);
    similarity_type_t tversky          (IWDYFP &);
    similarity_type_t tversky          (IWDYFP &, const Tversky &);

    similarity_type_t optimistic_distance(IWDYFP &, const Tversky &);

    similarity_type_t fvb_modified_tanimoto(IWDYFP &);
    similarity_type_t russel_rao(IWDYFP &);
    similarity_type_t forbes_similarity(IWDYFP &);
    similarity_type_t simple_matching(IWDYFP &);
    similarity_type_t sorensendice(IWDYFP &);
    similarity_type_t overlap(IWDYFP &);

    void operator +=(const IWDYFP &);

    IWDYFP & operator =(const IWDYFP &);

    int fold(int);

    void * copy_to_contiguous_storage(void *) const; 
    void * copy_to_contiguous_storage_gpu(void *) const; 

    const void * build_from_contiguous_storage(const void * p, int);
};

#endif
