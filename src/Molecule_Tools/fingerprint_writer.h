#ifndef MOLECULE_TOOLS_FINGERPRINT_WRITER_H_
#define MOLECULE_TOOLS_FINGERPRINT_WRITER_H_

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwbits/fixed_bit_vector.h"
#include "Foundational/iwmisc/iwdigits.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwstring/iwstring.h"

namespace fingerprint_writer {

// Several programmes generate fingerprints. This is a common means
// of producing different kinds of ourput.
//  Sparse fingerprints
//  Fixed width fingerprints
//  Descriptor file output.

class FingerprintWriter {
  private: 
    enum class OutputType {
      kSparse,
      kFixed,
      kDescriptor
    };

    OutputType _output_type;

    // When generating either sparse or fixed width fingerprints.
    IWString _tag;

    // The number of bits/columns if producing either a fixed width fingerprint or
    // descriptor file output.
    int _nbits = 0;

    // Inter-column separator if generating a descriptor file
    IWString _output_column_separator = ' ';

    // Feature name prefix if generating a descriptor file.
    IWString _feature_prefix = "fp";

    // Used when forming an output record for a descriptor file.
    int * _output_vector = nullptr;

    // To speed writing numbers during descriptor file output.
    IWDigits _iwdigits;

    // When writing fingerprints, we need a string to hold the encoded
    // form of the fingerprint prior to writing.
    IWString _fp;

    // When writing fixed width fingerprints, a bitvector for holding
    // the fingerprint being formed.
    fixed_bit_vector::FixedBitVector _fixed;

  // private functions
    void DisplayUsage(char flag, std::ostream& output) const;

    // The functions that do the actual writing, depending on what is in _output_type.
    int WriteFixedFingerprint(const Sparse_Fingerprint_Creator& sfc, IWString_and_File_Descriptor& output);
    int WriteSparseFingerprint(const Sparse_Fingerprint_Creator& sfc, IWString_and_File_Descriptor& output);
    int WriteDescriptors(const IWString& mname, const Sparse_Fingerprint_Creator& sfc, IWString_and_File_Descriptor& output);

  public:
    FingerprintWriter();
    ~FingerprintWriter();

    int Initialise(Command_Line& cl, char flag, int verbose);

    int SetSparseOutput(const char* tag);

    const IWString& tag() const {
      return _tag;
    }

    int IsWritingDescriptors() const {
      return _output_type == OutputType::kDescriptor;
    }

    // If we are writing KDescriptor form, then we can write a header.
    int WriteHeaderIfNeeded(IWString_and_File_Descriptor& output) const;

    // Write a fingerprint. What gets written, depends on what is in _output_type.
    // `mname` is needed if we are writing descriptors.
    int WriteFingerprint(const IWString& mname, const Sparse_Fingerprint_Creator& sfc, IWString_and_File_Descriptor& output);
};

}  // namespace fingerprint_writer

#endif  // MOLECULE_TOOLS_FINGERPRINT_WRITER_H_
