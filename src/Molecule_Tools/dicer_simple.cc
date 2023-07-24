// Dicer via dicer_api

#include <algorithm>
#include <iostream>
#include <memory>

#define IWQSORT_IMPLEMENTATION
#define IWQSORT_FO_IMPLEMENTATION
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwqsort/iwqsort.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/rwsubstructure.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/target.h"

#include "Molecule_Tools/dicer_api.h"

namespace dicer_api {

using std::cerr;

void
Usage(int rc) {
  cerr << "Dicer variant using dicer_api\n";

  cerr << " -c                remove chirality\n";
  cerr << " -l                strip to largest fragment\n";
  cerr << " -v                verbose output\n";

  ::exit(rc);
}

class LocalOptions {
  private:
    int _verbose;

    int _molecules_read;

    int _reduce_to_largest_fragment;

    int _remove_chirality;

    Element_Transformations _element_transformations;

    Chemical_Standardisation _chemical_standardisation;

  public:
    LocalOptions();

    int Initialise(Command_Line& cl);

    int Preprocess(Molecule& m);

    int Report(std::ostream& output) const;
};

LocalOptions::LocalOptions() {
  _verbose = 0;
  _molecules_read = 0;
  _reduce_to_largest_fragment = 0;
  _remove_chirality = 1;
}

int
LocalOptions::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present('c')) {
    _remove_chirality = 1;
    if (_verbose) {
      cerr << "Will remove chirality from input molecules\n";
    }
  }

  if (cl.option_present('g')) {
    if (! _chemical_standardisation.construct_from_command_line(cl, _verbose > 1, 'g')) {
      cerr << "Cannot initialise chemical standardisation\n";
      return 0;
    }
  }

  if (cl.option_present('l')) {
    _reduce_to_largest_fragment = 1;
    if (_verbose) {
      cerr << "Will reduce molecules to largest fragment\n";
    }
  }

  if (cl.option_present('T')) {
    if (! _element_transformations.construct_from_command_line(cl, _verbose, 'T')) {
      cerr << "LocalOptions::initialise:cannot initialise element transformations (-T)\n";
      return 0;
    }
  }

  return 1;
}

int
LocalOptions::Preprocess(Molecule& m) {
  if (_reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }

  if (_chemical_standardisation.active()) {
    _chemical_standardisation.process(m);
  }

  if (_remove_chirality) {
    m.remove_all_chiral_centres();
  }

  if (_element_transformations.active()) {
    _element_transformations.process(m);
  }

  if (m.empty()) {
    return 0;
  }

  return 1;
}

int
LocalOptions::Report(std::ostream& output) const {
  output << "LocalOptions:read " << _molecules_read << " molecules\n";
  if (_molecules_read == 0) {
    return 1;
  }

  return 1;
}

struct BitCount {
  uint32_t bit;
  uint32_t count;

  BitCount() {
    bit = 0;
    count = 0;
  }
};

class BitCountComparitor {
  public:
    int operator()(const BitCount& bc1, const BitCount& bc2) const {
      if (bc1.count < bc2.count) {
        return 1;
      }
      if (bc1.count == bc2.count) {
        return 0;
      }

      return -1;
    }
};

class BitCounts {
  private:
    uint32_t _ndx;
    BitCount* _bit_count;
    uint32_t _nbits;

  public:
    BitCounts() {
      _ndx = 0;
      _bit_count = nullptr;
      _nbits = 0;
    }
    ~BitCounts();

    int Initialise(uint32_t x, const Sparse_Fingerprint_Creator& sfc);

    uint32_t nbits() const {
      return _nbits;
    }

    const BitCount* bit_count() {
      return _bit_count;
    }

};

int
BitCounts::Initialise(uint32_t x, const Sparse_Fingerprint_Creator& sfc) {
  _ndx = x;
  const auto& bits_found = sfc.bits_found();
  _nbits = bits_found.size();
  if (_nbits == 0) {
    return 1;
  }

  _bit_count = new BitCount[_nbits];
  uint32_t ndx = 0;
  for (const auto [b, c] : bits_found) {
    _bit_count[ndx].bit = b;
    _bit_count[ndx].count = c;
  }

  static BitCountComparitor cmp;

  if (_nbits > 1) {
    iwqsort(_bit_count, _nbits, cmp);
  }

  return 1;
}

BitCounts::~BitCounts() {
  if (_bit_count != nullptr) {
    delete [] _bit_count;
  }
}

int
MaxCommonSubstructure(const DicerApi& dicer_api,
                      const BitCounts& fp1,
                      const BitCounts& fp2,
                      IWString_and_File_Descriptor& output) {
  if (fp1.nbits() < fp2.nbits()) {
    return MaxCommonSubstructure(dicer_api, fp2, fp1, output);
  }

  uint32_t nb1 = fp1.nbits();
  uint32_t nb2 = fp2.nbits();
  uint32_t iptr = 0;
  uint32_t jptr = 0;
  while (fp1.bit_count[iptr] < fp1.bit_count[jptr]) {
  }

  return 1;
}

int
MaxCommonSubstructure(const DicerApi& dicer_api,
                      resizable_array_p<Sparse_Fingerprint_Creator>& fingerprint,
                      IWString_and_File_Descriptor& output) {
  uint32_t nbits = dicer_api.nbits();
  std::unique_ptr<BitCount[]> bit_set = std::make_unique<BitCount[]>(nbits);

  for (uint32_t i = 0; i < nbits; ++i) {
    bit_set[i].bit = i;
  }

  for (const Sparse_Fingerprint_Creator* fp : fingerprint) {
    for (auto [bit, _] : fp->bits_found()) {
      ++bit_set[bit].count;
    }
  }
  BitCountComparitor cmp;
#ifdef DEBUG_MAX_COMMON_SUBSTRUCTURE
  cerr << "For sorting\n";
  for (int i = 0; i < nbits; ++i) {
    cerr << " bit " << bit_set[i].bit << " count " << bit_set[i].count << '\n';
  }
#endif

  iwqsort(bit_set.get(), nbits, cmp);
  for (uint32_t i = 0; i < nbits; ++i) {
    uint32_t bit = bit_set[i].bit;
    cerr << "Bit " << bit << " count " << bit_set[i].count << ' ' << dicer_api.Smiles(bit) << '\n';
  }

  const uint32_t nmolecules = fingerprint.size();

  std::unique_ptr<BitCounts[]> bit_counts = std::make_unique<BitCounts[]>(nmolecules);
  for (uint32_t i = 0; i < nmolecules; ++i) {
    bit_counts[i].Initialise(i, *fingerprint[i]);
  }

  for (uint32_t i = 0; i < nmolecules; ++i) {
    for (uint32_t j = i + 1; j < nmolecules; ++j) {
      MaxCommonSubstructure(dicer_api, bit_counts[i], bit_counts[j], output);
    }
  }

  return 1;
}

int
MaxCommonSubstructure(DicerApi& dicer,
            Molecule& m,
            resizable_array_p<Sparse_Fingerprint_Creator>& fingerprint) {
  std::unique_ptr<Sparse_Fingerprint_Creator> sfc = std::make_unique<Sparse_Fingerprint_Creator>();
  dicer.Process(m, *sfc);
  fingerprint << sfc.release();

  return 1;
}

int
Dicer(DicerApi& dicer,
      Molecule& m,
      IWString_and_File_Descriptor& output) {
  std::cout << m.smiles() << ' ' << m.name() << '\n';
  Sparse_Fingerprint_Creator sfc;
  return dicer.Process(m, sfc);
}

int
MaxCommonSubstructure(LocalOptions& local_options,
      DicerApi& dicer,
      data_source_and_type<Molecule>& input,
      resizable_array_p<Sparse_Fingerprint_Creator>& fingerprint) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! local_options.Preprocess(*m)) {
      return 0;
    }

    MaxCommonSubstructure(dicer, *m, fingerprint);
  }

  return 1;
}

int
Dicer(LocalOptions& local_options,
      DicerApi& dicer,
      data_source_and_type<Molecule>& input,
      IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr) {
    std::unique_ptr<Molecule> free_m(m);

    if (! local_options.Preprocess(*m)) {
      return 0;
    }

    Dicer(dicer, *m, output);
  }

  return 1;
}

int
Dicer(LocalOptions& local_options, DicerApi& dicer,
            const char * fname,
            FileType input_type,
            IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return Dicer(local_options, dicer, input, output);
}

int
MaxCommonSubstructure(LocalOptions& local_options, DicerApi& dicer,
            const char * fname,
            FileType input_type,
            resizable_array_p<Sparse_Fingerprint_Creator>& fingerprint) {
  if (input_type == FILE_TYPE_INVALID) {
    input_type = discern_file_type_from_name(fname);
  }

  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.ok()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return MaxCommonSubstructure(local_options, dicer, input, fingerprint);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vE:A:i:g:lcT:D:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
    cerr << "Cannot process aromaticity options\n";
    return 1;
  }

  if (! process_elements(cl, verbose, 'E')) {
    cerr << "Cannot process standard elements options (-E)\n";
    return 1;
  }

  LocalOptions local_options;
  if (! local_options.Initialise(cl)) {
    cerr << "Cannot initialise local options\n";
    Usage(1);
  }

  DicerApi dicer;
  if (! dicer.Initialise(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  int do_mcs = 0;
  if (cl.option_present('D')) {
    do_mcs = 1;
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      cerr << "Cannot determine input type\n";
      Usage(1);
    }
  } else if (all_files_recognised_by_suffix(cl)) {
  } else {
    input_type = FILE_TYPE_SMI;
  }

  IWString_and_File_Descriptor output(1);
  if (do_mcs) {
    resizable_array_p<Sparse_Fingerprint_Creator> fingerprint;
    for (const char* fname : cl) {
      if (! MaxCommonSubstructure(local_options, dicer, fname, input_type, fingerprint)) {
        cerr << "Fatal error processing '" << fname << "'\n";
        return 1;
      }
      MaxCommonSubstructure(dicer, fingerprint, output);
    }
    if (verbose) {
      cerr << "Read " << fingerprint.size() << " fingerprints\n";
    }
  } else {
    for (const char* fname : cl) {
      if (! Dicer(local_options, dicer, fname, input_type, output)) {
        cerr << "Fatal error processing '" << fname << "'\n";
        return 1;
      }
    }
  }

  if (verbose) {
    dicer.Report(cerr);
  }

  return 0;
}

}  // namespace dicer_api

int
main(int argc, char** argv) {
  int rc = dicer_api::Main(argc, argv);

  return rc;
}
