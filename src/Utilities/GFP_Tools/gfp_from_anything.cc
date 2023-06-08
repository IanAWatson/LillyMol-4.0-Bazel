// Convert a set of tokens to gfp a fingerprint

#include <iostream>

#include "absl/container/flat_hash_map.h"

#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/sparse_fp_creator.h"
#include "Foundational/iwstring/absl_hash.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

namespace gfp_from_anything {

using std::cerr;

void
Usage() {
  cerr << "Converts tokens in a file to a gfp fingerprint\n";
  cerr << "The input must be sorted by the id of the item that generated the token(s)\n";
  cerr << "For each parent in the input, gather all the tokens belonging to it - from consecutive lines\n";
  cerr << "Convert those tokens to bit numbers and write a gfp non colliding fingerprint\n";
  cerr << " -idcol <col>            the identifier of the owning object in column <col>\n";
  cerr << " -item <col>             the item to be converted to a bit is in column <col>\n";
  cerr << " -J <tag>                create fingerprint with tag <tag>\n";
  cerr << " -S <fname>              a smiles file with the smiles of the parent\n";
  cerr << " -v                      verbose output\n";

  ::exit(0);
}

class GfpFromAnyThing {
  private:
    int _verbose;

    // The column containing the id of the object that produced the tokens.
    int _id_col;

    // The column in which the identifiers are found.
    int _token_col;

    // Keep track of the current identifier being processed.
    IWString _previous_id;

    IWString _smiles_tag;
    IWString _identifier_tag;

    IWString _tag;

    Sparse_Fingerprint_Creator _sfp;

    IW_STL_Hash_Map<IWString, uint32_t> _feature_to_bit;

    absl::flat_hash_map<IWString, IWString> _smiles;

    // We always write the cross reference to a space separated file.
    IWString_and_File_Descriptor _stream_for_cross_reference;

    char _sep;

  // Private functions.
    int ReadSmiles(IWString& fname);
    int ReadSmiles(iwstring_data_source& input);
    int ReadSmilesRecord(const const_IWSubstring& buffer);

    int InsertSmiles(const IWString& id,
                              IWString_and_File_Descriptor& output) const;

  public:
    GfpFromAnyThing();

    int Initialise(Command_Line_v2& cl);

    int Process(const char* fname, IWString_and_File_Descriptor& output);
    int Process(iwstring_data_source& input, IWString_and_File_Descriptor& output);
    int ProcessRecord(const const_IWSubstring& buffer,
                               IWString_and_File_Descriptor& output);

    uint32_t FeatureToBit(const IWString& value);
    uint32_t WriteCurrentFingerprint(const IWString& id,
                                         IWString_and_File_Descriptor& output);

    int Report(std::ostream& output) const;
};

GfpFromAnyThing::GfpFromAnyThing() {
  _verbose = 0;
  _id_col = 0;
  _token_col = 1;
  _tag = "NCX<";
  _smiles_tag = "$SMI<";
  _identifier_tag = "PCN<";
  _sep = ' ';
}

int
GfpFromAnyThing::Initialise(Command_Line_v2& cl) {
  _verbose = cl.option_count('v');

  if (cl.option_present("idcol")) {
    if (! cl.value("idcol", _id_col) || _id_col < 1) {
      cerr << "The identifier column must be a whole +ve number\n";
      Usage();
    }
    if (_verbose) {
      cerr << "Identifiers in column " << _id_col << '\n';
    }
    --_id_col;
  }

  if (cl.option_present("itemcol")) {
    if (! cl.value("itemcol", _token_col) || _token_col < 1) {
      cerr << "The item column must be a whole +ve number\n";
      Usage();
    }
    if (_verbose) {
      cerr << "Items in column " << _token_col  << '\n';
    }
    --_token_col;
  }

  if (cl.option_present('J')) {
    cl.value('J', _tag);
    if (! _tag.ends_with('<')) {
      _tag << '<';
    }
  }

  if (cl.option_present("sep")) {
    IWString s = cl.string_value("sep");
    char_name_to_char(s);
    _sep = s[0];
    if (_verbose) {
      cerr << "Input token separator set to '" << _sep << "'\n";
    }
  }

  if (! cl.option_present("XREF")) {
    cerr << "Must specify the name of the cross reference file with the -XREF option\n";
    return 0;
  }

  if (cl.option_present("XREF")) {
    IWString fname = cl.string_value("XREF");
    if (! _stream_for_cross_reference.open(fname.null_terminated_chars())) {
      cerr << "Cannot open cross reference stream '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Value -> bit number cross reference written to '" << fname << "'\n";
    }
  }

  if (cl.option_present("smi")) {
    IWString fname = cl.string_value("smi");
    if (! ReadSmiles(fname)) {
      cerr << "Cannot read smiles file '" << fname << "'\n";
      return 0;
    }

    if (_verbose) {
      cerr << "Read " << _smiles.size() << " identifier->smiles mappings from '" << fname << "'\n";
    }
  }

  return 1;
}

int
GfpFromAnyThing::ReadSmiles(IWString& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "GfpFromAnyThing::ReadSmiles:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadSmiles(input);
}

int
GfpFromAnyThing::ReadSmiles(iwstring_data_source& input) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! ReadSmilesRecord(buffer)) {
      cerr << "GfpFromAnyThing::ReadSmiles:error reading '" << buffer << "'\n";
      return 0;
    }
  }

  return _smiles.size();
}

int
GfpFromAnyThing::ReadSmilesRecord(const const_IWSubstring& buffer) {
  static const char kSep = ' ';

  IWString smiles, id;
  if (! buffer.split(smiles, kSep, id) || smiles.empty() || id.empty()) {
    return 0;
  }

  _smiles[id] = smiles;

  return 1;
}

int
GfpFromAnyThing::Process(const char* fname,
                         IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "GfpFromAnyThing::Process:cannot open '" << fname << "'\n";
    return 0;
  }

  return Process(input, output);
}

int
GfpFromAnyThing::Process(iwstring_data_source& input,
                         IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (! ProcessRecord(buffer, output)) {
      cerr << "GfpFromAnyThing::Process:cannot process '" << buffer << "'\n";
      return 0;
    }
  }

  return 1;
}

int
GfpFromAnyThing::ProcessRecord(const const_IWSubstring& buffer,
                               IWString_and_File_Descriptor& output) {
  int i = 0;
  IWString id, value;
  const_IWSubstring token;
  for (int col = 0; buffer.nextword_single_delimiter(token, i, _sep); ++col) {
    if (col == _id_col) {
      id = token;
    } else if (col == _token_col) {
      value = token;
    }
  }

  if (id.empty() || value.empty()) {
    cerr << "GfpFromAnyThing::ProcessRecord:id '" << id << " or value '" << value << "' empty\n";
    return 0;
  }

  uint32_t bit = FeatureToBit(value);

  if (id == _previous_id) {
    _sfp.hit_bit(bit);
    return 1;
  }

  WriteCurrentFingerprint(id, output);

  _previous_id = id;
  _sfp.clear();
  _sfp.hit_bit(bit);

  return 1;
}

// Return the corresponding bit number for `value`.
uint32_t
GfpFromAnyThing::FeatureToBit(const IWString& value) {
  const auto iter = _feature_to_bit.find(value);
  if (iter != _feature_to_bit.end()) {
    return iter->second;
  }

  const uint32_t result = _feature_to_bit.size();

  _feature_to_bit.emplace(std::make_pair(value, result));

  static constexpr char kSep = ' ';

  _stream_for_cross_reference << value << kSep << result << '\n';
  _stream_for_cross_reference.write_if_buffer_holds_more_than(32768);

  return result;
}

int
GfpFromAnyThing::InsertSmiles(const IWString& id,
                              IWString_and_File_Descriptor& output) const {
  const auto iter = _smiles.find(id);
  if (iter == _smiles.end()) {
    return 0;
  }

  output << _smiles_tag << iter->second << ">\n";

  return 1;
}

uint32_t
GfpFromAnyThing::WriteCurrentFingerprint(const IWString& id,
                                         IWString_and_File_Descriptor& output) {
  if (_smiles.empty()) {
  } else if (InsertSmiles(id, output)) {
  } else {
    cerr << "GfpFromAnyThing::WriteCurrentFingerprint:no smiles for '" << id << "'\n";
    return 0;
  }

  output << _identifier_tag << id << ">\n";
  _sfp.write_fingerprint(_tag, output);

  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

int
GfpFromAnyThing::Report(std::ostream& output) const {
  output << "Assigned " << _feature_to_bit.size() << " feature to bit mappings\n";

  return 1;
}

int
Main(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v-J=s-idcol=ipos-item=ipos-J=s-sep=s-XREF=s-smi=sfile");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage();
  }

  const int verbose = cl.option_count('v');

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage();
  }

  GfpFromAnyThing doit;
  if (! doit.Initialise(cl)) {
    cerr << "Cannot initialise GfpFromAnyThing\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);

  for (const char* fname : cl) {
    if (! doit.Process(fname, output)) {
      cerr << "Cannot process '" << fname << "'\n";
      return 1;
    }
  }

  if (verbose) {
    doit.Report(cerr);
  }

  return 0;
}

}  // namespace gfp_from_anything

int
main(int argc, char ** argv)
{
  int rc = gfp_from_anything::Main(argc, argv);

  return rc;
}
