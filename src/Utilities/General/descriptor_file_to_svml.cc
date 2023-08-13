// Convert a tabular descriptor file to svml form - for xgboost CLI
#include <cstdlib>
#include <iostream>

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwmisc/proto_support.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

#include "Utilities/General/descriptor_model.pb.h"


namespace descriptor_file_to_svml {

using std::cerr;

void
Usage(int rc) {
  cerr << "Converts a descriptor file to svml form\n";
  cerr << " -A <fname>          activity file - when processing a training set\n";
  cerr << " -C <fname>          write the feature cross reference to <fname> when processing training set\n";
  cerr << " -U <fname>          use a previously created (-C) file when processing a test set\n";
  cerr << " -v          verbose\n";

  ::exit(rc);
}

class Options {
  private:
    int _verbose;;

    int _columns_in_input;

    char _input_separator;

    // The activity might be in a column of the name.
    int _activity_column;

    IWString _missing_value;

    IW_STL_Hash_Map_float _id_to_activity;

    // This is used both when a new proto is being created, and when one is read from
    // a file.
    DescriptorModel::Model _proto;

    // when using an existing cross reference proto we need a cross reference from
    // columns in the current input to columns in the output.
    uint32_t* _column_xref;

  // private functions
    int ReadActivityData(const const_IWSubstring& fname);
    int ReadActivityData(iwstring_data_source& input);
    int ReadActivityDataRecord(const const_IWSubstring& buffer);
    int GetResponseName(const const_IWSubstring& header);
    int ReadFeatures(IWString& fname);

    int AnalyseHeaderCreateNew(const const_IWSubstring& header);
    int AnalyseHeaderUsePrevious(const const_IWSubstring& header);

  public:
    Options();
    ~Options();

    char input_separator() const {
      return _input_separator;
    }

    int Initialise(Command_Line& cl);

    int AnalyseHeader(const const_IWSubstring& header);

    int ProcessRecord(const const_IWSubstring& buffer,
                       IWString_and_File_Descriptor& output);

    int WriteProto(const IWString& fname);
};

Options::Options() {
  _verbose = 0;
  _columns_in_input = 0;
  _input_separator = ' ';
  _activity_column = 1;
  _missing_value = '.';
  _column_xref = nullptr;
}

Options::~Options() {
  if (_column_xref != nullptr) {
    delete [] _column_xref;
  }
}

int
Options::Initialise(Command_Line& cl) {
  _verbose = cl.option_count('v');

  if (! cl.option_present('A')) {
    cerr << "Activity file not specified (-A) assuming test set processing\n";
  }

  if (cl.option_present('A')) {
    const_IWSubstring a;
    for (int i = 0; cl.value('A', a, i); ++i) {
      if (a.starts_with("col=")) {
        a.remove_leading_chars(4);
        if (! a.numeric_value(_activity_column) || _activity_column < 1) {
          cerr << "The activity column, col= must be valid column number\n";
          return 0;
        }
        if (_verbose) {
          cerr << "Activity values extracted from " << _activity_column << " of name\n";
        }
        --_activity_column;
      } else if (! ReadActivityData(a)) {
        cerr << "Cannot read activity data '" << a << "'\n";
        return 0;
      }

      if (_verbose) {
        cerr << "Read " << _id_to_activity.size() << " activity/data values\n";
      }
    }
  } else if (cl.option_present('U')) {
    IWString fname = cl.string_value('U');
    if (! ReadFeatures(fname)) {
      cerr << "Cannot read features to process from '" << fname << "'\n";
      return 0;
    }
  } else {
    cerr << "Must specify either -A -C (create new proto) or -U (use existing)\n";
    return 0;
  }

  return 1;
}

int
Options::ReadActivityData(const const_IWSubstring& fname) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Options::ReadActivityData:cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadActivityData(input);
}

int
Options::ReadActivityData(iwstring_data_source& input) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "Options::ReadActivityData:cannot read header record\n";
    return 0;
  }

  if (! GetResponseName(buffer)) {
    return 0;
  }


  while (input.next_record(buffer)) {
    if (! ReadActivityDataRecord(buffer)) {
      cerr << "Options::ReadActivityData:invalid data '" << buffer << "'\n";
      return 0;
    }
  }

  return _id_to_activity.size();
}

int
Options::GetResponseName(const const_IWSubstring& header) {
  int i = 0;
  const_IWSubstring token;
  for (int col = 0; header.nextword(token, i, _input_separator); ++col) {
    if (col == _activity_column) {
      _proto.set_response(token.data(), token.length());
      return 1;
    }
  }

  cerr << "Options::GetResponseName:did not find col " << (_activity_column+1) << " in '" << header << "'\n";
  return 0;
}

int
Options::ReadActivityDataRecord(const const_IWSubstring& buffer) {
  int i = 0;
  const_IWSubstring token;
  IWString id, activity;

  static const char kSep = ' ';

  for (int col = 0; buffer.nextword(token, i, kSep); ++col) {
    if (col == 0) {
      id = token;
    } else if (col == _activity_column) {
      activity = token;
      if (! id.empty()) {
        break;
      }
    }
  }

  if (id.empty() || activity.empty()) {
    cerr << "Options::ReadActivityDataRecord:no id '" << id << "' or activity '" << activity << "'\n";
    return 0;
  }

  float a;
  if (! activity.numeric_value(a)) {
    cerr << "Options:ReadActivityDataRecord:invalid activity '" << activity << "'\n";
    return 0;
  }

  _id_to_activity[id] = a;

  return 1;
}

int
Options::ReadFeatures(IWString& fname) {
  if (! iwmisc::ReadBinaryProto(fname, _proto)) {
    cerr << "Options::ReadFeatures:cannot read features from '" << fname << "'\n";
    return 0;
  }

  return 1;
}

int
Options::AnalyseHeader(const const_IWSubstring& header) {
  if (_proto.xref().size() > 0) {
    return AnalyseHeaderUsePrevious(header);
  } else {
    return AnalyseHeaderCreateNew(header);
  }
}

// We are reading a training set, and are creating a new proto.
int
Options::AnalyseHeaderCreateNew(const const_IWSubstring& header) {
  int i = 0;
  const_IWSubstring token;
  if (! header.nextword(token, i, _input_separator)) {
    cerr << "Options::AnalyseHeader:empty header\n";
    return 0;
  }

  for (int col = 0; header.nextword(token, i, _input_separator); ++col) {
    auto* x = _proto.mutable_xref()->Add();
    x->set_bit(col);
    x->set_name(token.data(), token.length());
  }

  return _proto.xref().size();
}

// We are reading a test set and can use the data in _proto to set up the
// column cross reference.
int
Options::AnalyseHeaderUsePrevious(const const_IWSubstring& header) {
  return 1;
}

int 
Options::ProcessRecord(const const_IWSubstring& buffer,
                       IWString_and_File_Descriptor& output) {
  if (_proto.xref().size() > 0) {
    return ProcessRecordUseExisting(buffer, output);
  }

  const_IWSubstring token;
  int i = 0;
  IWString id;
  if (! buffer.nextword(id, i, _input_separator)) {
    cerr << "options::ProcessRecord:empty record\n";
    return 0;
  }

  if (_id_to_activity.empty()) {
    output << '0';
  } else if (auto iter = _id_to_activity.find(id); iter != _id_to_activity.end()) {
    output << iter->second;
  } else {
    cerr << "Options::ProcessRecord:no activity data for '" << id << "'\n";
    for (const auto& [k, v] : _id_to_activity) {
      cerr << k << ' ' << v << '\n';
    }
    return 0;
  }

  for (int col = 0; buffer.nextword(token, i, _input_separator); ++col) {
    if (token == _missing_value) {
      continue;
    }
    output << ' ' << col << ':' << token;
  }
  output << '\n';

  return 1;
}

int 
Options::ProcessRecordUseExisting(const const_IWSubstring& buffer,
                       IWString_and_File_Descriptor& output) {
  int i = 0;
  const_IWSubstring token;
  if (! buffer.nextword(token, i, _input_separator)) {
    cerr << "Options:ProcessRecordUseExisting:empty record\n";
    return 0;
  }

  output << "0 ";
  for (int col = 0; buffer.nextword(token, i, _input_separator); ++col) {
    auto iter = _proto.xref().find(col);
    if (iter == _proto.xref().end()) {
      continue;
    }
  }

  return 1;
}

int
Options::WriteProto(const IWString& fname) {
  IWString myfname;
  myfname <<fname << ".textproto";
  iwmisc::WriteProtoAsText<DescriptorModel::Model>(_proto, myfname);

  myfname.resize_keep_storage(0);
  myfname << fname << ".dat";
  return iwmisc::WriteBinaryProto<DescriptorModel::Model>(_proto, myfname);
}

int
DescriptorFileToSvmlLine(Options& options,
                     const const_IWSubstring& buffer,
                     IWString_and_File_Descriptor& output) {
  return options.ProcessRecord(buffer, output);
  const_IWSubstring token;
}

int
DescriptorFileToSvml(Options& options,
                     iwstring_data_source& input,
                     IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  if (! input.next_record(buffer)) {
    cerr << "Cannot retrieve header record\n";
    return 0;
  }

  if (! options.AnalyseHeader(buffer)) {
    cerr << "Cannot process header record '" << buffer << "'\n";
    return 0;
  }

  while (input.next_record(buffer)) {
    if (! DescriptorFileToSvmlLine(options, buffer, output)) {
      cerr << "Fatal error on line " << input.lines_read() << '\n';
      cerr << buffer << '\n';
      return 0;
    }
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
DescriptorFileToSvml(Options& options,
                     const char* fname,
                     IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return DescriptorFileToSvml(options, input, output);
}

int
Main(int argc, char** argv) {
  Command_Line cl(argc, argv, "vi:A:C:U:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  const int verbose = cl.option_count('v');

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  if (cl.size() > 1) {
    cerr << "Cannot handle multiple files\n";
    Usage(1);
  }

  Options options;
  if (! options.Initialise(cl)) {
    cerr << "Cannot initialise job options\n";
    return 1;
  }

  IWString_and_File_Descriptor output(1);

  for (const char* fname : cl) {
    if (! DescriptorFileToSvml(options, fname, output)) {
      cerr << "Error processing '" << fname << "'\n";
      return 1;
    }
  }

  output.flush();

  if (cl.option_present('C')) {
    const IWString fname = cl.string_value('C');
    options.WriteProto(fname);
  }

  return 0;
}

}  // namespace descriptor_file_to_svml

int
main(int argc, char** argv) {
  return descriptor_file_to_svml::Main(argc, argv);
}
