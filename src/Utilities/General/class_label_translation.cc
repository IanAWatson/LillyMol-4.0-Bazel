// Forward and backward class label translation.
// Creating a translation table and using an existing translation table
// are quite separate tasks, so it would be quite reasonable to make
// two programmes to do those tasks separately.

#include <tuple>

#include "google/protobuf/text_format.h"
#include "google/protobuf/io/zero_copy_stream.h"
#include "google/protobuf/io/zero_copy_stream_impl.h"

#include "Foundational/cmdline_v2/cmdline_v2.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"
#include "Foundational/iwmisc/iwdigits.h"

#include "Utilities/General/class_label_translation.pb.h"

namespace class_label_translation {
using std::cerr;

int verbose = 0;

int identifier_column = 0;
int activity_column = 1;
int header_records = 1;
char input_separator = ' ';
char output_separator = ' ';

int too_many_classes = 10;

IWDigits class_numbers;

void
Usage(int rc) {
  cerr << "Performs multi-directional class label translations\n";
  cerr << " -C <fnmame>    create a cross lab4el cross reference proto\n";
  cerr << " -H <fnmame>    use    a cross lab4el cross reference proto\n";
  cerr << " -desc          input is a descriptor file - activity in col 2\n";
  cerr << " -smi           input is a smiles file - activity in col 3\n";
  cerr << " -c <col>       activity values are in column <col>\n";
  cerr << " -hdr <n>       number of header records in the input\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

int
UseCrossReferenceRecord(const const_IWSubstring& buffer,
                   const ClassLabelTranslation::ClassLabelTranslation& mapping,
                   IWString_and_File_Descriptor& output) {
  int i = 0;
  const_IWSubstring token;
  for(int col = 0; buffer.nextword(token, i, input_separator); ++col) {
    if (col > 0) {
      output << output_separator;
    }
    if (col != activity_column) {
      output << token;
      continue;
    }

    const std::string s(token.data(), token.length());
    const auto iter = mapping.to_numeric().find(s);
    if (iter == mapping.to_numeric().end()) {
      cerr << "Unrecognised label '" << token << "'\n";
      return 0;
    }
    class_numbers.append_number(output, iter->second);
  }

  output << '\n';

  return 1;
}

int
UseCrossReference(iwstring_data_source& input,
                   const ClassLabelTranslation::ClassLabelTranslation& mapping,
                   IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  for (int i = 0; i < header_records; ++i) {
    if (! input.next_record(buffer)) {
      cerr << "Cannot head header\n";
      return 0;
    }
    output << buffer << '\n';
    output.write_if_buffer_holds_more_than(8192);
  }

  while (input.next_record(buffer)) {
    if (! UseCrossReferenceRecord(buffer, mapping, output)) {
      cerr << "Cannot process '" << buffer << "'\n";
      return 0;
    }
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

int
UseCrossReference(const char * fname,
                   const ClassLabelTranslation::ClassLabelTranslation& mapping,
                   IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return UseCrossReference(input, mapping, output);
 }

int
UseCrossReference(const Command_Line_v2& cl,
                   const ClassLabelTranslation::ClassLabelTranslation& mapping,
                   IWString_and_File_Descriptor& output) {
  for (const char * fname : cl) {
    if (! UseCrossReference(fname, mapping, output)) {
      cerr << "Fatal error processing '" << fname << "'\n";
      return 0;
    }
  }

  return 1;
}

int
UseCrossReference(const Command_Line_v2& cl,
                   IWString& proto_file_name,
                   IWString_and_File_Descriptor& output) {
  iwstring_data_source input(proto_file_name.null_terminated_chars());
  if (! input.good()) {
    cerr << "Cannot open proto file '" << proto_file_name << "'\n";
    return 0;
  }

  google::protobuf::io::FileInputStream file_input_stream(input.fd());

  ClassLabelTranslation::ClassLabelTranslation mapping;

  if (! google::protobuf::TextFormat::Parse(&file_input_stream, &mapping)) {
    cerr << "GfpToFeature::FromProto:cannot parse '" << proto_file_name << "'\n";
    return 0;
  }

  return UseCrossReference(cl, mapping, output);
}

// End of code for using existing mapping.

int
WriteMapping(const ClassLabelTranslation::ClassLabelTranslation& mapping,
             IWString& fname) {
  IWString_and_File_Descriptor output;
  if (! output.open(fname.null_terminated_chars())) {
    cerr << "WriteMapping:cannot open " << fname << '\n';
    return 0;
  }

  using google::protobuf::io::ZeroCopyOutputStream;
  using google::protobuf::io::FileOutputStream;
  std::unique_ptr<ZeroCopyOutputStream> zero_copy_output(new FileOutputStream(output.fd()));
  if (! google::protobuf::TextFormat::Print(mapping, zero_copy_output.get())) {
    cerr << "GfpToFeature::WriteProto:cannot write " << fname << "\n";
    return 0;
  }

  return 1;
}

int
ProfileLabelsLine(const const_IWSubstring& buffer,
                  IW_STL_Hash_Map_int& count) {
  int i = 0;
  const_IWSubstring token;
  for (int col = 0; buffer.nextword(token, i, input_separator); ++col) {
    if (col != activity_column) {
      continue;
    }
    auto iter = count.find(token);
    if (iter == count.end()) {
      count.emplace(token, 1);
    } else {
      iter->second++;
    }
    return 1;
  }

  cerr << "Did not find column " << activity_column << "  in '" << buffer << "'\n";
  return 0;
}

int
ProfileLabels(iwstring_data_source& input,
              IW_STL_Hash_Map_int& count) {
  const_IWSubstring buffer;
  for (int i = 0; i < header_records; ++i) {
    if (! input.next_record(buffer)) {
      cerr << "Cannot read header record\n";
      return 0;
    }
  }

  while (input.next_record(buffer)) {
    if (! ProfileLabelsLine(buffer, count)) {
      cerr << "Cannot process '" << buffer << "'\n";
      return 0;
    }
  }

  return count.size();
}

int
ProfileLabels(const char * fname,
              IW_STL_Hash_Map_int& count) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ProfileLabels(input, count);
}

std::optional<ClassLabelTranslation::ClassLabelTranslation>
CreateCrossReference(const Command_Line_v2& cl,
                     IWString& output_fname) {
  IW_STL_Hash_Map_int count;
  for (const char * fname : cl) {
    if (! ProfileLabels(fname, count)) {
      cerr << "Fatal error profiling '" << fname << "'\n";
      return std::nullopt;
    }
  }

  const int nclasses = count.size();

  if (verbose) {
    cerr << "Read " << nclasses << " class labels\n";
    if (verbose > 1) {
      for (const auto& [key, value]: count) {
        cerr << value << " instances of '" << key << "'\n";
      }
    }
  }

  if (nclasses == 0) {
    cerr << "No data\n";
    return std::nullopt;
  }

  if (nclasses > too_many_classes) {
    cerr << "Too many classes " << nclasses << '\n';
    return std::nullopt;
  }

  if (nclasses == 1) {
    cerr << "Only one class, cannot continue\n";
    for (const auto& [key, value] : count) {
      cerr << value << " instances of '" << key << "'\n";
    }
    return std::nullopt;
  }

//using std::tuple<IWString, int> = ClassCount;
  using ClassCount = std::tuple<IWString, int>;
  std::unique_ptr<ClassCount[]> label_count(new ClassCount[nclasses]);

  int ndx = 0;
  for (const auto& [key, value] : count) {
    label_count[ndx] = {key, value};
    ndx++;
  }
  std::sort(label_count.get(), label_count.get() + nclasses, [](const ClassCount& c1, const ClassCount& c2) {
    return std::get<1>(c1) < std::get<1>(c2);
  });

  int class_number = -1;
  int class_delta = 2;
  if (nclasses > 2) {
    class_number = 0;
    class_delta = 1;
  }

  ClassLabelTranslation::ClassLabelTranslation mapping;  // to be returned.

  for (int i = 0; i < nclasses; ++i, class_number += class_delta) {
    const auto [label, _] = label_count[i];
    const std::string tmp(label.data(), label.length());
    google::protobuf::MapPair<std::string, int32_t> to_insert(tmp, class_number);
    mapping.mutable_to_numeric()->insert(to_insert);
  }

  for (int i = 0; i < nclasses; ++i) {
    const auto [label, n] = label_count[i];
    const std::string tmp(label.data(), label.length());
    google::protobuf::MapPair<std::string, uint32_t> to_insert(tmp, n);
    mapping.mutable_class_count()->insert(to_insert);
  }

  if (! WriteMapping(mapping, output_fname)) {
    return std::nullopt;
  }

  return mapping;
}

int
ClassLabelTranslation(int argc, char** argv) {
  Command_Line_v2 cl(argc, argv, "-v-C=s-U=s-c=ipos-smi-desc-hdr=ipos");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present("smi")) {
    identifier_column = 1;
    activity_column = 2;
    header_records = 0;
  } else if (cl.option_present("desc")) {
    identifier_column = 0;
    activity_column = 1;
    header_records = 1;
  }

  if (cl.option_present('c')) {
    cl.value('c', activity_column);
    if (verbose) {
      cerr << "Activity in column " << activity_column << '\n';
    }
    activity_column--;
  }

  if (cl.option_present("hdr")) {
    cl.value("hdr", header_records);
    if (verbose)
      cerr << "Inpus has " << header_records << " header records\n";
  }

  class_numbers.initialise(too_many_classes);

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);
  int rc = 0;
  if (cl.option_present('C')) {
    IWString cfile = cl.string_value('C');
    std::optional<ClassLabelTranslation::ClassLabelTranslation> mapping = CreateCrossReference(cl, cfile);
    if (! mapping) {
      cerr << "Cannot create cross reference file '" << cfile << "'\n";
      return 1;
    }
    rc = UseCrossReference(cl, *mapping, output);
  } else if (cl.option_present('U')) {
    IWString ufile = cl.string_value('U');
    rc = UseCrossReference(cl, ufile, output);
  } else {
    cerr << "Must specify one of -C (create) or -U (use) options\n";
    Usage(1);
  }

  if (rc == 0) {
    return 1;
  }
  return 0;
}

}  // namespace class_label_translation

int
main(int argc, char ** argv) {
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  return class_label_translation::ClassLabelTranslation(argc, argv);
}
