// process the output of tp1_pipe.sh and create a formatted report

#include <iostream>
#include <memory>
#include <string>

#define RESIZABLE_ARRAY_IMPLEMENTATION
#define RESIZABLE_ARRAY_IWQSORT_IMPLEMENTATION

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/data_source/iwstring_data_source.h"
#include "Foundational/iwqsort/iwqsort.h"
#include "Foundational/iwmisc/iwre2.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iw_stl_hash_map.h"

const char *prog_name = NULL;

namespace mc_summarise {

using std::cerr;

void
Usage(int rc) {
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << '\n';
  cerr << "Converts the results of tp1_pipe.sh into a summary\n";
  cerr << " -r             include the rejection reason with the output\n";
  cerr << " -z             include zero demerit molecules in the output\n";
  cerr << " -h             include a header record\n";
  cerr << " -D             exclude the D prefix\n";
  cerr << " -B <stem>      bad file stem (default 'bad')\n";
  cerr << " -s <string>    separator between output fields\n";
  cerr << " -u             change spaces in reason fields to underscores\n";
  cerr << " -t             give report of which reasons hit\n";
  cerr << " -T <fname>     write report on rejection reasons to <fname>\n";
  cerr << " -b             in the -T file, produce two column format rejections | demerits\n";
  cerr << " -X             produce table in LaTex format\n";
  cerr << " -A             produce table in AsciiDoc format\n";
  cerr << " -m <n>         in the table file, discard any reason with <n> or fewer examples\n";
  cerr << " -c             produce a demerit based scale factor file\n";
  cerr << " -f <n>         numeric demerit value for rejections (default 100)\n";
  cerr << " -k             only process the survivors file (do not process bad0, bad1...)\n";
  cerr << " -v             verbose output\n";

  exit(rc);
}

constexpr int kRejected = 100;

class Reason_and_Count {
 private:
  const IWString _reason;
  const int _count;

 public:
  Reason_and_Count(const IWString &r, int c) : _reason(r), _count(c) {}

  const IWString &
  reason() const {
    return _reason;
  }
  int
  count() const {
    return _count;
  }

  template <typename T>
  int latex_table(T &) const;
  template <typename T>
  int asciidoc_table(T &os) const;
};

template <typename T>
T &
operator<<(T &os, const Reason_and_Count &rc) {
  os << rc.count() << ' ' << rc.reason();

  return os;
}

template <typename T>
int
Reason_and_Count::latex_table(T &os) const {
  os << _count << " & ";

  if (_reason.contains('_')) {
    IWString tmp(_reason);
    tmp.gsub("_", "\\_");
    os << tmp;
  } else
    os << _reason;

  os << "\\\\\n";

  return 1;
}

template <typename T>
int
Reason_and_Count::asciidoc_table(T &os) const {
  os << '|' << _count << " | " << _reason;

  return 1;
}

class Reason_and_Count_Comparator {
 private:
 public:
  int operator()(const Reason_and_Count *, const Reason_and_Count *) const;
};

int
Reason_and_Count_Comparator::operator()(const Reason_and_Count *rc1,
                                        const Reason_and_Count *rc2) const {
  if (rc1->count() < rc2->count())
    return 1;
  if (rc1->count() > rc2->count())
    return -1;

  return 0;
}

struct Parameters {
  int verbose = 0;

  IWString bad_file_stem = "bad";

  IWString bad_file_stem_array = "BQTP";

  int molecules_written = 0;

  int rejected_molecules = 0;

  int demerited_molecules = 0;

  int include_reason = 0;

  int include_zero_demerit_molecules = 0;

  int include_header_record = 0;

  IWString prepend_string;

  int prepend_d = 1;

  IW_STL_Hash_Map_int reason_bad;
  IW_STL_Hash_Map_int reason_survivor;

  int accumulate_reasons = 0;

  extending_resizable_array<int> demerits_per_molecule;

  char output_separator = ' ';

  int latex_table = 0;

  int asciidoc_table = 0;

  int gsub_spaces_in_reason_to_underscore = 0;

  int bad0_demerit = 200;
  int bad12_demerit = 100;

  int produce_demerit_scaling_file = 0;

  int suppress_normal_output = 0;

  int process_rejection_files = 1;

  int write_reasons = 0;

  int min_reasons_needed_for_output = 0;

  const char * fname_for_reason_summary = nullptr;

  // functions.
  private:
    void DefaultValues();
    int AllFilesPresent(const char *okfile, const IWString &bad_stem) const;
    int WriteDemeritValue(const const_IWSubstring &id,
                  const int demerit,
                  IWString_and_File_Descriptor &output);
    int ProcessSingleRun(const char *okfile,
                   const IWString &bad_stem,
                   IWString_and_File_Descriptor &output);
    int ProcessFromIwdemerit(const char *fname,
                       IWString_and_File_Descriptor& output);
    int ProcessFromIwdemerit(const const_IWSubstring &buffer,
                       int must_have_demerit,
                       IWString_and_File_Descriptor &output);
    int ProcessRejectedDemerits(const IWString& bad_stem,
                       IWString_and_File_Descriptor& output);
    int ProcessSurvivingDemerits(const char* fname,
                IWString_and_File_Descriptor& output);
    int ProcessFromIwdemerit(iwstring_data_source &input,
                       int must_have_demerit,
                       IWString_and_File_Descriptor &output);
    int ProcessBad0(iwstring_data_source& input,
             IW_STL_Hash_Map_int& reasons,
             IWString_and_File_Descriptor& output);
    int ProcessBad0(const const_IWSubstring& buffer,
             IW_STL_Hash_Map_int& reasons,
             IWString_and_File_Descriptor& output);
    int ProcessBad0(const IWString &bad_stem,
             IW_STL_Hash_Map_int &reasons,
             IWString_and_File_Descriptor &output);
    int ProcessBad01(const IWString &fname,
              IW_STL_Hash_Map_int &reasons,
              IWString_and_File_Descriptor &output);
    int ProcessBad12(const const_IWSubstring &buffer,
              IW_STL_Hash_Map_int &reasons,
              IWString_and_File_Descriptor &output);
    int ProcessBad12(iwstring_data_source &input,
              IW_STL_Hash_Map_int &reasons,
              IWString_and_File_Descriptor &output);
    int ProcessBad12(const IWString &fname,
              IW_STL_Hash_Map_int &reasons,
              IWString_and_File_Descriptor &output);
    int ProcessBad1(const IWString& bad_stem,
             IW_STL_Hash_Map_int& reason,
             IWString_and_File_Descriptor& output);
    int ProcessBad2(const IWString& bad_stem,
             IW_STL_Hash_Map_int& reasons,
             IWString_and_File_Descriptor& output);
    int MaybeWriteHeader(IWString_and_File_Descriptor& output) const;
    template <typename T> void
    WriteRejectionReason(const resizable_array_p<Reason_and_Count> &r,
                       const int ndx,
                       T &output) const;
    template <typename T> int
    WriteReasons(const IW_STL_Hash_Map_int &reason_bad,
              const IW_STL_Hash_Map_int &reason_survivor,
              const int min_reasons_needed_for_output,
              T &output) const;
    template <typename T> int
    WriteReasons(const IW_STL_Hash_Map_int &reasons,
              const int min_reasons_needed_for_output,
              T &output) const;
  public:
    Parameters();

    int Initialise(Command_Line& cl);

    int MaybeWriteReasons(std::ostream& output) const;

    int Report(std::ostream& output) const;

    int Process(const Command_Line& cl,
                IWString_and_File_Descriptor& output);
};

void
Parameters::DefaultValues() {
  verbose = 0;

  bad_file_stem = "bad";

  bad_file_stem_array = "BQTP";

  molecules_written = 0;

  rejected_molecules = 0;

  demerited_molecules = 0;

  include_reason = 0;

  include_zero_demerit_molecules = 0;

  include_header_record = 0;

  prepend_string.resize(0);

  prepend_d = 1;

  reason_bad.clear();
  reason_survivor.clear();

  accumulate_reasons = 0;

  demerits_per_molecule.resize(0);

  output_separator = ' ';

  latex_table = 0;

  asciidoc_table = 0;

  gsub_spaces_in_reason_to_underscore = 0;

  bad0_demerit = 200;
  bad12_demerit = 100;

  produce_demerit_scaling_file = 0;

  suppress_normal_output = 0;

  process_rejection_files = 1;

  write_reasons = 0;

  min_reasons_needed_for_output = 0;

  fname_for_reason_summary = nullptr;
}

Parameters::Parameters() {
  DefaultValues();
}

int
Parameters::Initialise(Command_Line& cl) {
  verbose = cl.option_count('v');

  if (cl.option_present('B')) {
    cl.value('B', bad_file_stem);

    if (verbose)
      cerr << "Bad files have stem '" << bad_file_stem << "'\n";

    bad_file_stem_array = bad_file_stem;
  }

  if (cl.option_present('r')) {
    include_reason = 1;

    if (verbose)
      cerr << "The reason for the rejection will be included\n";

    if (cl.option_present('u')) {
      gsub_spaces_in_reason_to_underscore = 1;

      if (verbose)
        cerr << "Spaces in reason converted to underscore\n";
    }
  }

  if (cl.option_present('h')) {
    include_header_record = 1;

    if (verbose)
      cerr << "Will write a header record\n";
  }

  if (cl.option_present('z')) {
    include_zero_demerit_molecules = 1;

    if (verbose)
      cerr << "Will include zero demerit molecules\n";
  }

  if (cl.option_present('D')) {
    prepend_d = 0;

    if (verbose)
      cerr << "No D prefix\n";
  }

  if (cl.option_present('P')) {
    cl.value('P', prepend_string);

    if (verbose)
      cerr << "Will prepend '" << prepend_string << "' to all output\n";
  }

  if (cl.option_present('t')) {
    accumulate_reasons = 1;
    include_reason = 1;

    if (verbose)
      cerr << "Will accumulate reasons for rejections\n";
  }

  if (cl.option_present('T')) {
    fname_for_reason_summary = cl.option_value('T');

    if (verbose)
      cerr << "Summary of rejection reasons written to '" << fname_for_reason_summary << "'\n";

    accumulate_reasons = 1;
    include_reason = 1;
  }

  if (cl.option_present('c')) {
    produce_demerit_scaling_file = 1;

    if (verbose)
      cerr << "Will produce a file where demerit values have been converted to scaling factors\n";
  }

  if (cl.option_present('n')) {
    suppress_normal_output = 1;

    if (verbose)
      cerr << "Will suppress normal output\n";
  }

  if (cl.option_present('k')) {
    process_rejection_files = 0;

    if (verbose)
      cerr << "Will only process survivor files (not bad0, bad1...)\n";
  }

  const_IWSubstring fj;

  if (cl.option_present('j'))
    cl.value('j', fj);
  else if (cl.option_present('f'))
    cl.value('f', fj);

  if (fj.length()) {
    int j;
    if (!fj.numeric_value(j) || j < 1) {
      cerr << "The rejection threshold value (-f) must be a whole +ve number\n";
      Usage(4);
    }

    bad0_demerit = j;
    bad12_demerit = j;

    if (verbose)
      cerr << "Demerit value assigned to hard rejections " << j << '\n';
  }

  if (cl.option_present('s')) {
    IWString s;
    cl.value('s', s);

    if (!char_name_to_char(s)) {
      cerr << "Unrecognised output separator '" << s << "'\n";
      return 0;
    }

    output_separator = s[0];
  }

  if (cl.option_present('X')) {
    latex_table = 1;
  }

  if (cl.option_present('A')) {
    asciidoc_table = 1;
  }

  if (cl.option_present('m')) {
    accumulate_reasons = 1;
    if (!cl.value('m', min_reasons_needed_for_output) || min_reasons_needed_for_output < 0) {
      cerr << "The minimum reasons for output (-m) option must be a whole +ve number\n";
      Usage(1);
    }

    if (verbose)
      cerr << "Will not write any rejection/demerit reason to the -T file with fewer than "
           << min_reasons_needed_for_output << '\n';
  }

  if (cl.option_present('b')) {
    write_reasons = 1;
  }

  return 1;
}

int
Parameters::MaybeWriteHeader(IWString_and_File_Descriptor& output) const {
  if (! include_header_record) {
    return 1;
  }

  output << "Name" << output_separator;
  if (produce_demerit_scaling_file)
    output << "scale\n";
  else
    output << "demerit\n";

  return 1;
}

int
Parameters::Report(std::ostream& output) const {
  output << "Wrote " << molecules_written << " molecules, ";
  output << rejected_molecules << " rejected, " << demerited_molecules << " demerited\n";

  for (int i = 0; i < demerits_per_molecule.number_elements(); i++) {
    if (demerits_per_molecule[i])
      output << demerits_per_molecule[i] << " molecules had " << i << " demerits\n";
  }

  return output.good();
}

int
Parameters::WriteDemeritValue(const const_IWSubstring &id,
                  const int demerit,
                  IWString_and_File_Descriptor &output) {
  molecules_written++;

  if (demerit >= kRejected)
    rejected_molecules++;
  else if (demerit > 0)
    demerited_molecules++;

  if (suppress_normal_output)
    return 1;

  output << id << output_separator;

  if (prepend_string.length())
    output << prepend_string;

  if (prepend_d)
    output << 'D';

  if (produce_demerit_scaling_file) {
    if (demerit >= kRejected)
      output << '0';
    else
      output << static_cast<float>(kRejected - demerit) * 0.01f;
  } else
    output << demerit;

  output.write_if_buffer_holds_more_than(32768);

  return output.good();
}

std::unique_ptr<RE2> d_parentheses;

// Parse a record that looks like
// SMILES ID ... : D(140) too_many_atoms:crown_2_2
int
Parameters::ProcessFromIwdemerit(const const_IWSubstring &buffer,
                       int must_have_demerit,
                       IWString_and_File_Descriptor &output) {
  if (d_parentheses == nullptr) {
    d_parentheses = std::make_unique<RE2>("^D\\([0-9]+\\)$");
  }

  int i = 0;
  const_IWSubstring token;
  buffer.nextword(token, i);  // Smiles

  const_IWSubstring id;
  buffer.nextword(id, i);

  // cerr << "Processing '" << id << "'\n";

  int previous_was_colon = 0;

  int demerit = -1;

  while (buffer.nextword(token, i)) {
    if (':' == token) {
      previous_was_colon = 1;
      continue;
    }

    if (!previous_was_colon)
      continue;

    re2::StringPiece tmp(token.data(), token.length());
    if (! RE2::FullMatch(tmp, *d_parentheses, &demerit)) {
      previous_was_colon = 0;
      continue;
    }

    break;
  }

  // cerr << " Demerit for '" << id << " is '" << dmrt << "'\n";

  if (demerit >= 0) {  // Got a valid value, good.
  } else if (must_have_demerit) {
    cerr << "NO demerit value found for '" << id << "'\n";
    return 0;
  } else {
    if (include_zero_demerit_molecules) {
      WriteDemeritValue(id, 0, output);
      if (!suppress_normal_output)
        output << '\n';
    }

    return output.good();
  }

  if (demerit > kRejected)
    demerit = kRejected;

  WriteDemeritValue(id, demerit, output);

  if (include_reason) {
    int demerit_reasons_this_molecule = 0;

    static IWString myreason;

    while (buffer.nextword(myreason, i, ':')) {
      myreason.strip_leading_blanks();

      if (!suppress_normal_output)
        output << output_separator << myreason;

      if (!accumulate_reasons)
        ;
      else if (demerit >= kRejected)
        reason_bad[myreason]++;
      else
        reason_survivor[myreason]++;

      demerit_reasons_this_molecule++;
    }

    demerits_per_molecule[demerit_reasons_this_molecule]++;
  }

  if (!suppress_normal_output)
    output << '\n';

  output.write_if_buffer_holds_more_than(23768);

  return output.good();
}

int
Parameters::ProcessFromIwdemerit(iwstring_data_source &input,
                       int must_have_demerit,
                       IWString_and_File_Descriptor &output) {
  const_IWSubstring buffer;

  while (input.next_record(buffer)) {
    if (!ProcessFromIwdemerit(buffer, must_have_demerit, output)) {
      cerr << "Invalid from iwdemerit record, line " << input.lines_read() << '\n';
      cerr << buffer << '\n';
      return 0;
    }
  }

  return output.good();
}

std::unique_ptr<RE2> matches_to;

// Process a record that looks like
// SMILES ID ... (1 matches to 'activated_phthalimide')
// Note that this will fail if all queries have been checked,
// in which case there may be multiple matches with the molecule.
// TODO:ianwatson: fix this
int
Parameters::ProcessBad12(const const_IWSubstring &buffer,
              IW_STL_Hash_Map_int &reasons,
              IWString_and_File_Descriptor &output) {
  if (! matches_to) {
    matches_to = std::make_unique<RE2>("^\\S+ +(\\S+)..* \\(([0-9]) matches to '([^']+)'\\)");
  }

  re2::StringPiece tmp(buffer.data(), buffer.length());
  std::string id;
  int nhits;
  std::string reason;
  if (! RE2::FullMatch(tmp, *matches_to, &id, &nhits, &reason)) {
    cerr << "Parameters::ProcessBad12:no matches to in '" << buffer << "'\n";
    return 0;
  }

  WriteDemeritValue(id, bad12_demerit, output);

  if (include_reason) {
    const IWString myreason(reason);

    if (!suppress_normal_output)
      output << output_separator << myreason;

    if (accumulate_reasons)
      reasons[myreason]++;
  }

  if (!suppress_normal_output)
    output << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return output.good();
}

int
Parameters::ProcessBad12(iwstring_data_source &input,
              IW_STL_Hash_Map_int &reasons,
              IWString_and_File_Descriptor &output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (!ProcessBad12(buffer, reasons, output)) {
      cerr << "Invalid bad12 record '" << buffer << "'\n";
      return 0;
    }
  }

  return output.good();
}

int
Parameters::ProcessBad12(const IWString &fname,
              IW_STL_Hash_Map_int &reasons,
              IWString_and_File_Descriptor &output) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ProcessBad12(input, reasons, output);
}

int
Parameters::ProcessBad1(const IWString &bad_stem,
             IW_STL_Hash_Map_int &reason,
             IWString_and_File_Descriptor &output) {
  IWString fname;

  fname = bad_stem;
  fname << "1.smi";

  return ProcessBad12(fname, reason, output);
}

int
Parameters::ProcessBad2(const IWString &bad_stem,
             IW_STL_Hash_Map_int &reasons,
             IWString_and_File_Descriptor &output) {
  IWString fname;

  fname = bad_stem;
  fname << "2.smi";

  return ProcessBad12(fname, reasons, output);
}

int
Parameters::ProcessFromIwdemerit(const char *fname,
                       IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ProcessFromIwdemerit(input, 0 /* must_have_demerit */, output);
}

// Process ${bad_stem}3.smi.
// In that file, all records must contain a demerit value.
int
Parameters::ProcessRejectedDemerits(const IWString& bad_stem,
                       IWString_and_File_Descriptor& output) {
  IWString fname;

  fname = bad_stem;
  fname << "3.smi";

  constexpr int kMustHaveDemerit = 1;

  return ProcessFromIwdemerit(fname.null_terminated_chars(),
                              kMustHaveDemerit,
                              output);
}
// Process the final output of the rules. Written by
// iwdemerit, but this time molecules may, or may not
// have demerits.
int
Parameters::ProcessSurvivingDemerits(const char* fname,
                IWString_and_File_Descriptor& output) {
  constexpr int kOptionalDemerit = 0;
  return ProcessFromIwdemerit(fname, kOptionalDemerit, output);
}

int
Parameters::ProcessBad0(const const_IWSubstring& buffer,
             IW_STL_Hash_Map_int& reasons,
             IWString_and_File_Descriptor& output) {
  int i = 0;
  const_IWSubstring token;

  buffer.nextword(token, i);

  const_IWSubstring id;
  buffer.nextword(id, i);

  int got_tp1 = 0;

  while (buffer.nextword(token, i)) {
    if ("TP1" != token)
      continue;

    got_tp1 = 1;
    break;
  }

  if (!got_tp1) {
    cerr << "No 'TP1' token\n";
    return 0;
  }

  WriteDemeritValue(id, bad0_demerit, output);

  if (include_reason) {
    static IWString myreason;

    myreason.resize_keep_storage(0);

    while (buffer.nextword(token, i)) {
      myreason.append_with_spacer(token);
    }

    if (gsub_spaces_in_reason_to_underscore)
      myreason.gsub(' ', '_');

    if (!suppress_normal_output)
      output << output_separator << myreason;

    if (accumulate_reasons)
      reasons[myreason]++;
  }

  if (!suppress_normal_output)
    output << '\n';

  output.write_if_buffer_holds_more_than(32768);

  return output.good();
}

int
Parameters::ProcessBad0(iwstring_data_source& input,
             IW_STL_Hash_Map_int& reasons,
             IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    if (!ProcessBad0(buffer, reasons, output)) {
      cerr << "Invalid bad0 record '" << buffer << "'\n";
      return 0;
    }
  }

  return output.good();
}

int
Parameters::ProcessBad01(const IWString &fname,
              IW_STL_Hash_Map_int &reasons,
              IWString_and_File_Descriptor &output) {
  iwstring_data_source input(fname);

  if (!input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return ProcessBad0(input, reasons, output);
}

int
Parameters::ProcessBad0(const IWString &bad_stem,
             IW_STL_Hash_Map_int &reasons,
             IWString_and_File_Descriptor &output) {
  IWString fname;

  fname = bad_stem;
  fname << "0.smi";

  return ProcessBad01(fname, reasons, output);
}

int
Parameters::AllFilesPresent(const char *okfile, const IWString &bad_stem) const {
  for (int i = 0; i < 3; i++) {
    IWString fname;
    fname << bad_stem << i << ".smi";

    if (verbose > 1)
      cerr << "Checking '" << fname << "', result " << dash_f(fname.null_terminated_chars())
           << '\n';

    if (!dash_f(fname.null_terminated_chars())) {
      if (verbose)
        cerr << "File '" << fname << "' not present\n";
      return 0;
    }
  }

  return dash_f(okfile);
}

int
Parameters::ProcessSingleRun(const char *okfile,
                   const IWString &bad_stem,
                   IWString_and_File_Descriptor &output) {
  if (process_rejection_files) {
    ProcessBad0(bad_stem, reason_bad, output);
    ProcessBad1(bad_stem, reason_bad, output);
    ProcessBad2(bad_stem, reason_bad, output);
  }
  ProcessRejectedDemerits(bad_stem, output);
  ProcessSurvivingDemerits(okfile, output);

  return output.good();
}

static int
last_item_meeting_support_requirement(const resizable_array_p<Reason_and_Count> &r,
                                      const int min_reasons_needed_for_output) {
  for (int i = r.number_elements() - 1; i >= 0; --i) {
    if (r[i]->count() >= min_reasons_needed_for_output)
      return i;
  }

  return -1;
}

static void
sorted_list_of_reasons(const IW_STL_Hash_Map_int &reason, resizable_array_p<Reason_and_Count> &r) {
  int n = reason.size();

  r.resize(n);

  for (IW_STL_Hash_Map_int::const_iterator i = reason.begin(); i != reason.end(); ++i) {
    Reason_and_Count *t = new Reason_and_Count((*i).first, (*i).second);
    r.add(t);
  }

  Reason_and_Count_Comparator rcc;

  r.iwqsort(rcc);

  return;
}

template <typename T>
void
Parameters::WriteRejectionReason(const resizable_array_p<Reason_and_Count> &r,
                       const int ndx,
                       T &output) const {
  if (ndx < r.number_elements()) {
    if (latex_table)
      r[ndx]->latex_table(output);
    else if (asciidoc_table)
      r[ndx]->asciidoc_table(output);
    else
      output << *(r[ndx]);
  } else {
    if (latex_table)
      output << " & &";
    else if (asciidoc_table)
      output << "|.|.";
  }

  return;
}

// Write a multi column file containing rejection reasons, and demerit reasons
template <typename T>
int
Parameters::WriteReasons(const IW_STL_Hash_Map_int &reason_bad,
              const IW_STL_Hash_Map_int &reason_survivor,
              const int min_reasons_needed_for_output,
              T &output) const {
  resizable_array_p<Reason_and_Count> rej, dem;
  sorted_list_of_reasons(reason_bad, rej);
  sorted_list_of_reasons(reason_survivor, dem);

  int dstop = last_item_meeting_support_requirement(dem, min_reasons_needed_for_output);
  int rstop = last_item_meeting_support_requirement(rej, min_reasons_needed_for_output);

  if (dstop < 0 && rstop < 0)
    return 1;

  if (asciidoc_table) {
    output << "[width=\"30%\", options=\"header\"]\n";
    output << "|===============\n";
    output << "|N |Rejected |N |Demerits\n";
  } else if (latex_table) {
    output << "\\begin{center}\n";
    output << "\\begin{tabular}{r c r c}\n";
    output << "Molecules & Reason Rej & Molecules & Reason Demerit \\\\\n";
    output << "\\hline\n";
  }

  int istop = dstop;
  if (istop < rstop)
    istop = rstop;

  for (int i = 0; i < istop; ++i) {
    WriteRejectionReason(rej, i, output);
    WriteRejectionReason(dem, i, output);
    output << '\n';
  }

  if (asciidoc_table)
    output << "|===============\n";
  else if (latex_table) {
    output << "\\hline\n";
    output << "\\end{tabular}\n";
    output << "\\end{center}\n";
  }

  return 1;
}

template <typename T>
int
Parameters::WriteReasons(const IW_STL_Hash_Map_int &reasons,
              const int min_reasons_needed_for_output,
              T &output) const {
  resizable_array_p<Reason_and_Count> r;

  sorted_list_of_reasons(reasons, r);

  if (asciidoc_table) {
    output << "[width=\"30%\"]\n";
    output << "|===============\n";
  } else if (latex_table) {
    output << "\\begin{center}\n";
    output << "\\begin{tabular}{r l}\n";
    output << "Molecules & Reason \\\\\n";
    output << "\\hline\n";
  }

  for (int i = 0; i < r.number_elements(); i++) {
    const Reason_and_Count *ri = r[i];

    if (ri->count() < min_reasons_needed_for_output)
      break;

    if (latex_table)
      ri->latex_table(output);
    else if (asciidoc_table)
      ri->asciidoc_table(output);
    else
      output << (*ri);

    output << '\n';
  }

  if (asciidoc_table) {
    output << "|===============\n";
  } else if (latex_table) {
    output << "\\hline\n";
    output << "\\end{tabular}\n";
    output << "\\end{center}\n";
  }

  return 1;
}

int
Parameters::MaybeWriteReasons(std::ostream& output) const {
  if (! accumulate_reasons) {
    return 1;
  }

  output << "Encountered " << (reason_bad.size() + reason_survivor.size())
         << " different reasons\n";

  if (fname_for_reason_summary == nullptr) {
    for (IW_STL_Hash_Map_int::const_iterator i = reason_bad.begin(); i != reason_bad.end(); ++i) {
      if (i->second > min_reasons_needed_for_output)
        output << (*i).second << " occurrences of '" << (*i).first << "'\n";
    }
    for (const auto& [reason, count] : reason_bad) {
      if (count > min_reasons_needed_for_output)
        output << count << " occurrences of '" << reason << "'\n";
    }
    for (const auto& [reason, count] : reason_survivor) {
      if (count > min_reasons_needed_for_output)
        output << count << " occurrences of '" << reason << "'\n";
    }
  } else {
    IWString_and_File_Descriptor tfile;
    if (!tfile.open(fname_for_reason_summary)) {
      cerr << "Cannot open reason summary file '" << fname_for_reason_summary << "'\n";
      return 5;
    }

    if (write_reasons) {
      WriteReasons(reason_bad, reason_survivor, min_reasons_needed_for_output, tfile);
    }
    else {
      if (verbose)
        output << reason_bad.size() << " reasons associated with rejected molecules, "
             << reason_survivor.size() << " with survivors\n";
      WriteReasons(reason_bad, min_reasons_needed_for_output, tfile);
      WriteReasons(reason_survivor, min_reasons_needed_for_output, tfile);
    }
  }

  return 1;
}


int
Parameters::Process(const Command_Line& cl,
                    IWString_and_File_Descriptor& output) {
  cerr << "Processing " << cl.number_elements() << " files\n";

  MaybeWriteHeader(output);

  if (cl.number_elements() == 1) {
    if (! ProcessSingleRun(cl[0], bad_file_stem, output)) {
      return 0;
    }
  } else {
    for (int i = 0; i < cl.number_elements(); i++) {
      IWString bstem = bad_file_stem_array;
      bstem << (i + 1) << '_';

      if (!AllFilesPresent(cl[i], bstem))
        return 0;

      if (verbose)
        cerr << "Processing okfile '" << cl[i] << "', bad stem '" << bstem << "'\n";

      if (!ProcessSingleRun(cl[i], bstem, output)) {
        cerr << "Fatal error processing files with stem '" << bstem << "'\n";
        return 0;
      }
    }
  }

  output.flush();

  MaybeWriteReasons(cerr);
  if (verbose) {
    Report(cerr);
  }

  return 1;
}

static int
MCSummarise(int argc, char **argv) {
  Command_Line cl(argc, argv, "vB:S:rhzDs:tT:uj:f:XAP:cnkbm:");

  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(2);
  }

  Parameters params;
  if (! params.Initialise(cl)) {
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  if (! params.Process(cl, output)) {
    return 1;
  }

  return 0;
}

}  // namespace mc_summarise

int
main(int argc, char **argv) {
  prog_name = argv[0];

  int rc = mc_summarise::MCSummarise(argc, argv);

  return rc;
}
