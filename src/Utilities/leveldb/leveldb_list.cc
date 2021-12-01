// List the contents of a leveldb database.

#include <iomanip>

#include "leveldb/db.h"
#include "re2/re2.h"

#include "Foundational/cmdline/cmdline.h"

namespace leveldb_list {

using std::cerr;

const char * prog_name;

int verbose = 0;

int keys_traversed = 0;

int suppress_outut = 0;

std::unique_ptr<RE2> filter_key;
int no_match_filter_key = 0;
std::unique_ptr<RE2> filter_value;
int no_match_filter_value = 0;

int keys_written = 0;

IWString output_separator(" ");

void usage(int rc) {
  cerr << "List the contents of a LevelDb database\n";
  cerr << " -K <rx>        only process keys   that match <rx>\n";
  cerr << " -V <rx>        only process values that match <rx>\n";
  cerr << " -s <sep>       output separator\n";
  cerr << " -n             suppress output\n";
  cerr << " -v             verbose output\n";
  exit(rc);
}

bool matches_filter(const leveldb::Slice& s,
                        std::unique_ptr<RE2>& regex,
                        int & counter) {
  if (! regex) {  // Not active, so all vales match.
    return true;
  }
  re2::StringPiece tmp(s.data(), s.size());
  if (RE2::PartialMatch(tmp, *regex)) {
    return true;
  }

  counter++;
  return false;
}

void Process(const leveldb::Iterator* iter,
            IWString_and_File_Descriptor& output) {
  const leveldb::Slice key = iter->key();
  if (!matches_filter(key, filter_key, no_match_filter_key)) {
   return;
  }
  const leveldb::Slice value = iter->value();
  if (!matches_filter(value, filter_value, no_match_filter_value)) {
    return;
  }

  keys_written++;  // Even if suppressed.

  if (suppress_outut) {
    return;
  }

  const_IWSubstring tmp(key.data(), key.size());
  output << tmp << output_separator;
  tmp.set(value.data(), value.size());
  output << tmp << '\n';
  output.write_if_buffer_holds_more_than(8192);
}

int LevelDbList(leveldb::Iterator* iter,
                IWString_and_File_Descriptor& output) {
  while (iter->Valid()) {
    keys_traversed++;
    Process(iter, output);
    iter->Next();
  }
  return 1;
}

int LevelDbList(leveldb::DB * database,
                IWString_and_File_Descriptor& output) {
  leveldb::ReadOptions read_options;
  leveldb::Iterator * iter = database->NewIterator(read_options);
  std::unique_ptr<leveldb::Iterator> free_iter(iter);
  iter->SeekToFirst();
  return LevelDbList(iter, output);
}

int LevelDbList(const char * dbname,
                IWString_and_File_Descriptor& output) {
  leveldb::DB *database;
  leveldb::Options leveldb_options;
  const leveldb::Status status = leveldb::DB::Open(leveldb_options, dbname,  &database);
  if (!status.ok()) {
    cerr << "Cannot open " << std::quoted(dbname) << " " << status.ToString() << "\n";
    return 1;
  }
  std::unique_ptr<leveldb::DB> free_database(database);

  return LevelDbList(database, output);
}

int LevelDbList(int argc, char ** argv) {
  Command_Line cl(argc, argv, "vnK:V:s:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('n')) {
    suppress_outut = 1;
    if (verbose)
      cerr << "Output suppressed\n";
  }

  if (cl.option_present('K')) {
    const std::string k = cl.std_string_value('K');
    filter_key.reset(new RE2(k));
    if (verbose) {
      cerr << "Only keys that match " << std::quoted(filter_key->pattern()) << " will be written\n";
    }
  }

  if (cl.option_present('V')) {
    const std::string v = cl.std_string_value('V');
    filter_value.reset(new RE2(v));
    if (verbose) {
      cerr << "Only values that match " << std::quoted(filter_value->pattern()) << " will be written\n";
    }
  }

  if (cl.option_present('s')) {
    cl.value('s', output_separator);
    if (verbose) {
      cerr << "Output separator '" << output_separator << "'\n";
    }
  }

  IWString_and_File_Descriptor output(1);

  if (cl.empty()) {
    cerr << "Must specify database name on command line\n";
    usage(1);
  }

  for (const auto fname : cl) {
    if (! LevelDbList(fname, output)) {
      cerr << "Fatal error processing " << std::quoted(fname) << "\n";
      return 1;
    }
  }
  output.flush();

  if (verbose) {
    cerr << "Traversed " << keys_traversed << " keys\n";
    if (filter_key) {
      cerr << no_match_filter_key << " keys did not match '" << filter_key->pattern() << "'\n";
    }
    if (filter_value) {
      cerr << no_match_filter_value << " values did not match '" << filter_value->pattern() << "'\n";
    }
    cerr << keys_written << " keys written\n";
  }


  return 0;
}
}  // leveldb_list

int
main (int argc, char **argv) {
  leveldb_list::prog_name = argv[0];

  int rc = leveldb_list::LevelDbList(argc, argv);

  return rc;
}
