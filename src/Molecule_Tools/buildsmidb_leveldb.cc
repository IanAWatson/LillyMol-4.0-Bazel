#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

#include "leveldb/db.h"

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/report_progress.h"

#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/mol2graph_proto.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/molecule_database_options.pb.h"
#include "Molecule_Tools/molecule_database_options.h"

namespace molecule_database {

using std::cerr;

void Usage(int rc) {
  cerr << " -Q <smi>           kind of smiles to store, 'usmi' or 'kusmi'\n";
  cerr << " -G ...             if storing a graph form, options for molecular graph formation\n";
  cerr << " -c                 remove chirality from input molecules\n";
  cerr << " -l                 strip input molecules to largest fragment\n";
  cerr << " -r <rpt>           report progress every <rpt> molecules processed\n";
  cerr << " -A ...             default aromaticity options, enter '-A help' for info\n";
  cerr << " -E ...             default element options, enter '-A help' for info\n";
  cerr << " -i <type>          input file type\n";
  cerr << " -v                 verbose output\n";
  exit(rc);
}

struct Options {
  int verbose = 0;

  int molecules_read = 0;
  int molecules_stored = 0;
  int appended_to_existing = 0;

  LLYMol::MoleculeDatabaseOptions database_structure_options;

  // By default, leveldb's Put method overwrites what is present.
  // We might want to concatenate all identifiers.
  bool append_to_existing_entries = false;

  // Maybe this should be an option.
  const char * kSpacer = ":";
               
  // If appending to existing data, we first need to lookup every
  // new item. Some ReadOptions to control that lookup.
  leveldb::ReadOptions read_options;
  // And write options to control writing.
  leveldb::WriteOptions write_options;

  Report_Progress report_progress;

  Chemical_Standardisation chemical_standardisation;

  // For efficiency, we store a Mol2Graph class at file scope,
  // so it can be re-used for each molecule.
  Mol2Graph mol2graph;

  // The database we are loading.
  std::unique_ptr<leveldb::DB> database;

  public:
    int Build (Command_Line& cl);
  
    int Report(std::ostream& output) const;

    int DoStore(const IWString& smiles, const IWString& name);

    int AppendToExistingData(const IWString& smiles, const IWString& name);

    int BuildSmiDb(const IWString& smiles, const IWString& name);

    int BuildSmiDb(Molecule& m);

    int BuildSmiDb(data_source_and_type<Molecule>& input);

    int BuildSmiDb(const char * fname,
               FileType input_type);
};

int
Options::Build(Command_Line& cl) {
  verbose = cl.option_count('v');

  if (cl.option_present('r')) {
    if (! report_progress.initialise(cl, 'r', verbose)) {
      cerr << "Cannot initialize progress reporting\n";
      return 1;
    }
  }

  if (cl.option_present('g')) {
    if (! chemical_standardisation.construct_from_command_line(cl, verbose > 1, 'g'))
    {
      cerr << "Cannot initialise chemical standardisation\n";
      Usage(6);
    }
  }

  std::optional<LLYMol::MoleculeDatabaseOptions> maybe_options = molecule_database::OptionsFromCmdline(cl, verbose);
  if (! maybe_options.has_value()) {
    cerr << "Cannot initialize common storage options\n";
    return 1;
  }

  database_structure_options = *maybe_options;

  // Set file scope mol2graph if needed.
  if (maybe_options->mol2graph().active()) {
    mol2graph = Mol2GraphFromProto(maybe_options->mol2graph());
  }

  if (cl.option_present('a')) {
    append_to_existing_entries = true;
    if (verbose)
      cerr << "Will append new data to old\n";
  }

  if (! cl.option_present('d')) {
    cerr << "Options::Build: must specify database via the -d option\n";
    return 0;
  }

  const std::string dbname = cl.std_string_value('d');

  leveldb::Options leveldb_options;
  leveldb_options.create_if_missing = true;

  leveldb::DB *db;
  const leveldb::Status status = leveldb::DB::Open(leveldb_options, dbname, &db);
  if (!status.ok()) {
    cerr << "Options::Build:cannot open " << std::quoted(dbname) << " " << status.ToString() << "\n";
    return 0;
  }

  database.reset(db);

  if (! molecule_database::StoreStructureTypeInformation(database_structure_options, write_options, database.get())) {
    cerr << "Cannot store structure specification\n";
    return 1;
  }

  return 1;
}

int
Options::Report(std::ostream& output) const {
  output << "Read " << molecules_read << " molecules, stored " << molecules_stored << '\n';
  if (append_to_existing_entries) {
    output << appended_to_existing << " items appended to existing keys\n";
  }
  return 1;
}

int verbose = 0;

const char * prog_name;

int molecules_read = 0;
int molecules_stored = 0;
int appended_to_existing = 0;

// By default, leveldb's Put method overwrites what is present.
// We might want to concatenate all identifiers.
bool append_to_existing_entries = false;

// Maybe this should be an option.
const char * kSpacer = ":";
               
// If appending to existing data, we first need to lookup every
// new item. Some ReadOptions to control that lookup.

leveldb::ReadOptions read_options;

Report_Progress report_progress;

Chemical_Standardisation chemical_standardisation;

// For efficiency, we store a Mol2Graph class at file scope,
// so it can be re-used for each molecule.
Mol2Graph mol2graph;

int
Options::DoStore(const IWString& smiles, const IWString& name) {
  leveldb::Slice key(smiles.data(), smiles.length());
  leveldb::Slice value(name.data(), name.length());

  const leveldb::Status status = database->Put(write_options, key, value);
  if (status.ok()) {
    molecules_stored++;
    return 1;
  }

  cerr << "Options::DoStore:Did not store " << status.ToString() << "\n";

  return 0;
}

int
Options::AppendToExistingData(const IWString& smiles, const IWString& name) {
  leveldb::Slice key(smiles.data(), smiles.length());

  std::string value;
  const leveldb::Status status = database->Get(read_options, key, &value);
  if (!status.ok()) {
    return DoStore(smiles, name);
  }

  IWString new_value(value);
  new_value.append_with_spacer(name, kSpacer);
  appended_to_existing++;
  return DoStore(smiles, new_value);
}

int
Options::BuildSmiDb(const IWString& smiles, const IWString& name) {
  if (append_to_existing_entries) {
    return AppendToExistingData(smiles, name);
  }

  return DoStore(smiles, name);
}

int
Options::BuildSmiDb(Molecule& m) {
  molecule_database::Preprocess(m, database_structure_options, chemical_standardisation);

  if (mol2graph.active()) {
    // TOTO(ianiwatson@gmail.com) Implement graph processing...
    cerr << "Implement graph storage sometime\n";
    return 0;
  }

  if (database_structure_options.unique_smiles()) {
    return BuildSmiDb(m.unique_smiles(), m.name());
  }

  if (database_structure_options.unique_kekule_smiles()) {
    return BuildSmiDb(m.UniqueKekuleSmiles(), m.name());
  }

  cerr << "Not sure what to store\n";
  return 0;
}

int
Options::BuildSmiDb(data_source_and_type<Molecule>& input) {
  Molecule * m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);
    molecules_read++;
    report_progress.report("", "\n", std::cerr);
    if (! BuildSmiDb(*m)) {
      cerr << "Options::BuildSmiDb:fatal error processing '" << m->name() << "'\n";
      return 0;
    }
  }

  return 1;
}

int
Options::BuildSmiDb(const char * fname,
               FileType input_type) {
  if (FILE_TYPE_INVALID == input_type)
  {
    input_type = discern_file_type_from_name(fname);
    assert (FILE_TYPE_INVALID != input_type);
  }

  data_source_and_type<Molecule> input(input_type, fname);

  if (! input.ok())
  {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return BuildSmiDb(input);
}

int BuildSmiDb(int argc, char ** argv) {
  Command_Line cl(argc, argv, "vA:g:clGi:d:r:aQ:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    Usage(1);
  }

  Options options;
  if (! options.Build(cl)) {
    cerr << "Cannot initialise options\n";
    Usage(1);
  }

  if (cl.option_present('A')) {
    if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
      cerr << "Cannot initialise aromaticity\n";
      return 1;
    }
  }

  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      Usage(6);
    }
  }

  if (FILE_TYPE_INVALID != input_type)   // great, explicitly specified
    ;
  else if (1 == cl.number_elements() && 0 == strcmp("-", cl[0]))   // reading from a pipe, assume smiles input
  {
    if (verbose)
      cerr << "Assuming smiles input from pipe read\n";
    input_type = FILE_TYPE_SMI;
  }
  else if (all_files_recognised_by_suffix(cl))
    ;
  else
  {
    cerr << "Cannot discern file types from names\n";
    return 4;
  }

  for (int i = 0; i < cl.number_elements(); ++i) {
    const char * fname = cl[i];
    if (! options.BuildSmiDb(fname, input_type)) {
      cerr << "Fatal error processing " << std::quoted(fname) << "\n";
      return i + 1;
    }
  }

  if (verbose) {
    options.Report(cerr);
  }


  return 0;
}

}  // namespace molecule_database

int
main (int argc, char **argv)
{
  molecule_database::prog_name = argv[0];

  int rc = molecule_database::BuildSmiDb(argc, argv);

  return rc;
}
