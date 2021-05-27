#include <iomanip>
#include <memory>

#include "leveldb/db.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/mol2graph_proto.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/molecule.h"

#include "Molecule_Tools/molecule_database_options.pb.h"
#include "Molecule_Tools/molecule_database_options.h"

namespace molecule_database {

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

void usage(int rc) {
  exit(rc);
}

int DoStore(const IWString& smiles, const IWString& name,
               const leveldb::WriteOptions& write_options,
               leveldb::DB* database) {
  leveldb::Slice key(smiles.data(), smiles.length());
  leveldb::Slice value(name.data(), name.length());

  const leveldb::Status status = database->Put(write_options, key, value);
  if (status.ok()) {
    molecules_stored++;
    return 1;
  }
  cerr << "Did not store " << status.ToString() << "\n";

  return 0;
}

int AppendToExistingData(const IWString& smiles, const IWString& name,
               const leveldb::WriteOptions& write_options,
               leveldb::DB* database) {
  leveldb::Slice key(smiles.data(), smiles.length());

  std::string value;
  const leveldb::Status status = database->Get(read_options, key, &value);
  if (!status.ok()) {
    return DoStore(smiles, name, write_options, database);
  }

  IWString new_value(value);
  new_value.append_with_spacer(name, kSpacer);
  appended_to_existing++;
  return DoStore(smiles, new_value, write_options, database);
}

int BuildSmiDb(const IWString& smiles, const IWString& name,
               const leveldb::WriteOptions& write_options,
               leveldb::DB* database) {
  if (append_to_existing_entries) {
    return AppendToExistingData(smiles, name, write_options, database);
  }

  return DoStore(smiles, name, write_options, database);
}

int BuildSmiDb(Molecule& m,
               const LLYMol::MoleculeDatabaseOptions& options,
               const leveldb::WriteOptions& write_options,
               leveldb::DB* database) {
  molecule_database::Preprocess(m, options, chemical_standardisation);

  if (! options.mol2graph().active()) {
    return BuildSmiDb(m.unique_smiles(), m.name(), write_options, database);
  }

  return 1;
}

int BuildSmiDb(data_source_and_type<Molecule>& input,
               const LLYMol::MoleculeDatabaseOptions& options,
               const leveldb::WriteOptions& write_options,
               leveldb::DB* database) {
  Molecule * m;
  while (NULL != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);
    molecules_read++;
    report_progress.report("", "\n", std::cerr);
    if (! BuildSmiDb(*m, options, write_options, database)) {
      cerr << "Fatal error processing '" << m->name() << "'\n";
      return 0;
    }
  }

  return 1;
}

int BuildSmiDb(const char * fname,
               FileType input_type,
               const LLYMol::MoleculeDatabaseOptions& options,
               const leveldb::WriteOptions& write_options,
               leveldb::DB* database) {
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

  return BuildSmiDb(input, options, write_options, database);
}

int BuildSmiDb(int argc, char ** argv) {
  Command_Line cl(argc, argv, "vA:g:clGi:d:r:a");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

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
      usage(6);
    }
  }

  std::optional<LLYMol::MoleculeDatabaseOptions> options = molecule_database::OptionsFromCmdline(cl, verbose);
  if (! options.has_value()) {
    cerr << "Cannot initialize common storage options\n";
    return 1;
  }

  // Set file scope mol2graph if needed.
  if (options.value().mol2graph().active()) {
    mol2graph = Mol2GraphFromProto(options.value().mol2graph());
  }

  if (! cl.option_present('d')) {
    cerr << "Must specify name of database via the -d option\n";
    usage(1);
  }

  if (cl.option_present('a')) {
    append_to_existing_entries = true;
    if (verbose)
      cerr << "Will append new data to old\n";
  }

  const std::string dbname = cl.std_string_value('d');

  leveldb::Options leveldb_options;
  leveldb_options.create_if_missing = true;

  leveldb::DB *database;
  const leveldb::Status status = leveldb::DB::Open(leveldb_options, dbname,  &database);
  if (!status.ok()) {
    cerr << "Cannot open " << std::quoted(dbname) << " " << status.ToString() << "\n";
    return 1;
  }
  std::unique_ptr<leveldb::DB> free_database(database);

  leveldb::WriteOptions write_options;

  if (! molecule_database::StoreStructureTypeInformation(options.value(), write_options, database)) {
    cerr << "Cannot store structure specification\n";
    return 1;
  }


  FileType input_type = FILE_TYPE_INVALID;
  if (cl.option_present('i'))
  {
    if (! process_input_type(cl, input_type))
    {
      cerr << "Cannot determine input type\n";
      usage(6);
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
    if (! BuildSmiDb(fname, input_type, options.value(), write_options, database)) {
      cerr << "Fatal error processing " << std::quoted(fname) << "\n";
      return i + 1;
    }
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules, stored " << molecules_stored << '\n';
    if (append_to_existing_entries) {
      cerr << appended_to_existing << " items appended to existing keys\n";
    }
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
