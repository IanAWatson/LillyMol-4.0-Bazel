// Look up molecules in a LevelDB database.
#include <iomanip>
#include <iostream>

#include "leveldb/db.h"
#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/report_progress.h"
#include "Molecule_Lib/aromatic.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/standardise.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/mol2graph_proto.h"
#include "Molecule_Lib/output.h"

#include "Molecule_Tools/molecule_database_options.h"
#include "Molecule_Tools/molecule_database_options.pb.h"

namespace molecule_database {

using std::cerr;

int verbose = 0;

const char * prog_name;

int molecules_read = 0;
int molecules_found = 0;

Report_Progress report_progress;

Chemical_Standardisation chemical_standardisation;

Molecule_Output_Object stream_for_found_structures;
Molecule_Output_Object stream_for_not_found_structures;

// For efficiency, we store a Mol2Graph class at file scope,
// so it can be re-used for each molecule.
Mol2Graph mol2graph;


void usage(int rc) {
  exit(rc);
}

int HandleFoundMolecule(Molecule& m, const std::string& value) {
  molecules_found++;

  if (stream_for_found_structures.active()) {
    IWString new_name(m.name());
    const const_IWSubstring svalue(value.data(), value.size());
    new_name.append_with_spacer(svalue);
    m.set_name(new_name);
    stream_for_found_structures.write(m);
  }
  return 1;
}

int DoLookup(Molecule& m, const IWString& smiles,
               const leveldb::ReadOptions& read_options,
               leveldb::DB* database) {
  leveldb::Slice key(smiles.data(), smiles.length());

  std::string value;
  const leveldb::Status status = database->Get(read_options, key, &value);
  if (status.ok()) {
    return HandleFoundMolecule(m, value);
  }

  if (stream_for_not_found_structures.active()) {
    m.invalidate_smiles();
    stream_for_not_found_structures.write(m);
  }

  return 1;
}

int InDatabase(Molecule& m,
               const LLYMol::MoleculeDatabaseOptions& options,
               const leveldb::ReadOptions& read_options,
               leveldb::DB* database) {
  molecule_database::Preprocess(m, options, chemical_standardisation);

  if (! options.mol2graph().active()) {
    return DoLookup(m, m.unique_smiles(), read_options, database);
  }

  // Note that we do not actually use options.mol2graph().
  Molecule mcopy(m);
  mcopy.change_to_graph_form(mol2graph);
  return DoLookup(m, mcopy.unique_smiles(), read_options, database);
}

int InDatabase(data_source_and_type<Molecule>& input,
               const LLYMol::MoleculeDatabaseOptions& options,
               const leveldb::ReadOptions& read_options,
               leveldb::DB* database) {
  Molecule * m;
  while (nullptr != (m = input.next_molecule())) {
    std::unique_ptr<Molecule> free_m(m);
    molecules_read++;
    report_progress.report("", "\n", std::cerr);
    if (! InDatabase(*m, options, read_options, database)) {
      cerr << "Fatal error processing '" << m->name() << "'\n";
      return 0;
    }
  }

  return 1;
}

int InDatabase(const char * fname,
               FileType input_type,
               const LLYMol::MoleculeDatabaseOptions& options,
               const leveldb::ReadOptions& read_options,
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

  return InDatabase(input, options, read_options, database);
}

int InDatabase(int argc, char ** argv) {
  Command_Line cl(argc, argv, "vA:g:clGi:d:r:F:U:o:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "Unrecognised options encountered\n";
    usage(1);
  }

  verbose = cl.option_count('v');

  if (cl.option_present('A')) {
    if (! process_standard_aromaticity_options(cl, verbose, 'A')) {
      cerr << "Cannot initialise aromaticity\n";
      return 1;
    }
  }

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

  if (! cl.option_present('d')) {
    cerr << "Must specify name of database via the -d option\n";
    usage(1);
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

  leveldb::ReadOptions read_options;

  std::optional<LLYMol::MoleculeDatabaseOptions> options = molecule_database::GetStructureTypeInformation(read_options, database);
  if (! options.has_value()) {
    cerr << "Cannot initialize common storage options\n";
    return 1;
  }

  // Set file scope mol2graph if needed.
  if (options.value().mol2graph().active()) {
    mol2graph = Mol2GraphFromProto(options.value().mol2graph());
  }

  if (cl.option_present('F')) {
    if (! stream_for_found_structures.determine_output_types(cl, 'o')) {
      cerr << "Cannot determine output type(s)\n";
      return 1;
    }
    const IWString f = cl.option_value('F');
    if (stream_for_found_structures.would_overwrite_input_files(cl, f)) {
      cerr << "Cannot overwrite input(s) '" << f << "'\n";
    }
    if (! stream_for_found_structures.new_stem(f)) {
      cerr << "Cannot open stream for found molecules '" << f << "'\n";
      return 1;
    }
    if (verbose)
      cerr << "Molecules found written to '" << f << "'\n";
  }

  if (cl.option_present('U')) {
    if (! stream_for_not_found_structures.determine_output_types(cl, 'o')) {
      cerr << "Cannot determine output type(s)\n";
      return 1;
    }
    const IWString u = cl.option_value('U');
    if (stream_for_not_found_structures.would_overwrite_input_files(cl, u)) {
      cerr << "Cannot overwrite input(s) '" << u << "'\n";
    }
    if (! stream_for_not_found_structures.new_stem(u)) {
      cerr << "Cannot open stream for not found molecules '" << u << "'\n";
      return 1;
    }
    if (verbose)
      cerr << "Molecules not found written to '" << u << "'\n";
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
    if (! InDatabase(fname, input_type, options.value(), read_options, database)) {
      cerr << "Fatal error processing " << std::quoted(fname) << "\n";
      return i + 1;
    }
  }

  if (verbose) {
    cerr << "Read " << molecules_read << " molecules, found " << molecules_found << '\n';
  }


  return 0;
}

}  // namespace molecule_database

int
main (int argc, char **argv)
{
  molecule_database::prog_name = argv[0];

  int rc = molecule_database::InDatabase(argc, argv);

  return rc;
}
