// Generate topological torsion fingerprints.

#include "Foundational/cmdline/cmdline.h"
#include "Foundational/iwmisc/misc.h"
#include "Foundational/iwstring/iwstring.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/standardise.h"

#include "Molecule_Tools/fingerprint_writer.h"

#include "topological_torsion.h"

namespace topological_torsion {

IWString smiles_tag("$SMI<");
IWString identifier_tag("PCN<");

struct JobOptions {
  TorsionOptions torsion_options;

  Atom_Typing_Specification atom_typing;

  bool verbose = false;

  int molecules_read = 0;

  Chemical_Standardisation chemical_standardisation;

  int reduce_to_largest_fragment = 0;

  fingerprint_writer::FingerprintWriter _fp_writer;

  IWString tag;

  int fixed_width = 0;

  bool work_as_filter = false;

  int ntest = 0;
};

using std::cerr;

void
Usage(int rc) {
  cerr << "Generates topological torsion fingerprints\n";
  cerr << " -J ...          fingerprint specification\n";
  cerr << " -f              function as TDT filter\n";
  cerr << " -P <atype>      atom typing specifications\n";
  cerr << " -g ...          chemical standardisation\n";
  cerr << " -l              reduce to largest fragment\n";
  cerr << " -v              verbose output\n";
  exit(rc);
}

int
Preprocess(JobOptions& job_options, Molecule& m) {
  if (job_options.reduce_to_largest_fragment) {
    m.reduce_to_largest_fragment_carefully();
  }
  if (job_options.chemical_standardisation.active()) {
    job_options.chemical_standardisation.process(m);
  }
  return 1;
}

int
RunTest(Molecule& m,
        JobOptions& job_options) {
  const int matoms = m.natoms();
  std::unique_ptr<atom_type_t[]> atom_type = std::make_unique<atom_type_t[]>(matoms);
  if (! job_options.atom_typing.assign_atom_types(m, atom_type.get())) {
    cerr << "TopologicalTorsion:cannot assign atom types\n";
    return 0;
  }
  std::unique_ptr<int[]> include_atom(new_int(m.natoms(), 1));

  IWString expected;

  for (int i = 0; i < job_options.ntest; ++i) {
    IWString rsmi = m.random_smiles();
    Molecule m2;
    if (! m2.build_from_smiles(rsmi)) {
      cerr << "Cannot parse random smiles '" << rsmi << "'\n";
      return 0;
    }
    Sparse_Fingerprint_Creator sfc = TopologicalTorsion(m, atom_type.get(), include_atom.get(), job_options.torsion_options);
    IWString ascii;
    sfc.daylight_ascii_form_with_counts_encoded(ascii);
    if (expected.empty()) {
      expected = ascii;
    } else if (expected != ascii) {
      cerr << "Inconsistent fingerprints " << m.smiles() << ' ' << rsmi << '\n';
      return 0;
    }
  }

  return 1;
}

int
TopologicalTorsion(Molecule& m,
                   JobOptions& job_options,
                   IWString_and_File_Descriptor& output) {
  if (job_options.ntest) {
    return RunTest(m, job_options);
  }

  const int matoms = m.natoms();
  std::unique_ptr<atom_type_t[]> atom_type = std::make_unique<atom_type_t[]>(matoms);
  if (! job_options.atom_typing.assign_atom_types(m, atom_type.get())) {
    cerr << "TopologicalTorsion:cannot assign atom types\n";
    return 0;
  }

  std::unique_ptr<int[]> include_atom(new_int(m.natoms(), 1));
  Sparse_Fingerprint_Creator sfc = TopologicalTorsion(m, atom_type.get(),
                                 include_atom.get(), job_options.torsion_options);

  if (! job_options.work_as_filter) {
    output << smiles_tag << m.smiles() << ">\n";
    output << identifier_tag << m.name() << ">\n";
  }

  job_options._fp_writer.WriteFingerprint(m.name(), sfc, output);

  if (! job_options.work_as_filter) {
    output << "|\n";
  }

  output.write_if_buffer_holds_more_than(8192);

  return 1;
}

int
TopologicalTorsionFilter(iwstring_data_source& input,
                         JobOptions& job_options,
                         IWString_and_File_Descriptor& output) {
  const_IWSubstring buffer;
  while (input.next_record(buffer)) {
    output << buffer << '\n';
    if (! buffer.starts_with(smiles_tag)) {
      continue;
    }
    buffer.remove_leading_chars(smiles_tag.length());
    buffer.chop();
    Molecule m;
    if (! m.build_from_smiles(buffer)) {
      cerr << "TopologicalTorsionFilter:invalid smiles '" << buffer << "'\n";
      return 0;
    }
    if (! TopologicalTorsion(m, job_options, output)) {
      return 0;
    }
  }
  return 1;
}

int
TopologicalTorsion(data_source_and_type<Molecule>& input,
                   JobOptions& job_options,
                   IWString_and_File_Descriptor& output) {
  Molecule * m;
  while ((m = input.next_molecule()) != nullptr ) {
    std::unique_ptr<Molecule> free_m(m);
    job_options.molecules_read++;
    Preprocess(job_options, *m);
    if (! TopologicalTorsion(*m, job_options, output)) {
      return 0;
    }
  }

  return 1;
}

int
TopologicalTorsionFilter(const char * fname,
                   JobOptions& job_options,
                   IWString_and_File_Descriptor& output) {
  iwstring_data_source input(fname);
  if (! input.good()) {
    cerr << "Cannot open '" << fname << "'\n";
    return 0;
  }

  return TopologicalTorsionFilter(input, job_options, output);
}

int
TopologicalTorsion(const char * fname,
                   FileType input_type,
                   JobOptions& job_options,
                   IWString_and_File_Descriptor& output) {
  if (input_type == FILE_TYPE_INVALID)
  {
    input_type = discern_file_type_from_name(fname);
    assert(input_type != FILE_TYPE_INVALID);
  }
  data_source_and_type<Molecule> input(input_type, fname);
  if (! input.good()) {
    cerr << "Cannot open " << fname << '\n';
    return 0;
  }

  return TopologicalTorsion(input, job_options, output);
}

int
TopologicalTorsion(int argc, char ** argv) {
  Command_Line cl(argc, argv, "vi:A:g:lJ:fP:w:T:y:");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  JobOptions job_options;

  job_options.verbose = cl.option_count('v');

  if (cl.option_present('l')) {
    job_options.reduce_to_largest_fragment = 1;
    if (job_options.verbose) {
      cerr << "Will strip to largest fragment\n";
    }
  }

  if (job_options.verbose) {
    cerr << "Tag '" << job_options.tag << "'\n";
  }

  if (cl.option_present('P')) {
    const IWString atype = cl.string_value('P');
    if (! job_options.atom_typing.build(atype)) {
      cerr << "Invalid atom type specification '" << atype << "'\n";
      return 1;
    }
  } else {
    job_options.atom_typing.build("TT");
  }

  if (cl.option_present('w')) {
    if (! cl.value('w', job_options.fixed_width)) {
      cerr << "Invalid fixed width specification (-w)\n";
      return 1;
    }
    if (job_options.verbose) {
      cerr << "Will produce fixed width fingerprints " << job_options.fixed_width << '\n';
    }
  }

  if (cl.option_present('J')) {
    job_options.tag = cl.string_value('J');
    if (! job_options.tag.ends_with('<')) {
      job_options.tag << '<';
    }
  } else if (job_options.fixed_width) {
    job_options.tag = "FPTT<";
  } else {
    job_options.tag = "NCTT<";
  }

  if (cl.option_present('J')) {
    if (! job_options._fp_writer.Initialise(cl, 'J', job_options.verbose)) {
      cerr << "Cannot initialise fingerprint writer (-J)\n";
      return 1;
    }
  } else {
    job_options._fp_writer.SetSparseOutput("NCTT<");
  }

  if (cl.option_present('y')) {
    int y;
    for (int i = 0; cl.value('y', y, i); ++i) {
      if (y == 3) {
        job_options.torsion_options.fingerprint_3_membered_rings = true;
        if (job_options.verbose) {
          cerr << "Will fingerprint 3 membered rings\n";
        }
      } else if (y == 4) {
        job_options.torsion_options.fingerprint_4_membered_rings = true;
        if (job_options.verbose) {
          cerr << "Will fingerprint 4 membered rings\n";
        }
      } else {
        cerr << "Unrecognised -y qualifier\n";
        return 1;
      }
    }
  }

  if (cl.option_present('T')) {
    if (! cl.value('T', job_options.ntest)) {
      cerr << "Invalid ntest specification (-T)\n";
      return 1;
    }
    if (job_options.verbose) {
      cerr << "Test mode, " << job_options.ntest << " tests per molecule\n";
    }
  }

  FileType input_type = FILE_TYPE_SMI;

  if (cl.option_present('i')) {
    if (! process_input_type(cl, input_type)) {
      return 1;
    }
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  IWString_and_File_Descriptor output(1);

  if (cl.option_present('f')) {
    job_options.work_as_filter = true;
    if (! TopologicalTorsionFilter(cl[0], job_options, output)) {
      cerr << "Work as filter failed\n";
      return 1;
    }
  } else {
    for (const char * fname : cl) {
      if (! TopologicalTorsion(fname, input_type, job_options, output)) {
        cerr << "Fatal error processing '" << fname << "'\n";
        return 1;
      }
    }
  }

  output.flush();

  if (job_options.verbose) {
    cerr << "Read " << job_options.molecules_read << " molecules\n";
  }

  return 0;
}

}  // namespace topological_torsion

int
main(int argc, char ** argv) {
  return topological_torsion::TopologicalTorsion(argc, argv);
}
