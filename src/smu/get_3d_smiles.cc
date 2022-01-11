// Parse SMU TFDataRecords to extract smiles

#include "Foundational/data_source/tfdatarecord.h"
#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/molecule.h"
#include "Molecule_Lib/smiles.h"

#include "smu/dataset.pb.h"

namespace get_3d_smiles {

using std::cerr;
struct JobOptions {
  int verbose = 0;
  // A list of the FATE values to process.
  resizable_array<int> fate;

  int items_read = 0;
  int items_written = 0;

  int no_bond_topologies = 0;

  int no_molecule_generated = 0;

  int no_optimized_geometry = 0;

  int atom_count_mismatch = 0;

  int starting_topology_not_found = 0;

  int write_smiles = 1;
  
  // The position at which the starting BT is found.
  extending_resizable_array<int> starting_topology_found;

  // Number of BondTopologies in the input.
  extending_resizable_array<int> bt_count;

  IWString_and_File_Descriptor stream_for_multiple_bt;
};

std::optional<Molecule>
MoleculeFromBondTopology(const BondTopology& bond_topology) {
  Molecule result;

  for (int i = 0; i < bond_topology.atoms().size(); ++i) {
    const BondTopology::AtomType atype = bond_topology.atoms(i);
    switch (atype) {
      case BondTopology::ATOM_H: {
          const Element * e = get_element_from_atomic_number(1);
          result.add(e);
        }
        break;
      case BondTopology::ATOM_C: {
        const Element * e = get_element_from_atomic_number(6);
        result.add(e);
      }
        break;
      case BondTopology::ATOM_N: {
        const Element * e = get_element_from_atomic_number(7);
        result.add(e);
      }
        break;
      case BondTopology::ATOM_NPOS: {
        const Element * e = get_element_from_atomic_number(7);
        result.add(e);
        result.set_formal_charge(result.natoms() - 1, 1);
      }
        break;
      case BondTopology::ATOM_O: {
        const Element * e = get_element_from_atomic_number(8);
        result.add(e);
      }
        break;
      case BondTopology::ATOM_ONEG: {
        const Element * e = get_element_from_atomic_number(8);
        result.add(e);
        result.set_formal_charge(result.natoms() - 1, -1);
      }
        break;
      case BondTopology::ATOM_F: {
        const Element * e = get_element_from_atomic_number(9);
        result.add(e);
      }
        break;
      case BondTopology::ATOM_UNDEFINED: {
        cerr << "Undefined atom\n";
        return std::nullopt;
      }
        break;
      default:
        cerr << "Unrecognised atom type " << atype << '\n';
        return std::nullopt;
    }
  }

  for (int i = 0; i < bond_topology.bonds().size(); ++i) {
    const BondTopology::Bond& bond = bond_topology.bonds(i);
    const atom_number_t a1 = bond.atom_a();
    const atom_number_t a2 = bond.atom_b();
    switch (bond.bond_type()) {
      case BondTopology::BOND_SINGLE:
        result.add_bond(a1, a2, SINGLE_BOND, 1);
        break;
      case BondTopology::BOND_DOUBLE:
        result.add_bond(a1, a2, DOUBLE_BOND, 1);
        break;
      case BondTopology::BOND_TRIPLE:
        result.add_bond(a1, a2, TRIPLE_BOND, 1);
        break;
      default:
        cerr << "Unrecognised bond type " << bond.bond_type() << '\n';
        return  std::nullopt;
    }
  }

  return result;
}

constexpr float kBohrToAngstrom = 0.529177f;

int
AddGeometry(const Geometry& geometry,
            JobOptions& options,
            Molecule& m) {
  const int matoms = m.natoms();
  if (matoms != geometry.atom_positions().size()) {
    options.atom_count_mismatch++;
    return 0;
  }

  for (int i = 0; i < matoms; ++i) {
    const Geometry::AtomPos& apos = geometry.atom_positions(i);
    coord_t x = apos.x() * kBohrToAngstrom;
    coord_t y = apos.y() * kBohrToAngstrom;
    coord_t z = apos.z() * kBohrToAngstrom;
    m.setxyz(i, x, y, z);
  }

  return 1;
}

void
Usage(int rc) {
  exit(rc);
}

int
MaybeWriteMultiBT(const Conformer& conformer,
                  JobOptions& options) {
  if (! options.stream_for_multiple_bt.is_open()) {
    return 1;
  }

  const int matoms = conformer.bond_topologies(0).atoms().size();
  int * ring_bond_count = new_int(matoms);
  std::unique_ptr<int[]> free_ring_bond_count(ring_bond_count);

  int ring_bond_count_differs = 0;
  int find_first_at = -1;

  const int number_bond_topologies = conformer.bond_topologies().size();

  // Compute the smiles once.
  IWString* smiles = new IWString[number_bond_topologies];
  std::unique_ptr<IWString[]> free_smiles(smiles);

  for (int i = 0; i < number_bond_topologies; ++i) {
    const BondTopology & bt = conformer.bond_topologies(i);
    std::optional<Molecule> maybe_mol = MoleculeFromBondTopology(bt);
    if (! maybe_mol) {
      cerr << "No Molecule from BT\n";
      continue;
    }
    AddGeometry(conformer.optimized_geometry(), options, *maybe_mol);
    smiles[i] = maybe_mol->smiles();
    if (bt.is_starting_topology()) {
      find_first_at = i;
    }

    for (int j = 0; j < matoms; ++j) {
      const int rbcj = maybe_mol->ring_bond_count(j);

      if (i == 0) {
        ring_bond_count[j] = rbcj;
      } else if (rbcj != ring_bond_count[j]) {
        ring_bond_count_differs = 1;
      }
    }
  }

  if (find_first_at < 0) {
    options.starting_topology_not_found++;
  } else {
    options.starting_topology_found[find_first_at]++;
  }

  if (! options.write_smiles) {
    return 1;
  }
  constexpr char sep = ' ';
  for (int i = 0; i < number_bond_topologies; ++i) {
    char is_first = (i == find_first_at ? '*' : '.');
    options.stream_for_multiple_bt << smiles[i] <<
         sep << conformer.conformer_id() <<
         sep << conformer.properties().optimized_geometry_energy().value() <<
         sep << conformer.fate() << 
         sep << is_first <<
         sep << find_first_at <<
         sep << ring_bond_count_differs <<
         sep << i << '\n';

    options.stream_for_multiple_bt.write_if_buffer_holds_more_than(8192);
  }

  return options.stream_for_multiple_bt.good();
}

int
Get3DSmiles(const Conformer& conformer,
            JobOptions& options,
            IWString_and_File_Descriptor& output) {
  options.bt_count[conformer.bond_topologies().size()]++;

  if (conformer.bond_topologies().size() == 0) {
    options.no_bond_topologies++;
    return 0;
  }

  if (conformer.bond_topologies().size() > 1) {
    MaybeWriteMultiBT(conformer, options);
  }

  if (! options.write_smiles) {
    return 1;
  }

  std::optional<Molecule> maybe_mol = MoleculeFromBondTopology(conformer.bond_topologies(0));
  if (! maybe_mol) {
    options.no_molecule_generated++;
    return 0;
  }

  if (conformer.optimized_geometry().atom_positions().size() == 0) {
    options.no_optimized_geometry++;
    return 1;
  }

  if (! AddGeometry(conformer.optimized_geometry(), options, *maybe_mol)) {
    return 0;
  }

  const char sep = ' ';

  output << maybe_mol->smiles() << sep << conformer.conformer_id() <<
         sep << conformer.properties().optimized_geometry_energy().value() <<
         sep << conformer.fate() << 
         sep << conformer.bond_topologies().size() << '\n';

  output.write_if_buffer_holds_more_than(8192);

  options.items_written++;

  return 1;
}

int
Get3DSmiles(const const_IWSubstring& data,
            JobOptions& options,
            IWString_and_File_Descriptor& output) {
  const std::string as_string(data.data(), data.length());
  Conformer conformer;
  if (! conformer.ParseFromString(as_string)) {
    cerr << "Cannot decode proto\n";
    return 0;
  }
  options.items_read++;

  return Get3DSmiles(conformer, options, output);
}

int
IndexLastNonZero(const extending_resizable_array<int>& values) {
  int rc = -1;
  for (int i = 0; i < values.number_elements(); ++i) {
    if (values[i] > 0) {
      rc = i;
    }
  }

  return rc;
}

int
Get3DSmiles(iw_tf_data_record::TFDataReader& reader,
            JobOptions& options,
            IWString_and_File_Descriptor& output) {
  while (reader.good() && ! reader.eof()) {
    std::optional<const_IWSubstring> data = reader.Next();
    if (! data) {
      return 1;
    }
    if (! Get3DSmiles(*data, options, output)) {
      cerr << "Cannot process item " << reader.items_read() << '\n';
      return 0;
    }
  }

  return 1;
}

int
Get3DSmiles(const char * fname,
            JobOptions& options,
            IWString_and_File_Descriptor& output) {
  iw_tf_data_record::TFDataReader reader(fname);
  if (! reader.good()) {
    cerr << "Cannot open " << fname << '\n';
    return 0;
  }

  return Get3DSmiles(reader, options, output);
}

int
Get3DSmiles(int argc, char** argv) {
  Command_Line cl(argc, argv, "vM:n");
  if (cl.unrecognised_options_encountered()) {
    cerr << "unrecognised_options_encountered\n";
    Usage(1);
  }

  if (cl.empty()) {
    cerr << "Insufficient arguments\n";
    Usage(1);
  }

  JobOptions options;
  options.verbose = cl.option_count('v');

  if (cl.option_present('n')) {
    options.write_smiles = 0;
    if (options.verbose)
      cerr << "Output suppressed\n";
  }

  if (cl.option_present('M')) {
    IWString mfile = cl.option_value('M');
    if (! mfile.ends_with(".smi")) {
      mfile << ".smi";
    }
    if (! options.stream_for_multiple_bt.open(mfile)) {
      cerr << "Cannot open stream for multiple BondTopology '" << mfile << "'\n";
      return 1;
    }

    if (options.verbose) {
      cerr << "Molecules with multiple BT written to " << mfile << '\n';
    }
  }

  set_append_coordinates_after_each_atom(1);

  IWString_and_File_Descriptor output(1);
  for (const char * fname : cl) {
    if (! Get3DSmiles(fname, options, output)) {
      cerr << "Fatal error processing " << fname << '\n';
      return 1;
    }
  }

  if (options.verbose) {
    cerr << "Read " << options.items_read << " wrote " << options.items_written << '\n';
    cerr << options.atom_count_mismatch << " atom count mismatch\n";
    cerr << options.no_bond_topologies << " no bond topologies\n";
    cerr << options.no_molecule_generated << " no_molecule_generated\n";
    cerr << options.no_optimized_geometry << " no_optimized_geometry\n";
    cerr << options.starting_topology_not_found << " starting topology not recovered\n";
    int nprint = IndexLastNonZero(options.bt_count);
    for (int i = 0; i <= nprint; ++i) {
      cerr << options.bt_count[i] << " molecules had " << i << " bond topologoes\n";
    }
    nprint = IndexLastNonZero(options.starting_topology_found);
    for (int i = 0; i <= nprint; ++i) {
      cerr << options.starting_topology_found[i] << " starting bond topologies found in position " << i << '\n';
    }
  }
  return 0;
}

}  // namespace get_3d_smiles


int
main(int argc, char ** argv)
{
  int rc = get_3d_smiles::Get3DSmiles(argc, argv);

  return rc;
}
