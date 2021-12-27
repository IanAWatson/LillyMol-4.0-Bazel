// Parse SMU TFDataRecords to extract smiles

#include "Foundational/data_source/tfdatarecord.h"

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
  
  extending_resizable_array<int> bt_count;
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
    coord_t x = apos.x();
    coord_t y = apos.y();
    coord_t z = apos.z();
    m.setxyz(i, x, y, z);
  }

  return 1;
}

void
Usage(int rc) {
  exit(rc);
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
  Command_Line cl(argc, argv, "v");
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
    for (int i = 0; i < options.bt_count.number_elements(); ++i) {
      cerr << options.bt_count[i] << " molecules had " << i << " bond topologoes\n";
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
