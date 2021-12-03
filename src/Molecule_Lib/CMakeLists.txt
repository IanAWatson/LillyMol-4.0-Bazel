project(iwmolecule)
find_package(Protobuf REQUIRED)
include_directories(SYSTEM ${Protobuf_INCLUDE_DIRS})

set(Protobuf_IMPORT_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

#set(proto_src
#   ${CMAKE_CURRENT_SOURCE_DIR}/geometric_constraints.proto
#   ${CMAKE_CURRENT_SOURCE_DIR}/substructure.proto
#   ${CMAKE_CURRENT_SOURCE_DIR}/toggle_kekule_form.proto
#   ${CMAKE_CURRENT_SOURCE_DIR}/mol2graph.proto
#   ${PROTO_SRCS}
#   ${PROTO_HDRS}
#)
#protobuf_generate_cpp(PROTO_SRCS PROTO_HDRS ${proto_src})

file(GLOB proto_files ${CMAKE_CURRENT_SOURCE_DIR}/*.proto)
foreach(fname ${proto_files})
  get_filename_component(dir ${fname} DIRECTORY)
  get_filename_component(stem ${fname} NAME_WE)
  set(PROTO_HDRS "${PROTO_HDRS};${dir}/${stem}.pb.h")
  set(PROTO_SRCS "${PROTO_HDRS};${dir}/${stem}.pb.cc")
endforeach()

#set(proto_files
#   ${CMAKE_CURRENT_SOURCE_DIR}/geometric_constraints.proto
#   ${CMAKE_CURRENT_SOURCE_DIR}/substructure.proto
#   ${CMAKE_CURRENT_SOURCE_DIR}/toggle_kekule_form.proto
#   ${CMAKE_CURRENT_SOURCE_DIR}/mol2graph.proto
#)

add_custom_command(
  OUTPUT ${PROTO_HDRS} ${PROTO_SRCS}
  COMMAND ${Protobuf_PROTOC_EXECUTABLE}
  ARGS   -I ${CMAKE_CURRENT_SOURCE_DIR}
         --cpp_out ${CMAKE_CURRENT_SOURCE_DIR}
         ${proto_files}
  COMMENT "Running protoc on ${proto_files}"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
)

set(srcs
  ${CMAKE_CURRENT_SOURCE_DIR}/allowed_elements.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/aromatic.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/atom.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/atom_alias.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/atom_typing.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/bond.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/bond_list.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/careful_frag.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/charge_assigner.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/charge_calculation.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/charmm.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/chiral_centre.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/cif.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/cis_trans_bond.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/coordinates.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/coordinate_box.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/count_atoms_in_smiles.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/csubstructure.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/dihedral.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/donor_acceptor.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/element.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/element_hits_needed.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ematch.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/etrans.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/frag.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/geometric_constraints.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ISIS_Atom_List.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/is_actually_chiral.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwrcb.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwqry_wstats.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwrnm.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/iwsubstructure.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mdl_atom_record.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mdl_file_data.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mdl_molecule.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/misc2.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mol2graph.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mol2graph_proto.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/molecule.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/molecule_to_query.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/moleculeb.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/moleculed.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/moleculeh.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/moleculer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/molecule_smarts.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/mpr.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/numass.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ostream_and_type.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/output.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/parse_smarts_tmp.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/path.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/path_around_ring.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/path_scoring.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/pearlman.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ring_bond_iterator.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/rmele.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/rotbond_common.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/rwsubstructure.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/rwsubstructure_proto.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/set_of_atoms.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/smiles.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/smiles_support.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ss_atom_env.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ss_bonds.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ss_ring.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ss_ring_base.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ss_ring_sys.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/substructure_chiral.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/substructure_results.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/symm_class_can_rank.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/standardise.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/substructure_a.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/substructure_nmab.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/substructure_spec.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/substructure_env.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/symmetry.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/target.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/temp_detach_atoms.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/toggle_kekule_form.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/tokenise_atomic_smarts.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/unique.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/u3b.c
  ${CMAKE_CURRENT_SOURCE_DIR}/_substructure.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/_qry_wstats.cc
  ${PROTO_SRCS}
)

add_library(iwmolecule ${srcs})

enable_testing()
