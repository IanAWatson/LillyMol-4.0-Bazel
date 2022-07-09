#ifndef MOLECULE_LIB_MOLECULEIO_H_
#define MOLECULE_LIB_MOLECULEIO_H_

// Functions controlling aspects of reading writing and building molecules.

#include "Foundational/iwstring/iwstring.h"

namespace moleculeio {
int set_store_pdb_atom_information(int s);
int set_store_pdb_atom_information(int);
int set_use_stored_atom_information_when_writing_pdb_files(int);

const IWString & newline_string();

char input_file_delimiter();
void set_record_delimiter(char s);

int input_is_dos_mode();
void set_dos_mode(int);
void set_write_DOS_records(int);
int write_DOS_records();

int seek_to_from_command_line();
int max_offset_from_command_line();
int number_connection_table_errors_to_skip();
int skip_first_molecules();
int do_only_n_molecules();

//int valid_file_type(FileType);
////int suffix_for_file_type(FileType);

int discern_file_type_from_name(IWString const&);

int set_seek_to(long);

void set_discern_cis_trans_bonds(int);
int discern_cis_trans_bonds();
void set_discern_chirality_from_3d_coordinates(int);
int  discern_chirality_from_3d_coordinates();


void set_ignore_all_chiral_information_on_input(int s);
void set_ignore_incorrect_chiral_input(int s);
int ignore_incorrect_chiral_input();
int convert_from_mdl_charge(int);
int int3d(const_IWSubstring const&, int&, int&, int*);
int ignore_all_chiral_information_on_input();
int read_extra_text_info();
void set_read_extra_text_info(int s);
//int create_file_with_appropriate_name(const_IWSubstring const&, IWString&, FileType, int);
void set_flush_files_after_writing_each_molecule(int);
int flush_files_after_writing_each_molecule();
void set_write_extra_text_info(int);
int write_extra_text_info();
int string_to_file_type(const_IWSubstring const&);
int set_mol2_write_formal_charge_as_partial_charge(int);
int set_pdb_number_within_sequence(int);
int set_pdb_number_by_element_count(int);
int set_write_pdb_files_in_fragment_order(int);

//void set_tdt_identifier_dataitem(const const_IWSubstring &);
//void set_tdt_append_dataitem(const const_IWSubstring &);

void ResetFileScopeVariables();
}  // namespace moleculeio
#endif
