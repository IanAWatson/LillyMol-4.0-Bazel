// Holds the global settings defined in moleculeio.h

#include "moleculeio.h"

namespace moleculeio {

int _ignore_all_chiral_information_on_input = 0;

void
set_ignore_all_chiral_information_on_input(int s) {
  _ignore_all_chiral_information_on_input = s;

  return;
}

int
ignore_all_chiral_information_on_input() {
  return _ignore_all_chiral_information_on_input;
}

// Another way for fileconv to function is to use chiral information
// only if it is correct

int _ignore_incorrect_chiral_input = 0;

void
set_ignore_incorrect_chiral_input(int i) {
  _ignore_incorrect_chiral_input = i;

  return;
}

int
ignore_incorrect_chiral_input() {
  return _ignore_incorrect_chiral_input;
}

int _flush_files_after_writing_each_molecule = 0;

int
flush_files_after_writing_each_molecule() {
  return _flush_files_after_writing_each_molecule;
}

void
set_flush_files_after_writing_each_molecule(int i) {
  _flush_files_after_writing_each_molecule = i;
}

// When inputting MOLFILE's or TDT's we can optionally save all
// the non-connection-table records in the molecule's _extra_info
// array.

static int _read_extra_text_info = 0;

void
set_read_extra_text_info(int r)
{
  _read_extra_text_info = r;
}

int
read_extra_text_info() {
  return _read_extra_text_info;
}

// Does the extra text info get written or not

static int _write_extra_text_info = 0;

void set_write_extra_text_info(int w) {
  _write_extra_text_info = w;
}

int
write_extra_text_info() {
  return _write_extra_text_info;
}

// Especially with 3rd party molecules we can get files without a newline

char record_delimiter = '\n';

char 
input_file_delimiter() {
  return record_delimiter;
}

void
set_record_delimiter(char s) {
  record_delimiter = s;
}

int dos_mode = 1;    // Mar 2005. Change to default

int 
input_is_dos_mode() {
  return dos_mode;
}

void
set_dos_mode(int s) {
  dos_mode = s;
}

IWString file_scope_newline_string('\n');

void
generate_newline_string(IWString & newline_string)
{
  if (write_DOS_records()) {
    newline_string.resize_keep_storage(0);
    newline_string << static_cast<char>(13) << '\n';   // cannot put newline in src
  }
  else 
    newline_string = '\n';

  return;
}

const IWString & 
newline_string() {
  return file_scope_newline_string;
}

static int _write_DOS_records = 0;

void 
set_write_DOS_records(int s) {
  _write_DOS_records = s;

  generate_newline_string(file_scope_newline_string);
}

int
write_DOS_records() {
  return _write_DOS_records;
}

/*
  When reading MDL files we can look at the coordinates to perceive
  cis-trans bonds
*/

int _discern_cis_trans_bonds = 0;

int
discern_cis_trans_bonds()
{
  return _discern_cis_trans_bonds;
}

void
set_discern_cis_trans_bonds(int s)
{
  _discern_cis_trans_bonds = s;
}

int _discern_chirality_from_3d_coordinates = 0;

void
set_discern_chirality_from_3d_coordinates(int s)
{
  _discern_chirality_from_3d_coordinates = s;
}

int
discern_chirality_from_3d_coordinates()
{
  return _discern_chirality_from_3d_coordinates;
}

void
ResetFileScopeVariables() {
  _ignore_all_chiral_information_on_input = 0;
  _ignore_incorrect_chiral_input = 0;
  _flush_files_after_writing_each_molecule = 0;
  _read_extra_text_info = 0;
  _write_extra_text_info = 0;
  record_delimiter = '\n';
  dos_mode = 1;
  _write_DOS_records = 0;
  _discern_cis_trans_bonds = 0;
  _discern_chirality_from_3d_coordinates = 0;
  file_scope_newline_string = '\n';

  _discern_cis_trans_bonds = 0;
}

}  // namespace moleculeio
