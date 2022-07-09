// Reading and writing ThorDataTree files
#include <iostream>

#include "Foundational/data_source/iwstring_data_source.h"

#include "molecule.h"
#include "moleculeio.h"
#include "rwmolecule.h"

using std::cerr;

// We have the ability to append one or more dataitems to the name field

static resizable_array_p<IWString> dataitems_to_append;

void
set_tdt_append_dataitem (const const_IWSubstring & a)
{
  IWString * tmp = new IWString(a);

  if (! tmp->ends_with('<'))
    *tmp += '<';

  dataitems_to_append.add(tmp);

  return;
}

static IWString identifier_dataitem ("PCN<");

void
set_tdt_identifier_dataitem(const const_IWSubstring & d)
{
  identifier_dataitem = d;

  if (! identifier_dataitem.ends_with('<'))
    identifier_dataitem += '<';

  return;
}
static int ignore_tdts_with_no_smiles = 0;

void
set_ignore_tdts_with_no_smiles (int s)
{
  ignore_tdts_with_no_smiles = s;
}

static IWString smiles_tag = "$SMI<";

void
set_smiles_tag (const const_IWSubstring & tag)
{
  smiles_tag = tag;

  if (! smiles_tag.ends_with('<'))
    smiles_tag += '<';

  return;
}

// When we are concatenating things onto the name, we can add either the
// contents of the tdt record, or the whole record itself

static int append_dataitem_content = 1;

void
set_tdt_append_dataitem_content (int s)
{
  append_dataitem_content = s;
}

static int
check_for_append(const const_IWSubstring & buffer,    // as read in
                 IWString & add_to_name)              // what will finally be appended
{
  int na = dataitems_to_append.number_elements();

  for (int i = 0; i < na; i++)
  {
    if (buffer.starts_with(*(dataitems_to_append[i])))
    {
      if (add_to_name.nchars())
        add_to_name += ' ';

      if (append_dataitem_content)    // need to remove tag and <>
      {
        add_to_name += buffer.substr(dataitems_to_append[i]->length());
        add_to_name.chop();    // remove trailing >
      }
      else
        add_to_name += buffer;
    }
  }

  return 1;
}

int
Molecule::read_molecule_tdt_ds (iwstring_data_source & input)
{
  assert(input.good());

  if (input.eof())
    return 0;

  input.set_skip_blank_lines(1);

  const_IWSubstring buffer;

  int got_structure = 0;
  int got_identifier = 0;

  IWString add_to_name;

  while (1)
  {
    EXTRA_STRING_RECORD(input, buffer, "read mol tdt");

    if (buffer.starts_with(smiles_tag) && 0 == got_structure)    // we only parse the first smiles in a TDT
    {
      buffer += smiles_tag.length();     // skip over '$SMI<'

      int close_angle_bracket = buffer.index('>');
      assert(close_angle_bracket >= 1);

      buffer.chop();
      if (! build_from_smiles(buffer.rawchars(), close_angle_bracket))
      {
        cerr << "Molecule::read_molecule_tdt_ds: Cannot interpret smiles\n";
        return 0;
      }

      got_structure = 1;

//    If this is a tdt in dump form, look for a PCN<> attribute.

      int i = buffer.find(">PCN<");
      if (i < 0)
        continue;

      const_IWSubstring zrest = buffer.substr(i + 5);

//    cerr << "The rest is '" << zrest << "'\n";

      if (zrest.nchars())
      {
        int j = zrest.find('>');
        zrest.iwtruncate(j);
        set_name(zrest);
      }

      return 1;     // molecule and name is all we can do in this case
    }

    if (buffer.starts_with(identifier_dataitem) && 0 == got_identifier) // PCN follows $SMI
    {
      IWString tmp = substr(buffer, identifier_dataitem.nchars());
      tmp.chop();
      set_name(tmp);

      got_identifier = 1;

      continue;
    }

    if ('|' == buffer && got_structure)
    {
      if (add_to_name.length())
        _molecule_name << ' ' << add_to_name;

      return 1;
    }

    if ('|' == buffer)
    {
      cerr << "Molecule::read_molecule_tdt_ds: no structure in TDT\n";
      if (0 == ignore_tdts_with_no_smiles)
        return 0;
      if (got_structure || got_identifier)
        return 0;

      continue;       // let's skip this TDT and look at the next one
    }

    if (moleculeio::read_extra_text_info())
    {
      IWString * tmp = new IWString(buffer);
      _text_info.add(tmp);
    }

    if (got_structure) {
      check_for_append(buffer, add_to_name);
    }
  }
}

/*
  We need to be able to specify the datatype name used when writing
  structures.
  This was first implemented to allow creation of the '$GRF' datatype
*/

static IWString datatype_name_for_structure_in_tdt_files = "$SMI<";

int 
set_datatype_name_for_structure_in_tdt_files(const char * new_dtname)
{
  assert(strlen(new_dtname));

  datatype_name_for_structure_in_tdt_files = new_dtname;

  if (! datatype_name_for_structure_in_tdt_files.ends_with('<')) {
    datatype_name_for_structure_in_tdt_files += '<';
  }

  return 1;
}

/*
  Writes all the data items except $SMI for a tdt.
*/

int
Molecule::_write_molecule_tdt_pcn(std::ostream & os,
                                  const IWString & comment) const
{
  if (_molecule_name.length()) {
    if (_molecule_name.contains('>') || _molecule_name.contains('|'))
      os << "PCN<" << '"' << _molecule_name << "\">" << moleculeio::newline_string();
    else
      os << "PCN<" << _molecule_name << '>' << moleculeio::newline_string();
  }

  if (comment.length() && comment != _molecule_name) {
    os << "REM<" << comment << '>' << moleculeio::newline_string();
  }

  int nt = _text_info.number_elements();
  for (int i = 0; i < nt; i++) {
    os << (*_text_info[i]) << moleculeio::newline_string();
  }

  os << '|' << moleculeio::newline_string();

  if (moleculeio::flush_files_after_writing_each_molecule()) {
    os.flush();
  }

  return os.good();
}

int
Molecule::write_molecule_tdt(std::ostream & os, const IWString & comment)
{
  assert(ok());
  assert(os.good());

  os << datatype_name_for_structure_in_tdt_files << smiles() << '>' << moleculeio::newline_string();

  return _write_molecule_tdt_pcn(os, comment);
}

int
Molecule::write_molecule_tdt_unique(std::ostream & os, const IWString & comment)
{
  assert(ok());
  assert(os.good());

  os << datatype_name_for_structure_in_tdt_files << unique_smiles() << '>' << moleculeio::newline_string();

  return _write_molecule_tdt_pcn(os, comment);
}

int
Molecule::write_molecule_tdt_nausmi(std::ostream & os, const IWString & comment)
{
  assert(ok());
  assert(os.good());

  os << datatype_name_for_structure_in_tdt_files << non_aromatic_unique_smiles() << '>' << moleculeio::newline_string();

  return _write_molecule_tdt_pcn(os, comment);
}
void
reset_tdt_file_scope_variables() {
  smiles_tag = "$SMI<";
  identifier_dataitem  = "PCN<";
  dataitems_to_append.resize(0);
  datatype_name_for_structure_in_tdt_files = "$SMI<";
  ignore_tdts_with_no_smiles = 0;
  append_dataitem_content = 1;
}
