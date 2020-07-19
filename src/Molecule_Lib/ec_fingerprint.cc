#include <algorithm>
#include <limits>
#include <memory>

#include "Foundational/data_source/iwstring_data_source.h"

#include "ec_fingerprint.h"

namespace ec_fingerprint {

ShellInfo::ShellInfo(const Molecule& m, const int* include_atom, const atom_type_t* atom_type) :
       _m(m),
       _matoms(m.natoms()),
       _atom_type(atom_type)
{
  _status = new int[_matoms];
  if (NULL == include_atom) {
    std::fill_n(_status, _matoms, NOT_PROCESSED);
  } else {
    for (int i = 0; i < _matoms; ++i) {
      if (include_atom[i])
        _status[i] = NOT_PROCESSED;
      else
        _status[i] = EXCLUDED;
    }
  }

  _next_shell.make_room_for_extra_items(_matoms);

  return;
}

ShellInfo::~ShellInfo()
{
  delete [] _status;
}

void
ShellInfo::SetCentreAtom(const atom_number_t a)
{
  _a0 = a;

  for (int i = 0; i < _matoms; ++i) {
    if (EXCLUDED != _status[i])
      _status[i] = NOT_PROCESSED;
  }

  _status[a] = PROCESSING_COMPLETE;

  _next_shell.resize_keep_storage(0);

  return;
}

void
ShellInfo::GrabNextShell(Set_of_Atoms& next_shell) {

  // What used to be the next shall has been processed.
  _next_shell.set_vector(_status, PROCESSING_COMPLETE);

  next_shell.set_vector(_status, CURRENT_SHELL);

  _next_shell = std::move(next_shell);

  return;
}

ECFingerprint::ECFingerprint() {
  _min_radius = 0;
  _max_radius = 3;

  _bond_magic1 = 145;
  _bond_magic2 = 7001;
  _bond_magic3 = 23112;
  _bond_magic4 = 720114;

  _magic1 = 28231;
  _magic2 = 10101;
  _magic3 = 521;
  _magic4 = 9218;
  _magic5 = 2013;

  _precise = true;

  _additive_shell_formation = true;
}

atom_type_t
ECFingerprint::_BondConstant(const Bond& b) const {
#ifdef DEBUG_EC_FINGERPRINT
  cerr << "_BondConstant from " << b.a1() << " to " << b.a2() << " aromatic " << b.is_aromatic() << " single " << b.is_single_bond() << endl;
#endif
  if (b.is_aromatic())
    return _bond_magic4;

  if (b.is_single_bond())
    return _bond_magic1;

  if (b.is_double_bond())
    return _bond_magic2;

  if (b.is_triple_bond())
    return _bond_magic3;

  cerr << "EC_Fingerprint_:_BondConstant:what kind of bond ";
  b.debug_print(cerr);

  return 1;
}

void
ECFingerprint::_AddToRunningSum(ShellInfo& shell_info,
                       const atom_number_t a1,
                       const atom_type_t bond_constant,
                       const atom_number_t a2,
                       atom_type_t& running_sum) const {
#ifdef DEBUG_EC_FINGERPRINT
  cerr << "_AddToRunningSum from atom " << a1 << " type " << shell_info.atom_type(a1) << 
          " to " << a2 << " type " << shell_info.atom_type(a2) << endl;
#endif

  atom_type_t b;
  if (_precise) {
    b = shell_info.atom_type(a1) * _magic3 + _magic4 * bond_constant *
                       (_magic5 + shell_info.atom_type(a2));
//  cerr << "  _precise bond from " << a1 << " to " << a2 << " bond_constant " << bond_constant << " bit " << b << endl;
  } else {
    b = _magic4 * bond_constant + _magic5 * shell_info.atom_type(a2);
  }

  if (_additive_shell_formation)
    running_sum += b;
  else
    running_sum *= b;

#ifdef DEBUG_EC_FINGERPRINT
  cerr << "  bond from " << a1 << " to " << a2 << " bond_constant " << bond_constant << " bit " << b << " sum " << running_sum << endl;
#endif
}

int
ECBaseWithOutput::Open(IWString& fname)
{
  assert(! _output.is_open());

  if (! _output.open(fname))
  {
    cerr << "ECBaseWithOutput::Open:cannot open '" << fname << "'\n";
    return 0;
  }

  return 1;
}

JobParameters::JobParameters()
{
  function_as_tdt_filter = false;

  produce_output = true;
}

int
ProduceFingerprint::PrepareToProcess(Molecule& m)
{
  _sfc.clear();

  return 1;
}

int
ProduceFingerprint::DoAnyOutput(Molecule& m, const JobParameters& job_parameters,
                                IWString_and_File_Descriptor& output)
{
  if (! job_parameters.produce_output)
    return 1;

  if (! job_parameters.function_as_tdt_filter) {
    output << job_parameters.smiles_tag << m.smiles() << ">\n";
    output << job_parameters.identifier_tag << m.name() << ">\n";
  }

  IWString tmp;
  _sfc.daylight_ascii_form_with_counts_encoded(job_parameters.fingerprint_tag, tmp);
  output << tmp << "\n";

//if (! job_parameters.function_as_tdt_filter) {
//  output << "|\n";
//}

  return 1;
}

AtomMapCoverage::AtomMapCoverage()
{
  _matoms = 0;
  _coverage = nullptr;
}

AtomMapCoverage::~AtomMapCoverage()
{
  if (nullptr != _coverage)
    delete [] _coverage;
}

int
AtomMapCoverage::PrepareToProcess(Molecule& m)
{
  if (m.natoms() > _matoms)
  {
    if (nullptr != _coverage)
      delete [] _coverage;
    _matoms = m.natoms();
    _coverage = new int[_matoms];
  }

  std::fill_n(_coverage, m.natoms(), 0);

  return 1;
}

int
AtomMapCoverage::FingerprintingComplete(Molecule & m)
{
  assert (nullptr != _coverage);

  const int matoms = m.natoms();

  int min_coverage = _coverage[0];
  int max_coverage = _coverage[0];
  m.set_atom_map_number(0, _coverage[0]);

  for (int i = 1; i < matoms; ++i)
  {
    const int c = _coverage[i];
    m.set_atom_map_number(i, c);

    if (c < min_coverage)
      min_coverage = c;
    else if (c > max_coverage)
      max_coverage = c;
  }

  IWString tmp = m.name();
  tmp << " min " << min_coverage << " max " << max_coverage << " range " << (max_coverage - min_coverage);
  m.set_name(tmp);

  return 1;
}

void
WriteAllBits::Bit(const ShellInfo& shell_info, const atom_type_t running_sum, const int radius)
{
  const atom_number_t a0 = shell_info.a0();

  Molecule mcopy(shell_info.m());
  mcopy.set_isotope(a0, radius);
  
  _output << mcopy.smiles() << ' ' << shell_info.name() <<
         " bit " << running_sum << " atom type" << 
         shell_info.atom_type(a0) << " radius " << radius << "\n";

  _output.write_if_buffer_holds_more_than(8192);

  return;
}

void
ECCheckCollisions::Bit(const ShellInfo& shell_info,
                       const atom_type_t running_sum,
                       const int radius)
{
  const atom_number_t a0 = shell_info.a0();

  // Would be more efficient to use the Molecule inside shell_info, but this is safer.
  // This slows things down a LOT, even if we never look in the hash.
  Molecule mcopy(shell_info.m());
  mcopy.set_atom_map_number(a0, radius + 1);  // Because radius 0 will not get marked.

  std::tuple<IWString, atom_type_t, int> tmp{mcopy.smiles(), shell_info.atom_type(a0), radius};

  auto f = _bit_to_description.find(running_sum);
  if (f == _bit_to_description.end()) {
    _bit_to_description.emplace(running_sum, std::move(tmp));
    return;
  }

  // If same centre atom type and same radius, good enough.
  if (std::get<1>(tmp) == std::get<1>(f->second) &&
      std::get<2>(tmp) == std::get<2>(f->second)) {
    return;
  }

  _collisions_found++;

  _output << "Collision on bit " << running_sum << "\n";
  _output << mcopy.smiles() << ' ' << shell_info.name() << " atom type " << shell_info.atom_type(a0) << " radius " <<
  _output << "centre atom types " << std::get<1>(tmp) << ' ' << std::get<1>(f->second) << '\n';
  _output << "radii             " << std::get<2>(tmp) << ' ' << std::get<2>(f->second) << '\n';

  _output.write_if_buffer_holds_more_than(8192);

  return;
}

ECBuildPrecedent::ECBuildPrecedent()
{
}

ECBuildPrecedent::~ECBuildPrecedent()
{
}

void
ECBuildPrecedent::Bit(const ShellInfo& shell_info,
                      const atom_type_t running_sum,
                      const int radius)
{
  auto f = _bit_count.find(running_sum);
  
  if (f == _bit_count.end())  // Newly discovered bit.
  {
    Molecule mcopy(shell_info.m());
    const atom_number_t a0 = shell_info.a0();
    mcopy.set_atom_map_number(a0, radius + 1);
    bit_count_tuple_type tmp(mcopy.smiles(), shell_info.atom_type(a0),
                            radius, static_cast<count_type_t>(1));

    _bit_count.emplace(running_sum, std::move(tmp));
  } else {   
    std::get<3>(f->second)++;
  }

  return;
}

int
ECBuildPrecedent::WritePrecedentData(const char sep, const JobParameters& job_parameters, IWString& fname) const {
  IWString_and_File_Descriptor output;
  if (! output.open(fname))
  {
    cerr << "ECBuildPrecedent::WritePrecedentData:cannot open '" << fname << "'\n";
    return 0;
  }

  output << "#Atom type: \"" << job_parameters.atom_type_string << "\"\n";
  output << "#Bit Count Smiles Atype Radius\n";

  for (const auto iter : _bit_count) 
  {
    output << iter.first << sep <<                       // The bit
              std::get<3>(iter.second) << sep <<   // count
              std::get<IWString>(iter.second) << sep <<  // get<0>, smiles of first exemplar
              std::get<1>(iter.second) << sep <<         // atom type of centre atom
              std::get<2>(iter.second) <<                // radius
              "\n";
    output.write_if_buffer_holds_more_than(8192);
  }

  return 1;
}

ECUsePrecedent::ECUsePrecedent(const int max_radius) : _max_radius(max_radius)
{
  _count = new count_type_t[max_radius + 1];
  _atom = new atom_number_t[max_radius + 1];

  _missing_at_radius = new count_type_t[max_radius + 1];
  std::fill_n(_missing_at_radius, max_radius + 1, 0);

  return;
}

ECUsePrecedent::~ECUsePrecedent() {
  delete [] _count;
  delete [] _atom;
  delete [] _missing_at_radius;
}

int
ECUsePrecedent::PrepareToProcess(Molecule& m)
{
  std::fill_n(_count, _max_radius + 1, std::numeric_limits<count_type_t>::max());
  std::fill_n(_atom, _max_radius + 1, INVALID_ATOM_NUMBER);

  return 1;
}

//#define DEBUG_USE_PRECEDENT

void
ECUsePrecedent::Bit(const ShellInfo& shell_info,
                    const atom_type_t running_sum,
                    const int radius)
{
  const auto f = _precedent.find(running_sum);
  if (f == _precedent.end())
  {
#ifdef DEBUG_USE_PRECEDENT
    cerr << "Bit:no match for " << running_sum << ", radius " << radius << endl;
#endif
    _count[radius] = 0;
    _atom[radius] = shell_info.a0();
    _missing_at_radius[radius]++;
    return;
  }

#ifdef DEBUG_USE_PRECEDENT
  cerr << "Bit: we have info on " << running_sum << " radius " << radius << endl;
#endif
  const int db_radius = std::get<1>(f->second);
  if (db_radius != radius)  // collision
    return;

  const count_type_t db_count = std::get<0>(f->second);

#ifdef DEBUG_USE_PRECEDENT
  cerr << "No collision, count is " << db_count << endl;
#endif

  if (db_count < _count[radius])
  {
    _count[radius] = static_cast<int>(db_count);
    _atom[radius] = shell_info.a0();
  }

  return;
}

int
ECUsePrecedent::FingerprintingComplete(Molecule& m)
{
  IWString result;
  result << ' ';

  atom_number_t rarest_atom = INVALID_ATOM_NUMBER;
  count_type_t lowest_count = std::numeric_limits<count_type_t>::max();

//cerr << "ECUsePrecedent::FingerprintingComplete:radius " << _max_radius << endl;

  count_type_t previous_count = std::numeric_limits<count_type_t>::max();

  for (int r = 0; r <= _max_radius; ++r)
  {
    if (r > 0)
      result << ',';

#ifdef DEBUG_USE_PRECEDENT
    cerr << "Radius " << r << " count " << _count[r] << " ationm " << _atom[r] << endl;
#endif
    result << r << ',';
    if (std::numeric_limits<count_type_t>::max() == _count[r])
      result << '0';
    else
    {
      if (_count[r] < lowest_count)
      {
        lowest_count = _count[r];
        rarest_atom = _atom[r];
      }
      result << _count[r];
      if (_count[r] > previous_count)
      {
        cerr << "ECUsePrecedent::FingerprintingComplete:unexpected count increase, rad " <<
                (r-1) << " count " << previous_count <<
                " radius " << r << " count " << _count[r] << endl;
      }
      previous_count = _count[r];
    }
  }

  assert(rarest_atom != INVALID_ATOM_NUMBER);

#ifdef DEBUG_USE_PRECEDENT
  cerr << "lowest_count " << lowest_count << " rarest_atom " << rarest_atom << endl;
#endif

  m.set_atom_map_number(rarest_atom, lowest_count + 1);

  m.append_to_name(result);

  return 1;
}

int
ECUsePrecedent::DoAnyOutput(Molecule& m, const JobParameters& job_parameters,
                            IWString_and_File_Descriptor& output)
{
  if (! job_parameters.produce_output)  // No fingerprint output, must output smiles
  {
    output << m.smiles() << ' ' << m.name() << '\n';
  }
  else
  {
    output << job_parameters.smiles_tag << m.smiles() << ">\n";
    output << job_parameters.identifier_tag << m.name() << ">\n";
  }

  output.write_if_buffer_holds_more_than(8192);

  return output.good();
}

int
ECUsePrecedent::Report(std::ostream& output) const
{
  for (int i = 0; i <= _max_radius; ++i)
  {
    output << " radius " << i << " missing " << _missing_at_radius[i] << "\n";
  }

  return output.good();
}

int
ECUsePrecedent::ReadPrecedentData(IWString& fname)
{
  cerr << "ReadPrecedentData from '" << fname << "'\n";
  iwstring_data_source input;
  if (! input.open(fname))
  {
    cerr << "ECFingerprint::ReadPrecedentData:cannot open '" << fname << "'\n";
    return 0;
  }
    
  const_IWSubstring line;
  while (input.next_record(line))
  {
//  cerr << "Examining '" << line << "'\n";

    if (line.starts_with('#'))
      continue;

    if (! _ParsePrecedentRecord(line))
    {
      cerr << "ECFingerprint::ReadPrecedentData:invalid input '" << line << "'\n";
      return 0;
    }
  }
 
  return _precedent.size() > 0;
}

int
ECUsePrecedent::_ParsePrecedentRecord(const const_IWSubstring& line)
{
  int i = 0;
  const_IWSubstring token;

  if (! line.nextword(token, i))
  {
    cerr << "ECFingerprint::_ParsePrecedentRecord:cannot extract first token\n";
    return 0;
  }

  atom_type_t b;
  if (! token.numeric_value(b))
  {
    cerr << "efficient::_ParsePrecedentRecord:invalid bit number\n";
    return 0;
  }

  if (! line.nextword(token, i))
  {
    cerr << "ECFingerprint::_ParsePrecedentRecord:cannot extract second token\n";
    return 0;
  }

  count_type_t count;
  if (! token.numeric_value(count))
  {
    cerr << "ECFingerprint::_ParsePrecedentRecord:invalid count\n";
    return 0;
  }

  // Skip over smiles and central atom type to get radius token
  if (! line.nextword(token, i) || ! line.nextword(token, i) || ! line.nextword(token, i))
  {
    cerr << "ECFingerprint::_ParsePrecedentRecord:cannot extract tokens 3 and 4\n";
    return 0;
  }

  int radius;
  if (! token.numeric_value(radius))
  {
    cerr << "ECFingerprint::_ParsePrecedentRecord:Invalid radius\n";
    return 0;
  }

  std::tuple<count_type_t, int> tmp{count, radius};

  _precedent.emplace(b, std::move(tmp));

  return 1;
}


}  // namespace ec_fingerprint
