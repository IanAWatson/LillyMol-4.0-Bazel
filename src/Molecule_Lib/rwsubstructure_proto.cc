#include <cstdint>

#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <string>
using std::cerr;
using std::endl;

/*
  Functions associated with reading and writing substructure query objects and
  such.
*/

#include <google/protobuf/message.h>
#include <google/protobuf/text_format.h>

#include "Foundational/iwmisc/misc.h"

#include "Molecule_Lib/atom_typing.h"
#include "Molecule_Lib/istream_and_type.h"
#include "Molecule_Lib/mdl_molecule.h"
#include "Molecule_Lib/misc2.h"
#include "Molecule_Lib/molecule_to_query.h"
#include "Molecule_Lib/parse_smarts_tmp.h"
#include "Molecule_Lib/path.h"
#include "Molecule_Lib/set_of_atoms.h"
#include "Molecule_Lib/smiles.h"
#include "Molecule_Lib/substructure.h"
#include "Molecule_Lib/substructure.pb.h"
#include "Molecule_Lib/target.h"

constexpr uint32_t no_limit = std::numeric_limits<uint32_t>::max();

// Having unbalanced braces in the code messes up matching in the editor.
constexpr char open_brace = '{';

using google::protobuf::Descriptor;
using google::protobuf::FieldDescriptor;
using google::protobuf::Reflection;

template <typename P>
void
SetValue(uint32_t value,
         const char * min_max,
         const char * name_stem,
         const Descriptor* descriptor,
         const Reflection* reflection,
         P& proto) {
  IWString name;
  name << min_max << name_stem;
  const FieldDescriptor* fd = descriptor->FindFieldByName(name.null_terminated_chars());
  reflection->SetUInt32(&proto, fd, value);
}

// Convert `values` to `proto`.
// The name for the proto field is `name_stem`, with
// min_name_stem and max_name_stem possibly being populated.
template <typename P>
int
SetProtoValues(const Min_Max_Specifier<int> & values,
               const char * name_stem,
               P& proto) {
  const Descriptor* descriptor = proto.GetDescriptor();
  const Reflection * reflection = proto.GetReflection();
  if (values.size() > 0) {
    const FieldDescriptor * fd = descriptor->FindFieldByName(name_stem);
    for (int v : values) {
      reflection->AddUInt32(&proto, fd, static_cast<uint32_t>(v));
    }
  }

  int v;
  if (values.min(v)) {
    SetValue<P>(v, "min_", name_stem, descriptor, reflection, proto);
  }
  if (values.max(v)) {
    SetValue<P>(v, "max_", name_stem, descriptor, reflection, proto);
  }

  return 1;
}
template <typename P>
int
SetProtoValues(const iwmatcher::Matcher<int> & matcher,
               const char * name_stem,
               P& proto) {
  const Descriptor* descriptor = proto.GetDescriptor();
  const Reflection * reflection = proto.GetReflection();
  if (matcher.number_elements() > 0) {
    const FieldDescriptor * fd = descriptor->FindFieldByName(name_stem);
    const resizable_array<int> values = matcher.ValuesMatched();
    for (int v : values) {
      reflection->AddUInt32(&proto, fd, static_cast<uint32_t>(v));
    }
  }

  int v;
  if (matcher.min(v)) {
    SetValue<P>(v, "min_", name_stem, descriptor, reflection, proto);
  }
  if (matcher.max(v)) {
    SetValue<P>(v, "max_", name_stem, descriptor, reflection, proto);
  }

  return 1;
}

template <typename P>
int
SetProtoIfSet(int value,
              int minval,
              const char * name,
              P& proto) {
  if (value < minval) {
    return 0;
  }

  const Descriptor* descriptor = proto.GetDescriptor();
  const Reflection * reflection = proto.GetReflection();
  const FieldDescriptor* fd = descriptor->FindFieldByName(name);
  reflection->SetUInt32(&proto, fd, value);

  return 1;
}

// Convert between the enumeration in the proto and the operators
// needed by IW_Logical_Expression.
int
AddOperator(const SubstructureSearch::Operator op,
            IW_Logical_Expression& destination)
{
  switch(op)
  {
    case SubstructureSearch::SS_OR:
      destination.add_operator(IW_LOGEXP_OR);
      break;
    case SubstructureSearch::SS_AND:
      destination.add_operator(IW_LOGEXP_AND);
      break;
    case SubstructureSearch::SS_XOR:
      destination.add_operator(IW_LOGEXP_XOR);
      break;
    case SubstructureSearch::SS_LP_AND:
      destination.add_operator(IW_LOGEXP_LOW_PRIORITY_AND);
      break;
    default:
      cerr << "AddOperator:unrecognized operator " << op << endl;
      return 0;
  }

  return 1;
}

//static int set_element_hits_needed_during_molecule_to_query = 1;

/*
  In several places we need to quickly ascertain whether or not an msi_object
  specifies a rejection or not
*/


// Elements can be specified as either numeric (atomic number) or string (atomic symbol)
// forms, or both.
static int
fetch_elements(const SubstructureSearch::SubstructureAtomSpecifier & proto,
               resizable_array<const Element*>& ele,
               resizable_array<int>& element_unique_id)
{
  if (proto.atomic_symbol_size() > 0) {
    for (const auto& s : proto.atomic_symbol()) 
    {
      const_IWSubstring atomic_symbol(s);
      const Element * e = get_element_from_symbol_no_case_conversion(atomic_symbol);
      if (nullptr == e && auto_create_new_elements())
        e = create_element_with_symbol(atomic_symbol);

      if (nullptr == e)
      {
        cerr << "fetch_elements:no element for symbol '" << atomic_symbol << "'\n";
        return 0;
      }
      
      ele.add_if_not_already_present(e);
      element_unique_id.add_if_not_already_present(e->unique_id());
    }
  }

  if (proto.atomic_number_size() > 0) {
    for (const auto z : proto.atomic_number())
    {
      const Element * e = get_element_from_atomic_number(z);
      if (nullptr == e)
      {
        cerr << "fetch_elements:no element for atomic number " << z << endl;
        return 0;
      }
      
      ele.add_if_not_already_present(e);
      element_unique_id.add_if_not_already_present(e->unique_id());
    }
  }

  assert (ele.ok());

  return 1;
}


#define MIN_NOT_SPECIFIED 9393
#define MAX_NOT_SPECIFIED 14163

/*
  Process all the ATTRIBUTE_NAME attributes in an msi_object.
  Append resulting int values to SPECIFIER.
*/

template <typename T>
int
AppendIntValue(const T value,
               T min_val, T max_val,
               Min_Max_Specifier<T> & specifier)
{
  assert (specifier.ok());

  if (value < min_val || value > max_val)
  {
    cerr << "AppendIntValue:out of range " << value << " must bt btw " << min_val << " and " << max_val << endl;
    return 0;
  }

  specifier.add(value);

  return 1;
}

template <typename T, typename M>
bool
ReallyGruesome(const ::google::protobuf::RepeatedField<T> values,
               const std::optional<T>(minval), const std::optional<T>(maxval),
               const T min_allowed, const T max_allowed,
               const char * attribute,
               M& target)
{
  if (minval.has_value())
    target.set_min(minval.value());

  if (maxval.has_value())
    target.set_max(maxval.value());

  for (const auto v : values) {
    target.add(v);
  }

  return true;
}

#define MAYBEMIN(p, attribute, T) (p.has_min_ ## attribute() ? std::optional<T>(p.min_ ## attribute()) : std::optional<T>())
#define MAYBEMAX(p, attribute, T) (p.has_max_ ## attribute() ? std::optional<T>(p.max_ ## attribute()) : std::optional<T>())

#define GETVALUES(p, attribute, lowest_allowed, max_allowed)    ReallyGruesome<uint32_t>(p.attribute(), MAYBEMIN(p, attribute, uint32_t), MAYBEMAX(p, attribute, uint32_t), lowest_allowed, max_allowed, #attribute, _ ## attribute)
#define GETVALUESINT(p, attribute, lowest_allowed, max_allowed) ReallyGruesome<int32_t>(p.attribute(),  MAYBEMIN(p, attribute, int32_t), MAYBEMAX(p, attribute, int32_t), lowest_allowed, max_allowed, #attribute, _ ## attribute)
#define GETFLOATVALUES(p, attribute, lowest_allowed, max_allowed) ReallyGruesome<float>(p.attribute(),  MAYBEMIN(p, attribute, float32), MAYBEMAX(p, attribute, float), lowest_allowed, max_allowed, #attribute, _ ## attribute)

using ::google::protobuf::RepeatedField;

template <typename T>
bool
really_gruesome(const RepeatedField<T> values, const T minval, const T maxval,
                const T min_value_allowed, const T max_value_allowed,
                Min_Max_Specifier<T>& destination, const char * name)
{
  if (values.size() > 0)
  {
    for (const T v : values)
    {
      if (v < min_value_allowed || v > max_value_allowed)
      {
        cerr << "really_gruesome::value of " << name << " " << v << " out of range, must be btw " << min_value_allowed << " and " << max_value_allowed << "\n";
        return false;
      }

      destination.add_if_not_already_present(v);
    }

    return true;
  }

  if (minval < min_value_allowed)
  {
    cerr << "really_gruesome::value of " << name << " " << minval << " out of range, must be less than " << min_value_allowed << "\n";
    return false;
  }

  if (minval > max_value_allowed)
  {
    cerr << "really_gruesome::value of " << name << " " << maxval << " out of range, must be more than " << max_value_allowed << "\n";
    return false;
  }

  destination.set_min(minval);
  destination.set_min(maxval);

  return true;
}

/*
  Process a MinMaxSpecifierUint or MinMaxSpecifierint proto (type P).
*/



int
Substructure_Atom_Specifier::construct_from_proto(const SubstructureSearch::SubstructureAtomSpecifier & proto)
{
  assert (ok());

//cerr << "SubstructureAtomSpecifier::construct_from_proto:from " << proto.ShortDebugString() << endl;
  if (! fetch_elements(proto, _element, _element_unique_id))
    return 0;

  if (proto.has_preference_value())
    _preference_value = proto.preference_value();

  // Should relate these to the min and max reasonable formal charge values.
  if (!GETVALUESINT(proto, formal_charge, -12, 12))
    return 0;

  if (! GETVALUES(proto, ncon, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, nbonds, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, nrings, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, ring_bond_count, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, ncon2, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, hcount, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, ring_size, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, attached_heteroatom_count, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, lone_pair_count, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, unsaturation, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, daylight_x, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, isotope, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, aryl, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, vinyl, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, fused_system_size, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, symmetry_degree, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, fused_system_size, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, scaffold_bonds_attached_to_ring, 0, no_limit))
    return 0;

  if (proto.has_symmetry_group())
  {
    _symmetry_group = proto.symmetry_group();
    if (_symmetry_group == 0) {
      cerr << "SubstructureAtomSpecifier::construct_from_proto:_symmetry_group must be > 0\n";
      return 0;
    }
  }

  if (proto.has_aromatic())
  {
    if (proto.aromatic())
      _aromaticity = AROMATIC;
    else
      _aromaticity = NOT_AROMATIC;
  }

  if (proto.has_match_spinach_only())
    _match_spinach_only = proto.match_spinach_only();

  if (proto.has_chirality())
    _chirality = proto.chirality();

  if (!GETVALUES(proto, heteroatoms_in_ring, 0, no_limit))
    return 0;

  // Min should be related to min aromatic ring size.
  if (!GETVALUES(proto, aromatic_ring_size, 4, no_limit))
    return 0;

  if (!GETVALUES(proto, aliphatic_ring_size, 3, no_limit))
    return 0;

  if (proto.has_all_rings_kekule()) 
    _all_rings_kekule = proto.all_rings_kekule();

  if (proto.has_user_atom_type())
    _userAtomType = proto.user_atom_type();

  if (proto.has_atom_type())
    _atom_type = proto.atom_type();
 
  (void) count_attributes_specified();

  assert (ok());

  return 1;
}

int
Substructure_Atom::_create_preference_from_proto(const SubstructureSearch::SubstructureAtomSpecifier& proto)
{
  if (! proto.has_preference_value())
  {
    cerr << "Substructure_Atom_Specifier::_create_preference_from_proto:no preference value\n";
    return 0;
  }
//cerr << "_create_preference_from_proto " << proto.ShortDebugString() << endl;

  Substructure_Atom_Specifier * a = new Substructure_Atom_Specifier;

  if (! a->construct_from_proto(proto))
  {
    cerr << "Substructure_Atom_Specifier::_create_preference_from_proto:invalid preference specification\n";
    delete a;
    return 0;
  }

  _preferences.add(a);

  return 1;
}

/*
  There is a common function of examining an object and extracting its
  OPERATOR attribute and then adding the appropriate operator to a
  logical_expression
*/

static int
ExtractOperator(const google::protobuf::RepeatedField<int>& ops,
                const int ops_needed,
                IW_Logical_Expression& logexp,
                int default_operator,
                const char * caller)
{
  if (ops.empty())
  {
    for (int i = 0; i < ops_needed; ++i) {
      logexp.add_operator(default_operator);
    }

    return 1;
  }

  if (ops.size() != ops_needed)
  {
    cerr << "extract_operator::operator count mismatch, needed " << ops_needed << " got " << ops.size() << "\n";
    return 0;
  }

  for (const auto op : ops)
  {
    switch (op) {
      case SubstructureSearch::Operator::SS_OR:
        logexp.add_operator(IW_LOGEXP_OR);
        break;
      case SubstructureSearch::Operator::SS_AND:
        logexp.add_operator(IW_LOGEXP_AND);
        break;
      case SubstructureSearch::Operator::SS_XOR:
        logexp.add_operator(IW_LOGEXP_XOR);
        break;
      case SubstructureSearch::Operator::SS_LP_AND:
        logexp.add_operator(IW_LOGEXP_LOW_PRIORITY_AND);
        break;
      default:
        cerr << "extract_operator:unrecognized operator " << op << "\n";
        return 0;
    };
  }

  return 1;
}


/*int
Substructure_Atom::_add_bond (const msi_object & msi,
                              extending_resizable_array<Substructure_Atom *> & completed)
{
  Substructure_Bond * b = new Substructure_Bond;
  if (! b->construct_from_msi_object (msi, completed))
  {
    delete b;
    return 0;
  }

  if (this == b->a())
  {
    cerr << "Substructure_Atom::_add_bond: atom bonded to itself\n";
    cerr << (*msi);
    return 0;
  }

  return _add_bond (b);
}*/

/*
  Parse input that looks like:

    - 1
    -!@ 1 2 3

    The bond, then a set of atoms.
*/

static int
parse_smarts_bond_attribute(const IWString& attribute,
                            const atom_number_t must_not_be,
                            Set_of_Atoms & other_end,
                            IWString & bond_smarts)
{
  const int nw = attribute.nwords();

  if (nw < 2)
  {
    cerr << "parse_smarts_bond_attribute: bond smarts must have at least two tokens '" << attribute << "'\n";
    return 0;
  }

  const int natoms = nw-1;  // first token is the bond

  other_end.resize_keep_storage(natoms);

  int i = 0;
  attribute.nextword(bond_smarts, i);

  const_IWSubstring s;
  while (attribute.nextword(s, i))
  {
    atom_number_t zatom;
    if (! s.numeric_value(zatom) || zatom < 0)
    {
      cerr << "parse_smarts_bond_attribute: invalid atom number '" << s << "'\n";
      return 0;
    }

    if (zatom == must_not_be) 
    {
      cerr << "parse_smarts_bond_attribute:self bonds not allowed\n";
      return 0;
    }

    other_end.add_if_not_already_present(zatom);
  }

  return bond_smarts.length();
}

bond_type_t
BondTypeFromSSBond(const google::protobuf::RepeatedField<int>& bond_type)
{
  bond_type_t btype = 0;

  for (const auto b : bond_type)
  {
    switch (b) {
      case SubstructureSearch::BondType::SS_SINGLE_BOND:
        btype |= SINGLE_BOND;
        break;
      case SubstructureSearch::BondType::SS_DOUBLE_BOND:
        btype |= DOUBLE_BOND;
        break;
      case SubstructureSearch::BondType::SS_TRIPLE_BOND:
        btype |= TRIPLE_BOND;
        break;
      case SubstructureSearch::BondType::SS_AROMATIC_BOND:
        btype |= AROMATIC_BOND;
        break;
      default:
        cerr << "Substructure_Environment::_process_attachment_bonds:unrecognised bond type " << b << "\n";
        return 0;
    };
  }

  return btype;
}

int
Substructure_Atom::_process_substructure_bond(const SubstructureSearch::SubstructureBond& bond,
                                        extending_resizable_array<Substructure_Atom *> & completed)
{
  if (0 == bond.bond_type_size())
  {
    cerr << "Substructure_Atom::_process_substructure_bond:no bonds\n";
    return 0;
  }

  if (! bond.has_other_end())
  {
    cerr << "Substructure_Atom::_process_substructure_bond:no other end\n";
    return 0;
  }

  const atom_number_t other_end = bond.other_end();
  if (nullptr == completed[other_end])
  {
    cerr << "Substructure_Atom::_process_substructure_bond:other end atom " << other_end << " not defined\n";
    return 0;
  }

  const bond_type_t btype = BondTypeFromSSBond(bond.bond_type());

  if (0 == btype)
  {
    cerr << "Substructure_Atom::_process_substructure_bond:invalid bond type specification\n";
    cerr << bond.ShortDebugString() << endl;
    return 0;
  }

  // cerr << "Substructure_Atom::_process_substructure_bond:adding bond from " << _unique_id << " to " << other_end << " completed " << completed[other_end]->unique_id() << '\n';

  Substructure_Bond * b = new Substructure_Bond;

  b->set_bond_type(btype);

  b->set_atom(completed[other_end]);

  return _add_bond(b);
}

int
Substructure_Atom::_process_smarts_bond(const IWString& input,
                                        extending_resizable_array<Substructure_Atom *> & completed)
{
//cerr << "Parsing smarts bond '" << input << "'\n";
  Set_of_Atoms other_end;
  IWString bond_smarts;

  if (! parse_smarts_bond_attribute(input, _unique_id, other_end, bond_smarts))
  {
    cerr << "Substructure_Atom::_process_attribute_smarts_bond:invalid bond smarts " << input << endl;
    return 0;
  }

  for (const atom_number_t o : other_end)
  {
    if (nullptr == completed[o])
    {
      cerr << "Substructure_Atom::_process_attribute_smarts_bond: atom " << o << " has not been defined\n";
      return 0;
    }

    if (this == completed[o])
    {
      cerr << "Substructure_Atom::_process_attribute_smarts_bond: atom bonded to itself\n";
      return 0;
    }

    Substructure_Bond * b = new Substructure_Bond;

    int characters_processed;
    if (! b->construct_from_smarts(bond_smarts.rawchars(), bond_smarts.length(), characters_processed))
    {
      cerr << "parse_smarts_bond_attribute:cannot parse bond smarts " << input << "\n";
      delete b;
      return 0;
    }

    if (characters_processed != bond_smarts.length())
    {
      cerr << "parse_smarts_bond_attribute:extra junk at end of bond smarts " << input << "\n";
      cerr << characters_processed << " characters processed\n";
      delete b;
      return 0;
    }

    b->set_atom(completed[o]);

    _add_bond(b);
  }

  return 1;
}

/*
  An attribute bond will have the other end as the first int value
*/

// substructure_bond: "1 2 3 -!@"

int
Substructure_Environment::_process_attachment_via_substructure_bond(const SubstructureSearch::EnvironmentAttachment& proto,
                                             extending_resizable_array<Substructure_Atom *> & completed)
{
  Set_of_Atoms attachments;
  IWString bond_smarts;
  if (!parse_smarts_bond_attribute(proto.substructure_bond(), -1, attachments, bond_smarts)) {
    cerr << "Substructure_Environment::_process_attachment_via_substructure_bond:invalid substructure bond " << proto.substructure_bond() << "'\n";
    return 0;
  }

  for (const auto zatom : attachments) {
    if (nullptr == completed[zatom]) {
      cerr << "Substructure_Environment::_process_attachment_via_substructure_bond:invalid atom " << zatom <<"\n";
      return 0;
    }
    _possible_parents.add(completed[zatom]);
  }

  int characters_processed;
  if (! _bond.construct_from_smarts(bond_smarts.rawchars(), bond_smarts.length(), characters_processed))
  {
    cerr << "Substructure_Environment::_process_attribute_bond:invalid bond smarts " << bond_smarts << endl;
    return 0;
  }

  if (characters_processed != bond_smarts.length()) {
    cerr << "Substructure_Environment::_process_attachment_via_substructure_bond:invalid bond smarts '" << bond_smarts << "'\n";
    return 0;
  }

  return 1;
}

/*
  The environment has just one bond. If you want the environment
  multiply attached, add a ring bond to the first atom in the env
*/

int
Substructure_Environment::_process_how_connected(const SubstructureSearch::SubstructureEnvironment& proto,
                                             extending_resizable_array<Substructure_Atom *> & completed)
{
  if (!proto.has_attachment()) {
    cerr << "Substructure_Environment::_process_how_connected:no attachment\n";
    return 0;
  }

  const SubstructureSearch::EnvironmentAttachment& attachment = proto.attachment();

  if (attachment.attachment_point_size() > 0)
    return _process_attachment_bonds(attachment, completed);
  else if (attachment.has_substructure_bond())
    return _process_attachment_via_substructure_bond(attachment, completed);

  cerr << "Substructure_Atom::_process_how_connected:environment must have attachment\n";
  return 0;
}

int 
Substructure_Environment::_process_attachment_bonds(const SubstructureSearch::EnvironmentAttachment& proto,
                                             extending_resizable_array<Substructure_Atom *> & completed)
{
  if (0 == proto.attachment_point_size() || 0 == proto.btype_size())
  {
    cerr << "Substructure_Environment::_process_attachment_bonds:must have both attachment and bond attributes\n";
    cerr << proto.ShortDebugString() << endl;
    return 0;
  }

  Set_of_Atoms attachments;
  for (const auto a : proto.attachment_point())
  {
    if (attachments.add_if_not_already_present(a)) {
      if (nullptr == completed[a]) {
        cerr << "Substructure_Environment::_process_attachment_bonds:unrecognized atom number " << a << "\n";
        return 0;
      }
      _possible_parents.add(completed[a]);
    }
  }

  const bond_type_t btype = BondTypeFromSSBond(proto.btype());

  _bond.set_bond_type(btype);

  return 1;
}

/*
  When creating an environment, the environment components will
  register themselves as children. They are really environment components,
  so after they are done, we move any extra children to components.
*/

int
Substructure_Atom::_create_environment_from_proto(const SubstructureSearch::SubstructureAtomEnvironment& proto,
                                           extending_resizable_array<Substructure_Atom *> & completed)
{
  assert (_environment.empty());

  int initial_children = _children.number_elements();

  for (const auto& substructure_atom : proto.substructure_atom())
  {
    Substructure_Atom * a = new Substructure_Atom;
    if (! a->construct_from_proto(substructure_atom, completed))
    {
      cerr << "_create_environment_from_proto: cannot create environment component " << substructure_atom.ShortDebugString() << endl;
      delete a;
      return 0;
    }
  }

  int final_children = _children.number_elements();
  if (initial_children == final_children)
  {
    cerr << "Substructure_Atom::create_env_from_proto: children not attached to root!\n";
    cerr << proto.ShortDebugString();
    return 0;
  }

  while (_children.number_elements() > initial_children)
    _environment.transfer_in(_children, initial_children);

  return _environment.number_elements();
}


template <typename T, typename F>
bool
AssignValue(T proposed, T & target, F fn, const char * variable_name)
{
  if (! fn(proposed))
  {
    cerr << "AssignValue:invalid " << variable_name  << ' ' << proposed << "\n";
    return false;
  }

  target = proposed;

  return true;
}
 
/*
  The sub-objects of a Substructure_Atom will be either bonds, or
  the Substructure_Atoms which are its components.
  Note that components are not allowed to have bonds.
*/

int
Substructure_Atom::construct_from_proto(const SubstructureSearch::SubstructureAtom& proto,
                                        extending_resizable_array<Substructure_Atom *> & completed)
{
//#define DEBUG_COMPLTED
#ifdef DEBUG_COMPLTED
  cerr << "Building Substructure_Atom " << proto.ShortDebugString() << endl;
  for (int i = 0; i < 6; i++)
  {
    cerr << "Completed[" << i << "] = " << completed[i] << endl;
  }
#endif

  if (!proto.has_id())
  {
    cerr << "Substructure_Atom::_construct_from_proto:no id " << proto.ShortDebugString() << "\n";
    return 0;
  }

  _unique_id = proto.id();

  if (proto.has_match_as_match())
    _match_as_match_or_rejection = proto.match_as_match();

  if (nullptr != completed[_unique_id])
  {
    cerr << "Substructure_Atom::construct_from_msi_object: Yipes, atom " << _unique_id << " already allocated\n";
    cerr << proto.ShortDebugString() << endl;
    return 0;
  }

  if (proto.has_initial_atom_number())
    _initial_atom_number = proto.initial_atom_number();
  else
    _initial_atom_number = _unique_id;   // Not sure if this is the right thing to do or not.

#ifdef OR_ID_NO_LONGER_PROCESSED
  if (proto.has_or_id())
  {
    _or_id = proto.or_id();

    if (0 == _or_id)
    {
      cerr << "Substructure_Atom::_construct_from_proto: or id must be a positive number\n";
      return 0;
    }
  }
#endif

  if (proto.has_ring_id() && ! AssignValue(proto.ring_id(), _ring_id, [](int r) { return r > 0;}, "ring_id"))
    return 0;

  if (proto.has_fragment_id() && ! AssignValue(proto.fragment_id(), _fragment_id, [](int f) { return f > 0;}, "fragment_id"))
    return 0;

  if (proto.has_global_match_id()) {
    _global_match_id = proto.global_match_id();
  }

  if (proto.has_fused_system_id())
    _fused_system_id = proto.fused_system_id();

  if (proto.has_text_identifier())
    _text_identifier = proto.text_identifier();

  if (proto.has_numeric_value())
    _numeric_value = proto.numeric_value();

  if (proto.has_include_in_embedding())
    _include_in_embedding = proto.include_in_embedding();

  if (proto.has_sum_all_preference_hits())
    _sum_all_preference_hits = proto.sum_all_preference_hits();

  if (proto.atom_properties_size() == 1)
  {
    const auto& spec = proto.atom_properties(0);
    Substructure_Atom_Specifier* tmp = this;
    if (! tmp->construct_from_proto(spec))
    {
      cerr << "Substructure_Atom::_construct_from_proto:cannot parse Substructure_Atom_Specifier\n";
      cerr << spec.ShortDebugString() << "\n";
      return 0;
    }
  }
  else
  {
    for (int i = 0; i < proto.atom_properties_size(); ++i)
    {
      const auto& spec = proto.atom_properties(i);

      if (i > 0 && ! spec.has_logical_operator()) {
        cerr << "Substructure_Atom::_construct_from_proto:second Substructure_Atom_Specifier's must have operator\n";
        cerr << spec.ShortDebugString();
        return 0;
      }

      std::unique_ptr<Substructure_Atom_Specifier> tmp = std::make_unique<Substructure_Atom_Specifier>();
      if (!tmp->construct_from_proto(spec))
      {
        cerr << "Substructure_Atom::_construct_from_proto:cannot parse Substructure_Atom_Specifier\n";
        cerr << spec.ShortDebugString() << "\n";
        return 0;
      }

      _components.add(tmp.release());
      if (i > 0)
        AddOperator(spec.logical_operator(), _operator);
    }
  }

// For simplicity, one can specify only one of atom smarts, smiles and smarts

  const int smiles_and_smarts = proto.has_atom_smarts() +
                                proto.has_smarts() +
                                proto.has_smiles();

  if (_components.number_elements() > 0 && smiles_and_smarts)
  {
    cerr << "Substructure_Atom::_construct_from_proto:cannot mix Substructure_Atom_Specifier with smiles/smarts specifications\n";
    cerr << proto.ShortDebugString() << endl;
    return 0;
  }

  if (smiles_and_smarts > 1)
  {
    cerr << "Substructure_Atom::_construct_from_proto:can have only one of 'atom_smarts', 'smarts' or 'smiles'\n";
    cerr << proto.ShortDebugString() << endl;
    return 0;
  }

  if (proto.has_atom_smarts())
  {
    IWString csmarts = proto.atom_smarts();

    if (! (csmarts.starts_with('[') && csmarts.ends_with(']')))
    {
      cerr << "Substructure_Atom::_construct_from_proto: atomic smarts must be within [] '" << csmarts << "'\n";
      return 0;
    }

    if (! construct_from_smarts_token(csmarts))
    {
      cerr << "Substructure_Atom::_construct_from_proto: cannot interpret smarts '" << csmarts << "'\n";
      return 0;
    }
  }

  if (proto.has_smiles())
  {
    const IWString smiles = proto.smiles();
    if (! parse_smiles_specifier(smiles))
    {
      cerr << "Substructure_Atom::_construct_from_proto:invalid smiles '" << proto.smiles() << "'\n";
      return 0;
    }
  }

  if (proto.has_smarts())
  {
    const_IWSubstring smarts = proto.smarts();
    if (! parse_smarts_specifier(smarts))
    {
      cerr << "Substructure_Atom::_construct_from_proto:invalid smarts '" << proto.smarts() << "'\n";
      return 0;
    }
  }

  completed[_unique_id] = this;

  for (const auto& env : proto.environment())
  {
    if (! _create_environment_from_proto(env, completed))
    {
      cerr << "Substructure_Atom::_construct_from_proto:invalid environment\n";
      cerr << env.ShortDebugString() << "\n";
      return 0;
    }
  }

  for (const auto & pref : proto.preference())
  {
    if (! _create_preference_from_proto(pref))
    {
      cerr << "Substructure_Atom::_construct_from_proto:invalid preference\n";
      cerr << pref.ShortDebugString() << "\n";
      return 0;
    }
  }

  if (proto.single_bond_size() > 0)
  {
    if (! _add_bonds(proto.single_bond(), SINGLE_BOND, completed))
    {
      cerr << "Substructure_Atom::construct_from_proto:invalid single bond\n";
      cerr << proto.ShortDebugString() << endl;
      return 0;
    }
  }

  if (proto.double_bond_size() > 0)
  {
    if (! _add_bonds(proto.double_bond(), DOUBLE_BOND, completed))
    {
      cerr << "Substructure_Atom::construct_from_proto:invalid double bond\n";
      cerr << proto.ShortDebugString() << endl;
      return 0;
    }
  }

  if (proto.triple_bond_size() > 0)
  {
    if (! _add_bonds(proto.triple_bond(), TRIPLE_BOND, completed))
    {
      cerr << "Substructure_Atom::construct_from_proto:invalid triple bond\n";
      cerr << proto.ShortDebugString() << endl;
      return 0;
    }
  }

  if (proto.aromatic_bond_size() > 0)
  {
    if (! _add_bonds(proto.aromatic_bond(), AROMATIC_BOND, completed))
    {
      cerr << "Substructure_Atom::construct_from_proto:invalid aromatic bond\n";
      cerr << proto.ShortDebugString() << endl;
      return 0;
    }
  }

  if (proto.has_bond_smarts())
  {
    if (! _process_smarts_bond(proto.bond_smarts(), completed))
    {
      cerr << "Substructure_Atom::construct_from_proto:invalid bond smarts bond\n";
      cerr << proto.ShortDebugString() << endl;
      return 0;
    }
  }

#ifdef QUERY_BONDS_PROCESSED_GREEDY
  if (proto.query_bond_size() > 0)
  {
    for (const auto & bond : proto.query_bond())
    {
      if (! _process_substructure_bond(bond, completed)) 
      {
        cerr << "Substructure_Atom::construct_from_proto:invalid substructure bond\n";
        cerr << proto.ShortDebugString() << endl;
        return 0;
      }
    }
  }
#endif

  if (!GETVALUES(proto, unmatched_atoms_attached, 0, no_limit))
    return 0;

  if (proto.has_atom_type_group())
    _atom_type_group = proto.atom_type_group();

  assert (ok());
  return 1;
}

// Bonds are discerned after all the atoms have been constructed.
int
Substructure_Atom::FormBonds(const SubstructureSearch::SubstructureAtom& proto,
                             extending_resizable_array<Substructure_Atom *> & completed)
{
  if (proto.query_bond_size() == 0) {
    return 1;
  }

#ifdef DEBUG_FORM_BONDS
  cerr << "Substructure_Atom::FormBonds:adding " << proto.query_bond_size() << " query bonds\n";
#endif
  for (const auto & bond : proto.query_bond())
  {
    if (! _process_substructure_bond(bond, completed)) 
    {
      cerr << "Substructure_Atom::construct_from_proto:invalid substructure bond\n";
      cerr << proto.ShortDebugString() << endl;
      return 0;
    }
  }

  return 1;
}


int
Substructure_Atom::_add_bonds(const google::protobuf::RepeatedField<uint32_t>& atoms,
    const bond_type_t btype,
    extending_resizable_array<Substructure_Atom *> & completed)
{
  for (const auto a : atoms)
  {
    if (nullptr == completed[a])
    {
      cerr << "Substructure_Atom::_add_bonds:non existent atom " << a << endl;
      return 0;
    }

    if (static_cast<int>(a) == _unique_id)
    {
      cerr << "Substructure_Atom::_add_bonds:self bonds not allowed, atom " << a << endl;
      return 0;
    }

    Substructure_Bond * b = new Substructure_Bond;
    b->set_atom(completed[a]);
    b->set_type(btype);
    _add_bond(b);
  }

  return 1;
}

int
Single_Substructure_Query::WriteProto(const char * fname)
{
  assert (nullptr != fname);

  std::ofstream os(fname, std::ios::out);
  if (! os.good())
  {
    cerr << "Single_Substructure_Query::write: cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteProto(os);
}

int
Single_Substructure_Query::WriteProto(std::ostream& output) 
{
  cerr << "Single_Substructure_Query::WriteProto:not implemented yet\n";
  return 0;
}

// for all the individual values in Matcher `p`, add them as
// repeated fields of the same to proto `p`.
#define ADDREPEATEDFIELDMATCHER(p, field) { \
  const resizable_array<int> values = _ ##field.ValuesMatched(); \
  for (auto v : values) { \
    proto.add_ ## field(v); \
  }  \
}

#define ADDREPEATEDFIELD(p, field) { \
  for (auto value : _ ## field) { \
    proto.add_ ## field(value); \
  } \
}

// Populates a Min_Max_Specifier.
#define SETPROTOVALUES(p, attribute, type) { \
  ADDREPEATEDFIELD(p, attribute) \
  type tmp; \
  if (_ ## attribute.min(tmp)) \
    p.set_min_ ## attribute(tmp); \
  if (_ ## attribute.max(tmp)) \
    p.set_max_ ## attribute(tmp); \
}

// Populates a Matcher.
#define SETPROTOVALUESMATCHER(p, attribute, type) { \
  ADDREPEATEDFIELDMATCHER(p, attribute) \
  type tmp; \
  if (_ ## attribute.min(tmp)) \
    p.set_min_ ## attribute(tmp); \
  if (_ ## attribute.max(tmp)) \
    p.set_max_ ## attribute(tmp); \
}

int
Single_Substructure_Query::BuildProto(SubstructureSearch::SingleSubstructureQuery& proto) const
{
  if (_comment.length())
    proto.set_comment(_comment.rawchars(), _comment.length());

  proto.set_one_embedding_per_start_atom(_find_one_embedding_per_start_atom);
  proto.set_normalise_rc_per_hits_needed(_normalise_rc_per_hits_needed);
  proto.set_subtract_from_rc(_subtract_from_rc);
  proto.set_max_matches_to_find(_max_matches_to_find);
  proto.set_save_matched_atoms(_save_matched_atoms);
  proto.set_ncon_ignore_singly_connected(_ncon_ignore_singly_connected);
  if (_do_not_perceive_symmetry_equivalent_matches) {
    proto.set_perceive_symmetric_equivalents(false);
  } else {
    proto.set_perceive_symmetric_equivalents(true);
  }
  if (_implicit_ring_condition >= 0) {
    proto.set_implicit_ring_condition(_implicit_ring_condition);
  }
  proto.set_all_hits_in_same_fragment(_all_hits_in_same_fragment);
  proto.set_only_match_largest_fragment(_only_keep_matches_in_largest_fragment);
  proto.set_embeddings_do_not_overlap(_embeddings_do_not_overlap);
  proto.set_sort_by_preference_value(_sort_by_preference_value);

  ADDREPEATEDFIELD(proto, numeric_value);

  // TODO: implement this
#ifdef IMPLEMENT_NO_MATCHED_ATOMS_BETWEEN
  for (const auto * nmb : _no_matched_atoms_between) {
    SubstructureSearch::NoMatchedAtomsBetween * p = proto.add_no_matched_atoms_between();
    nmb->BuildProtoNoBond(*p);
  }
#endif
  proto.set_no_matched_atoms_between_exhaustive(_no_matched_atoms_between_exhaustive);

  for (const auto * lnk : _link_atom) {
    SubstructureSearch::LinkAtoms* p = proto.add_link_atoms();
    lnk->BuildProto(*p);
  }
  proto.set_fail_if_embeddings_too_close(_fail_if_embeddings_too_close);

  if (_matched_atoms_to_check_for_hits_too_close > 0)
    proto.set_distance_between_hits_ncheck(_matched_atoms_to_check_for_hits_too_close);

  SetProtoValues(_attached_heteroatom_count, "attached_heteroatom_count", proto);
  SetProtoValues(_hits_needed, "hits_needed", proto);
  SetProtoValues(_ring_atoms_matched, "ring_atoms_matched", proto);
  SetProtoValues(_heteroatoms_matched, "heteroatoms_matched", proto);
  SetProtoValues(_heteroatoms_in_molecule, "heteroatoms_in_molecule", proto);
  SetProtoValues(_natoms, "natoms", proto);
  SetProtoValues(_nrings, "nrings", proto);
  SetProtoValues(_ncon, "ncon", proto);
  SetProtoValues(_fused_rings, "fused_rings", proto);
  SetProtoValues(_strongly_fused_rings, "strongly_fused_rings", proto);
  SetProtoValues(_isolated_rings, "isolated_rings", proto);
  SetProtoValues(_isolated_ring_objects, "isolated_ring_objects", proto);
  SetProtoValues(_aromatic_rings, "aromatic_rings", proto);
  SetProtoValues(_non_aromatic_rings, "non_aromatic_rings", proto);
  SetProtoValues(_distance_between_hits, "distance_between_hits", proto);
  SetProtoValues(_number_isotopic_atoms, "number_isotopic_atoms", proto);
  SetProtoValues(_number_fragments, "number_fragments", proto);
  SetProtoValues(_distance_between_root_atoms, "distance_between_root_atoms", proto);
  SetProtoValues(_atoms_in_spinach, "atoms_in_spinach", proto);
  SetProtoValues(_inter_ring_atoms, "inter_ring_atoms", proto);
  SetProtoValues(_unmatched_atoms, "unmatched_atoms", proto);
  SetProtoValues(_net_formal_charge, "net_formal_charge", proto);
  proto.set_environment_must_match_unmatched_atoms(_environment_must_match_unmatched_atoms);

  if (_min_fraction_atoms_matched > 0.0f) {
    proto.set_min_fraction_atoms_matched(_min_fraction_atoms_matched);
  }
  if (_max_fraction_atoms_matched > 0.0f) {
    proto.set_max_fraction_atoms_matched(_max_fraction_atoms_matched);
  }

  // TODO: Do something with _environments_can_share_attachment_points
  for (const auto* env : _environment) {
    env->BuildProto(*proto.add_environment());
  }
  for (const auto* env : _environment_rejections) {
    env->BuildProto(*proto.add_environment_no_match());
  }

  cerr << "Processing " << _ring_system_specification.size() << " ring system specifiers\n";
  for (const auto * ring_sys : _ring_system_specification) {
    SubstructureSearch::SubstructureRingSystemSpecification * p = proto.add_ring_system_specifier();
    ring_sys->BuildProto(*p);
  }

  for (const auto * ring : _ring_specification) {
    SubstructureSearch::SubstructureRingSpecification* p = proto.add_ring_specifier();
    ring->BuildProto(*p);
  }

  for (const auto * ele : _elements_needed) {
    ele->BuildProto(*proto.add_elements_needed());
  }
  SetProtoValues(_aromatic_atoms, "aromatic_atoms", proto);
  proto.set_unique_embeddings_only(_find_unique_embeddings_only);
  ADDREPEATEDFIELD(proto, heteroatoms);
  proto.set_respect_initial_atom_numbering(_respect_initial_atom_numbering);
  proto.set_compress_embeddings(_compress_embeddings);

  for (const auto * c: _chirality) {
    c->BuildProto(*proto.add_chiral_centre());
  }

  // atom type

  for (const auto* geom : _geometric_constraints) {
    geom->BuildProto(*proto.add_geometric_constraints());
  }

  for (const auto* sep : _separated_atoms) {
    sep->BuildProto(*proto.add_separated_atoms());
  }

  extending_resizable_array<Substructure_Atom*> atoms;
  for (Substructure_Atom * a : _root_atoms) {
    a->collect_all_atoms(atoms);
  }

  for (const Substructure_Atom * atom : atoms) {
    if (atom != nullptr) {
      atom->BuildProto(*proto.add_query_atom());
    }
  }

  return 1;
}

// Many operations involve checking if a value is > `minval` and if so
// then set the proto attribute.
#define SET_PROTO_IF_SET(p, attribute, minval) { \
  if (_ ## attribute > minval) { \
    p.set_ ## attribute(_ ## attribute); \
  } \
}

int
Substructure_Atom::BuildProto(SubstructureSearch::SubstructureAtom& proto) const {
  proto.set_id(_unique_id);  // 1

  if (! _match_as_match_or_rejection) {  // Default is true.
    proto.set_match_as_match(_match_as_match_or_rejection);  // 2
  }

  if (_text_identifier.length() > 0) {
    proto.set_text_identifier(_text_identifier.data(), _text_identifier.length());  // 3
  }

  SetProtoIfSet(_atom_map_number, -1, "atom_map_number", proto);  // 4
  SetProtoIfSet(_initial_atom_number, -1, "initial_atom_number", proto);  // 4
#ifdef OR_ID_NO_LONGER_PROCESSED
  SET_PROTO_IF_SET(proto, or_id, 0);  // 6
#endif
  SET_PROTO_IF_SET(proto, ring_id, 0);  // 9
  SET_PROTO_IF_SET(proto, fused_system_id, 0);  // 10
  SET_PROTO_IF_SET(proto, fragment_id, 0);  // 11
  SetProtoIfSet(_global_match_id, 0, "global_match_id", proto);  // 34

  double nv;
  if (_numeric_value.value(nv)) {  // 12
    proto.set_numeric_value(nv);
  }
  if (! _include_in_embedding) {
    proto.set_include_in_embedding(false);  // 13
  }

  //  smarts // 14
  //  atom_smarts // 15
  //  smiles // 16

  // environment looks hard.  // 17

  for (const Substructure_Bond* b : _bonds) {  // 21
    b->BuildProto(*proto.add_query_bond());
  }

  // bond smarts 22

  // single, double triple aromatic, bond 25-29

  for (const Substructure_Atom_Specifier* a : _preferences) {  // 23
    a->BuildProto(*proto.add_preference());
  }

  if (_sum_all_preference_hits) {
    proto.set_sum_all_preference_hits(true);  // 24
  }

  if (_unmatched_atoms_attached.is_set()) {
    SETPROTOVALUESMATCHER(proto, unmatched_atoms_attached, int);
  }

//SET_PROTO_IF_SET(proto, atom_type_group, 0);
  SetProtoIfSet(_atom_map_number, 0, "atom_type_group", proto);

  if (_components.number_elements() == 0) {
    const Substructure_Atom_Specifier* me = this;
    me->BuildProto(*proto.add_atom_properties());
  } else {
    for (const Substructure_Atom_Specifier* c : _components) {
      c->BuildProto(*proto.add_atom_properties());
    }
  }

  if (_components.number_elements() > 1) {
    for (int i = 0; i < _components.number_elements(); ++i) {
      SubstructureSearch::SubstructureAtomSpecifier* atom_prop = proto.add_atom_properties();
      _components[i]->BuildProto(*atom_prop);
      if (i == 0) {
        continue;
      }
      const int op = _operator.op(i - 1);
      if (op == IW_LOGEXP_AND) {
        atom_prop->set_logical_operator(SubstructureSearch::SS_AND);
      } else if (op == IW_LOGEXP_OR) {
        atom_prop->set_logical_operator(SubstructureSearch::SS_OR);
      } else if (op == IW_LOGEXP_XOR) {
        atom_prop->set_logical_operator(SubstructureSearch::SS_XOR);
      } else if (op == IW_LOGEXP_LOW_PRIORITY_AND) {
        atom_prop->set_logical_operator(SubstructureSearch::SS_LP_AND);
      } else {
        cerr << "Substructure_Atom::BuildProto:unrecognized operator type " << op << '\n';
        return 0;
      }
    }
  }

  return 1;
}

int
Substructure_Bond::BuildProto(SubstructureSearch::SubstructureBond& proto) const {
  if (_bond_types & SINGLE_BOND) {
    proto.add_bond_type(SubstructureSearch::SS_SINGLE_BOND);
  }
  if (_bond_types & DOUBLE_BOND) {
    proto.add_bond_type(SubstructureSearch::SS_DOUBLE_BOND);
  }
  if (_bond_types & TRIPLE_BOND) {
    proto.add_bond_type(SubstructureSearch::SS_TRIPLE_BOND);
  }
  if (_bond_types & AROMATIC_BOND) {
    proto.add_bond_type(SubstructureSearch::SS_AROMATIC_BOND);
  }

  proto.set_other_end(_a1->unique_id());

  if (_b == nullptr) {
    return 1;
  }

  cerr << "Substructure_Bond::BuildProto:unhandled bond conditions\n";
  return 0;
}

int
Substructure_Atom_Specifier::BuildProto(SubstructureSearch::SubstructureAtomSpecifier& proto) const {
  for (const Element* e : _element) {
    if (e->is_in_periodic_table()) {
      proto.add_atomic_number(e->atomic_number());
    } else {
      proto.add_atomic_symbol(e->symbol().data(), e->symbol().length());
    }
  }

  SetProtoValues(_ncon, "ncon", proto);  // 3
  SetProtoValues(_ncon2, "ncon2", proto);  // 6
  SetProtoValues(_nbonds, "nbonds", proto);  // 9
  SetProtoValues(_formal_charge, "formal_charge", proto);  // 12
  SetProtoValues(_nrings, "nrings", proto);  // 15
  SetProtoValues(_ring_bond_count, "ring_bond_count", proto);  // 18
  SetProtoValues(_ring_size, "ring_size", proto);  // 21
  SetProtoValues(_hcount, "hcount", proto);  // 24
  if (_aromaticity == SUBSTRUCTURE_NOT_SPECIFIED) {
  }  else if (_aromaticity == AROMATIC) {  // 27
    proto.set_aromatic(true);
  }  else if (_aromaticity == NOT_AROMATIC) {
    proto.set_aromatic(false);
  }
  if (_chirality != SUBSTRUCTURE_NOT_SPECIFIED)
    proto.set_chirality(true);  // 28
  SetProtoValues(_aromatic_ring_size, "aromatic_ring_size", proto);  // 30
  SetProtoValues(_aliphatic_ring_size, "aliphatic_ring_size", proto);  // 33
  SetProtoValues(_attached_heteroatom_count, "attached_heteroatom_count", proto);  // 36
  SetProtoValues(_lone_pair_count, "lone_pair_count", proto);  // 39
  SetProtoValues(_unsaturation, "unsaturation", proto);  // 42
  SetProtoValues(_daylight_x, "daylight_x", proto);  // 45
  SetProtoValues(_isotope, "isotope", proto);  // 48
  SetProtoValues(_aryl, "aryl", proto);  // 51
  SetProtoValues(_vinyl, "vinyl", proto);  // 54
  SetProtoValues(_fused_system_size, "fused_system_size", proto);  // 57
  if (_all_rings_kekule != SUBSTRUCTURE_NOT_SPECIFIED) {
    proto.set_all_rings_kekule(true);  // 60
  }
  SetProtoValues(_heteroatoms_in_ring, "heteroatoms_in_ring", proto);  // 61
  SET_PROTO_IF_SET(proto, match_spinach_only, -1);  // 64
  SetProtoValues(_scaffold_bonds_attached_to_ring, "scaffold_bonds_attached_to_ring", proto);  // 65
  SET_PROTO_IF_SET(proto, preference_value, 0);  // 68
  SetProtoValues(_symmetry_degree, "symmetry_degree", proto);  // 69
  SET_PROTO_IF_SET(proto, symmetry_group, 0);  // 72

  // Not sure what to do with operator...
  // atom typing is not implemented, seems ambiguous.

  return 1;
}

int
Single_Substructure_Query::_parse_ring_specifier(const SubstructureSearch::SingleSubstructureQuery& proto)
{
  for (const auto& ring_spec : proto.ring_specifier())
  {
    Substructure_Ring_Specification * r = new Substructure_Ring_Specification();
    if (! r->ConstructFromProto(ring_spec))
    {
      delete r;

      cerr << "SingleSubstructureQuery::_parse_ring_specifier:could not create ring specifier object from msi object\n";
      cerr << ring_spec.ShortDebugString();
      return 0;
    }

    _ring_specification.add(r);
  }

  const int ring_specification_count = _ring_specification.number_elements();

  if (1 == ring_specification_count)    // no operators with first object
    return 1;

  if (0 == proto.ring_specification_logexp_size())
  {
    for (int i = 0; i < ring_specification_count - 1; ++i)
    {
      _ring_specification_logexp.add_operator(IW_LOGEXP_AND);
    }

    return 1;
  }

  if (! ExtractOperator(proto.ring_specification_logexp(), ring_specification_count - 1,
          _ring_specification_logexp, IW_LOGEXP_AND,
          "Single_Substructure_Query::_parse_ring_specifier_object"))
    return 0;

  return 1;
}

int
Single_Substructure_Query::_parse_ring_system_specifier(const SubstructureSearch::SingleSubstructureQuery& proto)
{
  if (0 == proto.ring_system_specifier_size())
    return 1;

  for (const auto& spec : proto.ring_system_specifier())
  {
    Substructure_Ring_System_Specification * r = new Substructure_Ring_System_Specification();
    if (! r->ConstructFromProto(spec))
    {
      delete r;

      cerr << "SingleSubstructureQuery::_parse_ring_system_specifier:could not create ring system specifier object from proto\n";
      cerr << spec.ShortDebugString() << "\n";
      return 0;
    }

    _ring_system_specification.add(r);
  }

  const int ring_system_spec_count = _ring_system_specification.number_elements();

  if (1 == ring_system_spec_count)  // No operators to worry about.
    return 1;

  if (proto.ring_system_specifier_logexp().empty())
//if (0 == proto.ring_system_specification_logexp_size())
  {
    for (int i = 0; i < ring_system_spec_count - 1; ++i) {
      _ring_system_specification_logexp.add_operator(IW_LOGEXP_AND);
    }

    return 1;
  }

  if (! ExtractOperator(proto.ring_system_specifier_logexp(), ring_system_spec_count - 1,
        _ring_system_specification_logexp, IW_LOGEXP_AND, "Single_Substructure_Query::_parse_ring_system_specifier_object"))
    return 0;

  return 1;
}

/*
  Is an msi object a root atom or not.  Basically, if it has an
  attribute which ends in "bond" (so as to not also match "nbonds"),
  then it is not a root
*/

static bool
IsRootSubstructureAtom(const SubstructureSearch::SubstructureAtom& proto)
{
  if (proto.query_bond_size() > 0)
    return false;

  if (proto.has_bond_smarts())
    return false;

  if (proto.single_bond_size() > 0)
    return false;

  if (proto.double_bond_size() > 0)
    return false;

  if (proto.triple_bond_size() > 0)
    return false;

  if (proto.aromatic_bond_size() > 0)
    return false;

  if (proto.bond_size() > 0)
    return false;

  return true;
}

/*
  The environment object 
  It can specify query_atoms, smiles or smarts
*/

int
Single_Substructure_Query::_construct_environment_from_proto(
    const google::protobuf::RepeatedPtrField<SubstructureSearch::SubstructureEnvironment>& env,
    extending_resizable_array<Substructure_Atom *> & completed,
    resizable_array_p<Substructure_Environment> & destination)
{
//cerr << "Constructing env from " << env.size() << " components\n";

  for (const auto& e : env)
  {
    Substructure_Environment * a = new Substructure_Environment();

    _collect_all_atoms(completed);

    int rc = a->construct_from_proto(e, completed);

    if (0 == rc)
    {
      cerr << "Single_Substructure_Query::_construct_environment_from_proto:invalid environment " << e.ShortDebugString() << "\n";
      delete a;
      return 0;
    }

//  cerr << "Just processed env " << e.ShortDebugString() << endl;
    destination.add(a);
  }

  return 1;
}

int
Single_Substructure_Query::_construct_matched_atoms_match(
    const google::protobuf::RepeatedPtrField<SubstructureSearch::MatchedAtomMatch>& matched_atoms_match) {
  for (const auto & mam : matched_atoms_match) {

    std::unique_ptr<MatchedAtomMatch> m = std::make_unique<MatchedAtomMatch>();
    if (! m->ConstructFromProto(mam)) {
      cerr << "Single_Substructure_Query::_construct_matched_atoms_match:invalid " << mam.ShortDebugString() << '\n';
      return 0;
    }
    if (mam.has_logexp()) {
      AddOperator(mam.logexp(), _matched_atom_match_operator);
    } else if (_matched_atom_match.size() > 0) {
      AddOperator(SubstructureSearch::SS_AND, _matched_atom_match_operator);
    }
    _matched_atom_match << m.release();
  }

  return 1;
}

#ifdef PROBABLY_DUP
int
Substructure_Environment::_add_possible_parent (atom_number_t possible_parent,
                                                bond_type_t possible_parent_bond_type,
                                                extending_resizable_array<Substructure_Atom *> & completed)
{
  if (0 != possible_parent_bond_type)     // only add if bond_type has defined a particular type
    _bond.set_type(possible_parent_bond_type);
  else
    _bond.set_match_any();

  _possible_parents.add(completed[possible_parent]);

  return 1;
}
#endif

int
Substructure_Environment::BuildProto(SubstructureSearch::SubstructureEnvironment& proto) const
{
  proto.set_id(_unique_id);

  SET_PROTO_IF_SET(proto, or_id, 0);  // 8
  SET_PROTO_IF_SET(proto, and_id, 0);  // 9

  SETPROTOVALUESMATCHER(proto, hits_needed, int);  // 10
  if (_no_other_substituents_allowed) {
    proto.set_no_other_substituents_allowed(true);  // 13
  }
  if (_environments_can_share_attachment_points) {
    proto.set_env_matches_can_share_attachment_points(true);  // 15
  }
  SET_PROTO_IF_SET(proto, max_matches_to_find, 0);  // 16
  if (_hydrogen_ok_as_environment_match) {
    proto.set_hydrogen_ok(true);  // 17
  }
  if (_max_environment_matches_per_attachment_point > 0) {
    proto.set_max_env_matches_per_anchor(_max_environment_matches_per_attachment_point);
  }

  // The rest of this is too hard.
  for (const Substructure_Atom* a : _possible_parents) {
    a->BuildProto(*proto.add_query_atom());
  }
  cerr << "Substructure_Environment::BuildProto:implement this sometime\n";

  return 1;
}

int
Substructure_Environment::construct_from_proto(const SubstructureSearch::SubstructureEnvironment& proto,
                               extending_resizable_array<Substructure_Atom *> & completed,
                               atom_number_t possible_parent,
                               bond_type_t possible_parent_bond_type)
{
// Process any bonds which have been specified as attributes

  if (! _process_how_connected(proto, completed))
  {
    cerr << "Substructure_Environment::construct_from_proto:cannot connect " << proto.ShortDebugString() << "\n";
    return 0;
  }

  if (INVALID_ATOM_NUMBER != possible_parent)
    _add_possible_parent(possible_parent, possible_parent_bond_type, completed);

  const int np = _possible_parents.number_elements();

  if (0 == np)
  {
    cerr << "Substructure_Environment::construct_from_proto: environment not connected\n";
    return 0;
  }

  if (proto.hits_needed_size() > 0) {
    for (const auto x : proto.hits_needed()) {
        _hits_needed.add_if_not_already_present(x);
    }
  }
  if (proto.has_min_hits_needed()) {
    _hits_needed.set_min(proto.min_hits_needed());
  }
  if (proto.has_max_hits_needed()) {
    _hits_needed.set_max(proto.max_hits_needed());
  }

  if (!_hits_needed.ok()) {
    cerr << "Substructure_Environment::construct_from_proto:invalid hits needed\n";
    _hits_needed.debug_print(cerr);
  }

  if (proto.has_no_other_substituents_allowed())
    _no_other_substituents_allowed = proto.no_other_substituents_allowed();

  if (proto.has_hydrogen_ok())
    _hydrogen_ok_as_environment_match = proto.hydrogen_ok();

  if (proto.has_max_env_matches_per_anchor())
    _max_environment_matches_per_attachment_point = proto.max_env_matches_per_anchor();

  if (proto.has_env_matches_can_share_attachment_points())
    _environments_can_share_attachment_points = proto.env_matches_can_share_attachment_points();

  if (proto.has_max_matches_to_find())
    _max_matches_to_find = proto.max_matches_to_find();

  if (proto.has_or_id())
  {
    _or_id = proto.or_id();

    if (_or_id <= 0) 
    {
      cerr << "Substructure_Environment::construct_from_proto: or_id values must be whole +ve numbers " << _or_id << " invalid\n";
      return 0;
    }
  }

  if (proto.has_and_id()) 
  {
    _and_id = proto.and_id();
    if (_and_id <= 0)
    {
      cerr << "Substructure_Environment::construct_from_proto: and_id values must be whole +ve numbers " << _and_id << " invalid\n";
      return 0;
    }
  }

  if (_and_id && _or_id)
  {
    cerr << "Substructure_Environment::construct_from_proto:cannot have both OR and AND specifications\n";
    return 0;
  }

// Now construct the structural specification. There can be any number

  for (const std::string& smarts: proto.smarts()) {
    Substructure_Atom * a = new Substructure_Atom;
//  a->set_unique_id(msi.object_id());

    const_IWSubstring x = smarts;
    const char * s = smarts.data();

    if (x.length() && (isdigit(s[0]) || '>' == s[0] || '<' == s[0] || s[0] == open_brace)) {
      int chars_consumed = substructure_spec::SmartsNumericQualifier(s, smarts.length(), _hits_needed);
      if (chars_consumed == 0) {
        cerr << "Substructure_Environment::construct_from_proto:invalid numeric qualifier '" << x << "'\n";
        delete a;
        return 0;
      }
#ifdef OLD_VERSION
      int value, qualifier;
      int chars_consumed = substructure_spec::SmartsFetchNumeric(s, value, qualifier);
      if (0 == chars_consumed)
      {
        cerr << "Substructure_Environment::construct_from_proto:invalid numeric qualifier '" << x << "'\n";
        delete a;
        return 0;
      }
      if (qualifier > 0)
        _hits_needed.set_min(value+1);
      else if (qualifier < 0)
        _hits_needed.set_max(value-1);
      else
        _hits_needed.add(value);
#endif

      x += chars_consumed;
      if (! a->parse_smarts_specifier(x)) {
        cerr << "Substructure_Environment::construct_from_proto:invalid smarts '" << x << "'\n";
        return 0;
      }
    }

    else if (! a->parse_smarts_specifier(x))
    {
      delete a;
      return 0;
    }

//  cerr << "Build query from smarts " << (*att) << endl;
    add(a);
  }

  for (const auto& smiles : proto.smiles()) {
    Substructure_Atom * a = new Substructure_Atom;
//  a->set_unique_id(msi.object_id());

    if (! a->parse_smiles_specifier(smiles))
    {
      delete a;
      return 0;
    }

//  cerr << "Build query from smiles " << (*att) << endl;
    add(a);
  }

  for (const auto& query_atom: proto.query_atom()) {
    std::unique_ptr<Substructure_Atom> a = std::make_unique<Substructure_Atom>();
    if (! a->construct_from_proto(query_atom, completed))
      return 0;

    if (IsRootSubstructureAtom(query_atom))
      cerr << "Root atom " << query_atom.ShortDebugString() << endl;
    if (IsRootSubstructureAtom(query_atom))
      add(a.release());
    else
    {
      cerr << "Substructure_Environment::construct_from_proto:non root atom\n";
      cerr << query_atom.ShortDebugString() << endl;
      return 0;
    }
  }

  return 1;
}


template <typename T>
void
transfer_to_our_array(resizable_array_p<T> & to,
                      const resizable_array<T *> & from)
{
  for (T* x : from)
  {
    to.add(x);
  }

  return;
}

template void transfer_to_our_array (resizable_array_p<Substructure_Atom> &, const resizable_array<Substructure_Atom *> &);
template void transfer_to_our_array (resizable_array_p<Bond> &, const resizable_array<Bond *> &);
template void transfer_to_our_array (resizable_array_p<Link_Atom> &, const resizable_array<Link_Atom *> &);

SeparatedAtoms::SeparatedAtoms() {
  _a1 = INVALID_ATOM_NUMBER;
  _a2 = INVALID_ATOM_NUMBER;
}

int
SeparatedAtoms::Build(const SubstructureSearch::SeparatedAtoms& proto) {
  _a1 = proto.a1();
  _a2 = proto.a2();
  if (_a1 == _a2) {
    cerr << "SeparatedAtoms::Build:separated atoms must be distinct\n";
    return 0;
  }
  if (proto.bonds_between().size() > 0) {
    for (auto s : proto.bonds_between()) {
      _separation.add_if_not_already_present(s);
    }
  }
  if (proto.min_bonds_between() > 0) {
    _separation.set_min(proto.min_bonds_between());
  }
  if (proto.max_bonds_between() > 0) {
    if (! _separation.set_max(proto.max_bonds_between())) {
      cerr << "SeparatedAtoms::Build:invalid max bonds between\n";
      return 0;
    }
  }

  return 1;
}

//#define DEBUG_SSSQ_PARSE_SMARTS_SPECIFIER

/*
  Parsing a smarts is complicated by the possible presence of the '.' character,
  and the need to process the '...' construct (no_matched_atoms_between)

  Valid inputs are

  'C-O-C'
  'C-O-C.[NH2]'
  'Br...C(=O)Cl'

*/




template <typename FROM, typename TO>
int
FetchRepeatedField(const google::protobuf::RepeatedField<FROM> & values,
                   resizable_array<TO> & destination)
{
  for (const FROM z : values)
  {
    destination.add(static_cast<TO>(z));
  }

  return destination.number_elements();
}

int
Single_Substructure_Query::_construct_from_proto(const SubstructureSearch::SingleSubstructureQuery & proto)
{
  assert (ok());
  assert (_root_atoms.empty());

  _numeric_value.resize(0);

  int environments_can_share_attachment_points = -1;    // indicates not set

  if (proto.has_comment())
    _comment = proto.comment();
  else if (proto.has_label())
    _comment = proto.label();

  if (proto.has_one_embedding_per_start_atom())
    set_find_one_embedding_per_atom(proto.one_embedding_per_start_atom());

  if (proto.has_unique_embeddings_only())
    set_find_unique_embeddings_only(proto.unique_embeddings_only());

  if (proto.has_normalise_rc_per_hits_needed())
    set_normalise_rc_per_hits_needed(proto.normalise_rc_per_hits_needed());

  if (proto.has_subtract_from_rc())
    _subtract_from_rc = proto.subtract_from_rc();

  if (proto.has_respect_initial_atom_numbering())
    _respect_initial_atom_numbering = proto.respect_initial_atom_numbering();

  if (proto.has_max_matches_to_find())
    _max_matches_to_find = proto.max_matches_to_find();

  if (proto.has_save_matched_atoms())
    _save_matched_atoms = proto.save_matched_atoms();

  if (proto.has_ncon_ignore_singly_connected())
    _ncon_ignore_singly_connected = proto.ncon_ignore_singly_connected();

  if (proto.has_perceive_symmetric_equivalents())
    _do_not_perceive_symmetry_equivalent_matches = !proto.perceive_symmetric_equivalents();

  if (proto.has_embeddings_do_not_overlap())
    _embeddings_do_not_overlap = proto.embeddings_do_not_overlap();

  if (proto.has_implicit_ring_condition())
    _implicit_ring_condition = proto.implicit_ring_condition();

  if (proto.has_all_hits_in_same_fragment())
    _all_hits_in_same_fragment = proto.all_hits_in_same_fragment();

  if (proto.has_only_match_largest_fragment())
    _only_keep_matches_in_largest_fragment = proto.only_match_largest_fragment();

  if (proto.has_one_embedding_per_start_atom())
    _find_one_embedding_per_start_atom = proto.one_embedding_per_start_atom();

  if (proto.has_sort_by_preference_value())
    _sort_by_preference_value = proto.sort_by_preference_value();

  if (proto.numeric_value_size() > 0)
    FetchRepeatedField( proto.numeric_value(), _numeric_value);

  if (proto.has_fail_if_embeddings_too_close())
    _fail_if_embeddings_too_close = proto.fail_if_embeddings_too_close();

  if (proto.has_distance_between_hits_ncheck())
  {
    _matched_atoms_to_check_for_hits_too_close = proto.distance_between_hits_ncheck();
    if (_matched_atoms_to_check_for_hits_too_close == 0)
    {
      cerr << "Single_Substructure_Query::_construct_from_proto:invalid distance_between_hits_ncheck\n";
      return 0;
    }
  }

  if (proto.has_environment_must_match_unmatched_atoms())
    _environment_must_match_unmatched_atoms = proto.environment_must_match_unmatched_atoms();

  if (proto.has_environments_can_share_attachment_points())
    environments_can_share_attachment_points = proto.environments_can_share_attachment_points();

  if (proto.heteroatoms_size() > 0)
  {
    if (!FetchRepeatedField(proto.heteroatoms(), _heteroatoms))
    {
      cerr << "SingleSubstructureQuery::_construct_from_proto:invalid heteroatoms specification\n";
      cerr << proto.ShortDebugString() << endl;
      return 0;
    }
  }

  // The smiles and smarts directives scramble the _unique_id.
  if (_respect_initial_atom_numbering &&
      (proto.has_smiles() || proto.has_smarts()))
  {
    cerr << "SingleSubstructureQuery::construct_from_proto:cannot use respect_initial_atom_numbering with smiles or smarts directives\n";
    return 0;
  }

  if (proto.has_smiles())
  {
    const IWString smiles = proto.smiles();
    Substructure_Atom * a = new Substructure_Atom;
    if (! a->parse_smiles_specifier(smiles))
    {
      delete a;
      cerr << "Single_Substructure_Query::_construct_from_proto: cannot parse '" << smiles << "'\n";
      return 0;
    }

    _root_atoms.add(a);
  }

  if (proto.has_smarts())
  {
    const_IWSubstring smarts = proto.smarts();

    if (! create_from_smarts(smarts))
    {
      cerr << "Single_Substructure_Query::_construct_from_proto: cannot parse '" << smarts << "'\n";
      return 0;
    }
  }

  if (proto.no_matched_atoms_between_size() > 0)
  {
    for (const auto& nma : proto.no_matched_atoms_between())
    {
      if (nma.a1() == nma.a2())
      {
        cerr << "Single_Substructure_Query::_construct_from_proto:invalid no matched atoms between " << nma.ShortDebugString() << "'\n";
        return 0;
      }

      Bond * b = new Bond(nma.a1(), nma.a2(), SINGLE_BOND);

      _no_matched_atoms_between.add(b);
    }
  }

  if (proto.link_atoms_size() > 0)
  {
    for (const auto& link : proto.link_atoms())
    {
      Link_Atom * l = new Link_Atom;

      if (! l->ConstructFromProto(link))
      {
        cerr << "Single_Substructure_Query::_construct_from_proto:invalid link atom specification\n";
        cerr << link.ShortDebugString() << endl;
        delete l;
        return 0;
      }

      _link_atom.add(l);
    }
  }

  if (proto.has_sort_matches() > 0)
  {
    _sort_matches_by = 0;

    const IWString s(proto.sort_matches());

    int j = 0;
    const_IWSubstring token;
    while (s.nextword(token, j))
    {
      int direction = 1;
      if (token.ends_with('+'))
        token.chop();
      else if (token.ends_with('-'))
      {
        token.chop();
        direction = -1;
      }
      if ("ncon" == token)
      {
        if (1 == direction)
          _sort_matches_by |= SORT_MATCHES_BY_INCREASING_NCON;
        else
          _sort_matches_by |= SORT_MATCHES_BY_DECREASING_NCON;
      }
      else if ("maxd" == token)
      {
        if (1 == direction)
          _sort_matches_by |= SORT_MATCHES_BY_INCREASING_MAXD;
        else
          _sort_matches_by |= SORT_MATCHES_BY_DECREASING_MAXD;
      }
      else
      {
        cerr << "Single_Substructure_Query:construct_from_proto:unrecognised sort attribute '" << s << "'\n";
        return 0;
      }
    }
  }

  if (_only_keep_matches_in_largest_fragment && _all_hits_in_same_fragment)
  {
    cerr << "Single_Substructure_Query::_construct_from_proto:the _only_keep_matches_in_largest_fragment and _all_hits_in_same_fragment attributes are mutually inconsistent\n";
    return 0;
  }

  constexpr uint32_t no_limit = std::numeric_limits<uint32_t>::max();

  if (!GETVALUES(proto, attached_heteroatom_count, 0, no_limit))
    return 0;

  if (_hits_needed.is_set())   // from a numeric qualifier on a smarts
    ;
  else if (! GETVALUES(proto, hits_needed, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, ring_atoms_matched, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, heteroatoms_matched, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, heteroatoms_in_molecule, 0, no_limit))
    return 0;

  if (_heteroatoms.number_elements() && 
      (! _attached_heteroatom_count.is_set() && ! _heteroatoms_in_molecule.is_set()))
  {
    cerr << "Single_Substructure_Query::_construct_from_proto: " << _heteroatoms.number_elements() <<
            " heteroatoms defined, but no attached_heteroatom_count\n";
    return 0;
  }

  if (! GETVALUES(proto, natoms, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, nrings, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, ncon, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, fused_rings, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, strongly_fused_rings, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, isolated_rings, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, isolated_ring_objects, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, aromatic_atoms, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, aromatic_rings, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, non_aromatic_rings, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, distance_between_hits, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, number_isotopic_atoms, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, number_fragments, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, distance_between_root_atoms, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, atoms_in_spinach, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, inter_ring_atoms, 0, no_limit))
    return 0;

  if (! GETVALUES(proto, unmatched_atoms, 0, no_limit))
    return 0;

  // Need to sync up with reasonable formal charge
  if (! GETVALUESINT(proto, net_formal_charge, -12, 12))
    return 0;

  if (proto.has_min_fraction_atoms_matched())
    _min_fraction_atoms_matched = proto.min_fraction_atoms_matched();

  if (proto.has_max_fraction_atoms_matched())
    _max_fraction_atoms_matched = proto.max_fraction_atoms_matched();

  if (proto.has_min_fraction_atoms_matched() && !proto.has_max_fraction_atoms_matched())
    _max_fraction_atoms_matched = 1.0f;
  else if (_min_fraction_atoms_matched > _max_fraction_atoms_matched)
  {
    cerr << "Single_Substructure_Query::_construct_from_proto:invalid min/max fraction atoms matched\n";
    cerr << proto.ShortDebugString() << "\n";
    return 0;
  }

  if (proto.ring_specifier_size() > 0)
  {
    if ( ! _parse_ring_specifier(proto))
      return 0;
  }

  if (proto.ring_system_specifier_size() > 0)
  {
    if (! _parse_ring_system_specifier(proto))
      return 0;
  }

  if (proto.element_hits_needed_size() > 0)
  {
    if (! _parse_element_hits_needed(proto))
      return 0;
  }

  if (proto.elements_needed_size() > 0)
  {
    if (! _parse_elements_needed(proto))
      return 0;
  }

  if (proto.has_compress_embeddings())
    _compress_embeddings = proto.compress_embeddings();

  extending_resizable_array<Substructure_Atom *> completed;

  // Form root atoms first.
  for (const auto& atom : proto.query_atom())
  {
    if (IsRootSubstructureAtom(atom))
    {
      Substructure_Atom * r = new Substructure_Atom;

      if (! r->construct_from_proto(atom, completed))
      {
        delete r;
        return 0;
      }

      _root_atoms.add(r);
    }
  }

  // And now non-root atoms.
  for (const auto & atom : proto.query_atom())
  {
    if (IsRootSubstructureAtom(atom))
      continue;

    Substructure_Atom* a = new Substructure_Atom();
    if (! a->construct_from_proto(atom, completed))
    {
      cerr << "Single_Substructure_Query::_construct_from_proto:invalid query atom " << atom.ShortDebugString() << "\n";
      delete a;
      return 0;
    }
  }

  // Now that all atoms are available in the completed array, form bonds.
  // Note that this only handles query_bond specifications, not single_bond, double_bond...
  for (const auto& atom : proto.query_atom()) {
    Substructure_Atom * a = completed[atom.id()];
    if (! a->FormBonds(atom, completed)) {
      cerr << "Single_Substructure_Query::_construct_from_proto:FormBonds failed\n";
      return 0;
    }
    if (!IsRootSubstructureAtom(atom) && a->nbonds() == 0) {
      cerr << "Non root Substructure_Atom not bonded, id " << a->unique_id() << "\n";
      cerr << atom.ShortDebugString() << "\n";
      return 0;
    }
  }

// OK if the root is not built, but let's make sure there are some global
// attributes available in that case

  if (_root_atoms.empty())
  {
//  cerr << "Single_Substructure_Query::_construct_from_proto: Root atom not found\n";

//  apr 99. There are some global conditions which don't make sense with no root atom

//  TODO add check for atom types here.
    if (_element_hits_needed.number_elements() || _distance_between_root_atoms.is_set() || 
        _no_matched_atoms_between.number_elements() || _attached_heteroatom_count.is_set() ||
        _ncon.is_set() || _heteroatoms_matched.is_set() || _ring_atoms_matched.is_set())
    {
      cerr << "Single_Substructure_Query::_construct_from_proto: global conditions incompatible with no root atom\n";
      return 0;
    }

//  Should list all the possible match conditions

    const int nat = _compute_attribute_counts();

    if (nat > 0) {
      // Great, got things to match.
    } else if (_natoms.is_set() || _nrings.is_set() || _aromatic_rings.is_set() ||
             _atoms_in_spinach.is_set() || _inter_ring_atoms.is_set() ||
             _net_formal_charge.is_set() || 
             _ring_specification.size() > 0 ||
             _ring_system_specification.size() > 0) {
    } else {
      cerr << "No root atoms and no global attributes specified, rejected\n";
      return 0;
    }
  }

  if (! _construct_environment_from_proto(proto.environment(), completed, _environment))
  {
    cerr << "Single_Substructure_Query::_construct_from_proto: environment interpretation failed\n";
    cerr << proto.ShortDebugString() << endl;
    return 0;
  }

  if (! _construct_environment_from_proto(proto.environment_no_match(), completed, _environment_rejections))
  {
    cerr << "Single_Substructure_Query::_construct_from_proto: environment interpretation failed\n";
    cerr << proto.ShortDebugString() << endl;
    return 0;
  }

  if (environments_can_share_attachment_points >= 0)
  {
    for (int i = 0; i < _environment.number_elements(); ++i)
    {
      _environment[i]->set_environments_can_share_attachment_points(environments_can_share_attachment_points);
    }
    for (int i = 0; i < _environment_rejections.number_elements(); ++i)
    {
      _environment_rejections[i]->set_environments_can_share_attachment_points(environments_can_share_attachment_points);
    }
  }

  if (proto.has_atom_type())
  {
    const auto t1 = proto.atom_type();
    const_IWSubstring t2(t1.data(), t1.size());
    _atom_typing = new Atom_Typing_Specification();
    if (! _atom_typing->build(t2))
    {
      cerr << "Single_Substructure_Query::_construct_from_proto:invalid atom type specification " << t2 << "'\n";
      return 0;
    }
  }

  if (proto.no_matched_atoms_between().size() > 0) {
  }

  if (proto.geometric_constraints().size() > 0) {
    for (const auto& constraint : proto.geometric_constraints()) {
      std::unique_ptr<geometric_constraints::SetOfGeometricConstraints> c =
              std::make_unique<geometric_constraints::SetOfGeometricConstraints>();
      if (! c->BuildFromProto(constraint)) {
        cerr << "Single_Substructure_Query::_construct_from_proto:invalid distance constraint " << constraint.ShortDebugString() << '\n';
        return 0;
      }
      if (! _atom_numbers_in_geometric_constraints_ok(c.get())) {
        cerr << "Single_Substructure_Query::_construct_from_proto:invalid atom numbers in geometric constraint " << constraint.ShortDebugString() << '\n';
        return 0;
      }
      _geometric_constraints.add(c.release());
    }
  }

  if (proto.separated_atoms().size() > 0) {
    for (const auto& separated_atoms : proto.separated_atoms()) {
      std::unique_ptr<SeparatedAtoms> s(new SeparatedAtoms);
      if (! s->Build(separated_atoms)) {
        cerr << "Single_Substructure_Query::_construct_from_proto:invalid separated atoms proto " << separated_atoms.ShortDebugString() << '\n';
        return 0;
      }
      _separated_atoms.add(s.release());
    }
  }

  if (!proto.matched_atom_must_be().empty()) {
    if (! _construct_matched_atoms_match(proto.matched_atom_must_be())) {
      cerr << "Single_Substructure_Query::_construct_from_proto:invalid matched_atom_must_be\n";
      return 0;
    }
  }

//if (_root_atoms.number_elements())
//  _adjust_for_internal_consistency();

  return 1;
}

int
Single_Substructure_Query::_atom_numbers_in_geometric_constraints_ok(const geometric_constraints::SetOfGeometricConstraints* constraints) const {
  const resizable_array<int> atom_numbers_present = constraints->AtomNumbersPresent();
  for (int a : atom_numbers_present) {
    const Substructure_Atom* q = query_atom_with_initial_atom_number(a);
    if (q == nullptr) {
      cerr << "Single_Substructure_Query::_atom_numbers_in_geometric_constraints_ok:no query atom " << a << '\n';
      return 0;
    }
  }

  return 1;
}

// Note that this silently does the wrong thing if called multiple times.
int
Single_Substructure_Query::ConstructFromProto(const SubstructureSearch::SingleSubstructureQuery& proto)
{
  assert (ok());

  if (! _construct_from_proto(proto))
  {
    cerr << "SingleSubstructureQuery::ConstructFromProto:cannot build from proto\n";
    cerr << proto.ShortDebugString() << "\n";
    return 0;
  }

// Now that all the atoms have been defined, we can read in any chirality

  for (const auto& chiral : proto.chiral_centre())
  {
    if (! _build_chirality_specification_from_proto(chiral))
    {
      cerr << "Single_Substructure_Query::ConstructFromProto:invalid chirality '" << chiral.ShortDebugString() << "'\n";
      return 0;
    }
  }

  _preferences_present = 0;
  for (const auto* root_atom : _root_atoms) {
    if (root_atom->preferences_present()) {
      _preferences_present = 1;
      break;
    }
  }

  _fragment_ids_present = 0;
  for (const auto * root_atom : _root_atoms)
  {
    if (root_atom->fragment_ids_present())
    {
      _fragment_ids_present = 1;
      break;
    }
  }

  assert (ok());

  (void) min_atoms_in_query();

  for (int i = 0; i < _no_matched_atoms_between.number_elements(); i++)
  {
    const Bond * b = _no_matched_atoms_between[i];
    if (b->a1() >= _min_atoms_in_query || b->a2() >= _min_atoms_in_query)
    {
      cerr << "Single_Substructure_Query::ConstructFromProto: illegal no_matched_atoms_between specifier\n";
      cerr << "There are as few as " << _min_atoms_in_query << " query atoms\n";
      cerr << "No matched atoms between " << i << " specifies atoms " << b->a1() << " and " << b->a2() << endl;
      cerr << "this is impossible\n";
      return 0;
    }
  }

  return 1;
}

int
Single_Substructure_Query::ReadProto(iwstring_data_source & input)
{
  assert (input.ok());

  input.set_ignore_pattern("^#");
  input.set_skip_blank_lines(1);

  IWString contents;
  const_IWSubstring line;
  while (input.next_record(line))
  {
    if (contents.length() > 0)
      contents << "\n";
    contents << line;
    if (line.starts_with('}'))
      break;
  }

  const std::string s(contents.rawchars(), contents.length());

  SubstructureSearch::SingleSubstructureQuery proto;

  if (!google::protobuf::TextFormat::ParseFromString(s, &proto))
  {
    cerr << "SingleSubstructureQuery::ReadProto:cannot parse proto data\n";
    cerr << s << endl;
    return 0;
  }

  return ConstructFromProto(proto);
}

int
Single_Substructure_Query::ReadProto(const char * fname)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Single_Substructure_Query::read: cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadProto(input);
}

int
Single_Substructure_Query::ReadProto(const IWString & fname)
{
  iwstring_data_source input(fname);

  if (! input.ok())
  {
    cerr << "Single_Substructure_Query::read: cannot open '" << fname << "'\n";
    return 0;
  }

  return ReadProto(input);
}

SubstructureSearch::SubstructureQuery
Substructure_Query::BuildProto() const
{
  SubstructureSearch::SubstructureQuery to_be_returned;

  for (int i = 0; i < _number_elements; i++)
  {
    _things[i]->assign_unique_numbers();
  }

  to_be_returned.set_comment(_comment.rawchars(), _comment.length());

  for (const auto* query : *this)
  {
    query->BuildProto(*to_be_returned.add_query());
  }

  if (_each_component_search)
    to_be_returned.set_match_each_component(true);

  return to_be_returned;
}

int
Substructure_Query::WriteProto(IWString & fname) const
{
  IWString_and_File_Descriptor output;
  if (! output.open(fname.null_terminated_chars()))
  {
    cerr << "Single_Substructure_Query::write: cannot open '" << fname << "'\n";
    return 0;
  }

  return WriteProto(output);
}

int
Substructure_Query::WriteProto(IWString_and_File_Descriptor& output) const
{
  const SubstructureSearch::SubstructureQuery proto = BuildProto();

  std::string s;
  if (! google::protobuf::TextFormat::PrintToString(proto, &s))
  {
    cerr << "Substructure_Query::WriteProto:PrintToString failed\n";
    return 0;
  }

  return output.write(s.data(), s.length());
}

int
Substructure_Ring_Base::ConstructFromProto(const SubstructureSearch::SubstructureRingBase& proto)
{
  if (proto.has_match_as_match())
    _match_as_match_or_rejection = ! proto.match_as_match();

  if (proto.has_all_hits_in_same_fragment())
    _all_hits_in_same_fragment = proto.all_hits_in_same_fragment();
    
  if (proto.has_environment())
  {
    if (_environment_atom.number_elements())
    {
      cerr << "Substructure_Ring_Base::ConstructFromProto:environment already specified '" << proto.ShortDebugString() << "'\n";
      return 0;
    }

    const IWString env = proto.environment();
    
    if (! _construct_environment(env))
    {
      cerr << "Substructure_Ring_Base::ConstructFromProto:invalid environment '" << env << "'\n";
      return 0;
    }
  }

  if (proto.has_environment_can_match_in_ring_atoms())
    _environment_can_match_in_ring_atoms = proto.environment_can_match_in_ring_atoms();

  if (proto.has_set_global_id()) {
    _set_global_id = proto.set_global_id();
  }

  static constexpr uint32_t no_limit = std::numeric_limits<uint32_t>::max();

  if (!GETVALUES(proto, hits_needed, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, ncon, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, heteroatom_count, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, attached_heteroatom_count, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, within_ring_unsaturation, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, atoms_with_pi_electrons, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, largest_number_of_bonds_shared_with_another_ring, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, strongly_fused_ring_neighbours, 0, no_limit))
    return 0;

  return 1;
}


int
Substructure_Ring_Specification::ConstructFromProto(SubstructureSearch::SubstructureRingSpecification const& proto)
{
  static constexpr uint32_t no_limit = std::numeric_limits<uint32_t>::max();

  if (! Substructure_Ring_Base::ConstructFromProto(proto.base()))
    return 0;

  if (!GETVALUES(proto, ring_size, 3, no_limit))
    return 0;

  if (!GETVALUES(proto, fused, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, fused_aromatic_neighbours, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, fused_non_aromatic_neighbours, 0, no_limit))
    return 0;

  if (proto.has_aromatic())
  {
    if (proto.aromatic())
      _aromatic = AROMATIC;
    else
      _aromatic = NOT_AROMATIC;
  }

  return 1;
}

int
Substructure_Ring_System_Specification::ConstructFromProto(const SubstructureSearch::SubstructureRingSystemSpecification& proto)
{
  if (! Substructure_Ring_Base::ConstructFromProto(proto.base()))
    return 0;

  if (_atoms_with_pi_electrons.is_set())
    _need_per_atom_array = 1;

  static constexpr uint32_t no_limit = std::numeric_limits<uint32_t>::max();

  if (!GETVALUES(proto, rings_in_system, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, ring_sizes, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, aromatic_ring_count, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, non_aromatic_ring_count, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, degree_of_fusion, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, atoms_in_system, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, strongly_fused_ring_count, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, number_spinach_groups, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, number_non_spinach_groups, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, atoms_in_spinach_group, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, length_of_spinach_group, 0, no_limit))
    return 0;

  if (!GETVALUES(proto, distance_to_another_ring, 0, no_limit))
    return 0;

  extending_resizable_array<int> ring_sizes_encountered;

  for (const auto& ring_size_count : proto.ring_size_count()) {
    RingSizeCount * r = new RingSizeCount();
    if (! r->ConstructFromProto(ring_size_count))
    {
      cerr << "Substructure_Ring_System_Specification::ConstructFromProto:invalid ring size\n";
      cerr << ring_size_count.ShortDebugString() << endl;
      delete r;
      return 0;
    }

    if (ring_sizes_encountered[ring_size_count.ring_size()])
    {
      cerr << "Substructure_Ring_System_Specification::ConstructFromProto:duplicate ring size specificiation in ring size count\n";
      cerr << proto.ShortDebugString() << endl;
      delete r;
      return 0;
    }

    ring_sizes_encountered[ring_size_count.ring_size()] = 1;

    _ring_size_count.add(r);
  }

  return 1;
}

RingSizeCount::RingSizeCount() {
  _ring_size = 0;
}

int
RingSizeCount::ConstructFromProto(const SubstructureSearch::RingSizeRequirement& proto) 
{
  if (!proto.has_ring_size())
  {
    cerr << "RingSizeCount::ConstructFromProto:no ring_size " << proto.ShortDebugString() << endl;
    return 0;
  }

  static constexpr uint32_t no_limit = std::numeric_limits<uint32_t>::max();

  _ring_size = proto.ring_size();

  if (!GETVALUES(proto, count, 0, no_limit))
    return 0;

  return 1;
}

int
Single_Substructure_Query::_parse_element_hits_needed(const SubstructureSearch::SingleSubstructureQuery& proto)
{
  for (const auto& spec : proto.element_hits_needed())
  {
    Elements_Needed * e = new Elements_Needed();
    if (! e->ConstructFromProto(spec))
    {
      cerr << "SingleSubstructureQuery::_parse_element_hits_needed:invalid input " << spec.ShortDebugString() << "\n";
      delete e;
      return 0;
    }

    _element_hits_needed.add(e);
  }

  return 1;
}

int
Single_Substructure_Query::_parse_elements_needed(const SubstructureSearch::SingleSubstructureQuery& proto)
{
  for (const auto & spec : proto.elements_needed())
  {
    Elements_Needed * e = new Elements_Needed();
    if (! e->ConstructFromProto(spec))
    {
      cerr << "Single_Substructure_Query::_parse_elements_needed:invalid " << spec.ShortDebugString() << "\n";
      delete e;
      return 0;
    }

    _elements_needed.add(e);
  }

  return 1;
}

int
Elements_Needed::ConstructFromProto(const SubstructureSearch::ElementsNeeded& proto)
{
  if (proto.has_atomic_number())
  {
    _z = proto.atomic_number();
    if (! REASONABLE_ATOMIC_NUMBER(_z))
    {
      cerr << "Elements_Needed::construct_from_proto:invalid atomic number " << proto.ShortDebugString() << "\n";
      return 0;
    }
  }
  else if (proto.has_atomic_symbol())
  {
    const IWString s = proto.atomic_symbol();
    int notused;
    const Element * e = get_element_from_symbol(s, notused);
    if (nullptr == e)
    {
      cerr << "Elements_Needed::construct_from_proto:invalid atomic symbol " << proto.ShortDebugString() << "\n";
      return 0;
    }
    _z = e->atomic_number();
  }

  static constexpr uint32_t no_limit = std::numeric_limits<uint32_t>::max();

  if (! GETVALUES(proto, hits_needed, 0, no_limit))
    return 0;

  return 1;
}

int
Link_Atom::ConstructFromProto(const SubstructureSearch::LinkAtoms & proto)
{
  if (!proto.has_a1() || !proto.has_a2())
  {
    cerr << "LinkAtoms::ConstructFromProto:incomplete atom spec " << proto.ShortDebugString() << "\n";
    return 0;
  }

  if (proto.distance_size() > 0)
    ;
  else if (proto.has_min_distance() || proto.has_max_distance())
    ;
  else
  {
    cerr << "LinkAtoms::ConstructFromProto:incomplete distance spec " << proto.ShortDebugString() << "\n";
    return 0;
  }

  _a1 = proto.a1();
  _a2 = proto.a2();
  if (_a1 == _a2)
  {
    cerr << "Link_Atom::ConstructFromProto:atoms must be distinct " << proto.ShortDebugString() << "\n";
    return 0;
  }

  static constexpr uint32_t no_limit = std::numeric_limits<uint32_t>::max();

  if (!GETVALUES(proto, distance, 0, no_limit))
    return 0;

  return 1;
}


#ifdef IMPLEMENT_THIS
TODO
int
Single_Substructure_Query::_add_component(const SubstructureSearch::SubstructureChiralCenter::AtomNumberHydrogenLonePair&& atom_or,
    void (Substructure_Chiral_Centre::*)(atom_number_t) setter,
    resizable_array<int>& numbers_encountered)
{
}
#endif

int
Single_Substructure_Query::_build_chirality_component(const SubstructureSearch::SubstructureChiralCenter& proto,
                                                      void (Substructure_Chiral_Centre::* ptrmfn)(const Substructure_Atom*),
                                                      const SubstructureSearch::AtomNumberHydrogenLonePair& atom_or,
                                                      Substructure_Chiral_Centre& c,
                                                      resizable_array<int>& numbers_encountered)
{
  static constexpr int hydrogen = -9;

  if (atom_or.has_atom_number())
  {
    const int uid = atom_or.atom_number();
    const Substructure_Atom * a = query_atom_with_initial_atom_number(uid);
    if (nullptr == a)
    {
      cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:no top front atom '" << proto.ShortDebugString() << "'\n";
      return 0;
    }

    if (!numbers_encountered.add_if_not_already_present(uid))
    {
      cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:duplicate atom number " << proto.ShortDebugString() << endl;
      return 0;
    }
    (c.*ptrmfn)(a);
  }
  else  // hydrogen or lone pair 
  { 
    const auto h_or_lp = atom_or.h_or_lp();
    if (h_or_lp == SubstructureSearch::AtomNumberHydrogenLonePair::HYDROGEN)
      ;
    else if (h_or_lp == SubstructureSearch::AtomNumberHydrogenLonePair::LONEPAIR) 
    {
      cerr <<"Single_Substructure_Query::_build_chirality_specification_from_proto:lp directive not supported yet " << proto.ShortDebugString() << endl;
      return 0;
    }
    else
    {
      cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:unrecognised chirality component " << proto.ShortDebugString() << endl;
      return 0;
    }

    if (!numbers_encountered.add_if_not_already_present(hydrogen))
    {
      cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:multiple H " << proto.ShortDebugString() << endl;
      return 0;
    }
  }

  return 1;
}

int
Single_Substructure_Query::_build_chirality_specification_from_proto(const SubstructureSearch::SubstructureChiralCenter& proto)
{
  if (! proto.has_center())
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:chirality with no center atom " << proto.ShortDebugString() << endl;
    return 0;
  }

  if (! proto.has_top_front() || ! proto.has_top_back() ||
      ! proto.has_left_down() || ! proto.has_right_down())
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:all components of a chiral center must be specified\n";
    cerr << proto.ShortDebugString() << endl;
    return 9;
  }

  Substructure_Chiral_Centre * c = new Substructure_Chiral_Centre;

  const int uid = proto.center();

  const Substructure_Atom * a = query_atom_with_initial_atom_number(uid);
  if (nullptr == a)
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:no centre atom id " <<
            uid << " in '" << proto.ShortDebugString() << "'\n";
    return 0;
  }

  c->set_centre(a);

  resizable_array<int> numbers_encountered;   // make sure no duplicates

  numbers_encountered.add(uid);

  // Top front.

  auto atom_or = proto.top_front();

  if (! _build_chirality_component(proto, &Substructure_Chiral_Centre::set_top_front, proto.top_front(), *c, numbers_encountered))
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:invalid top front " << proto.ShortDebugString() << endl;
    return 0;
  }

  if (! _build_chirality_component(proto, &Substructure_Chiral_Centre::set_top_back, proto.top_back(), *c, numbers_encountered))
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:invalid top back " << proto.ShortDebugString() << endl;
    return 0;
  }

  if (! _build_chirality_component(proto, &Substructure_Chiral_Centre::set_left_down, proto.left_down(), *c, numbers_encountered))
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:invalid left down back " << proto.ShortDebugString() << endl;
    return 0;
  }

  if (! _build_chirality_component(proto, &Substructure_Chiral_Centre::set_right_down, proto.right_down(), *c, numbers_encountered))
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_proto:invalid right down back " << proto.ShortDebugString() << endl;
    return 0;
  }

  if (5 != numbers_encountered.number_elements())
  {
    cerr << "Single_Substructure_Query::_build_chirality_specification_from_msi_attribute:duplicate atom numbers '" << proto.ShortDebugString() << "'\n";
    return 0;
  }

  _chirality.add(c);

  return 1;
}

int
Substructure_Query::ConstructFromProto(const SubstructureSearch::SubstructureQuery& proto)
{
  if (proto.has_comment())
    _comment = proto.comment();

  if (0 == proto.query_size())
  {
    cerr << "Substructure_Query::ConstructFromProto:no components in query\n";
    return 0;
  }

  for (const auto& query : proto.query())
  {
    Single_Substructure_Query * q = new Single_Substructure_Query();

    if (! q->ConstructFromProto(query))
    {
      cerr << "Substructure_Query::ConstructFromProto:cannot parse component\n";
      cerr << query.ShortDebugString() << endl;
      delete q;
      return 0;
    }

    add(q);
  }

  if (proto.has_match_each_component())
    _each_component_search = proto.match_each_component();

  if (1 == _number_elements)
  {
    if (0 == proto.logexp_size())  // Great.
      return 1;

    cerr << "Substructure_Query::ConstructFromProto:operator(s) with one component\n";
    cerr << proto.ShortDebugString() << endl;
    return 0;
  }

  if (0 == proto.logexp_size())
    return _number_elements;

  if (_each_component_search && proto.logexp_size() > 0)
  {
    cerr << "Substructure_Query::ConstructFromProto:cannot combine match_each_component with operators\n";
    cerr << proto.ShortDebugString() << endl;
    return 0;
  }

  _operator.RemoveAllOperators();

  if (! ExtractOperator(proto.logexp(), _number_elements - 1, _operator, IW_LOGEXP_OR, "Substructure_Query::ConstructFromProto"))
  {
    cerr << "Substructure_Query::ConstructFromProto:cannot process operators\n";
    cerr << proto.ShortDebugString() << endl;
    return 0;
  }

  return _number_elements;
}

int
Substructure_Ring_Base::BuildProto(SubstructureSearch::SubstructureRingBase& proto) const {
  if (! _match_as_match_or_rejection) {  // Default is true.
    proto.set_match_as_match(_match_as_match_or_rejection);  // 1
  }

  SETPROTOVALUES(proto, hits_needed, int);  // 2
  SETPROTOVALUES(proto, attached_heteroatom_count, int);  // 5
  SETPROTOVALUES(proto, heteroatom_count, int);  // 8
  SETPROTOVALUES(proto, ncon, int);  // 11
  if (_all_hits_in_same_fragment) {  // 14
    proto.set_all_hits_in_same_fragment(true);
  }
  SETPROTOVALUES(proto, within_ring_unsaturation, int);  // 16
  SETPROTOVALUES(proto, largest_number_of_bonds_shared_with_another_ring, int);  // 19
  SETPROTOVALUES(proto, atoms_with_pi_electrons, int);  // 26
  SETPROTOVALUES(proto, strongly_fused_ring_neighbours, int);  // 29
  // environment as string....
  if (_environment_can_match_in_ring_atoms) {
    proto.set_environment_can_match_in_ring_atoms(true);  // 23
  }
  if (_set_global_id) {
    proto.set_set_global_id(_set_global_id);
  }

  return 1;
}

int
Substructure_Ring_Specification::BuildProto(SubstructureSearch::SubstructureRingSpecification& proto) const {

  ::Substructure_Ring_Base::BuildProto(*proto.mutable_base());

  SETPROTOVALUES(proto, ring_size, int);  // 2
  if (_aromatic) {
    proto.set_aromatic(true);   // 5
  }
  SETPROTOVALUES(proto, fused, int);   // 6
  SETPROTOVALUES(proto, fused_aromatic_neighbours, int);   // 9
  SETPROTOVALUES(proto, fused_non_aromatic_neighbours, int);   // 12

  return 1;
}

int
Substructure_Ring_System_Specification::BuildProto(SubstructureSearch::SubstructureRingSystemSpecification& proto) const {
  SETPROTOVALUES(proto, rings_in_system, int);  // 2
  SETPROTOVALUES(proto, ring_sizes, int);  // 5

  for (const RingSizeCount * rsc : _ring_size_count) {
    SubstructureSearch::RingSizeRequirement* s = proto.add_ring_size_count();
    rsc->BuildProto(*s);
  }
  SETPROTOVALUES(proto, aromatic_ring_count, int);  // 11
  SETPROTOVALUES(proto, non_aromatic_ring_count, int);  // 14
  SETPROTOVALUES(proto, degree_of_fusion, int);  // 17
  SETPROTOVALUES(proto, atoms_in_system, int);  // 20
  SETPROTOVALUES(proto, number_spinach_groups, int);  // 23
  SETPROTOVALUES(proto, number_non_spinach_groups, int);  // 26
  SETPROTOVALUES(proto, atoms_in_spinach_group, int);  // 29
  SETPROTOVALUES(proto, length_of_spinach_group, int);  // 32
  SETPROTOVALUES(proto, distance_to_another_ring, int);  // 35
  SETPROTOVALUES(proto, strongly_fused_ring_count, int);  // 38

  return 1;
}

int
RingSizeCount::BuildProto(SubstructureSearch::RingSizeRequirement & proto) const {
  proto.set_ring_size(_ring_size);
  SETPROTOVALUES(proto, count, int);
  return 1;
}

int
Elements_Needed::BuildProto(SubstructureSearch::ElementsNeeded& proto) const {
  proto.set_atomic_number(_z);
  SETPROTOVALUESMATCHER(proto, hits_needed, int);
  return 1;
}

int
Link_Atom::BuildProto(SubstructureSearch::LinkAtoms& proto) const {
  proto.set_a1(_a1);
  proto.set_a2(_a2);
  SETPROTOVALUESMATCHER(proto, distance, int);
  SetProtoValues(_distance, "distance", proto);
  return 1;
}

int
Substructure_Chiral_Centre::BuildProto(SubstructureSearch::SubstructureChiralCenter& proto) const {
  proto.set_center(_numeric->a());
  const auto tf = _numeric->top_front();
  if (tf == CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN) {
    proto.mutable_top_front()->set_h_or_lp(SubstructureSearch::AtomNumberHydrogenLonePair::HYDROGEN);
  } else if (tf == CHIRAL_CONNECTION_IS_LONE_PAIR) {
    proto.mutable_top_front()->set_h_or_lp(SubstructureSearch::AtomNumberHydrogenLonePair::LONEPAIR);
  } else {
    proto.mutable_top_front()->set_atom_number(tf);
  }

  const auto tb = _numeric->top_back();
  if (tb == CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN) {
    proto.mutable_top_back()->set_h_or_lp(SubstructureSearch::AtomNumberHydrogenLonePair::HYDROGEN);
  } else if (tb == CHIRAL_CONNECTION_IS_LONE_PAIR) {
    proto.mutable_top_back()->set_h_or_lp(SubstructureSearch::AtomNumberHydrogenLonePair::LONEPAIR);
  } else {
    proto.mutable_top_back()->set_atom_number(tb);
  }

  const auto ld = _numeric->left_down();
  if (ld == CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN) {
    proto.mutable_left_down()->set_h_or_lp(SubstructureSearch::AtomNumberHydrogenLonePair::HYDROGEN);
  } else if (ld == CHIRAL_CONNECTION_IS_LONE_PAIR) {
    proto.mutable_left_down()->set_h_or_lp(SubstructureSearch::AtomNumberHydrogenLonePair::LONEPAIR);
  } else {
    proto.mutable_left_down()->set_atom_number(ld);
  }

  const auto rd = _numeric->right_down();
  if (rd == CHIRAL_CONNECTION_IS_IMPLICIT_HYDROGEN) {
    proto.mutable_right_down()->set_h_or_lp(SubstructureSearch::AtomNumberHydrogenLonePair::HYDROGEN);
  } else if (rd == CHIRAL_CONNECTION_IS_LONE_PAIR) {
    proto.mutable_right_down()->set_h_or_lp(SubstructureSearch::AtomNumberHydrogenLonePair::LONEPAIR);
  } else {
    proto.mutable_right_down()->set_atom_number(rd);
  }
  return 1;
}

int
SeparatedAtoms::BuildProto(SubstructureSearch::SeparatedAtoms& proto) const {
  proto.set_a1(_a1);
  proto.set_a2(_a2);
  SetProtoValues(_separation, "bonds_between", proto);
  return 1;
}

int
MatchedAtomMatch::ConstructFromProto(const SubstructureSearch::MatchedAtomMatch& proto) {
  if (proto.atom_size() == 0) {
    cerr << "MatchedAtomMatch::ConstructFromProto:no atom " << proto.ShortDebugString() << '\n';
    return 0;
  }

  if (proto.smarts_size() == 0) {
    cerr << "MatchedAtomMatch::ConstructFromProto:no smarts " << proto.ShortDebugString() << '\n';
    return 0;
  }

  for (int atom : proto.atom()) {
    _atoms << atom;
  }

  for (const std::string smarts : proto.smarts()) {
    if (smarts.empty()) {
      cerr << "MatchedAtomMatch::ConstructFromProto:skipping empty smarts\n";
      continue;
    }
    bool match_as_match;
    const_IWSubstring s;
    if (smarts[0] == '!') {
      match_as_match = false;
      s = smarts;
      s += 1;
    } else {
      match_as_match = true;
      s = smarts;
    }

    std::unique_ptr<Single_Substructure_Query> a = std::make_unique<Single_Substructure_Query>();
    if (! a->create_from_smarts(s)) {
      cerr << "MatchedAtomMatch::ConstructFromProto:invalid smarts '" << smarts << "'\n";
      cerr << proto.ShortDebugString() << '\n';
      return 0;
    }

    if (match_as_match) {
      _positive_matches << a.release();
    } else {
      _negative_matches << a.release();
    }
  }

  return 1;
}

int
MatchedAtomMatch::BuildProto(SubstructureSearch::MatchedAtomMatch& proto) const {
  for (int atom : _atoms) {
    proto.add_atom(atom);
  }

  cerr << "MatchedAtomMatch::BuildProto:cannot emit queries\n";

  return 1;
}
