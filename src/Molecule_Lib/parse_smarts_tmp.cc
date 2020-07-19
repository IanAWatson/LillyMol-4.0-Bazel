#include <stdlib.h>

#include "substructure.h"
#include "parse_smarts_tmp.h"

Parse_Smarts_Tmp::Parse_Smarts_Tmp ()
{
  _last_query_atom_created = -1;

  return;
}

Parse_Smarts_Tmp::~Parse_Smarts_Tmp ()
{
  return;
}

int
Parse_Smarts_Tmp::set_natoms (int n)
{
  assert (n > 0);

//cerr << "Parse_Smarts_Tmp:set_natoms: natoms " << n << endl;

  return 1;
}

int
ThreeDots::BuildFromProto(const SubstructureSearch::NoMatchedAtomsBetween& proto)
{
  if (! proto.has_a1() || ! proto.has_a2()) {
    cerr << "ThreeDots::BuildFromProto:proto missing atom(s) " << proto.ShortDebugString() << endl;
    return 0;
  }

  _matched_atom_1 = proto.a1();
  _matched_atom_2 = proto.a2();

  if (_matched_atom_1 == _matched_atom_2) {
    cerr << "ThreeDots::BuildFromProto:Invalid matched atoms " << proto.ShortDebugString() << endl;
    return 0;
  }

  if (proto.has_qualifier()) {
    _qualifier = proto.qualifier();
  }

  return 1;
}
