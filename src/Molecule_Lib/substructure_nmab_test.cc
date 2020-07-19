// Tester for the No_Matched_Atoms

#include "substructure.h"

#include "googlemock/include/gmock/gmock.h"
#include "googletest/include/gtest/gtest.h"
#include "google/protobuf/text_format.h"

namespace {

using testing::Eq;
using testing::ElementsAre;
using testing::IsEmpty;
using testing::Property;
using testing::ResultOf;

class TestNMABToken : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    const_IWSubstring s;

    resizable_array_p<NMAB_Token> tokens;
};

void
TestNMABToken::SetUp() 
{
}

TEST_F(TestNMABToken, Empty) {
  s = "{}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, Incomplete1) {
  s = "{";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, Incomplete2) {
  s = "}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, Incomplete3) {
  s = "a}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, Incomplete4) {
  s = "{bb";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, SingleNumber) {
  s = "{3}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 1);

  EXPECT_EQ(tokens[0]->numbers()[0], 3);
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::numbers, ElementsAre(3)));
}

TEST_F(TestNMABToken, InvalidRange1) {
  s = "{3-}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, InvalidRange2) {
  s = "{3-2}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, InvalidRange3) {
  s = "{-2}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, Range1) {
  s = "{3-3}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 1);
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::numbers, ElementsAre(3)));
}

TEST_F(TestNMABToken, Range35) {
  s = "{3-5}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 1);
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::numbers, ElementsAre(3, 4, 5)));
}

TEST_F(TestNMABToken, LessThan) {
  s = "{<3}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 1);
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::relational, -1));
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::numbers, ElementsAre(3)));
}

TEST_F(TestNMABToken, GreaterThan) {
  s = "{<0}";

  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 1);
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::relational, -1));
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::numbers, ElementsAre(0)));
//EXPECT_THAT(*tokens[0], ResultOf([](const NMAB_Token& nmbt) { return nmbt.number1();}, Eq(0)));
}

TEST_F(TestNMABToken, JustOperator) {
  s = "{&}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);

  s = "{,}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);

  s = "{^}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);

  s = "{;}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, JustRelational) {
  s = "{<}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);

  s = "{>}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, JustOpAndRelational) {
  s = "{&<}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);

  s = "{;>}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, EmptySmarts) {
  s = "{[]}";
  cerr << "first\n";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
  cerr << "second\n";
  s = "{1[]}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
  cerr << "third\n";
  s = "{>1[]}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
  s = "{,>1[]}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
  s = "{,4-5[]}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, ValidSmarts1) {
  s = "{[c]}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 1);
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::op, Eq(0)));
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::relational, 0));
// Not sure why these fail. I want to test negative, not a specific -ve value.
//EXPECT_THAT(*tokens[0], Property(&NMAB_Token::number1, testing::Lt(0)));
//EXPECT_THAT(*tokens[0], Property(&NMAB_Token::number2, testing::Lt(0)));
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::numbers, IsEmpty()));

  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::smarts, Eq("c")));
}

TEST_F(TestNMABToken, ValidSmarts2) {
  s = "{[cc]}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 1);
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::op, Eq(0)));
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::relational, 0));
// Not sure why these fail. I want to test negative, not a specific -ve value.
//EXPECT_THAT(*tokens[0], Property(&NMAB_Token::number1, testing::Lt(0)));
//EXPECT_THAT(*tokens[0], Property(&NMAB_Token::number2, testing::Lt(0)));

  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::smarts, Eq("cc")));
}

TEST_F(TestNMABToken, NoFirstOp) {
  s = "{;3,4}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, NoOperator) {
  s = "{[c]3}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 0);
}

TEST_F(TestNMABToken, Combinations1) {
  s = "{<2[cc],3-5[n]}";
  EXPECT_EQ(TokeniseNMABSpecification(s, tokens), 2);

  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::op, Eq(0)));
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::relational, -1));
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::numbers, ElementsAre(2)));
  EXPECT_THAT(*tokens[0], Property(&NMAB_Token::smarts, Eq("cc")));

  EXPECT_THAT(*tokens[1], Property(&NMAB_Token::op, Eq(2)));
  EXPECT_THAT(*tokens[1], Property(&NMAB_Token::relational, 0));
  EXPECT_THAT(*tokens[1], Property(&NMAB_Token::numbers, ElementsAre(3, 4, 5)));
  EXPECT_THAT(*tokens[1], Property(&NMAB_Token::smarts, Eq("n")));
}

class TestNMAB : public testing::Test
{
  protected:
    void SetUp() override;

  protected:
    Single_Substructure_Query _query;
    Molecule _m;
    Substructure_Results _sresults;
};

void
TestNMAB::SetUp()
{
}

TEST_F(TestNMAB, TestBadSmarts1) {
  EXPECT_EQ(_query.create_from_smarts("C..."), false);
}

TEST_F(TestNMAB, TestBadSmarts2) {
  EXPECT_EQ(_query.create_from_smarts("..."), false);
}

TEST_F(TestNMAB, TestBadSmarts3) {
  EXPECT_EQ(_query.create_from_smarts("...C"), false);
}

TEST_F(TestNMAB, EstersMatch) {
  ASSERT_TRUE(_query.create_from_smarts("COC(=O)...C(=O)OC"));
  ASSERT_TRUE(_m.build_from_smiles("COC(=O)CCCCCC(=O)OC"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, EstersNoMatch) {
  ASSERT_TRUE(_query.create_from_smarts("COC(=O)...C(=O)OC"));
  ASSERT_TRUE(_m.build_from_smiles("CC(=O)OCCCCCC(=O)OC"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 0);
}

TEST_F(TestNMAB, Match1) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cc(O)ccc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1Distance) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{1}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cc(O)ccc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1NoMatchDistance) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{2}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cc(O)ccc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 0);
}

TEST_F(TestNMAB, Match1NoMatchDistanceLongPath) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{3}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cc(O)ccc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 0);
}

TEST_F(TestNMAB, Match1MatchSmartsAll) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{[c]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cc(O)ccc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1MatchSmartsCount) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{1[c]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cc(O)ccc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1MatchSmartsGeater) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{>0[c]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cc(O)ccc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1MatchSmartsLess) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{<3[c]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cc(O)ccc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1MatchSmartsNotAllSame1) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{>0[c]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cnc(O)cc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1MatchSmartsNotAllSame2) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{>0[n]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cnc(O)cc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
  ASSERT_TRUE(_m.build_from_smiles("Oc1ccc(O)nc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1MatchSmartsNotAllSame3) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{>0[n],>0[c]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cnc(O)cc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
  ASSERT_TRUE(_m.build_from_smiles("Oc1ccc(O)nc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1MatchSmartsExactCount) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{1[n];1[c]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cnc(O)cc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
  ASSERT_TRUE(_m.build_from_smiles("Oc1ccc(O)nc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, Match1MatchSmartsOr) {
  ASSERT_TRUE(_query.create_from_smarts("[OH]-c...{[n,c]}c-[OH]"));
  ASSERT_TRUE(_m.build_from_smiles("Oc1cnc(O)cc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
  ASSERT_TRUE(_m.build_from_smiles("Oc1ccc(O)nc1"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, AtomsOrSeparation) {
  ASSERT_TRUE(_query.create_from_smarts("O...{>5,>2[C]}O"));
  ASSERT_TRUE(_m.build_from_smiles("ONNNNNNO"));
  _query.set_find_unique_embeddings_only(1);
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 1);
  ASSERT_TRUE(_m.build_from_smiles("OCNCCO"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 1);
}

TEST_F(TestNMAB, AllRingAtoms) {
  ASSERT_TRUE(_query.create_from_smarts("O...{[R]}O"));
  ASSERT_TRUE(_m.build_from_smiles("ONNNNNNO"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 0);
  ASSERT_TRUE(_m.build_from_smiles("OC1C2C(O)C12"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}

TEST_F(TestNMAB, AllRingAtomsButNotBonds) {
  ASSERT_TRUE(_query.create_from_smarts("O...{[R]}O"));
  ASSERT_TRUE(_m.build_from_smiles("OC1CC1C2CC2O"));
  EXPECT_EQ(_query.substructure_search(&_m, _sresults), 2);
}


}  // namespace
