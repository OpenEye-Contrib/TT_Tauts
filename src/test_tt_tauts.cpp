//
// file tet_tt_tauts.cpp
// David Cosgrove
// CozChemIx Ltd
// 3rd April 2017.
//
// Tests for the tt_tauts code.
// Used the Catch testing library in catch.hpp.
// https://github.com/philsquared/Catch

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <limits>
#include <string>
#include <vector>

unsigned int make_taut_skeleton( const std::string &in_smi ,
                                 const std::string &mol_name ,
                                 float max_time ,
                                 std::string &t_skel_smi,
                                 bool &timed_out );
void make_taut_skeleton_and_tauts(const std::string &in_smi,
                                  const std::string &mol_name ,
                                  std::string &t_skel_smi ,
                                  std::vector<std::string> &taut_smis,
                                  bool &timed_out,
                                  float max_time = std::numeric_limits<float>::max() );
// take a SMILES string and a molecule name, returns the SMILES for
// the t_skel. Max_time is in seconds.
std::string test_make_tautomer_skeleton(const std::string &in_smi,
                                        const std::string &mol_name,
                                        float max_time = std::numeric_limits<float>::max() ) {
  std::string t_skel_smi;
  bool timed_out = false;
  unsigned int num_tauts = make_taut_skeleton( in_smi, mol_name,
                                               max_time, t_skel_smi,
                                               timed_out );
  return t_skel_smi;
}

// take a SMILES string, genereate t_skel and all tautomers, make sure the
// tautomers all give the same t_skel - the "round-tripping test".
bool test_round_trips(const std::string &in_smi,
                      const std::string &mol_name,
                      float max_time = std::numeric_limits<float>::max() ) {
  std::string orig_t_skel_smi;
  std::vector<std::string> taut_smis;
  bool orig_timed_out = false;
  make_taut_skeleton_and_tauts(in_smi, mol_name, orig_t_skel_smi,
                               taut_smis, orig_timed_out, max_time );
  if(orig_timed_out) {
    std::cout << "Warning : " << mol_name << " timed out generatnig t_skel."
             << std::endl;
  }
#ifdef NOTYET
  std::cout << "Original t_skel_smi : " << orig_t_skel_smi << std::endl;
#endif
  for(size_t i = 0, is = taut_smis.size(); i < is; ++i){
    std::string t_skel_smi = test_make_tautomer_skeleton(taut_smis[i], mol_name);
    if( t_skel_smi != orig_t_skel_smi) {
      std::cout << taut_smis[i] << " gives " << t_skel_smi
                << " not " << orig_t_skel_smi << std::endl;
      return false;
    }
  }
  return true;
}

TEST_CASE( "Simple Test Cases" , "[test_tt_tauts]" ) {
  CHECK( test_make_tautomer_skeleton( "CCC", "simple_test1") == "CCC" );
  CHECK( test_make_tautomer_skeleton( "CC(O)C", "simple_test2") == "CC(C)O" );
  CHECK( test_make_tautomer_skeleton( "CC(=O)C", "simple_test3") == "[CH2][C]([CH2])[O]" );
  CHECK( test_make_tautomer_skeleton( "CC(=C)O", "simple_test4") == "[CH2][C]([CH2])[O]" );
  CHECK( test_make_tautomer_skeleton( "CC(=O)CC(=O)C", "simple_test5") == "[CH2][C]([CH][C]([CH2])[O])[O]" );
  CHECK( test_make_tautomer_skeleton( "CC(=N)NC", "simple_test6") == "C[N][C]([CH2])[NH]" );
  CHECK( test_make_tautomer_skeleton( "CC(=NC)N", "simple_test6") == "C[N][C]([CH2])[NH]" );
  CHECK( test_make_tautomer_skeleton( "CNC(=C)N", "simple_test6") == "C[N][C]([CH2])[NH]" );
}

TEST_CASE( "Historical Test Cases" , "[test_tt_tauts]" ) {
  // these are tests from before there was this automated testing.
  // Mostly taken from Chembl structures, some of them edited for size.
  // Test name may have Chembl number in it.
  CHECK( test_make_tautomer_skeleton( "c1cc(ccc1C/C(=N/c2ccc(cc2)O)/c3ccc(c(c3O)O)O)Cl", "hist_test_1082532") ==
          "c1cc(ccc1[CH][C]([C]2[CH][CH][C]([C]([C]2[O])[O])[O])[N][C]3[CH][CH][C]([CH][CH]3)[O])Cl" );
  CHECK( test_make_tautomer_skeleton( "CC(=O)CC(c1ccccc1)c2c(c3ccccc3oc2=O)O", "hist_test_1464") ==
          "[CH2][C]([CH]C(c1ccccc1)[C]2[C](c3ccccc3O[C]2[O])[O])[O]" );
  CHECK( test_make_tautomer_skeleton( "Cc1cc(cn1C)c2csc(n2)N=C(N)N", "hist_test_153534a") == "Cc1cc(cn1C)C2=CS[C]([N]2)[N][C]([NH])[NH]" );
  CHECK( test_make_tautomer_skeleton( "Cc1cc(cn1C)c2csc(=NC(=N)N)[nH]2", "hist_test_153534b") == "Cc1cc(cn1C)C2=CS[C]([N]2)[N][C]([NH])[NH]" );
  CHECK( test_make_tautomer_skeleton( "c1cc(ccc1C(=O)NO)I", "hist_test_155287a") == "c1cc(ccc1[C]([N][O])[O])I" );
  CHECK( test_make_tautomer_skeleton( "c1cc(ccc1C(=NO)O)I", "hist_test_155287b") == "c1cc(ccc1[C]([N][O])[O])I" );
  CHECK( test_make_tautomer_skeleton( "c1cc(ccc1C(N=O)O)I", "hist_test_155287c") == "c1cc(ccc1[C]([N][O])[O])I" );
  CHECK( test_make_tautomer_skeleton( "CC1=NNC(=C)NN1", "hist_test_18048a") == "[CH2][C]1[N][N][C]([N][N]1)[CH2]" );
  CHECK( test_make_tautomer_skeleton( "CC1=NNC(=C)NN1", "hist_test_18048b") == "[CH2][C]1[N][N][C]([N][N]1)[CH2]" );
  CHECK( test_make_tautomer_skeleton( "CC1NN=C(N=N1)C", "hist_test_18048c") == "[CH2][C]1[N][N][C]([N][N]1)[CH2]" );
  CHECK( test_make_tautomer_skeleton( "CC1N=NC(N=N1)C", "hist_test_18048d") == "[CH2][C]1[N][N][C]([N][N]1)[CH2]" );
  CHECK( test_make_tautomer_skeleton( "C=c1[nH][nH]c(=C)[nH][nH]1", "hist_test_18048e") == "[CH2][C]1[N][N][C]([N][N]1)[CH2]" );
  CHECK( test_make_tautomer_skeleton( "C[C@H](C(C)(C)C)/N=c\\1/c(=N/c2ccc(cc2O)C#N)/c(c1O)O", "hist_test_19253a") == "CC(C)(C)[C]([CH2])[N][C]1[C]([C]([C]1[O])[O])[N][C]2[CH][CH][C]([CH][C]2[O])[C]=[N]" );
  CHECK( test_make_tautomer_skeleton( "CC(=NC1=C(C(=C1NC2C=CC(=CC2=O)C#N)O)O)C(C)(C)C", "hist_test_19253b") == "CC(C)(C)[C]([CH2])[N][C]1[C]([C]([C]1[O])[O])[N][C]2[CH][CH][C]([CH][C]2[O])[C]=[N]" );
  CHECK( test_make_tautomer_skeleton( "C[C@H](C(C)(C)C)/N=c\\1/c(=N/c2ccc(cc2O)C#N)/c(c1O)O", "hist_test_19253c") == "CC(C)(C)[C]([CH2])[N][C]1[C]([C]([C]1[O])[O])[N][C]2[CH][CH][C]([CH][C]2[O])[C]=[N]" );
  CHECK( test_make_tautomer_skeleton( "C[C@H](C(C)(C)C)/N=c\\1/c(c(c1=O)O)Nc2ccc(cc2O)C#N", "hist_test_19253d") == "CC(C)(C)[C]([CH2])[N][C]1[C]([C]([C]1[O])[O])[N][C]2[CH][CH][C]([CH][C]2[O])[C]=[N]" );
  CHECK( test_make_tautomer_skeleton( "CC(C)(C)C(=C)Nc1c(c(=O)c1=O)NC2C=CC(C=C2O)C#N", "hist_test_19253e") == "CC(C)(C)[C]([CH2])[N][C]1[C]([C]([C]1[O])[O])[N][C]2[CH][CH][C]([CH][C]2[O])[C]=[N]" );
  CHECK( test_make_tautomer_skeleton( "C[C@@H](Cc1c2ccoc2c(c3c1occ3)Br)N.Cl", "hist_test_263083") == "C[C@@H](Cc1c2ccoc2c(c3c1occ3)Br)N.Cl" );
  CHECK( test_make_tautomer_skeleton( "CN(C)S(=O)(=O)CCNC(=O)C(CCCl)N=O", "hist_test_31034a") == "CN(C)S(=O)(=O)CCNC(=O)[C]([CH]CCl)[N][O]" );
  CHECK( test_make_tautomer_skeleton( "CN(C)S(=O)(=O)CCNC(=O)C(=NO)CCCl", "hist_test_31034b") == "CN(C)S(=O)(=O)CCNC(=O)[C]([CH]CCl)[N][O]" );
  CHECK( test_make_tautomer_skeleton( "CN(C)S(=O)(=O)CCNC(=O)C(CCCl)N=O", "hist_test_31034c") == "CN(C)S(=O)(=O)CCNC(=O)[C]([CH]CCl)[N][O]" );
  CHECK( test_make_tautomer_skeleton( "CCOc1ccccc1c2nc-3ccc[nH]c3n2", "hist_test_31359a") == "CCOc1ccccc1[C]2[N][C]3[CH][CH][CH][N][C]3[N]2" );
  CHECK( test_make_tautomer_skeleton( "CCOc1ccccc1c2[nH]c3c(n2)cccn3", "hist_test_31359b") == "CCOc1ccccc1[C]2[N][C]3[CH][CH][CH][N][C]3[N]2" );
  CHECK( test_make_tautomer_skeleton( "c1cc[n+](cc1)CCCCC[n+]2cccc(c2)NC(=O)C=NO", "hist_test_3306810a") == "c1cc[n+](cc1)CCCCC[n+]2cccc(c2)NC(=O)[CH][N][O]" );
  CHECK( test_make_tautomer_skeleton( "c1cc[n+](cc1)CCCCC[n+]2cccc(c2)NC(=O)CN=O", "hist_test_3306810b") == "c1cc[n+](cc1)CCCCC[n+]2cccc(c2)NC(=O)[CH][N][O]" );
  CHECK( test_make_tautomer_skeleton( "CN1CNC=C(C1=O)F", "hist_test_34961a") == "CN1C[N][CH][C](C1=O)F" );
  CHECK( test_make_tautomer_skeleton( "CN1CNC=C(C1=O)F", "hist_test_34961b") == "CN1C[N][CH][C](C1=O)F" );
  CHECK( test_make_tautomer_skeleton( "c1ccc-2c(c1)CCCc3c2[nH]c(n3)c4cccnc4", "hist_test_414205") == "c1ccc2c(c1)CCCC3=C2[N][C]([N]3)[C]4[CH][CH][CH][N][CH]4" );
  CHECK( test_make_tautomer_skeleton( "CC(C)[C@H]1C(=O)N[C@@H](CSSC[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N1)CCCCN)Cc2c[nH]c3c2cccc3)Cc4c[nH]cn4)NC(=O)[C@H](Cc5ccc6ccccc6c5)N)C(=O)N[C@@H](Cc7cccc(c7)c8ccccc8)C(=O)N", "hist_test_440067a") == "CC(C)[C@H]1C(=O)N[C@@H](CSSC[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N1)CCCCN)Cc2c[nH]c3c2cccc3)CC4=C[N][CH][N]4)NC(=O)[C@H](Cc5ccc6ccccc6c5)N)C(=O)N[C@@H](Cc7cccc(c7)c8ccccc8)C(=O)N" );
  CHECK( test_make_tautomer_skeleton( "CC(C)[C@H]1C(=O)N[C@@H](CSSC[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N1)CCCCN)Cc2c[nH]c3c2cccc3)Cc4cnc[nH]4)NC(=O)[C@H](Cc5ccc6ccccc6c5)N)C(=O)N[C@@H](Cc7cccc(c7)c8ccccc8)C(=O)N", "hist_test_440067b") == "CC(C)[C@H]1C(=O)N[C@@H](CSSC[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N1)CCCCN)Cc2c[nH]c3c2cccc3)CC4=C[N][CH][N]4)NC(=O)[C@H](Cc5ccc6ccccc6c5)N)C(=O)N[C@@H](Cc7cccc(c7)c8ccccc8)C(=O)N" );
  CHECK( test_make_tautomer_skeleton( "COc1ccccc1NC2=C(C(=O)CO2)C(=O)Nc3ccccc3OC", "hist_test_440484a") == "COc1ccccc1NC(=O)[C]2[C]([CH]O[C]2[N]c3ccccc3OC)[O]" );
  CHECK( test_make_tautomer_skeleton( "COc1ccccc1Nc2c(c(co2)O)C(=O)Nc3ccccc3OC", "hist_test_440484b") == "COc1ccccc1NC(=O)[C]2[C]([CH]O[C]2[N]c3ccccc3OC)[O]" );
  CHECK( test_make_tautomer_skeleton( "COc1ccccc1NC(=O)C2C(=O)COC2=Nc3ccccc3OC", "hist_test_440484c") == "COc1ccccc1NC(=O)[C]2[C]([CH]O[C]2[N]c3ccccc3OC)[O]" );
  CHECK( test_make_tautomer_skeleton( "CC(C(=O)c1nc2ccccc2o1)NC=O", "hist_test_500474a") == "[CH2][C]([C]([C]1[N]c2ccccc2O1)[O])[N][CH][O]" );
  CHECK( test_make_tautomer_skeleton( "CC(=NCO)C(=O)c1nc2ccccc2o1", "hist_test_500474b") == "[CH2][C]([C]([C]1[N]c2ccccc2O1)[O])[N][CH][O]" );
  CHECK( test_make_tautomer_skeleton( "CC1Cc2c(ccc(c2c3ccc(cc3)O)O)C(=C)N1", "hist_test_501944a") == "CC1C[C]2[C]([CH][CH][C]([C]2c3ccc(cc3)O)[O])[C]([N]1)[CH2]" );
  CHECK( test_make_tautomer_skeleton( "CC1CC2=C(C(=O)C=CC2C(=C)N1)c3ccc(cc3)O", "hist_test_501944b") == "CC1C[C]2[C]([CH][CH][C]([C]2c3ccc(cc3)O)[O])[C]([N]1)[CH2]" );
  CHECK( test_make_tautomer_skeleton( "CC1CC2=C(C(=O)C=CC2C(=N1)C)c3ccc(cc3)O", "hist_test_501944c") == "CC1C[C]2[C]([CH][CH][C]([C]2c3ccc(cc3)O)[O])[C]([N]1)[CH2]" );
  CHECK( test_make_tautomer_skeleton( "CC1Cc2c(ccc(c2c3ccc(cc3)O)O)C(=C)N1", "hist_test_501944d") == "CC1C[C]2[C]([CH][CH][C]([C]2c3ccc(cc3)O)[O])[C]([N]1)[CH2]" );
  CHECK( test_make_tautomer_skeleton( "CC1Cc2c(ccc(c2c3ccc(cc3)O)O)C(=N1)C", "hist_test_501944e") == "CC1C[C]2[C]([CH][CH][C]([C]2c3ccc(cc3)O)[O])[C]([N]1)[CH2]" );
  CHECK( test_make_tautomer_skeleton( "C1CNC(=N)NC1C(C(=O)N)NC=O", "hist_test_503551a") == "C1C[N][C]([N]C1C(C(=O)N)NC=O)[NH]" );
  CHECK( test_make_tautomer_skeleton( "C1CN=C(NC1C(C(=O)N)NC=O)N", "hist_test_503551b") == "C1C[N][C]([N]C1C(C(=O)N)NC=O)[NH]" );
  CHECK( test_make_tautomer_skeleton( "C1CNC(=N)NC1C(C(=O)N)NC=O", "hist_test_503551c") == "C1C[N][C]([N]C1C(C(=O)N)NC=O)[NH]" );
  CHECK( test_make_tautomer_skeleton( "C1CNC(=NC1C(C(=O)N)NC=O)N", "hist_test_503551d") == "C1C[N][C]([N]C1C(C(=O)N)NC=O)[NH]" );
  CHECK( test_make_tautomer_skeleton( "COc1cc2c(cc1OC)nc(nc2N)N3CCN(CC3)/C(=N/c4ccc(cc4)N=[N+]=[N-])/S", "hist_test_6313a") == "COc1cc2c(cc1OC)nc(nc2N)N3CCN(CC3)[C]([N]c4ccc(cc4)N=[N+]=[N-])[S]" );
  CHECK( test_make_tautomer_skeleton( "COc1cc2c(cc1OC)nc(nc2N)N3CCN(CC3)C(=S)Nc4ccc(cc4)N=[N+]=[N-]", "hist_test_6313b") == "COc1cc2c(cc1OC)nc(nc2N)N3CCN(CC3)[C]([N]c4ccc(cc4)N=[N+]=[N-])[S]" );
  CHECK( test_make_tautomer_skeleton( "COc1cc2c(cc1OC)nc(nc2N)N3CCN(CC3)C(=O)/C=C/C(=O)c4ccccc4", "hist_test_6321") == "COc1cc2c(cc1OC)nc(nc2N)N3CCN(CC3)C(=O)/C=C/C(=O)c4ccccc4" );
  CHECK( test_make_tautomer_skeleton( "c1cnccc1C(=O)NN", "hist_test_64a") == "[CH]1[CH][N][CH][CH][C]1[C]([N][NH])[O]" );
  CHECK( test_make_tautomer_skeleton( "C1=CNC=CC1C(=O)N=N", "hist_test_64b") == "[CH]1[CH][N][CH][CH][C]1[C]([N][NH])[O]" );
  CHECK( test_make_tautomer_skeleton( "C1C=NC=CC1=C(N=N)O", "hist_test_64c") == "[CH]1[CH][N][CH][CH][C]1[C]([N][NH])[O]" );
  CHECK( test_make_tautomer_skeleton( "C1C=NC=CC1C(=O)N=N", "hist_test_64d") == "[CH]1[CH][N][CH][CH][C]1[C]([N][NH])[O]" );
  CHECK( test_make_tautomer_skeleton( "c1cnccc1C(=NN)O", "hist_test_64e") == "[CH]1[CH][N][CH][CH][C]1[C]([N][NH])[O]" );
  CHECK( test_make_tautomer_skeleton( "c1cnccc1C(=O)NN", "hist_test_64f") == "[CH]1[CH][N][CH][CH][C]1[C]([N][NH])[O]" );
  CHECK( test_make_tautomer_skeleton( "c1cnccc1C(N=N)O", "hist_test_64g") == "[CH]1[CH][N][CH][CH][C]1[C]([N][NH])[O]" );
  CHECK( test_make_tautomer_skeleton( "CCC1C(Cc2cc3c(cc(=O)oc3cc2N1)C(F)(F)F)C", "hist_test_6654a") == "CCC1C(C[C]2[CH][C]3[C]([CH][C](O[C]3[CH][C]2[N]1)[O])C(F)(F)F)C" );
  CHECK( test_make_tautomer_skeleton( "CCC1C(Cc2cc3=C(CC(=O)Oc3cc2=N1)C(F)(F)F)C", "hist_test_6654b") == "CCC1C(C[C]2[CH][C]3[C]([CH][C](O[C]3[CH][C]2[N]1)[O])C(F)(F)F)C" );
  CHECK( test_make_tautomer_skeleton( "CCC1C(Cc2cc3c(cc(=O)oc3cc2N1)C(F)(F)F)C", "hist_test_6654c") == "CCC1C(C[C]2[CH][C]3[C]([CH][C](O[C]3[CH][C]2[N]1)[O])C(F)(F)F)C" );
  CHECK( test_make_tautomer_skeleton( "c1cc2c(cc1O)c-3cn[nH]c(c3n2)NN", "hist_test_8387a") == "[CH]1[CH][C]2[C]([CH][C]1[O])[C]3[CH][N][N][C]([C]3[N]2)[N][NH]" );
  CHECK( test_make_tautomer_skeleton( "c1cc2c(cc1O)c3cnnc(c3[nH]2)NN", "hist_test_8387b") == "[CH]1[CH][C]2[C]([CH][C]1[O])[C]3[CH][N][N][C]([C]3[N]2)[N][NH]" );
  CHECK( test_make_tautomer_skeleton( "c1cc-2nc-3c([nH][nH]cc3c2cc1=O)NN", "hist_test_8387c") == "[CH]1[CH][C]2[C]([CH][C]1[O])[C]3[CH][N][N][C]([C]3[N]2)[N][NH]" );
  // CHECK( test_make_tautomer_skeleton( "", "hist_test_") == "" );
}

TEST_CASE( "Extended Rule 5 Cases" , "[test_tt_tauts]" ) {
  // Based on substructure of Chembl6306
  // No tautomers - amide rule over-rides
  CHECK( test_make_tautomer_skeleton( "COC(=O)CCC=CC(=O)N1CCNCC1", "ext_rule_5_test1" ) == "COC(=O)CCC=CC(=O)N1CCNCC1");
  // extended rule 5 applies
  CHECK( test_make_tautomer_skeleton( "COC(=O)CCC=CC(=O)C1CCNCC1", "ext_rule_5_test2" ) == "CO[C]([CH][CH][CH][CH][C]([C]1CCNCC1)[O])[O]");
  // extended rule 5 applies, because cis double bond
  CHECK( test_make_tautomer_skeleton( "COC(=O)CC/C=C\\C(=O)C1CCNCC1", "ext_rule_5_test3" ) == "CO[C]([CH][CH][CH][CH][C]([C]1CCNCC1)[O])[O]");
  // extended rule 5 not applicable, because trans double bond
  CHECK( test_make_tautomer_skeleton( "COC(=O)CC/C=C/C(=O)C1CCNCC1", "ext_rule_5_test4" ) == "COC(=O)CC/C=C/[C]([C]1CCNCC1)[O]");
}

TEST_CASE( "Adding H to N already with H" , "[test_tt_tauts]" ) {
  // Chembl13554
  std::string taut_res = "CC(=O)O[CH][C]1CSC2[C@@H](C(=O)N2[C]1[C]([O])OC(c3ccccc3)c4ccccc4)[N][C]5[N][C]6[CH][CH][C]([CH][C]6[N]5)[CH][C]=[N]";
  CHECK( test_make_tautomer_skeleton( "CC(=O)OCC1=C(N2C([C@@H](C2=O)Nc3[nH]c4ccc(cc4n3)C=C=N)SC1)C(=O)OC(c5ccccc5)c6ccccc6", "test_13554a") == taut_res );
  CHECK( test_make_tautomer_skeleton( "CC(=O)OCC1=C(N2C([C@@H](C2=O)N=c3[nH]c4ccc(cc4[nH]3)CC#N)SC1)C(=O)OC(c5ccccc5)c6ccccc6", "test_13554b") == taut_res );
  CHECK( test_make_tautomer_skeleton( "CC(=O)OCC1=C(N2C([C@@H](C2=O)Nc3[nH]c4cc(ccc4n3)C=C=N)SC1)C(=O)OC(c5ccccc5)c6ccccc6", "test_13554c") == taut_res );
  CHECK( test_make_tautomer_skeleton( "CC(=O)OCC1=C(N2C([C@@H](C2=O)Nc3[nH]c4cc(ccc4n3)CC#N)SC1)C(=O)OC(c5ccccc5)c6ccccc6", "test_13554d") == taut_res );
  CHECK( test_make_tautomer_skeleton( "CC(=O)OCC1=C(N2C([C@@H](C2=O)Nc3[nH]c4ccc(cc4n3)C=C=N)SC1)C(=O)OC(c5ccccc5)c6ccccc6", "test_13554e") == taut_res );
  CHECK( test_make_tautomer_skeleton( "CC(=O)OCC1=C(N2C([C@@H](C2=O)Nc3[nH]c4ccc(cc4n3)CC#N)SC1)C(=O)OC(c5ccccc5)c6ccccc6", "test_13554f") == taut_res );
}

TEST_CASE( "Simple acids ignored", "[test_tt_tauts]" ) {
  CHECK( test_make_tautomer_skeleton( "CC1COc2c3n1cc(c(=O)c3cc(c2N4CCN(CC4)C)F)C(=O)O", "acid_test_chembl_4") == "CC1COc2c3n1cc(c(=O)c3cc(c2N4CCN(CC4)C)F)C(=O)O" );
}

TEST_CASE( "Kekule form problem", "[test_tt_tauts]") {
  // this structure failed the round-trip test because the 2 tautomers gave different
  // kekule forms for the anilino-pyridine ring, one of which failed rule 3, the other
  // didn't.
  std::string taut_res = "CN(C)c1cccc2c1cccc2S(=O)(=O)[N][C]3[CH][CH][C]([C]([N]3)[CH2])Br";
  CHECK( test_make_tautomer_skeleton( "CN(C)c1cccc2c(cccc12)S(=O)(=O)Nc3ccc(Br)c(C)n3", "kekule_6566") == taut_res );
  CHECK( test_make_tautomer_skeleton( "Cc1c(ccc(=NS(=O)(=O)c2cccc3c2cccc3N(C)C)[nH]1)Br", "kekule_6566") == taut_res );
}

TEST_CASE( "Round Trips", "[test_tt_tauts]") {
  // As above, numbers refer to the CHEMBL structure

  CHECK( test_round_trips("c12ccccc2nc(o1)C(=O)C(C)NC=O", "500474 bit"));
  CHECK( test_round_trips("NC(=O)C(NC=O)C1NC(=N)NCC1", "503551"));
  CHECK( test_round_trips("c1cc2c(cc1O)c-3cn[nH]c(c3n2)NN", "8387"));
  CHECK( test_round_trips("CN(C)c1cccc2c(cccc12)S(=O)(=O)Nc3ccc(Br)c(C)n3", "6566"));
  CHECK( test_round_trips("CC\\1=C(c2/cc\\3/c(c(/c(/[nH]3)c/c4n/c(c(\\c5c(c(c([nH]5)/cc1\\n2)C)C(=O)OC)/C(=O)C(=O)OC)/[C@H]([C@@H]4C)CCC(=O)OC)C)C=C)C", "443682"));
  CHECK( test_round_trips("CCC1Nc2cc3OC(=O)C=C(c3cc2CC1C)C(F)(F)F", "6654"));
  CHECK( test_round_trips("C1CNC(=N)NC1C(C(=O)N)NC=O", "503551"));
  // CHECK( test_round_trips("", ""));
}
