//
// file TTTautsBatchSettings.cpp
// David Cosgrove
// CozChemIx Limited
//

#include "TTTautsBatchSettings.h"

#include <iostream>

using namespace std;
namespace po = boost::program_options;

// ****************************************************************************
TTTautsBatchSettings::TTTautsBatchSettings() :
  TTTautsSettings() , start_mol_(0), stop_mol_(-1) {

  build_program_options();

}

// ****************************************************************************
void TTTautsBatchSettings::build_program_options() {

  desc_.add_options()
      ("start-mol", po::value<int>(&start_mol_)->default_value(0),
       "First molecule to process.")
      ("stop-mol", po::value<int>(&stop_mol_)->default_value(-1),
       "Last molecule to process. Default means go to end.")
      ("max-atoms", po::value<int>(&max_atoms_)->default_value(-1),
       "Maximum number of heavy atoms for molecule. Default means no maximum.");

}
