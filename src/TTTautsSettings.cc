//
// file TTTautsSettings.cc
// David Cosgrove
// AstraZeneca
// 8th June 2015
//

#include "TTTautsSettings.H"

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <iostream>

using namespace std;
namespace po = boost::program_options;

// *****************************************************************************
TTTautsSettings::TTTautsSettings() :
  desc_( po::options_description("Allowed Options") ),
  dont_standardise_mols_(true), max_tauts_(2500) {

  build_program_options();

}

// ***************************************************************************
void TTTautsSettings::parse_options(int argc, char **argv) {

  po::variables_map vm;
  po::store( po::parse_command_line( argc , argv , desc_ ) , vm );
  po::notify( vm );

  if( vm.count( "help" ) ) {
    cout << desc_ << endl;
    exit( 1 );
  }

}

// ***************************************************************************
void TTTautsSettings::print_usage( ostream &os ) const {

  if(usage_text_.empty()) {
    ostringstream oss;
    oss << desc_;
    usage_text_ = oss.str();
  }
  os << usage_text_ << endl;

}

// **************************************************************************
void TTTautsSettings::build_program_options() {

  desc_.add_options()
    ( "help" , "Produce this help text" )
      ( "in-mol-file,I" , po::value<string>( &in_mol_file_ ) ,
        "Input molecule filename." )
      ( "dont-standardise-input-molecules",
        po::value<bool>( &dont_standardise_mols_ )->default_value(false, "false")->zero_tokens() ,
        "Turn off application of standardisation SMIRKS to molecules on input.")
      ( "max-time", po::value<float>( &max_time_ )->default_value(numeric_limits<float>::max(), "No limit"),
        "Maximum time for each t_skel generation, in CPU seconds.")
      ( "max-tauts", po::value<int>( &max_tauts_ )->default_value(2500, "2500"),
        "Maximum number of tautomers for each t_skel generation.");

}
