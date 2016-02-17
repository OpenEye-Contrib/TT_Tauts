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
TTTautsSettings::TTTautsSettings( int argc , char **argv ) {

  po::options_description desc( "Allowed Options" );
  build_program_options( desc );

  po::variables_map vm;
  po::store( po::parse_command_line( argc , argv , desc ) , vm );
  po::notify( vm );

  if( vm.count( "help" ) ) {
    cout << desc << endl;
    exit( 1 );
  }

  ostringstream oss;
  oss << desc;
  usage_text_ = oss.str();

}

// ***************************************************************************
void TTTautsSettings::print_usage( ostream &os ) const {

  os << usage_text_ << endl;

}

// **************************************************************************
void TTTautsSettings::build_program_options( po::options_description &desc ) {

  desc.add_options()
    ( "help" , "Produce this help text" )
      ( "in-mol-file,I" , po::value<string>( &in_mol_file_ ) ,
        "Input molecule filename." )
      ( "out-mol-file,O" , po::value<string>( &out_mol_file_ ) ,
        "Output molecule filename." );

}

