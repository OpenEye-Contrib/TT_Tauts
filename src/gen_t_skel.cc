//
// file gen_t_skel.cc
// David Cosgrove
// AstraZeneca
// 21st August 2015
//
// This is a command-line program that just reads a molecule file and writes the
// tautomer skeleton SMILES to another.

#include "FileExceptions.H"

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <oechem.h>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

using namespace boost;
using namespace std;
using namespace OEChem;
using namespace OESystem;

extern string BUILD_TIME; // in build_time.cc

unsigned int make_taut_skeleton( const string &in_smi , const string &mol_name ,
                                 float max_time , string &t_skel_smi );

// ****************************************************************************
void read_molecule_file( const string &filename ,
                         vector<string> &in_smiles ,
                         vector<string> &mol_names ) {

  oemolistream ims;
  if( !ims.open( filename ) ) {
    throw DACLIB::FileReadOpenError( filename.c_str() );
  }

  OEMol mol;
  while( ims >> mol ) {
    string smi;
    OECreateSmiString( smi , mol , OESMILESFlag::AtomStereo | OESMILESFlag::BondStereo | OESMILESFlag::Canonical );
    in_smiles.push_back( smi );
    mol_names.push_back( mol.GetTitle() );
    if( mol_names.back().empty() ) {
      mol_names.back() = string( "Str_" ) + boost::lexical_cast<string>( mol_names.size() );
    }
  }

}

// *****************************************************************************
int main( int argc , char **argv ) {

  cout << "gen_t_skel" << endl
       << "Built " << BUILD_TIME << " using OEToolkits version "
       << OEChem::OEChemGetRelease() << "." << endl << endl;
  if( argc < 2 ) {
    cout << "Need the name of a SMILES file for input." << endl;
    exit( 1 );
  }
  float max_time = numeric_limits<float>::max();
  if( argc > 2 ) {
    max_time = lexical_cast<float>( argv[2] );
  }

  vector<string> in_smiles , mol_names;
  read_molecule_file( string( argv[1] ) , in_smiles , mol_names );

  cout << "read " << in_smiles.size() << " molecules." << endl;
  OEStopwatch watch;
  watch.Start();

  float worst_time = 0.0F;
  string worst_name , worst_smiles;

  for( size_t i = 0 , is = in_smiles.size() ; i < is ; ++i ) {
#ifdef NOTYET
    cout << "Doing " << mol_names[i] << " : " << in_smiles[i] << endl;
#endif
    string t_skel_smi;
    float beg = watch.Elapsed();
    unsigned int num_tauts = make_taut_skeleton( in_smiles[i] , mol_names[i] ,
                                                 max_time , t_skel_smi );
    float fin = watch.Elapsed();
    float dur = fin - beg;
    if( dur > worst_time ) {
      worst_time = dur;
      worst_name = mol_names[i];
      worst_smiles = in_smiles[i];
    }
    cout << "T_Skel for " << mol_names[i] << " : " << in_smiles[i] << " : " << t_skel_smi
         << " num_tauts = " << num_tauts << " time = " << dur << endl;
    if( i && !( i % 10 ) ) {
      cout << "Avge per mol : " << watch.Elapsed() / float( i ) << " for " << i
           << " mols. Worst time for 1 molecule = " << worst_time
           << " (" << worst_name << ")" << endl;
    }
  }

  cout << "Avge per mol : " << watch.Elapsed() / float( in_smiles.size() )
       << " for " << in_smiles.size() << " mols." << endl;
  cout << "Worst time for 1 molecule : " << worst_time << " for "
       << worst_name << " : " << worst_smiles << endl;

}
