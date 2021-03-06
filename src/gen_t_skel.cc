//
// file gen_t_skel.cc
// David Cosgrove
// AstraZeneca
// 21st August 2015
//
// This is a command-line program that just reads a molecule file and writes the
// tautomer skeleton SMILES to another.

#include "FileExceptions.H"
#include "TTTautsBatchSettings.h"

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
                                 float max_time , int max_tauts,
                                 bool standardise_mol,
                                 string &t_skel_smi ,bool &timed_out );

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

  TTTautsBatchSettings ttts;
  if( argc < 2 ) {
    ttts.print_usage(cout);
    exit( 1 );
  }
  ttts.parse_options(argc, argv);

  vector<string> in_smiles , mol_names;
  read_molecule_file( ttts.in_mol_file() , in_smiles , mol_names );
  if(in_smiles.empty()) {
    cout << "No molecules read, so nothing to do." << endl;
    exit(0);
  }

  cout << "read " << in_smiles.size() << " molecules." << endl;
  OEStopwatch watch;
  watch.Start();

  float worst_time = 0.0F;
  string worst_name , worst_smiles;

  size_t last_mol = ttts.stop_mol() == -1 ? in_smiles.size() - 1 : ttts.stop_mol();
  last_mol = last_mol >= in_smiles.size() ? in_smiles.size() - 1 : last_mol;
  int num_done = 0;
  for( size_t i = ttts.start_mol() ; i <= last_mol ; ++i ) {
#ifdef NOTYET
    cout << "Doing " << i << " : " << mol_names[i] << " : " << in_smiles[i] << endl;
#endif
    string t_skel_smi;
    float beg = watch.Elapsed();
    bool timed_out = false;

    if(-1 != ttts.max_atoms()) {
      OEGraphMol mol;
      OESmilesToMol(mol, in_smiles[i]);
      if(static_cast<int>(mol.NumAtoms()) > ttts.max_atoms()) {
        cout << "Skipping " << mol.GetTitle() << " : " << in_smiles[i]
             << " which has " << mol.NumAtoms() << " atoms." << endl;
        continue;
      }
    }

    unsigned int num_tauts = make_taut_skeleton( in_smiles[i] , mol_names[i] ,
                                                 ttts.max_time(), ttts.max_tauts(),
                                                 ttts.standardise_mols(),
                                                 t_skel_smi , timed_out );
    ++num_done;
    float fin = watch.Elapsed();
    float dur = fin - beg;
    if( dur > worst_time ) {
      worst_time = dur;
      worst_name = mol_names[i];
      worst_smiles = in_smiles[i];
    }
    cout << "T_Skel for " << mol_names[i] << " : " << in_smiles[i] << " : " << t_skel_smi
         << " num_tauts = " << num_tauts << " time = " << dur;
    if( timed_out ) {
      cout << " BUT timed out";
    }
    cout << "." << endl;
    if( num_done && !( num_done % 10 ) ) {
      cout << "Avge per mol : " << watch.Elapsed() / float( num_done ) << " for " << num_done
           << " mols. Worst time for 1 molecule = " << worst_time
           << " (" << worst_name << ")" << endl;
    }
  }

  cout << "Avge per mol : " << watch.Elapsed() / float( num_done )
       << " for " << num_done << " mols." << endl;
  cout << "Worst time for 1 molecule : " << worst_time << " for "
       << worst_name << " : " << worst_smiles << endl;

}
