//
// File tt_tauts_batch.cc
// David Cosgrove
// AstraZeneca
// 27th July 2015
//
// This is the main line for the program tt_tauts_batch.cc.
// The program implements the tautomer enumeration algoirthm of Thalheim et al.,
// 'A Branch-and-Bound Approach for Tautomer Enumeration',
// Molecular Informatics, 2015 - At present I have the ASAP version, without a full
// reference.  DOI: 10.1002/minf.201400128.
// It's behind a paywall.
//
// The purpose of this program is for testing - it reads a molecule file, generates
// the tautomer skeleton and all tautomers of each molecule, does the same with all
// those tautomers, and checks that everything matches up.
// It requires 1 commnd-line argument, the name of the input file.  All output
// goes to stdout.

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

// in make_taut_skeleton.cc
void make_taut_skeleton_and_tauts( const string &in_smi , const string &mol_name ,
                                   string &t_skel_smi ,
                                   vector<string> &taut_smis,
                                   bool &timed_out,
                                   bool standardise_mols,
                                   float max_time,
                                   int max_tauts);

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

  cout << "tt_tauts_batch" << endl
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
  map<string,vector<string> > t_skels;

  cout << "read " << in_smiles.size() << " molecules." << endl;
  size_t last_mol = ttts.stop_mol() == -1 ? in_smiles.size() - 1 : ttts.stop_mol();
  last_mol = last_mol >= in_smiles.size() ? in_smiles.size() - 1 : last_mol;
  for( size_t i = ttts.start_mol() ; i <= last_mol ; ++i ) {
    cout << "Doing " << i << " : " << mol_names[i] << " : " << in_smiles[i] << endl;
    string t_skel_smi;
    vector<string> taut_smis;
    bool orig_timed_out = false;
    if(-1 != ttts.max_atoms()) {
      OEGraphMol mol;
      OESmilesToMol(mol, in_smiles[i]);
      if(static_cast<int>(mol.NumAtoms()) > ttts.max_atoms()) {
        cout << "Skipping " << mol.GetTitle() << " : " << in_smiles[i]
             << " which has " << mol.NumAtoms() << " atoms." << endl;
        continue;
      }
    }

    make_taut_skeleton_and_tauts( in_smiles[i] , mol_names[i] , t_skel_smi ,
                                  taut_smis , orig_timed_out,
                                  ttts.standardise_mols(), ttts.max_time(),
                                  ttts.max_tauts());
    cout << mol_names[i] << " : " << in_smiles[i] << " : " << t_skel_smi << endl;
    size_t max_tauts = taut_smis.size() > 100 ? 100 : taut_smis.size();
    for( size_t j = 0 ; j < max_tauts ; ++j ) {
      cout << taut_smis[j] << " T_" << j << endl;
      string new_t_skel_smi;
      vector<string> new_taut_smis;
      bool timed_out = false;
      make_taut_skeleton_and_tauts( taut_smis[j] , mol_names[i] , new_t_skel_smi ,
                                    new_taut_smis , timed_out ,
                                    ttts.standardise_mols(),
                                    ttts.max_time(), ttts.max_tauts());
      if( t_skel_smi != new_t_skel_smi ) {
        cout << "AWOOGA : different tautomer skeleton for " << mol_names[i]
                << " : " << j << " : " << taut_smis[j] << " : " << new_t_skel_smi
                << " : " << t_skel_smi << " in file " << ttts.in_mol_file();
        if(orig_timed_out) {
          cout << " BUT original t_skel generation timed out";
        }
        if(timed_out) {
          cout << " BUT this t_skel generation timed out";
        }
        cout << "." << endl;
      }
      if( new_taut_smis.empty() ) {
        cout << "AWOOGA : no tautomers for " << mol_names[i]
                << " : " << j << " : " << taut_smis[j] << " : "
                << " : " << t_skel_smi << " in file " << argv[1] << endl;
      }
    }
    map<string,vector<string> >::iterator p = t_skels.find( t_skel_smi );
    if( p != t_skels.end() ) {
      cout << "WAIT : t_skel for " << mol_names[i] << " (" << t_skel_smi << ") is same as for";
      BOOST_FOREACH( string tss , p->second ) {
        cout << " " << tss;
      }
      cout << " : " << t_skel_smi << endl;
      p->second.push_back( mol_names[i] );
    } else {
      t_skels.insert( make_pair( t_skel_smi , vector<string>( 1 , mol_names[i] ) ) );
    }
  }

}
