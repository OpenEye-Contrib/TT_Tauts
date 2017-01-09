//
// File TautomerGenerator.cc
// David Cosgrove
// AstraZeneca
// 16th December 2016
//

#include "TautomerGenerator.H"
#include "DACOEMolAtomIndex.H"
#include "DACOEMolBondIndex.H"

#include <oechem.h>

#include <boost/bind.hpp>
#include <boost/current_function.hpp>

using namespace boost;
using namespace std;
using namespace OEChem;
using namespace OESystem;

namespace DACLIB {
// in eponymous file
void apply_daylight_aromatic_model( OEMolBase &mol );
string create_cansmi( const OEMolBase &mol );
string create_noncansmi( const OEMolBase &mol );
}

// ****************************************************************************
TautomerGenerator::TautomerGenerator( pOEMolBase &input_mol ) :
  mol_( input_mol ){

  can_smi_ = DACLIB::create_cansmi( *mol_ );
  global_t_skel_mol_ = mol_;
  global_t_skel_smi_ = can_smi_;

}

// ****************************************************************************
TautomerGenerator::TautomerGenerator( pOEMolBase &input_mol ,
                                      vector<vector<int> > &ts_bnds_to_1 ,
                                      vector<vector<vector<int> > > &bnds_to_1 ,
                                      vector<vector<pair<unsigned int,unsigned int> > > &h_mov ,
                                      vector<vector<vector<unsigned int> > > &us_bonds ) :
mol_( input_mol ) , t_skel_bonds_to_1_( ts_bnds_to_1 ) ,
bonds_to_1_( bnds_to_1 ) , h_moves_( h_mov ) , unsat_bond_idxs_( us_bonds ) {

  can_smi_ = DACLIB::create_cansmi( *mol_ );
  generate_t_skels();

}

// ****************************************************************************
// Generate a set of tautomers, applying each mobile_h_ etc. in turn.
vector<pOEMolBase> TautomerGenerator::generate_conn_set_tauts() const {

  vector<pOEMolBase> ret_tauts;

  for( size_t i = 0 , is = h_moves_.size() ; i < is ; ++i ) {
    vector<pOEMolBase> these_tauts;
    generate_tautomers( *mol_ , h_moves_[i] ,
                        unsat_bond_idxs_[i] , bonds_to_1_[i] , these_tauts );
    ret_tauts.insert( ret_tauts.end() , these_tauts.begin() ,
                      these_tauts.end() );
  }

  return ret_tauts;

}

// ****************************************************************************
vector<string> TautomerGenerator::generate_conn_set_taut_smiles() const {

  if( conn_set_taut_smis_.empty() ) {
    vector<pOEMolBase> ret_tauts = generate_conn_set_tauts();

    conn_set_taut_smis_.reserve( ret_tauts.size() );
    for( size_t i = 0 , is = ret_tauts.size() ; i < is ; ++i ) {
      conn_set_taut_smis_.push_back( DACLIB::create_cansmi( *ret_tauts[i] ) );
    }

    sort( conn_set_taut_smis_.begin() , conn_set_taut_smis_.end() );
    conn_set_taut_smis_.erase( unique( conn_set_taut_smis_.begin() ,
                                       conn_set_taut_smis_.end() ) ,
                               conn_set_taut_smis_.end() );
  }

  return conn_set_taut_smis_;

}

// ****************************************************************************
unsigned int TautomerGenerator::num_conn_set_tauts() const {

  unsigned int ret_val = 0;

  for( size_t i = 0 , is = h_moves_.size() ; i < is ; ++i ) {
    ret_val += static_cast<unsigned int>( h_moves_[i].size() );
  }

  return ret_val;

}
// ****************************************************************************
vector<pOEMolBase> TautomerGenerator::generate_all_tautomers() const {

#ifdef NOTYET
  cout << "entering " << BOOST_CURRENT_FUNCTION << endl;
#endif

  vector<pOEMolBase> ret_tauts;
  ret_tauts.push_back( pOEMolBase( OENewMolBase( *mol_ , OEMolBaseType::OEDefault ) ) );

  unsigned int cs = 0;
  while( cs < h_moves_.size() ) {
    vector<pOEMolBase> next_tauts;
    for( size_t i = 0 , is = ret_tauts.size() ; i < is ; ++i ) {
      vector<pOEMolBase> these_tauts;
      generate_tautomers( *ret_tauts[i] , h_moves_[cs] , unsat_bond_idxs_[cs] ,
                          bonds_to_1_[cs] , these_tauts );
      next_tauts.insert( next_tauts.end() , these_tauts.begin() ,
                         these_tauts.end() );
    }
    ret_tauts = next_tauts;
    ++cs;
  }

#ifdef NOTYET
  cout << "returning " << ret_tauts.size() << " tautomers from " << BOOST_CURRENT_FUNCTION << endl;
#endif

  return ret_tauts;

}

// ****************************************************************************
vector<string> TautomerGenerator::generate_all_tautomer_smiles() const {

  vector<pOEMolBase> all_tauts = generate_all_tautomers();

  vector<string> ret_smis;
  ret_smis.reserve( all_tauts.size() );
  for( size_t i = 0 , is = all_tauts.size() ; i < is ; ++i ) {
    ret_smis.push_back( DACLIB::create_cansmi( *all_tauts[i] ) );
  }

  sort( ret_smis.begin() , ret_smis.end() );
  ret_smis.erase( unique( ret_smis.begin() , ret_smis.end() ) , ret_smis.end() );

  return ret_smis;

}

// ****************************************************************************
unsigned int TautomerGenerator::num_all_tautomers() const {

  unsigned int ret_val = 1;

  for( size_t i = 0 , is = h_moves_.size() ; i < is ; ++i ) {
    ret_val *= static_cast<unsigned int>( h_moves_[i].size() );
  }

  return ret_val;

}

// ****************************************************************************
// return accumulation of everything in bonds_to_1_
vector<vector<int> > TautomerGenerator::get_all_bonds_to_1() const {

  vector<vector<int> > ret_val;
  for( size_t i = 0 , is = bonds_to_1_.size() ; i < is ; ++i ) {
    ret_val.insert( ret_val.end() , bonds_to_1_[i].begin() ,
                    bonds_to_1_[i].end() );
  }

  return ret_val;

}

// ****************************************************************************
// return accumulation of everything in unsat_bond_idxs_
vector<vector<unsigned int> > TautomerGenerator::get_all_unsat_bond_idxs() const {

  vector<vector<unsigned int> > ret_val;
  for( size_t i = 0 , is = unsat_bond_idxs_.size() ; i < is ; ++i ) {
    ret_val.insert( ret_val.end() , unsat_bond_idxs_[i].begin() ,
                    unsat_bond_idxs_[i].end() );
  }

  return ret_val;

}

// ****************************************************************************
bool TautomerGenerator::operator==( const TautomerGenerator &rhs ) {

  if( global_t_skel_smi_ != rhs.global_t_skel_smi_ ) {
    return false;
  }
  if( num_conn_set_tauts() != rhs.num_conn_set_tauts() ) {
    return false;
  }
  if( num_all_tautomers() != rhs.num_all_tautomers() ) {
    return false;
  }
  if( generate_conn_set_taut_smiles() != rhs.generate_conn_set_taut_smiles() ) {
    return false;
  }

  return true;

}

// ****************************************************************************
bool TautomerGenerator::operator!=( const TautomerGenerator &rhs ) {

  return !( *this == rhs );

}

// ****************************************************************************
// take the TautomerGenerator old_one passed in, and take out any
// transformations in this one that produces a tautomer that's already
// in old_one
void TautomerGenerator::prune( const TautomerGenerator &old_one ) {

  // old_smis will come back sorted as part of the uniquification process
  vector<string> old_smis = old_one.generate_conn_set_taut_smiles();

  pair<unsigned int,unsigned int> bad_pair( numeric_limits<unsigned int>::max() ,
                                            numeric_limits<unsigned int>::max() );
  for( size_t i = 0 , is = h_moves_.size() ; i < is ; ++i ) {
    vector<pOEMolBase> these_tauts;
    generate_tautomers( *mol_ , h_moves_[i] , unsat_bond_idxs_[i] ,
                        bonds_to_1_[i] , these_tauts );
    for( size_t j = 0 , js = these_tauts.size() ; j < js ; ++j ) {
      string ts = DACLIB::create_cansmi( *these_tauts[j] );
      if( binary_search( old_smis.begin() , old_smis.end() , ts ) ) {
        bonds_to_1_[i][j].clear();
        h_moves_[i][j] = bad_pair;
        unsat_bond_idxs_[i][j].clear();
      } else {
        conn_set_taut_smis_.push_back( ts );
      }
    }
    bonds_to_1_[i].erase( remove_if( bonds_to_1_[i].begin() , bonds_to_1_[i].end() ,
                                    boost::bind( &vector<int>::empty , _1 ) ) ,
                          bonds_to_1_[i].end() );
    h_moves_[i].erase( remove( h_moves_[i].begin() , h_moves_[i].end() , bad_pair ) ,
                       h_moves_[i].end() );
    unsat_bond_idxs_[i].erase( remove_if( unsat_bond_idxs_[i].begin() , unsat_bond_idxs_[i].end() ,
                                          boost::bind( &vector<unsigned int>::empty , _1 ) ) ,
                               unsat_bond_idxs_[i].end() );
    if( bonds_to_1_[i].empty() ) {
      h_moves_[i].clear();
      t_skel_bonds_to_1_[i].clear();
      t_skel_mols_[i].reset();
      t_skel_smis_[i] = string( "" );
    }
  }

  h_moves_.erase( remove_if( h_moves_.begin() , h_moves_.end() ,
                             boost::bind( &vector<pair<unsigned int,unsigned int> >::empty , _1 ) ) ,
                  h_moves_.end() );
  t_skel_bonds_to_1_.erase( remove_if( t_skel_bonds_to_1_.begin() , t_skel_bonds_to_1_.end() ,
                                       boost::bind( &vector<int>::empty , _1 ) ) ,
                            t_skel_bonds_to_1_.end() );
  t_skel_mols_.erase( remove_if( t_skel_mols_.begin() , t_skel_mols_.end() ,
                                 !boost::bind( &pOEMolBase::get , _1 ) ) ,
                      t_skel_mols_.end() );
  t_skel_smis_.erase( remove_if( t_skel_smis_.begin() , t_skel_smis_.end() ,
                                 boost::bind( &string::empty , _1 ) ) ,
                      t_skel_smis_.end() );

  global_t_skel_mol_.reset();
  global_t_skel_smi_ = string( "" );
  if( !h_moves_.empty() ) {
    generate_t_skels();
    make_global_t_skel();
  }

}

// ****************************************************************************
void TautomerGenerator::generate_t_skels() {

  for( size_t i = 0 , is = h_moves_.size() ; i < is ; ++i ) {

    pOEMolBase this_t_skel_mol( OENewMolBase( *mol_ , OEMolBaseType::OEDefault ) );
    build_t_skel_mol( h_moves_[i] , t_skel_bonds_to_1_[i] , *this_t_skel_mol );
    // for the t_skel_smi, don't use a canonical SMILES because different t_skel_mol
    // might give the same canonical SMILES but with different atom numberings,
    // so when we check to see if we've seen this t_skel_mol before, the numberings
    // of atoms in the duplicate might be wrong, which will give rubbish tautomers.
    // This was shown by different input tautomers of CHEMBL18048
    string this_t_skel_smi = DACLIB::create_noncansmi( *this_t_skel_mol );

    t_skel_mols_.push_back( this_t_skel_mol );
    t_skel_smis_.push_back( this_t_skel_smi );
  }

  make_global_t_skel();

}

// ****************************************************************************
void TautomerGenerator::make_global_t_skel() {

  vector<int> glob_mob_h( DACLIB::max_atom_index( *mol_ ) , 0 );
  vector<int> glob_bnds_to_1( DACLIB::max_bond_index( *mol_ ) , 0 );

  global_t_skel_mol_ = pOEMolBase( OENewMolBase( *mol_ , OEMolBaseType::OEDefault ) );
  for( size_t i = 0 , is = h_moves_.size() ; i < is ; ++i ) {
    for( size_t j = 0 , js = h_moves_[i].size() ; j < js ; ++j ) {
      glob_mob_h[h_moves_[i][j].first] = 1;
    }
  }

  for( size_t i = 0 , is = t_skel_bonds_to_1_.size() ; i < is ; ++i ) {
    for( size_t j = 0 , js = t_skel_bonds_to_1_[i].size() ; j < js ; ++j ) {
      if( t_skel_bonds_to_1_[i][j] ) {
        glob_bnds_to_1[j] = 1;
      }
    }
  }

  build_t_skel_mol( glob_mob_h , glob_bnds_to_1 , *global_t_skel_mol_ );
  global_t_skel_smi_ = DACLIB::create_cansmi( *global_t_skel_mol_ );

}

// ****************************************************************************
void build_t_skel_mol( unsigned int h_from ,
                       const vector<int> &bonds_to_1 ,
                       OEMolBase &t_skel_mol ) {

#ifdef NOTYET
  cout << "build_t_skel_mol" << endl;
#endif

  remove_h_from_t_skel( h_from , t_skel_mol );
  set_bonds_to_1( bonds_to_1 , t_skel_mol );

  DACLIB::apply_daylight_aromatic_model( t_skel_mol );


}

// ****************************************************************************
void build_t_skel_mol( const vector<int> &mobile_h ,
                       const vector<int> &bonds_to_1 ,
                       OEMolBase &t_skel_mol ) {

#ifdef NOTYET
  cout << "build_t_skel_mol" << endl;
#endif

  for( size_t i = 0 , is = mobile_h.size() ; i < is; ++i ) {
    if( mobile_h[i] ) {
      remove_h_from_t_skel( static_cast<unsigned int>( i ) , t_skel_mol );
    }
  }
  set_bonds_to_1( bonds_to_1 , t_skel_mol );

  DACLIB::apply_daylight_aromatic_model( t_skel_mol );


}

// ****************************************************************************
void build_t_skel_mol( const vector<pair<unsigned int,unsigned> > &h_moves ,
                       const vector<int> &bonds_to_1 ,
                       OEMolBase &t_skel_mol ) {

#ifdef NOTYET
  cout << "build_t_skel_mol" << endl;
#endif

  for( size_t i = 0 , is = h_moves.size() ; i < is; ++i ) {
    remove_h_from_t_skel( h_moves[i].first , t_skel_mol );
  }
  set_bonds_to_1( bonds_to_1 , t_skel_mol );

  DACLIB::apply_daylight_aromatic_model( t_skel_mol );


}

// ****************************************************************************
void remove_h_from_t_skel( unsigned int h_from , OEMolBase &t_skel_mol ) {

  OEAtomBase *atom = t_skel_mol.GetAtom( DACLIB::HasAtomIndex( h_from ) );
  if( atom->GetImplicitHCount() ) {
    atom->SetImplicitHCount( atom->GetImplicitHCount() - 1 );
#ifdef NOTYET
    cout << "Remove H from " << h_from + 1 << endl;
#endif
      }
  atom->SetStereo( vector<OEAtomBase *>() , OEAtomStereo::Tetra ,
                   OEAtomStereo::Undefined );
}

// ****************************************************************************
void set_bonds_to_1( const vector<int> &bonds_to_1 ,
                     OEMolBase &mol ) {

  for( size_t i = 0 , is = bonds_to_1.size() ; i < is ; ++i ) {
    if( bonds_to_1[i] ) {
      OEBondBase *bond = mol.GetBond( DACLIB::HasBondIndex( static_cast<unsigned int>( i ) ) );
      if( bond ) {
        int bo = bond->GetOrder();
        if( bo > 1 ) {
          bond->SetOrder( bo - bonds_to_1[i] );
          bond->GetBgn()->SetStereo( vector<OEAtomBase *>() , OEAtomStereo::Tetra ,
                                     OEAtomStereo::Undefined );
          bond->GetEnd()->SetStereo( vector<OEAtomBase *>() , OEAtomStereo::Tetra ,
                                     OEAtomStereo::Undefined );
#ifdef NOTYET
    cout << "set bond order for " << DACLIB::atom_index( *bond->GetBgn() ) + 1
         << " to " << DACLIB::atom_index( *bond->GetEnd() ) + 1
         << " to " << bond->GetOrder() << endl;
#endif
        }
      }
    }
  }

}

// ****************************************************************************
void generate_tautomers( const OEMolBase &master_mol ,
                         const vector<pair<unsigned int,unsigned int> > &h_moves ,
                         const vector<vector<unsigned int> > &unsat_bond_idxs ,
                         const vector<vector<int> > &bonds_to_1 ,
                         vector<pOEMolBase > &tauts ) {

  for( size_t i = 0 , is = h_moves.size() ; i < is ; ++i ) {

#ifdef NOTYET
    cout << "Generating new tautomer : " << i << " of " << atoms_for_hs.size() << endl;
    cout << "Master mol : " << DACLIB::create_noncansmi( master_mol ) << endl;
    cout << "h move : " << h_moves[i].first + 1 << " -> " << h_moves[i].second + 1 << endl;
    cout << "   unsat_bond_idxs :";
    for( size_t ii = 0 , iis = unsat_bond_idxs[i].size() ; ii < iis ; ++ii ) {
      OEBondBase *b = master_mol.GetBond( DACLIB::HasBondIndex( unsat_bond_idxs[i][ii] ) );
      cout << " (" << unsat_bond_idxs[i][ii] << ") " << DACLIB::atom_index( *b->GetBgn() ) + 1 << "->"
           << DACLIB::atom_index( *b->GetEnd() ) + 1;
    }
    cout << endl;
    cout << "bonds_to_1 :";
    for( size_t ii = 0 , iis = bonds_to_1[i].size() ; ii < iis ; ++ii ) {
      if( bonds_to_1[i][ii] ) {
        OEBondBase *b = master_mol.GetBond( DACLIB::HasBondIndex( static_cast<unsigned int>( ii ) ) );
        cout << " (" << ii << ") " << DACLIB::atom_index( *b->GetBgn() ) + 1 << "->"
             << DACLIB::atom_index( *b->GetEnd() ) + 1;
      }
    }
    cout << endl;
    cout << "master_mol smi : " << DACLIB::create_noncansmi( master_mol ) << endl;
#endif

    pOEMolBase taut_mol( OENewMolBase( master_mol , OEMolBaseType::OEDefault ) );
    build_t_skel_mol( h_moves[i].first , bonds_to_1[i] , *taut_mol );
    OEAtomBase *atom = taut_mol->GetAtom( DACLIB::HasAtomIndex( h_moves[i].second ) );
    atom->SetImplicitHCount( atom->GetImplicitHCount() + 1 );

    for( size_t j = 0 , js = unsat_bond_idxs[i].size() ; j < js ; ++j ) {
      OEBondBase *b = taut_mol->GetBond( DACLIB::HasBondIndex( unsat_bond_idxs[i][j] ) );
      b->SetOrder( b->GetOrder() + 1 );
    }

    DACLIB::apply_daylight_aromatic_model( *taut_mol );

    string taut_smi( DACLIB::create_cansmi( *taut_mol ) );
    if( string::npos != taut_smi.find( "[CH" ) || string::npos != taut_smi.find( "[C]" ) ||
        string::npos != taut_smi.find( '$' ) ) {
      cout << "dodgy SMILES " << taut_smi << endl;
//      return;
      exit( 1 );
    } else {
      tauts.push_back( taut_mol );
    }

  }

}
