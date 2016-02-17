//
// file make_taut_skeleton.cc
// David Cosgrove
// AstraZeneca
// 15th June 2015
//
// This file contains functions to take a SMILES string and return the SMILES
// of the corresponding tautomer skeleton, as defined by Thalheim et al.

#include "stddefs.H"
#include "Combinator.H"
#include "DACOEMolAtomIndex.H"

#include <iostream>
#include <limits>
#include <list>
#include <numeric>
#include <string>
#include <vector>

#include <oechem.h>

#include <boost/bind.hpp>

using namespace boost;
using namespace std;
using namespace OEChem;
using namespace OESystem;

namespace DACLIB {
// in eponymous file
void apply_daylight_aromatic_model( OEMolBase &mol );
}

// ****************************************************************************
bool atom_has_multiple_bond( OEAtomBase *atom ) {

  OEIter<OEBondBase> bond = atom->GetBonds( OENot<OEBondBase>( OEHasOrder( 1 ) ) );
  if( bond ) {
    return true;
  } else {
    return false;
  }

}

// ****************************************************************************
int count_multiple_bonds( OEAtomBase *atom ) {

  OEIter<OEBondBase> bond = atom->GetBonds( OENot<OEBondBase>( OEHasOrder( 1 ) ) );
  int num_mb = 0;
  for( ; bond ; ++bond ) {
    ++num_mb;
  }

  return num_mb;

}

// ****************************************************************************
void extend_bond_path( unsigned int finish_atom ,
                       vector<OEAtomBase *> &hads ,
                       vector<int> &is_had , OEMolBase &mol ,
                       vector<OEAtomBase *> &curr_path ,
                       vector<vector<OEAtomBase *> > &bond_paths ) {

#ifdef NOTYET
  cout << "Extending bond path :";
  for( unsigned int i = 0 , is = curr_path.size() ; i < is ; ++i ) {
    cout << " " << DACLIB::atom_index( *curr_path[i] ) + 1;
  }
  cout << endl;
#endif

  for( OEIter<OEAtomBase> nb = curr_path.back()->GetAtoms( OEIsHeavy() ) ; nb ; ++nb ) {
#ifdef NOTYET
    cout << "from " << DACLIB::atom_index( *curr_path.back() ) + 1 << " to "
         << DACLIB::atom_index( *nb ) + 1 << endl;
#endif
    if( curr_path.end() != find( curr_path.begin() , curr_path.end() , nb ) ) {
#ifdef NOTYET
      cout << "already in path" << endl;
#endif
      continue;
    }
    // if it's an atom like an ether oxygen or an N with 3 heavy atoms
    // attached, it can't be in a tautomeric bond path
    bool accept_nb = false;
    if( !nb->GetImplicitHCount() && !atom_has_multiple_bond( nb ) ) {
      continue;
    }
    if( 1 == curr_path.size() ) {
      // just take this atom as it's the first extension
      accept_nb = true;
    } else {
      // it's a bit more complicated, as whether nb is accepted depends on the
      // previous bond in the path
      OEAtomBase *last_at = curr_path.back();
      OEAtomBase *prev_at = curr_path[curr_path.size() - 2];
#ifdef NOTYET
      cout << "  last_at idx = " << DACLIB::atom_index( *last_at ) + 1
           << " and prev_at = " << DACLIB::atom_index( *prev_at ) + 1 << endl;
      cout << "  last_at at num = " << last_at->GetAtomicNum() << endl;
      cout << "  bond last_at -> prev_at : " << mol.GetBond( last_at , prev_at )->GetOrder() << endl;
      cout << "  multiple bond at prev_at : " << atom_has_multiple_bond( prev_at ) << endl;
      cout << "  multiple bond at nb : " << atom_has_multiple_bond( nb ) << endl;
      cout << "prev, last, next had status : " << is_had[DACLIB::atom_index( *prev_at )]
           << " , " << is_had[DACLIB::atom_index( *last_at )] << " , " << is_had[DACLIB::atom_index( *nb )] << endl;
#endif
      if( 1 == mol.GetBond( last_at , prev_at )->GetOrder() &&
          1 < mol.GetBond( last_at , nb )->GetOrder() ) {
        // last one single, this one > single
#ifdef NOTYET
        cout << "Accepted reason 1" << endl;
#endif
        accept_nb = true;
      }
      if( !accept_nb && 1 < mol.GetBond( last_at , prev_at )->GetOrder() ) {
        // last one > single, this one can be anything (for allenes, for example)#ifdef NOTYET
#ifdef NOTYET
        cout << "Accepted reason 2" << endl;
#endif
        accept_nb = true;
      }
      if( !accept_nb && 1 == mol.GetBond( last_at , prev_at )->GetOrder() &&
          ( atom_has_multiple_bond( prev_at ) || is_had[DACLIB::atom_index( nb )] ) &&
          ( OEElemNo::C != prev_at->GetAtomicNum() ||
            OEElemNo::C != last_at->GetAtomicNum() ||
            OEElemNo::C != nb->GetAtomicNum() ) ) {
        // Deal with the case where the last atom in the current path is C,
        // for example, and we're adding X to make X=NCX, the path NCX needs
        // to be all HAD (especially the first X) because a valid tautomer
        // might be X-N=CX from which one can make X-N-C=X by another 1,3
        // shift. The rule is that if the last atom in the path is C, and the
        // bond to the previous one is single, and the previous atom has a
        // multiple bond, accept this atom. Also apply this if prev_at is
        // already established as a HAD, as it could have a double bond going
        // to it in a previous round of tautomerism.
        // BUT don't accept 3 C atoms in a row under this rule (sort of
        // Rule 5)
        // This is an extension to TT rules
#ifdef NOTYET
        cout << "Accepted reason 3" << endl;
#endif
        accept_nb = true;
      }
      if( !accept_nb && 1 == mol.GetBond( last_at , prev_at )->GetOrder() &&
          1 == mol.GetBond( last_at , nb )->GetOrder() &&
          ( atom_has_multiple_bond( nb ) || is_had[DACLIB::atom_index( nb )] ) &&
          ( OEElemNo::C != prev_at->GetAtomicNum() ||
            OEElemNo::C != last_at->GetAtomicNum() ||
            OEElemNo::C != nb->GetAtomicNum() ) ) {
        // and the reverse of the above, building the path from the other end
        // This is an extension to TT rules
#ifdef NOTYET
        cout << "Accepted reason 4" << endl;
#endif
        accept_nb = true;
      }
      if( !accept_nb ) {
        // the final option is that all the atoms bar the first are sp2 hybridised
        // and nb is too, or it's the finish atom
        unsigned int sp2_cnt = 1; // so we don't have to subtract 1 from curr_path.size()
        for( int i = 1 , is = curr_path.size() ; i < is ; ++i ) {
          if( atom_has_multiple_bond( curr_path[i] ) ) {
            ++sp2_cnt;
          } else {
            break;
          }
        }
        if( sp2_cnt == curr_path.size() ) {
          if( atom_has_multiple_bond( nb ) || nb == hads[finish_atom] ) {
#ifdef NOTYET
            cout << "Accepted reason 5" << endl;
#endif
            accept_nb = true;
          }
        }
      }
    }

    if( accept_nb ) {
      curr_path.push_back( nb );
      if( nb == hads[finish_atom] ) {
        // to be acceptable, the bond path must be at least 3 atoms (2 bonds)
        // and an even number of bonds
        unsigned int num_bonds = curr_path.size() - 1;
        if( num_bonds >= 2 && !( num_bonds % 2 ) ) {
          bond_paths.push_back( curr_path );
        }
      } else {
        extend_bond_path( finish_atom , hads , is_had , mol , curr_path , bond_paths );
      }
      curr_path.pop_back();
    }
  }

}

// ****************************************************************************
// find all bond paths between the start and finish atoms.  A bond path is
// a conjugated path, i.e.
// alternating single and double bonds OR
// with the exception of the donor, conjugation of sp2 hybridised atoms OR
// an n >= 3 ordered bond instead of a double bond
void find_bond_paths( vector<OEAtomBase *> &hads ,
                      vector<int> &is_had ,
                      unsigned int start_atom ,
                      unsigned int finish_atom , OEMolBase &mol ,
                      vector<vector<OEAtomBase *> > &bond_paths ) {

#ifdef NOTYET
  cout << "find_bond_paths from " << DACLIB::atom_index( *hads[start_atom] ) + 1 << " to "
       << DACLIB::atom_index( *hads[finish_atom] ) + 1 << endl;
#endif

  bond_paths.clear();
  vector<OEAtomBase *> curr_path = ( vector<OEAtomBase *>( 1 , hads[start_atom] ) );

  extend_bond_path( finish_atom , hads , is_had , mol , curr_path , bond_paths );

#ifdef NOTYET
  cout << "XXXX Bond paths from " << DACLIB::atom_index( *hads[start_atom] ) + 1 << " to "
       << DACLIB::atom_index( *hads[finish_atom] ) + 1 << endl;
  for( unsigned int j = 0 , js = bond_paths.size() ; j < js ; ++j ) {
    cout << "Path " << j << " :: ";
    for( unsigned int i = 0 , is = bond_paths[j].size() ; i < is ; ++i ) {
      cout << " " << DACLIB::atom_index( *bond_paths[j][i] ) + 1;
    }
    cout << endl;
  }
#endif

}

// ****************************************************************************
// rule 3 : a HAD without an H-atom (acceptor) must have a bond-path to
// another HAD (donor) that can release an H-atom
void apply_rule_3( vector<vector<OEAtomBase *> > &bond_paths ,
                   vector<OEAtomBase *> &hads ,
                   vector<int> &is_had ) {

  for( int i = 0 , is = hads.size() ; i < is ; ++i ) {
    if( !hads[i]->GetTotalHCount() ) {
      bool ok_to_keep = false;
      for( int j = 0 , js = bond_paths.size() ; j < js ; ++j ) {
        if( bond_paths[j].front() == hads[i] ) {
          if( bond_paths[j].back()->GetTotalHCount() ) {
            ok_to_keep = true;
            break;
          }
        } else if( bond_paths[j].back() == hads[i] ) {
          if( bond_paths[j].front()->GetTotalHCount() ) {
            ok_to_keep = true;
            break;
          }
        }
      }
      if( !ok_to_keep ) {
        cout << "Removing " << DACLIB::atom_index( *hads[i] ) + 1 << " from HAD list due to rule 3" << endl;
#ifdef NOTYET
#endif
        is_had[DACLIB::atom_index( *hads[i] )] = 0;
        hads[i] = 0;
      }
    }
  }

  hads.erase( remove( hads.begin() , hads.end() ,
                      static_cast<OEAtomBase *>( 0 ) ) , hads.end() );

}

// ****************************************************************************
// rule 5 : a carbon atom HAD must have a bondpath of length 2 to at least one
// heteroatom HAD
// OR (this is an extension to the TT rules) if both bonds are single, and
// there's a heteroatom in the path and the first or last atom has a multiple
// bond. In this case, there's a possible tautomer that has a multiple bond
// in the path, so rule 5 is resurrected.
// CHEMBL19253 gives rise to a more general rule - only apply to paths of
// 3 atoms. Assume longer ones are ok.
void apply_rule_5( vector<vector<OEAtomBase *> > &bond_paths ,
                   vector<OEAtomBase *> &hads ,
                   vector<int> &is_had ) {

  for( int i = 0 , is = hads.size() ; i < is ; ++i ) {
    if( OEElemNo::C == hads[i]->GetAtomicNum() ) {
      bool ok_to_keep = false;
      for( int j = 0 , js = bond_paths.size() ; j < js ; ++j ) {
        if( 3 != bond_paths[j].size() ) {
          // 3 atoms === 2 bonds
          continue;
        }
        if( hads[i] != bond_paths[j].front() ) {
          continue;
        }
        cout << "Rule 5 on Path " << j << " :: HAD " << DACLIB::atom_index( *hads[i] ) + 1
             << " :: ";
        for( unsigned int ii = 0 , iis = bond_paths[j].size() ; ii < iis ; ++ii ) {
          cout << " " << DACLIB::atom_index( *bond_paths[j][ii] ) + 1;
        }
        cout << endl;
#ifdef NOTYET
#endif
        // this is original rule 5
        if( bond_paths[j][0] == hads[i] ) {
          if( OEElemNo::C != bond_paths[j][2]->GetAtomicNum() ) {
            ok_to_keep = true;
            cout << "ok_to_keep 1" << endl;
            break;
          }
        } else if( bond_paths[j][2] == hads[i] ) {
          if( OEElemNo::C != bond_paths[j][0]->GetAtomicNum() ) {
            ok_to_keep = true;
            cout << "ok_to_keep 2" << endl;
            break;
          }
        }
#ifdef NOTYET
        // this is the extension
        if( hads[i] == bond_paths[j][0] || hads[i] == bond_paths[j][2] ) {
          if( ( OEElemNo::C != bond_paths[j][0]->GetAtomicNum() ||
                OEElemNo::C != bond_paths[j][1]->GetAtomicNum() ||
                OEElemNo::C != bond_paths[j][2]->GetAtomicNum() ) &&
              ( atom_has_multiple_bond( bond_paths[j][0] ) ||
                atom_has_multiple_bond( bond_paths[j][2] ) ) ) {
            ok_to_keep = true;
            cout << "ok_to_keep 3" << endl;
            break;
          }
        }
#endif
      }
      if( !ok_to_keep ) {
        cout << "Removing " << DACLIB::atom_index( *hads[i] ) + 1 << " from HAD list due to rule 5" << endl;
#ifdef NOTYET
#endif
        is_had[DACLIB::atom_index( *hads[i] )] = 0;
        // also remove any other paths with this had at the start
        // as they can't be valid, either
        for( int j = 0 , js = bond_paths.size() ; j < js ; ++j ) {
          if( hads[i] == bond_paths[j].front() ) {
            bond_paths[j].clear();
          }
        }
        hads[i] = 0;
      }
    }
  }

  hads.erase( remove( hads.begin() , hads.end() ,
                      static_cast<OEAtomBase *>( 0 ) ) , hads.end() );
  bond_paths.erase( remove_if( bond_paths.begin() , bond_paths.end() ,
                               bind( &vector<OEAtomBase *>::empty , _1 ) ) ,
                    bond_paths.end() );

}

// ****************************************************************************
// rule 6 : a heteroatom HAD must have a bondpath of even length to at least
// one heteroatom HAD of a bondpath of length 2 to at least one carbon atom HAD
void apply_rule_6( vector<vector<OEAtomBase *> > &bond_paths ,
                   vector<OEAtomBase *> &hads ,
                   vector<int> &is_had ) {

  for( int i = 0 , is = hads.size() ; i < is ; ++i ) {
    if( OEElemNo::C != hads[i]->GetAtomicNum() ) {
      bool ok_to_keep = false;
      for( int j = 0 , js = bond_paths.size() ; j < js ; ++j ) {
        // bondpaths are already forced to be of even length, so no need to
        // check that
        if( bond_paths[j].front() == hads[i] ) {
          if( OEElemNo::C != bond_paths[j].back()->GetAtomicNum() ) {
            ok_to_keep = true;
            break;
          } else if( OEElemNo::C == bond_paths[j].back()->GetAtomicNum() &&
                     3 == bond_paths[j].size() ) {
            ok_to_keep = true;
            break;
          }
        } else if( bond_paths[j].back() == hads[i] ) {
          if( OEElemNo::C != bond_paths[j].front()->GetAtomicNum() ) {
            ok_to_keep = true;
            break;
          } else if( OEElemNo::C == bond_paths[j].front()->GetAtomicNum() &&
                     3 == bond_paths[j].size() ) {
            ok_to_keep = true;
            break;
          }
        }
      }
      if( !ok_to_keep ) {
        cout << "Removing " << DACLIB::atom_index( *hads[i] ) + 1 << " from HAD list due to rule 6" << endl;
#ifdef NOTYET
#endif
        is_had[DACLIB::atom_index( *hads[i] )] = 0;
        hads[i] = 0;
      }
    }
  }

  hads.erase( remove( hads.begin() , hads.end() ,
                      static_cast<OEAtomBase *>( 0 ) ) , hads.end() );

}

// ****************************************************************************
void remove_hads_not_in_paths( vector<vector<OEAtomBase *> > &bond_paths ,
                               std::vector<OEAtomBase *> &hads ,
                               vector<int> &is_had ) {

#ifdef NOTYET
  cout << "remove_hads_not_in_paths" << endl;
  for( unsigned int j = 0 , js = bond_paths.size() ; j < js ; ++j ) {
    cout << "Path " << j << " :: ";
    for( unsigned int i = 0 , is = bond_paths[j].size() ; i < is ; ++i ) {
      cout << " " << DACLIB::atom_index( *bond_paths[j][i] ) + 1;
    }
    cout << endl;
  }
  cout << "is_had : ";
  copy( is_had.begin() , is_had.end() , intOut );
  cout << endl;
  cout << "hads :";
  for( int i = 0 , is = hads.size() ; i < is ; ++i ) {
    cout << " " << DACLIB::atom_index( *hads[i] ) + 1;
  }
  cout << endl;
#endif

  // remove any hads that are no longer in a path
  vector<int> new_is_had = vector<int>( is_had.size() , 0 );
  for( int i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    for( int j = 0 , js = bond_paths[i].size() ; j < js ; ++j ) {
      if( is_had[DACLIB::atom_index( *bond_paths[i][j] )] ) {
        new_is_had[DACLIB::atom_index( *bond_paths[i][j] )] = 1;
      }
    }
  }
  is_had = new_is_had;

  for( int i = 0 , is = hads.size() ; i < is ; ++i ) {
    if( !is_had[DACLIB::atom_index( *hads[i] )] ) {
      hads[i] = 0;
    }
  }

  hads.erase( remove( hads.begin() , hads.end() ,
                      static_cast<OEAtomBase *>( 0 ) ) , hads.end() );

}

// ****************************************************************************
// remove paths that are fully unsaturated and not part of a wider network
// of paths i.e. all atoms are only in that path.  The actual test is that only
// 1 atom isn't also in another path.
// e.g. the CCO chain in C1NC=C(C(=O)N1COCCO)F CHEMBL34961 needs to be removed
// but the C number 11 of CN(C)S(=O)(=CC=NC(=O)C(CCCl)NO)O T_50 (a tautomer of
// CHEMBL31034) stays in.
void apply_saturated_rule( vector<vector<OEAtomBase *> > &bond_paths ,
                           vector<OEAtomBase *> &hads ,
                           vector<int> &is_had ) {

#ifdef NOTYET
  cout << "apply_saturated_rule" << endl;
#endif
  vector<int> in_path_cnt( is_had.size() , 0 );
  for( int i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    for( int j = 0 , js = bond_paths[i].size() ; j < js ; ++j ) {
      ++in_path_cnt[DACLIB::atom_index( *bond_paths[i][j] )];
    }
  }

  for( int i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    bool has_unsatd = false;
    int num_in_1_path = 0;
    for( int j = 0 , js = bond_paths[i].size() ; j < js ; ++j ) {
      if( atom_has_multiple_bond( bond_paths[i][j] ) ) {
        has_unsatd = true;
        break;
      }
      if( 1 == in_path_cnt[DACLIB::atom_index( *bond_paths[i][j] )] ) {
        ++num_in_1_path;
      }
    }
#ifdef NOTYET
    cout << "Path " << i << " :: ";
    for( unsigned int j = 0 , js = bond_paths[i].size() ; j < js ; ++j ) {
      cout << " " << DACLIB::atom_index( *bond_paths[i][j] ) + 1;
    }
    cout << " : " << has_unsatd << " : " << num_in_1_path << endl;
#endif
    if( !has_unsatd && num_in_1_path > 1 ) {
#ifdef NOTYET
      cout << "removing path ";
      for( unsigned int j = 0 , js = bond_paths[i].size() ; j < js ; ++j ) {
        cout << " " << DACLIB::atom_index( *bond_paths[i][j] ) + 1;
      }
      cout << " due to all unsatd rule." << endl;
#endif
      bond_paths[i].clear();
    }
  }

  bond_paths.erase( remove_if( bond_paths.begin() , bond_paths.end() ,
                               bind( &vector<OEAtomBase *>::empty , _1 ) ) ,
                    bond_paths.end() );

  remove_hads_not_in_paths( bond_paths , hads , is_had );

}

// ****************************************************************************
// another one of mine, due to over-zealous inclusion prior to this - we
// don't want to be shuffling double bonds between C atoms, so remove all
// carbon-only paths
// Example of where this matters is a bit of CHEMBL583299 :
// OCC(C)=C(O)C(C)=CO where carbon 7 needs to be removed as a HAD.
void apply_carbon_only_rule( vector<vector<OEAtomBase *> > &bond_paths ,
                             vector<OEAtomBase *> &hads ,
                             vector<int> &is_had ) {

  for( int i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    bool all_c = true;
    for( int j = 0 , js = bond_paths[i].size() ; j < js ; ++j ) {
      if( OEElemNo::C != bond_paths[i][j]->GetAtomicNum() ) {
        all_c = false;
        break;
      }
    }
    if( all_c ) {
#ifdef NOTYET
      cout << "removing path ";
      for( unsigned int j = 0 , js = bond_paths[i].size() ; j < js ; ++j ) {
        cout << " " << DACLIB::atom_index( *bond_paths[i][j] ) + 1;
      }
      cout << " due to all carbon rule." << endl;
#endif
      bond_paths[i].clear();
    }
  }

  bond_paths.erase( remove_if( bond_paths.begin() , bond_paths.end() ,
                               bind( &vector<OEAtomBase *>::empty , _1 ) ) ,
                    bond_paths.end() );

  remove_hads_not_in_paths( bond_paths , hads , is_had );

}

// ****************************************************************************
// return true is had_atom is bonded to only 1 other had
bool is_conn_to_1_had( OEAtomBase *had_atom ,
                       vector<int> &is_had ) {

  int nc = 0;
  for( OEIter<OEAtomBase> nb = had_atom->GetAtoms() ; nb ; ++nb ) {
    if( is_had[DACLIB::atom_index( *nb )] ) {
      ++nc;
    }
  }

  return( 1 == nc );

}

// ****************************************************************************
void build_initial_had_list( OEMolBase &mol , vector<OEAtomBase *> &hads ,
                             vector<int> &is_had ) {

  is_had = vector<int>( DACLIB::max_atom_index( mol ) , 0 );
  for( OEIter<OEAtomBase> atom = mol.GetAtoms( OEIsHeavy() ) ; atom ; ++atom ) {
    // rule 2
    if( !atom->GetTotalHCount() && !atom_has_multiple_bond( atom ) ) {
      continue;
    }
    // rule 4
    if( OEElemNo::C == atom->GetAtomicNum() && atom->IsAromatic() ) {
      continue;
    }
    cout << DACLIB::atom_index( atom ) + 1 << " is possible HAD" << endl;
#ifdef NOTYET
#endif
    hads.push_back( atom );
    is_had[DACLIB::atom_index( atom )] = 1;
  }

}

// ****************************************************************************
// HAD is an H-atom acceptor atom or H-atom donor atom. Not to be confused with
// hydrogen-bond donors and acceptors. See Thalheim et al p4.
void find_hads( OEMolBase &mol , vector<vector<OEAtomBase *> > &bond_paths ,
                vector<OEAtomBase *> &hads , vector<int> &is_had ) {

  build_initial_had_list( mol , hads , is_had );

  // Now check bond-paths between hads, removing any hads that no longer fit
  // the definition
  for( unsigned int i = 0 ; i < hads.size() ; ++i ) {
    unsigned int num_paths_from_i = 0;
    for( unsigned int j = 0 ; j < hads.size() ; ++j ) {
      if( i == j || !hads[j] ) {
        continue;
      }
      vector<vector<OEAtomBase *> > these_paths;
      find_bond_paths( hads , is_had , i , j , mol , these_paths );
      num_paths_from_i += these_paths.size();
      bond_paths.insert( bond_paths.end() , these_paths.begin() , these_paths.end() );
    }
    if( !num_paths_from_i ) {
#ifdef NOTYET
      cout << "AHAH - no paths from " << DACLIB::atom_index( *hads[i] ) + 1 << "  can it be a had?" << endl;
#endif
      is_had[DACLIB::atom_index( *hads[i] )] = 0;
      hads[i] = 0;
    }
  }

  hads.erase( remove( hads.begin() , hads.end() , static_cast<OEAtomBase *>( 0 ) ) ,
              hads.end() );
#ifdef NOTYET
  cout << "All bond paths" << endl;
  for( unsigned int j = 0 , js = bond_paths.size() ; j < js ; ++j ) {
    cout << "Path " << j << " :: ";
    for( unsigned int i = 0 , is = bond_paths[j].size() ; i < is ; ++i ) {
      cout << " " << DACLIB::atom_index( *bond_paths[j][i] ) + 1;
    }
    cout << endl;
  }
#endif

  // remove atoms that don't adhere to remainder of rules
  apply_rule_3( bond_paths , hads , is_had );
  apply_rule_5( bond_paths , hads , is_had );
  apply_rule_6( bond_paths , hads , is_had );

  // my own rule - if a path has all saturated atoms, it's not a valid
  // path.  This has snuck in somewhere as a possibility with some of
  // the exceptions I put in, e.g. the CCO chain in
  // C1NC=C(C(=O)N1COCCO)F CHEMBL34961
  apply_saturated_rule( bond_paths , hads , is_had );

  // another one of mine, due to over-zealous inclusion prior to this - we
  // don't want to be shuffling double bonds between C atoms, so remove all
  // carbon-only paths
  // Example of where this matters is a bit of CHEMBL583299 :
  // OCC(C)=C(O)C(C)=CO where carbon 7 needs to be removed as a HAD.
  apply_carbon_only_rule( bond_paths , hads , is_had );

}

// ****************************************************************************
void remove_mobile_h_atoms( OEGraphMol &mol , vector<OEAtomBase *> &hads ,
                            vector<int> &abes , vector<int> &mobile_h ) {

  cout << "remove_mobile_h_atoms" << endl;
#ifdef NOTYET
#endif
  for( int i = 0 , is = hads.size() ; i < is ; ++i ) {
    // if it's a carbon atom and it has a non-single bond, it's an H atom
    // acceptor but not a donor no matter how many H atoms it holds
    if( OEElemNo::C == hads[i]->GetAtomicNum() &&
        atom_has_multiple_bond( hads[i] ) ) {
      continue;
    }
    // likewise, if it's any atom and only connected to 1 other heavy atom,
    // and has a multiple bond, don't remove an H atom, because otherwise
    // when the multiple bond is removed we'll have taken out at least 2 which
    // is too many.
    // This test may fail if there are explicit H atoms as well.
    if( 1 == hads[i]->GetHvyDegree() && atom_has_multiple_bond( hads[i] ) ) {
#ifdef NOTYET
      OEIter<OEAtomBase> nb = hads[i]->GetAtoms(); // we know there's only 1
      if( count_multiple_bonds( nb ) < 2 ) {
        continue;
      }
#endif
    }

    if( hads[i]->GetTotalHCount() ) {
#ifdef NOTYET
      cout << "Removing H from " << DACLIB::atom_index( *hads[i] ) + 1 << endl;
#endif
      mobile_h[DACLIB::atom_index( *hads[i] )] = 1;
      if( hads[i]->GetImplicitHCount() ) {
        hads[i]->SetImplicitHCount( hads[i]->GetImplicitHCount() - 1 );
      } else if( hads[i]->GetExplicitHCount() ) {
        OEIter<OEAtomBase> hat = hads[i]->GetAtoms( OEHasAtomicNum( OEElemNo::H ) );
        if( hat ) {
          mol.DeleteAtom( hat );
        }
      }
      ++abes[DACLIB::atom_index( *hads[i] )];
      // deleting the H will have removed any stereochem from the
      // atom, but that needs to be done explicitly for tetrahedral, but not,
      // experiments suggest, for E/Z.
      if( hads[i]->HasStereoSpecified( OEAtomStereo::Tetra ) ) {
        hads[i]->SetStereo( vector<OEAtomBase *>() , OEAtomStereo::Tetra , OEAtomStereo::Undefined );
      }
    }
  }

}

// ****************************************************************************
// make sure all paths have at least 1 abe on each atom. Remove paths
// that don't, adjusting hads if they no longer occur in any paths.
void check_abe_counts( vector<vector<OEAtomBase *> > &bond_paths ,
                       vector<int> &abes ,
                       vector<OEAtomBase *> &hads ,
                       vector<int> &is_had ) {

  for( int i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    bool no_abe = false;
    for( int j = 0 , js = bond_paths[i].size() ; j < js ; ++j ) {
      if( !abes[DACLIB::atom_index( *bond_paths[i][j] )] ) {
        no_abe = true;
        cout << "Path " << j << " :: ";
        for( unsigned int ii = 0 , iis = bond_paths[i].size() ; ii < iis ; ++ii ) {
          cout << " " << DACLIB::atom_index( *bond_paths[i][ii] ) + 1;
        }
        cout << " no abe for " << DACLIB::atom_index( *bond_paths[i][j] ) + 1 << endl;
#ifdef NOTYET
#endif
        break;
      }
    }
    if( no_abe ) {
      bond_paths[i].clear();
    }
  }

  bond_paths.erase( remove_if( bond_paths.begin() , bond_paths.end() ,
                               bind( &vector<OEAtomBase *>::empty , _1 ) ) ,
                    bond_paths.end() );

  remove_hads_not_in_paths( bond_paths , hads , is_had );

}

// ****************************************************************************
void set_bondpath_bonds_to_1( OEGraphMol &mol ,
                              vector<vector<OEAtomBase *> > &bond_paths ,
                              vector<OEAtomBase *> &hads ,
                              vector<int> &is_had ,
                              vector<int> &abes ) {

  for( int i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    for( int j = 0 , js = bond_paths[i].size() - 1 ; j < js ; ++j ) {
      OEBondBase *bond = mol.GetBond( bond_paths[i][j] , bond_paths[i][j+1] );
      if( bond ) {
        int bo = bond->GetOrder();
        if( bo > 1 ) {
          abes[DACLIB::atom_index( *bond->GetBgn() )] += bo - 1;
          abes[DACLIB::atom_index( *bond->GetEnd() )] += bo - 1;
          bond->SetOrder( 1 );
        }
      }
#ifdef NOTYET
      // there are also so-called Tautomer Shift relevant Atoms (TSAs) that
      // aren't HADs but are also in the bond-path. They need to have an abe
      // assigned, by removing an H atom if appropriate
      if( !is_had[DACLIB::atom_index( *bond_paths[i][j] )] &&
          !atom_has_multiple_bond( bond_paths[i][j] ) &&
          bond_paths[i][j]->GetImplicitHCount() &&
          !abes[DACLIB::atom_index( *bond_paths[i][j] )] ) {
        bond_paths[i][j]->SetImplicitHCount( bond_paths[i][j]->GetImplicitHCount() - 1 );
        ++abes[DACLIB::atom_index( *bond_paths[i][j] )];
      }
#endif
    }
  }

}

// ****************************************************************************
int sum_had_h_counts( vector<OEAtomBase *> &hads ) {

  int num_h = 0;
  for( int i = 0 , is = hads.size() ; i < is ; ++i ) {
    num_h += hads[i]->GetImplicitHCount();
  }

  return num_h;

}

// ****************************************************************************
// an atom with a bonding electron (abes is non-zero)
// that doesn't have neighbour with a bonding electron is isolated, and needs
// to be removed from had list
void remove_isolated_atoms( vector<int> &abes ,
                            OEGraphMol &mol , vector<OEAtomBase *> &hads ,
                            vector<int> &is_had ,
                            vector<int> &mobile_h , int &hcount ) {

  for( unsigned int i = 0 , is = abes.size() ; i < is ; ++i ) {
    // if this atom has no bonding electrons, it doesn't matter
    if( !abes[i] ) {
      continue;
    }
    cout << "looking at atom " << i + 1 << " abe = " << abes[i] << endl;
#ifdef NOTYET
#endif
    bool abe_nb = false;
    OEAtomBase *atom = mol.GetAtom( DACLIB::HasAtomIndex( i ) );
    for( OEIter<OEAtomBase> nb = atom->GetAtoms() ; nb ; ++nb ) {
      if( abes[DACLIB::atom_index( *nb )] ) {
        abe_nb = true;
      }
    }
    if( !abe_nb ) {
      cout << "isolated" << endl;
#ifdef NOTYET
#endif
      is_had[i] = 0;
      abes[i] = 0;
      for( int j = 0 , js = hads.size() ; j < js ; ++j ) {
        if( hads[j] && i == DACLIB::atom_index( *hads[j] ) ) {
          if( mobile_h[i] ) {
            --hcount;
            hads[j]->SetImplicitHCount( hads[j]->GetImplicitHCount() + 1 );
          }
          hads[j] = 0;
          break;
        }
      }
    }
  }

  hads.erase( remove( hads.begin() , hads.end() ,
                      static_cast<OEAtomBase *>( 0 ) ) , hads.end() );

}

// ****************************************************************************
void build_abe_counts( OEGraphMol &mol , vector<OEAtomBase *> &hads ,
                       vector<int> &is_had ,
                       vector<vector<OEAtomBase *> > &bond_paths ,
                       vector<int> &abes , int &num_act_h ) {

#ifdef NOTYET
  cout << "Calc orig hcount" << endl;
#endif

  int orig_hcount = sum_had_h_counts( hads );

#ifdef NOTYET
  cout << "orig_ h count = " << orig_hcount << endl;
#endif

  vector<int> mobile_h( is_had.size() , 0 );
  remove_mobile_h_atoms( mol , hads , abes , mobile_h );

#ifdef NOTYET
  cout << "abes after mobile h removal : ";
  for( int i = 0 , is = abes.size() ; i < is ; ++i ) {
    cout << abes[i] << " ";
  }
  cout << endl;
#endif

#ifdef NOTYET
  cout << "calc final hcount" << endl;
#endif

  int final_hcount = sum_had_h_counts( hads );

#ifdef NOTYET
  cout << "final h count = " << final_hcount << endl;
#endif

  set_bondpath_bonds_to_1( mol , bond_paths , hads , is_had , abes );

#ifdef NOTYET
  remove_isolated_atoms( abes , mol , hads , is_had , mobile_h , orig_hcount );
#endif

  // final check that abes balance between neighbours.  If an atom has a
  // lower abe than its neighbour, and it has a hydrogen, remove the H and
  // balance the abe.  This is so, for example, N#C-C-* can form from
  // HN=C=C-*. CHEMBL19253 is an example
  for( OEIter<OEAtomBase> atom = mol.GetAtoms() ; atom ; ++atom ) {
    unsigned int atidx = DACLIB::atom_index( *atom );
    if( abes[atidx] ) {
      for( OEIter<OEAtomBase> nb = atom->GetAtoms() ; nb ; ++nb ) {
        unsigned int nbidx = DACLIB::atom_index( *nb );
        if( abes[atidx] > abes[nbidx] ) {
          cout << "Imbalance between " << atidx + 1 << " and " << nbidx + 1 << endl;
          if( nb->GetImplicitHCount() ) {
            nb->SetImplicitHCount( nb->GetImplicitHCount() - 1 );
            --final_hcount;
            abes[nbidx]++;
          }
        }
      }
    }
  }

  num_act_h = orig_hcount - final_hcount;

}

// ****************************************************************************
void build_t_skel_mol( OEMolBase &mol , const vector<int> &is_had ,
                       OEMolBase &t_skel_mol ) {

  // build a new molecule based on mol, but without all the aromatic flags
  // and stuff which can't be got rid of and upset the SMILES generation due
  // to whinges about kekulisation, for example.
  // This will destroy chirality, so needs re-visiting, probably by copying
  // the molecule, deleting the hads and re-adding them.  Bits of the molecule
  // that are invariant to tautomers will then preserve chirality. It'll be
  // awkward due to atom indices changing as they are deleted and added.
  vector<OEAtomBase *> new_atoms( DACLIB::max_atom_index( mol ) , static_cast<OEAtomBase *>( 0 ) );
  for( OEIter<OEAtomBase> atom = mol.GetAtoms() ; atom ; ++atom ) {
    OEAtomBase *new_at = t_skel_mol.NewAtom( atom->GetAtomicNum() );
    DACLIB::set_atom_index( *new_at , DACLIB::atom_index( *atom ) );
    new_at->SetImplicitHCount( atom->GetImplicitHCount() );
    new_atoms[DACLIB::atom_index( atom )] = new_at;
  }

  for( OEIter<OEBondBase> bond = mol.GetBonds() ; bond ; ++bond ) {
    if( is_had[DACLIB::atom_index( *bond->GetBgn() )] &&
        is_had[DACLIB::atom_index( *bond->GetEnd() )] ) {
      t_skel_mol.NewBond( new_atoms[DACLIB::atom_index( *bond->GetBgn() )] ,
          new_atoms[DACLIB::atom_index( *bond->GetEnd() )] , 1 );
    } else {
      t_skel_mol.NewBond( new_atoms[DACLIB::atom_index( *bond->GetBgn() )] ,
          new_atoms[DACLIB::atom_index( *bond->GetEnd() )] ,
          bond->GetOrder() );
    }
  }

}

// ****************************************************************************
void build_t_skel_smi( OEMolBase &mol , const vector<int> &is_had ,
                       string &t_skel_smi ) {

  OEGraphMol t_skel_mol;
  build_t_skel_mol( mol , is_had , t_skel_mol );
  OECreateSmiString( t_skel_smi , t_skel_mol ,
                     OESMILESFlag::AtomStereo | OESMILESFlag::BondStereo | OESMILESFlag::Canonical );

}

// ****************************************************************************
// get the TSAs, which are atoms in bond_paths and not in hads
void find_tsas( const vector<vector<OEAtomBase *> > &bond_paths ,
                const vector<int> &is_had ,
                vector<OEAtomBase *> &tsas ,
                vector<int> &is_tsa ) {

  is_tsa = vector<int>( is_had.size() , 0 );
  for( int i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    for( int j = 0 , js = bond_paths[i].size() ; j < js ; ++j ) {
      if( !is_had[DACLIB::atom_index( *bond_paths[i][j] )] &&
          !is_tsa[DACLIB::atom_index( *bond_paths[i][j] )] ) {
        tsas.push_back( bond_paths[i][j] );
        is_tsa[DACLIB::atom_index( *bond_paths[i][j] )] = 1;
      }
    }
  }

}

// ****************************************************************************
// find an atom which has an abe that has the minimum number of neighbours
// with abes
OEAtomBase *find_unsat_start_atom( const vector<int> &abes ,
                                   OEGraphMol &t_skel_mol ) {

  OEAtomBase *ret_val = static_cast<OEAtomBase *>( 0 );
  unsigned int min_count = numeric_limits<unsigned int>::max();

  for( int i = 0 , is = abes.size() ; i < is ; ++i ) {
    if( !abes[i] ) {
      continue;
    }
    OEAtomBase *atom = t_skel_mol.GetAtom( DACLIB::HasAtomIndex( i ) );
    unsigned int nb_abe_count = 0;
    for( OEIter<OEAtomBase> nb = atom->GetAtoms() ; nb ; ++nb ) {
      if( abes[DACLIB::atom_index( *nb )] ) {
        ++nb_abe_count;
      }
    }
    if( nb_abe_count < min_count ) {
      ret_val = atom;
      min_count = nb_abe_count;
    }
  }

  return ret_val;

}

// ****************************************************************************
void build_connect_sets( vector<int> &abes , OEGraphMol &t_skel_mol ,
                         vector<vector<OEAtomBase *> > &connect_sets ,
                         vector<vector<int> > &abe_sets ) {

  vector<int> done_ats = abes;
  while( 1 ) {
    OEAtomBase *start_atom = find_unsat_start_atom( done_ats , t_skel_mol );
    if( !start_atom ) {
      break;
    }
#ifdef NOTYET
    cout << "start_atom : " << DACLIB::atom_index( *start_atom ) + 1 << endl;
#endif
    connect_sets.push_back( vector<OEAtomBase *>() );
    abe_sets.push_back( vector<int>( abes.size() , 0 ) );

    list<OEAtomBase *> to_do( 1 , start_atom );
    while( !to_do.empty() ) {
#ifdef NOTYET
      cout << "to_do :";
      for( list<OEAtomBase *>::iterator ii = to_do.begin() ; ii != to_do.end() ; ++ii ) {
        cout << " " << DACLIB::atom_index( *(*ii) ) + 1;
      }
      cout << endl;
#endif
      OEAtomBase *next_at = to_do.front();
      to_do.pop_front();
      if( done_ats[DACLIB::atom_index( *next_at )] ) {
        abe_sets.back()[DACLIB::atom_index( *next_at )] = abes[DACLIB::atom_index( *next_at )];
        connect_sets.back().push_back( next_at );
        done_ats[DACLIB::atom_index( *next_at )] = 0;
        for( OEIter<OEAtomBase> nb = next_at->GetAtoms() ; nb ; ++nb ) {
          if( done_ats[DACLIB::atom_index( *nb )] ) {
            to_do.push_back( nb );
          }
        }
      }
    }
  }

#ifdef NOTYET
  cout << "connect sets" << endl;
  for( int i = 0 , is = connect_sets.size() ; i < is ; ++i ) {
    cout << i << " :";
    for( int j = 0 , js = connect_sets[i].size() ; j < js ; ++j ) {
      cout << " " << DACLIB::atom_index( *connect_sets[i][j] ) + 1;
    }
    cout << endl;
  }
#endif

}

// ****************************************************************************
void build_bond_atom_pairs( vector<OEAtomBase *> &connect_set ,
                            vector<int> &abes ,
                            vector<pair<OEAtomBase *,OEAtomBase *> > &bond_atom_pairs ) {

  for( int i = 0 , is = connect_set.size() ; i < is ; ++i ) {
    for( OEIter<OEAtomBase> nb = connect_set[i]->GetAtoms() ; nb ; ++nb ) {
      if( abes[DACLIB::atom_index( *nb)] ) {
        if( DACLIB::atom_index( *nb ) > DACLIB::atom_index( *connect_set[i] ) ) {
          bond_atom_pairs.push_back( make_pair( connect_set[i] , nb ) );
        }
      }
    }
  }

#ifdef NOTYET
  cout << "Bonds for this connect set :" << endl;
  for( int i = 0 , is = bond_atom_pairs.size() ; i < is ; ++i ) {
    cout << DACLIB::atom_index( *bond_atom_pairs[i].first ) + 1 << " to "
         << DACLIB::atom_index( *bond_atom_pairs[i].second ) + 1 << endl;
  }
#endif

}

// ****************************************************************************
// take the bonds and pull out the ones that have to be set to double
// to make the tautomer
bool build_double_bond_set( vector<pair<OEAtomBase *,OEAtomBase *> > bond_atom_pairs ,
                            vector<int> this_abes , int num_abes ,
                            vector<pair<OEAtomBase *,OEAtomBase *> > &double_bond_atom_pairs ) {

#ifdef NOTYET
  cout << "build_double_bond_set" << endl;
  for( int i = 0 , is = bond_atom_pairs.size() ; i < is ; ++i ) {
    cout << DACLIB::atom_index( *bond_atom_pairs[i].first ) + 1 << " to "
         << DACLIB::atom_index( *bond_atom_pairs[i].second ) + 1 << endl;
  }
  cout << "abes :";
  for( int i = 0 , is = this_abes.size() ; i < is ; ++i ) {
    cout << " " << this_abes[i];
  }
  cout << endl;

  cout << "double bonds so far" << endl;
  for( int i = 0 , is = double_bond_atom_pairs.size() ; i < is ; ++i ) {
    cout << DACLIB::atom_index( *double_bond_atom_pairs[i].first ) + 1 << " to "
         << DACLIB::atom_index( *double_bond_atom_pairs[i].second ) + 1 << endl;
  }
#endif

  for( int i = 0 , is = bond_atom_pairs.size() ; i < is ; ++i ) {
    // try the next bond
    double_bond_atom_pairs.push_back( bond_atom_pairs[i] );
    --this_abes[DACLIB::atom_index( *bond_atom_pairs[i].first )];
    --this_abes[DACLIB::atom_index( *bond_atom_pairs[i].second )];
    num_abes -= 2;
    if( !num_abes ) {
#ifdef NOTYET
      cout << "RESULT 1" << endl;
      for( int ii = 0 , iis = double_bond_atom_pairs.size() ; ii < iis ; ++ii ) {
        cout << DACLIB::atom_index( *double_bond_atom_pairs[ii].first ) + 1 << " to "
             << DACLIB::atom_index( *double_bond_atom_pairs[ii].second ) + 1 << endl;
      }
#endif
      return true;
    }
    // make a new vector of bond_atom_pairs without any bond that uses either of
    // the atoms we've just tried if they now have a 0 for abes.
    // NB The same bond may have to go in twice for a triple bond and this is
    // not yet dealt with.
    vector<pair<OEAtomBase *,OEAtomBase *> > next_bond_atom_pairs;
    for( int j = 0 , js = bond_atom_pairs.size() ; j < js ; ++j ) {
      if( this_abes[DACLIB::atom_index( *bond_atom_pairs[j].first )] &&
          this_abes[DACLIB::atom_index( *bond_atom_pairs[j].second )] ) {
        next_bond_atom_pairs.push_back( bond_atom_pairs[j] );
      }
    }
    if( !next_bond_atom_pairs.empty() ) {
      if( build_double_bond_set( next_bond_atom_pairs , this_abes , num_abes ,
                                 double_bond_atom_pairs ) ) {
        return true;
      }
    }
    // set it up for another try
    double_bond_atom_pairs.pop_back();
    ++this_abes[DACLIB::atom_index( *bond_atom_pairs[i].first )];
    ++this_abes[DACLIB::atom_index( *bond_atom_pairs[i].second )];
    num_abes += 2;
  }

  return false;

}

// ****************************************************************************
// take the t_skel_mol and the current abe list and try and make it a full
// molecule by adding unsaturated bonds
bool add_unsaturated_bonds( vector<int> &abes , OEGraphMol &t_skel_mol ) {

#ifdef NOTYET
  cout << "add_unsaturated_bonds, abes :";
  for( int i = 0 , is = abes.size() ; i < is ; ++i ) {
    cout << " " << abes[i];
  }
  cout << endl;
#endif

  // start by getting the connect sets together - these are contiguous bits
  // of abe atoms.  Plonking the H atoms in might have split up any previously
  // contiguous pieces
  vector<vector<OEAtomBase *> > connect_sets;
  vector<vector<int> > abe_sets;
  build_connect_sets( abes , t_skel_mol , connect_sets , abe_sets );

  vector<pair<OEAtomBase *,OEAtomBase *> > final_atom_pairs;
  unsigned int total_abes = 0;

  // now do each connect set in turn
  for( int i = 0 , is = connect_sets.size() ; i < is ; ++i ) {
    vector<pair<OEAtomBase *,OEAtomBase *> > bond_atom_pairs;
    // get all the bonds involved in this connect set
    build_bond_atom_pairs( connect_sets[i] , abe_sets[i] ,
                           bond_atom_pairs );

    vector<pair<OEAtomBase *,OEAtomBase *> > double_bond_atom_pairs;

    vector<int> &this_abes = abe_sets[i];
    int num_abes = accumulate( this_abes.begin() , this_abes.end() , 0 );
#ifdef NOTYET
    cout << "num_abes : " << num_abes << endl;
#endif
    total_abes += num_abes;
    if( build_double_bond_set( bond_atom_pairs , this_abes , num_abes ,
                               double_bond_atom_pairs ) ) {
      final_atom_pairs.insert( final_atom_pairs.end() ,
                               double_bond_atom_pairs.begin() ,
                               double_bond_atom_pairs.end() );
    }

  }

  if( total_abes == 2 * final_atom_pairs.size() ) {
#ifdef NOTYET
    cout << "ITS A GOOD ONE" << endl;
    cout << "final_atom_pairs" << endl;
#endif
    for( int i = 0 , is = final_atom_pairs.size() ; i < is ; ++i ) {
#ifdef NOTYET
      cout << DACLIB::atom_index( *final_atom_pairs[i].first ) + 1 << " to "
           << DACLIB::atom_index( *final_atom_pairs[i].second ) + 1 << endl;
#endif
      OEBondBase *bond = t_skel_mol.GetBond( t_skel_mol.GetAtom( DACLIB::HasAtomIndex( DACLIB::atom_index( *final_atom_pairs[i].first ) ) ) ,
                                             t_skel_mol.GetAtom( DACLIB::HasAtomIndex( DACLIB::atom_index( *final_atom_pairs[i].second ) ) ) );
      bond->SetOrder( bond->GetOrder() + 1 );
    }
    return true;
  } else {
#ifdef NOTYET
    cout << "ITS INCOMPLETE : " << total_abes << " vs " << final_atom_pairs.size() << endl;
#endif
    return false;
  }

}

// ****************************************************************************
// see if there is an atom with a bonding electron (abes is non-zero)
// that doesn't have neighbour with a bonding electron
bool are_there_isolated_atoms( const vector<int> &abes ,
                               OEGraphMol &mol ) {

  for( int i = 0 , is = abes.size() ; i < is ; ++i ) {
    // if this atom has no bonding electrons, it doesn't matter
    if( !abes[i] ) {
      continue;
    }
#ifdef NOTYET
    cout << "looking at atom " << i + 1 << " abe = " << abes[i] << endl;
#endif
    bool abe_nb = false;
    OEAtomBase *atom = mol.GetAtom( DACLIB::HasAtomIndex( i ) );
    for( OEIter<OEAtomBase> nb = atom->GetAtoms() ; nb ; ++nb ) {
      if( abes[DACLIB::atom_index( *nb )] ) {
        abe_nb = true;
      }
    }
    if( !abe_nb ) {
#ifdef NOTYET
      cout << "isolated" << endl;
#endif
      return true;
    }
  }

  return false;

}

// ****************************************************************************
// make all the possible tautomers.  This may identify atoms that can't be
// part of the T-Skeleton as it's not possible to generate a valid tautomer
// using them. They will be removed. is_had may change.
void generate_tautomers( OEGraphMol &inmol ,
                         vector<int> &abes , int num_act_h ,
                         vector<int> &is_had ,
                         vector<OEAtomBase *> &hads ,
                         vector<int> &is_tsa ,
                         vector<OEAtomBase *> &tsas ,
                         vector<string> &taut_smis ) {

  vector<unsigned int> had_idxs;
  for( int i = 0 , is = hads.size() ; i < is ; ++i ) {
    had_idxs.push_back( DACLIB::atom_index( *hads[i] ) );
  }

  // num_act_h is the number of tautomer active H atoms.  We want to distribute
  // them across the hads in all combinations and then try and put double bonds
  // to mend the valences and then see if we get a proper molecule out

  taut_smis.clear();

  cout << "Num act h : " << num_act_h << endl;
#ifdef NOTYET
#endif
  DACLIB::Combinator comb( num_act_h , hads.size() );
  do {
    const vector<int> &next_comb = comb.value();
    for( int i = 0 , is = hads.size() ; i < is ; ++i ) {
      cout << i << " : " << hads[i] << endl;
    }
    cout << "Next comb :";
    copy( next_comb.begin() , next_comb.end() , intOut );
    cout << endl;
    for( int i = 0 ; i < num_act_h ; ++i ) {
      cout << " " << DACLIB::atom_index( *hads[next_comb[i]] ) + 1;
    }
    cout << endl;
#ifdef NOTYET
#endif

    OEGraphMol t_skel_mol;
    build_t_skel_mol( inmol , is_had , t_skel_mol );
    for( int i = 0 ; i < num_act_h ; ++i ) {
      OEAtomBase *atom = t_skel_mol.GetAtom( DACLIB::HasAtomIndex( DACLIB::atom_index( *hads[next_comb[i]] ) ) );
      atom->SetImplicitHCount(( atom->GetImplicitHCount() + 1 ) );
    }
#ifdef NOTYET
    string taut_smi;
    OECreateSmiString( taut_smi , t_skel_mol ,
                       OESMILESFlag::AtomStereo | OESMILESFlag::BondStereo );
    cout << "taut_smi : " << taut_smi << endl;
#endif

    // check to see if this comb has produced isolated tsas. If so, this is
    // not one that can work
    vector<int> this_abes = abes;
    for( int i = 0 ; i < num_act_h ; ++i ) {
      --this_abes[had_idxs[next_comb[i]]];
    }

    if( are_there_isolated_atoms( this_abes , t_skel_mol ) ) {
#ifdef NOTYET
      cout << "  Bad comb - isolated atoms" << endl;
#endif
      continue;
    } else {
#ifdef NOTYET
      cout << "  Good comb";
      for( int i = 0 ; i < num_act_h ; ++i ) {
        cout << " " << DACLIB::atom_index( *hads[next_comb[i]] ) + 1;
      }
      cout << endl;
#endif
      if( add_unsaturated_bonds( this_abes , t_skel_mol ) ) {
        string taut_smi;
        OECreateSmiString( taut_smi , t_skel_mol ,
                           OESMILESFlag::AtomStereo | OESMILESFlag::BondStereo | OESMILESFlag::Canonical );
#ifdef NOTYET
        cout << "Got a Good Tautomer : " << taut_smi << endl;
#endif
        taut_smis.push_back( taut_smi );
      } else {
#ifdef NOTYET
        cout << "  Bad comb - unsaturated bonds" << endl;
#endif
      }
    }

  } while( comb.step() , !comb.at_end() );

  sort( taut_smis.begin() , taut_smis.end() );
  taut_smis.erase( unique( taut_smis.begin() , taut_smis.end() ) , taut_smis.end() );

}

// ****************************************************************************
void make_taut_skeleton( const string &in_smi , string &t_skel_smi ,
                         vector<string> &taut_smis ) {

  OEGraphMol mol;
  OEParseSmiles( mol , in_smi );
  DACLIB::apply_daylight_aromatic_model( mol );

#ifdef NOTYET
  for( OEIter<OEBondBase> bond = mol.GetBonds() ; bond ; ++bond ) {
    cout << DACLIB::atom_index( *bond->GetBgn() ) + 1 << " to "
         << DACLIB::atom_index( *bond->GetEnd() ) + 1 << " order = " << bond->GetOrder() << endl;
  }
#endif

  // HAD is an H-atom acceptor atom or H-atom donor atom.
  vector<OEAtomBase *> hads;
  vector<int> is_had;
  vector<vector<OEAtomBase *> > bond_paths;
  find_hads( mol , bond_paths , hads , is_had );

  // ABEs are Atom Bonding Electrons and are the free valences when the mobile H
  // atoms are removed and the bondpath orders are set to 1. Hads and is_hads
  // can be changed in this function.
  int num_act_h = 0;
  vector<int> abes( DACLIB::max_atom_index( mol ) , 0 );
  build_abe_counts( mol , hads , is_had , bond_paths , abes ,
                    num_act_h );

  // get the TSAs, which are atoms in bond_paths and not in hads
  vector<OEAtomBase *> tsas;
  vector<int> is_tsa;
  find_tsas( bond_paths , is_had , tsas , is_tsa );

  cout << "final HAD list :";
  for( int i = 0 , is = hads.size() ; i < is ; ++i ) {
    cout << " " << DACLIB::atom_index( *hads[i] ) + 1;
  }
  cout << endl;
  cout << "final TSA list :";
  for( int i = 0 , is = tsas.size() ; i < is ; ++i ) {
    cout << " " << DACLIB::atom_index( *tsas[i] ) + 1;
  }
  cout << endl;
#ifdef NOTYET
  cout << "Final bond paths" << endl;
  for( unsigned int j = 0 , js = bond_paths.size() ; j < js ; ++j ) {
    cout << "Path " << j << " :: ";
    for( unsigned int i = 0 , is = bond_paths[j].size() ; i < is ; ++i ) {
      cout << " " << DACLIB::atom_index( *bond_paths[j][i] ) + 1;
    }
    cout << endl;
  }
#endif

#ifdef NOTYET
  cout << "Num tautomer active H : " << num_act_h << endl;
  cout << "abes_counts :";
  for( int i = 0 , is = abes.size() ; i < is ; ++i ) {
    cout << " " << abes[i];
  }
  cout << endl;
#endif

  // make all the possible tautomers.  This may identify atoms that can't be
  // part of the T-Skeleton as it's not possible to generate a valid tautomer
  // using them. They will be removed. is_had may change.
  generate_tautomers( mol , abes , num_act_h , is_had , hads ,
                      is_tsa , tsas , taut_smis );

  build_t_skel_smi( mol , is_had , t_skel_smi );

#ifdef NOTYET
  cout << "t_skel_smi : " << t_skel_smi << endl;
  cout << "taut_smis : " << endl;
  for( int i = 0 , is = taut_smis.size() ; i < is ; ++i ) {
    cout << taut_smis[i] << endl;
  }
#endif

}
