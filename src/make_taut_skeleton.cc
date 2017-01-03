//
// file make_taut_skeleton.cc
// David Cosgrove
// AstraZeneca
// 15th June 2015
//
// This file contains functions to take a SMILES string and return the SMILES
// of the corresponding tautomer skeleton, as defined by Thalheim et al.
// 'A Branch-and-Bound Approach for Tautomer Enumeration',
// Molecular Informatics, 2015 - At present I have the ASAP version, without a full
// reference.  DOI: 10.1002/minf.201400128.
// It's behind a paywall.
//
// Abbreviations used:
// HAD - Hydrogen atom Acceptor or Donor. A sink or source of an H in a
// tautomer system. Not to be confused with an h-bond donor or acceptor.
//
// NUTPOACTSYO - Never Underestimate The Power Of A Chemist To Screw You Over.
// Cosgrove's 1st law of Cheminformatics.

#include "stddefs.H"
#include "DACOEMolAtomIndex.H"
#include "DACOEMolBondIndex.H"

#include "TautomerGenerator.H"

#include <algorithm>
#include <iostream>
#include <limits>
#include <list>
#include <numeric>
#include <string>
#include <vector>

#include <oechem.h>

#include <boost/bind.hpp>
#include <boost/current_function.hpp>
#include <boost/shared_ptr.hpp>

using namespace boost;
using namespace std;
using namespace OEChem;
using namespace OESystem;

namespace DACLIB {
// in eponymous file
void apply_daylight_aromatic_model( OEMolBase &mol );
string create_cansmi( const OEMolBase &mol );
string create_noncansmi( const OEMolBase &mol );
void radical_atoms( OEMolBase &mol , vector<OEAtomBase *> &rad_atoms );
}

typedef boost::shared_ptr<OEMolBase> pOEMolBase;
typedef boost::shared_ptr<TautomerGenerator> pTautGen;

// ****************************************************************************
OEAtomBase *atom_has_multiple_bond( OEAtomBase *atom ) {

  OEIter<OEBondBase> bond = atom->GetBonds( OENot<OEBondBase>( OEHasOrder( 1 ) ) );
  if( bond ) {
    if( bond->GetBgn() == atom ) {
      return bond->GetEnd();
    } else {
      return bond->GetBgn();
    }
  } else {
    return static_cast<OEAtomBase *>( 0 );
  }

}
// ****************************************************************************
OEAtomBase *atom_has_het_nbour( OEAtomBase *atom ) {

  OEIter<OEAtomBase> het_at = atom->GetAtoms( OEIsHetero() );
  if( het_at ) {
    return het_at;
  } else {
    return static_cast<OEAtomBase *>( 0 );
  }

}

// ****************************************************************************
// Rule 5 is that a C had must have a bond path of length 2 to a hetero atom
// had, but I've also extended it. So this pre-filter just checks that there's
// a hetero 2 or 4 bonds from the atom. The latter is for the possibility of
// a 1,5 shift.
bool rule_5_pre_filter( OEAtomBase *atom ) {

  for( OEIter<OEAtomBase> nb1 = atom->GetAtoms() ; nb1 ; ++nb1 ) {
    if( OEElemNo::C != nb1->GetAtomicNum() ) {
      return true;
    }
    for( OEIter<OEAtomBase> nb2 = nb1->GetAtoms() ; nb2 ; ++nb2 ) {
      if( OEElemNo::C != nb2->GetAtomicNum() ) {
        return true;
      }
      for( OEIter<OEAtomBase> nb3 = nb2->GetAtoms() ; nb3 ; ++nb3 ) {
        if( nb1->GetIdx() == nb3->GetIdx() ) {
          continue;
        }
        for( OEIter<OEAtomBase> nb4 = nb3->GetAtoms() ; nb4 ; ++nb4 ) {
          if( nb2->GetIdx() == nb4->GetIdx() ) {
            continue;
          }
          if( OEElemNo::C != nb4->GetAtomicNum() ) {
#ifdef NOTYET
            cout << "Atom " << DACLIB::atom_index( *atom ) + 1 << " passes 4 bond Rule 5 pre-filter" << endl;
#endif
            return true;
          }
        }
      }
    }
  }

#ifdef NOTYET
  cout << "Atom " << DACLIB::atom_index( *atom ) + 1 << " fails Rule 5 pre-filter" << endl;
#endif

  return false;

}

// ****************************************************************************
// rule 4 - extended because of CHEMBL440484. This contains O=C1COC=C1
// which can be converted to Oc1cocc1 but won't go back again due to the
// aromatic ring that's formed. Relax rule 4 to allow hads on aromatic atoms
// if the atom is adjacent to a ring atom containing a hetero with an H
// and it's in a five-membered ring (which, if aromatic, must be a hetero
// ring).
bool relaxed_rule_4( OEMolBase &mol , OEAtomBase *atom ) {

#ifdef NOTYET
  cout << "relaxed_rule_4 for " << DACLIB::atom_index( *atom ) + 1 << endl;
#endif
  if( OEAtomIsInRingSize( atom , 5 ) ) {
    for( OEIter<OEAtomBase> nb = atom->GetAtoms() ; nb ; ++nb ) {
      if( !OEAtomIsInRingSize( *nb , 5 ) || !nb->IsInRing() || 1 == mol.GetBond( atom , nb )->GetOrder() ) {
        continue;
      }
      for( OEIter<OEAtomBase> nb_nb = nb->GetAtoms() ; nb_nb ; ++nb_nb ) {
#ifdef NOTYET
        cout << "nbor_nbor : " << DACLIB::atom_index( *nb_nb ) + 1
             << " : " << nb_nb->GetTotalHCount() << " : "
             << nb_nb->GetAtomicNum() << " : "
             << nb_nb->IsInRing() << endl;
#endif
        if( nb_nb->GetTotalHCount() && OEElemNo::C != nb_nb->GetAtomicNum() &&
            !nb_nb->IsInRing() ) {
#ifdef NOTYET
          cout << DACLIB::atom_index( *atom ) + 1 << " passes relaxed rule 4" << endl;
#endif
          return true;
        }
      }
    }
  }

  return false;

}

// ****************************************************************************
void extend_ring_path( vector<OEAtomBase *> &path ,
                       const vector<unsigned int> &ring_systs ,
                       unsigned int ring_syst ,
                       vector<vector<OEAtomBase *> > &rings ) {

  for( OEIter<OEAtomBase> atom = path.back()->GetAtoms() ; atom ; ++ atom ) {
    if( ring_syst != ring_systs[atom->GetIdx()] ) {
      continue;
    }
    vector<OEAtomBase *>::iterator p = find( path.begin() , path.end() , atom );
    if( p == path.end() ) {
      path.push_back( atom );
      extend_ring_path( path , ring_systs , ring_syst , rings );
      path.pop_back();
    } else {
      // there's a ring
      if( distance( p , path.end() ) > 2 ) {
        rings.push_back( vector<OEAtomBase *>( p , path.end() ) );
      }
    }
  }

}

// ****************************************************************************
// find all the rings in the given ring system.
void find_rings_in_ring_syst( const vector<unsigned int> &ring_systs ,
                              unsigned int ring_syst , OEAtomBase *start_atom ,
                              vector<vector<OEAtomBase *> > &rings ) {

  vector<OEAtomBase *> path( 1 , start_atom );
  extend_ring_path( path , ring_systs , ring_syst , rings );

#ifdef NOTYET
  cout << "RINGS" << endl;
  for( size_t i = 0 , is = rings.size() ; i < is ; ++i ) {
    cout << i << " :";
    for( size_t j = 0 , js = rings[i].size() ; j < js ; ++j ) {
      cout << " " << DACLIB::atom_index( *rings[i][j] ) + 1;
    }
    cout << endl;
  }
#endif

}

// ****************************************************************************
// returns true if the n_atom is in the same 6-membered aromatic ring as
// another N.
bool six_mem_aromatic_ring_rule( OEAtomBase *n_atom ,
                                 OEMolBase &mol ) {
  if( !OEAtomIsInAromaticRingSize( n_atom , 6 ) ) {
    return false;
  }

  // get the aromatic ring systems.
  vector<unsigned int> arom_ring_systs( mol.GetMaxAtomIdx() , 0 );
  OEDetermineAromaticRingSystems( mol , &arom_ring_systs[0] );

  unsigned int ring_syst = arom_ring_systs[n_atom->GetIdx()];

  vector<vector<OEAtomBase *> > rings;
  // find other N atoms in this ring system
  for( OEIter<OEAtomBase> other_n = mol.GetAtoms( OEHasAtomicNum( OEElemNo::N ) ) ; other_n ; ++other_n ) {
    if( n_atom == other_n ) {
      continue;
    }
    if( arom_ring_systs[other_n->GetIdx()] != ring_syst ) {
      continue;
    }
    if( !OEAtomIsInAromaticRingSize( *other_n , 6 ) ) {
      return false; // can't be in the same ring
    }
    // if the shortest path between these 2 atoms is more than 3 bonds,
    // they can't be in the same 6-membered ring.
    unsigned int min_path_len = OEGetPathLength( n_atom , other_n , 3 );
    if( min_path_len > 3 ) {
      return false;
    }
    // Now it gets complicated.  We need all rings in this ring system.
    if( rings.empty() ) {
      find_rings_in_ring_syst( arom_ring_systs , ring_syst , n_atom , rings );
    }
    for( size_t i = 0 , is = rings.size() ; i < is ; ++i ) {
      if( 6 != rings[i].size() ) {
        continue;
      }
      vector<OEAtomBase *>::iterator p = find( rings[i].begin() , rings[i].end() , n_atom );
      if( p != rings[i].end() ) {
        vector<OEAtomBase *>::iterator q = find( rings[i].begin() , rings[i].end() , other_n );
        if( q != rings[i].end() ) {
          // n_atom and other_n in same  6-membered ring
#ifdef NOTYET
          cout << DACLIB::atom_index( *n_atom ) + 1 << " and "
               << DACLIB::atom_index( *other_n ) + 1 << " in same ring" << endl;
#endif
          return true;
        }
      }
    }
  }

  return false;

}

// ****************************************************************************
// returns true if the O atom is part of a nitro group
bool nitro_group( OEAtomBase *o_atom ) {

  if( o_atom->GetHvyDegree() > 1 ) {
    return false;
  }

  OEIter<OEAtomBase> n_atom = o_atom->GetAtoms( OEHasAtomicNum( OEElemNo::N ) );
  if( !n_atom ) {
    return false;
  }

  int o_count = 0;
  for( OEIter<OEAtomBase> oo_atom = n_atom->GetAtoms( OEHasAtomicNum( OEElemNo::O ) ) ; oo_atom ; ++oo_atom ) {
    ++o_count;
  }

  return ( 2 == o_count );

}

// ****************************************************************************
// returns true if the N atom is the central one of an azido group, [N-]-[N+]#N
// or N=[N+]=[N-], so 2-connected and positive and attached to 2 N atoms
bool azido_group( OEAtomBase *n_atom ) {

  if( n_atom->GetHvyDegree() != 2 ) {
    return false;
  }
  if( n_atom->GetFormalCharge() != +1 ) {
    return false;
  }

  int n_count = 0;
  for( OEIter<OEAtomBase> no_atom = n_atom->GetAtoms() ; no_atom ; ++no_atom ) {
    ++n_count;
  }

  return ( 2 == n_count );

}


// ****************************************************************************
void build_initial_had_list( OEMolBase &mol , vector<OEAtomBase *> &hads ,
                             vector<int> &is_had ) {

  is_had = vector<int>( DACLIB::max_atom_index( mol ) , 0 );
  bool is_het = false;
  for( OEIter<OEAtomBase> atom = mol.GetAtoms( OEIsHeavy() ) ; atom ; ++atom ) {
    // rule 2
    if( !atom->GetTotalHCount() && !atom_has_multiple_bond( atom ) ) {
      continue;
    }
    // rule 4 - now extended because of CHEMBL440484. This contains O=C1COC=C1
    // which can be converted to Oc1cocc1 but won't go back again due to the
    // aromatic ring that's formed. Relax rule 4 to allow hads on aromatic atoms
    // if the atom is adjacent to a ring atom containing a hetero with an H
    // and it's in a five-membered ring (which, if aromatic, must be a hetero
    // ring).
    if( OEElemNo::C == atom->GetAtomicNum() && atom->IsAromatic() &&
        !relaxed_rule_4( mol , atom ) ) {
      continue;
    }
    // Rule 5 states that a C atom must have a bond path of length 2 to at
    // least 1 het atom. For the pre-filter, just check the C atom has a het
    // within 2 bonds. Exact definition of a bond path comes later.
    // Also, I've extended it to allow for 1,5 shifts in non-cyclic systems,
    // as in CC=CC=O where the H on the first C forms a pseudo 6-membered
    // ring that can pass the H. The latter is in agl1.smi.
    if( OEElemNo::C == atom->GetAtomicNum() && !rule_5_pre_filter( atom ) ) {
      continue;
    }
    if( OEElemNo::C != atom->GetAtomicNum() ) {
      is_het = true;
    }
    // This is a new rule suggested by Martin Packer.  If a nitrogen is in a
    // 6-membered aromatic ring, it can't accept a proton if there's another N
    // in the same ring.
    if( OEElemNo::N == atom->GetAtomicNum() &&
        six_mem_aromatic_ring_rule( atom , mol  ) ) {
      continue;
    }
    // nitro groups don't look sensible, either
    if( OEElemNo::O == atom->GetAtomicNum() && nitro_group( atom ) ) {
      continue;
    }
    // and azido groups probably best ignored (c.f. CHEMBL6313)
    if( OEElemNo::N == atom->GetAtomicNum() && azido_group( atom ) ) {
      continue;
    }
    hads.push_back( atom );
    is_had[DACLIB::atom_index( atom )] = 1;
  }

  // As shown by CHEMBL439119, you can get buckminsterfullerenes where all the
  // HADs are carbon, and it's not really feasible to generate all the paths
  // to show that there won't be any tautomers. So check here and save a lot
  // of grief.
  if( !is_het ) {
    hads.clear();
    is_had = vector<int>( DACLIB::max_atom_index( mol ) , 0 );
  }

}

// ****************************************************************************
void extend_bond_path( vector<int> &is_had ,
                       vector<OEAtomBase *> &curr_path ,
                       vector<vector<OEAtomBase *> > &bond_paths ) {

#ifdef NOTYET
  cout << "Extending bond path :";
  for( size_t i = 0 , is = curr_path.size() ; i < is ; ++i ) {
    cout << " " << DACLIB::atom_index( *curr_path[i] ) + 1;
  }
  cout << endl;
#endif

  for( OEIter<OEAtomBase> nb = curr_path.back()->GetAtoms( OEIsHeavy() ) ; nb ; ++nb ) {
    if( curr_path.end() != find( curr_path.begin() , curr_path.end() , nb ) ) {
      continue; // already in
    }
#ifdef NOTYET
    cout << "from " << DACLIB::atom_index( *curr_path.back() ) + 1 << " to "
         << DACLIB::atom_index( *nb ) + 1 << endl;
#endif
    // if it's an atom like an ether oxygen or an N with 3 heavy atoms
    // attached, it can't be in a tautomeric bond path
    if( !nb->GetTotalHCount() && !atom_has_multiple_bond( nb ) ) {
      continue;
    }
    // don't want any S-Het bonds - we really don't want sulfonamides
    // included in tautomer systems, for example.
    if( ( OEElemNo::S == curr_path.back()->GetAtomicNum() &&
          curr_path.back()->GetDegree() > 2 &&
          OEElemNo::C != nb->GetAtomicNum() ) ||
        ( OEElemNo::C != curr_path.back()->GetAtomicNum() &&
          OEElemNo::S == nb->GetAtomicNum() &&
          nb->GetDegree() > 2 ) ) {
#ifdef NOTYET
      cout << "Ignoring S atom : " << DACLIB::atom_index( *curr_path.back() ) + 1 << endl;
#endif
      continue;
    }

    // Likewise for phosphates.
    if( ( OEElemNo::P == curr_path.back()->GetAtomicNum() &&
          curr_path.back()->GetDegree() > 2 &&
          OEElemNo::C != nb->GetAtomicNum() ) ||
        ( OEElemNo::C != curr_path.back()->GetAtomicNum() &&
          OEElemNo::P == nb->GetAtomicNum() &&
          nb->GetDegree() > 2 ) ) {
#ifdef NOTYET
      cout << "Ignoring P atom : " << DACLIB::atom_index( *curr_path.back() ) + 1 << endl;
#endif
      continue;
    }

    curr_path.push_back( nb );
    // to be acceptable, the bond path must be at least 3 atoms (2 bonds)
    // and an even number of bonds and one end must be unsaturated (to receive
    // the H) and the other end must have an H atom
    if( is_had[DACLIB::atom_index( *nb )] ) {
      size_t num_bonds = curr_path.size() - 1;
      if( num_bonds >= 2 && !( num_bonds % 2 ) ) {
        if( ( atom_has_multiple_bond( curr_path.front() ) &&
              curr_path.back()->GetTotalHCount() ) ||
            ( atom_has_multiple_bond( curr_path.back() ) &&
              curr_path.front()->GetTotalHCount() ) ) {
          bond_paths.push_back( curr_path );
        }
      }
    }
    extend_bond_path( is_had , curr_path , bond_paths );
    curr_path.pop_back();
  }

}

// ****************************************************************************
// build all bond paths between hads.
void build_bond_paths( vector<OEAtomBase *> &hads ,
                       vector<int> &is_had ,
                       vector<vector<OEAtomBase *> > &bond_paths ) {

  for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
#ifdef NOTYET
    cout << "Making paths from " << DACLIB::atom_index( *hads[i] ) + 1 << endl;
#endif
    vector<OEAtomBase *> curr_path = ( vector<OEAtomBase *>( 1 , hads[i] ) );

    extend_bond_path( is_had , curr_path , bond_paths );

#ifdef NOTYET
    cout << "Number of paths ; " << bond_paths.size() << endl;
#endif
  }

}

// ****************************************************************************
// take out any paths that start or finish with a non-had, which can happen
// if a had has been removed from the list
void remove_paths_of_non_hads( vector<int> &is_had ,
                               vector<vector<OEAtomBase *> > &bond_paths ) {

  for( size_t i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    if( !is_had[DACLIB::atom_index( *bond_paths[i].front() )] ||
        !is_had[DACLIB::atom_index( *bond_paths[i].back() )] ) {
#ifdef NOTYET
      cout << "removing bond_path :";
      for( size_t ii = 0 , iis = bond_paths[i].size() ; ii < iis ; ++ii ) {
        cout << " " << DACLIB::atom_index( *bond_paths[i][ii] ) + 1 ;
      }
      cout << " because " << DACLIB::atom_index( *bond_paths[i].front() ) + 1
           << " or " << DACLIB::atom_index( *bond_paths[i].back() ) + 1
           << " not a had any more" << endl;
#endif
      bond_paths[i].clear();
    }
  }

  bond_paths.erase( remove_if( bond_paths.begin() , bond_paths.end() ,
                               bind( &vector<OEAtomBase *>::empty , _1 ) ) ,
                               bond_paths.end() );

}

// ****************************************************************************
// take out any paths that start AND finish with a non-had, which can happen
// if a had has been removed from the list
void remove_paths_with_non_hads_both_ends( vector<int> &is_had ,
                                           vector<vector<OEAtomBase *> > &bond_paths ) {

  for( size_t i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    if( !is_had[DACLIB::atom_index( *bond_paths[i].front() )] &&
        !is_had[DACLIB::atom_index( *bond_paths[i].back() )] ) {
#ifdef NOTYET
      cout << "removing bond_path :";
      for( size_t ii = 0 , iis = bond_paths[i].size() ; ii < iis ; ++ii ) {
        cout << " " << DACLIB::atom_index( *bond_paths[i][ii] ) + 1 ;
      }
      cout << " because " << DACLIB::atom_index( *bond_paths[i].front() ) + 1
           << " or " << DACLIB::atom_index( *bond_paths[i].back() ) + 1
           << " not a had any more" << endl;
#endif
      bond_paths[i].clear();
    }
  }

  bond_paths.erase( remove_if( bond_paths.begin() , bond_paths.end() ,
                               bind( &vector<OEAtomBase *>::empty , _1 ) ) ,
                               bond_paths.end() );

}

// ****************************************************************************
// rule 3 - a had without an H atom (the acceptor) must have a bond-path to
// another had that can release an H (the donor).
void apply_rule_3( vector<vector<OEAtomBase *> > &bond_paths ,
                   vector<OEAtomBase *> &hads ,
                   vector<int> &is_had ) {

  // only need to consider first and last elements of each path. Some atoms,
  // such as imine nitrogens, can be both acceptor and donor.
  for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
#ifdef NOTYET
    cout << "rule 3 for " << DACLIB::atom_index( *hads[i] ) + 1 << endl;
#endif
    if( hads[i]->GetTotalHCount() ) {
#ifdef NOTYET
      cout << "Skipping " << DACLIB::atom_index( *hads[i] ) + 1 << " due to H" << endl;
#endif
      continue;
    }
    int ok_to_keep = -1;
    for( int j = 0 , js = static_cast<int>( bond_paths.size() ) ; j < js ; ++j ) {
#ifdef NOTYET
      cout << "Path " << j << " :: ";
      for( size_t ii = 0 , iis = bond_paths[j].size() ; ii < iis ; ++ii ) {
        cout << " " << DACLIB::atom_index( *bond_paths[j][ii] ) + 1;
      }
      cout << endl;
#endif
      if( ( hads[i] == bond_paths[j].front() &&
            ( bond_paths[j].back()->GetTotalHCount() ) ) ||
          ( hads[i] == bond_paths[j].back() &&
            ( bond_paths[j].front()->GetTotalHCount() ) ) ) {
        ok_to_keep = j;
        break;
      }
    }
    if( -1 == ok_to_keep ) {
#ifdef NOTYET
      cout << "Removing HAD " << DACLIB::atom_index( *hads[i] ) + 1
           << " due to rule 3" << endl;
#endif
      is_had[DACLIB::atom_index( *hads[i] )] = 0;
      hads[i] = 0;
    } else {
#ifdef NOTYET
      cout << "Keeping HAD " << DACLIB::atom_index( *hads[i] ) + 1
           << " due to rule 3 with "
           << "Path " << ok_to_keep << " ::";
      for( size_t ii = 0 , iis = bond_paths[ok_to_keep].size() ; ii < iis ; ++ii ) {
        cout << " " << DACLIB::atom_index( *bond_paths[ok_to_keep][ii] ) + 1;
      }
      cout << endl;
#endif
    }
  }

  hads.erase( remove( hads.begin() , hads.end() , static_cast<OEAtomBase *>( 0 ) ) ,
              hads.end() );

  // take out any paths that start or finish with a non-had, which can happen
  // if a had has been removed from the list
  remove_paths_of_non_hads( is_had , bond_paths );

}

// ****************************************************************************
// rule 5 is that you can't have bond paths of length 2 between carbon atoms.
// I had spent a lot of time worrying about the extra cases, as shown by
// CHEMBL19253 where if we have CC(=X)C where X is hetero, or the hidden
// form CC(XH)=C, then a hydrogen can be passed from one carbon to the other
// via the intermediate enol-like form. This clearly matters in asymmetric
// systems.  It is now being dealt with by iterating round all the intermediate
// forms and expanding the t_skel as appropriate.
// This function should only be called for bond_paths of size 3, both ends
// being C.
// Check CHEMBL501944_small also works.
bool rule_5_test( vector<OEAtomBase *> bond_path ) {

  if( 3 == bond_path.size() &&
      ( OEElemNo::C != bond_path.back()->GetAtomicNum() ||
        OEElemNo::C != bond_path.front()->GetAtomicNum() ) &&
      !bond_path[0]->IsAromatic() && !bond_path[1]->IsAromatic() &&
      !bond_path[2]->IsAromatic() ) {
    return true;
  }

  return false;

}

// ****************************************************************************
// This is an extension to rule 5 that allows for H atom shifts between C and
// non-C via a 6-membered pseudo-ring akin to an intramolecular H bond. This
// should only have been called if bond_path is size 5 (4 bonds)
bool rule_5_test_ext( const vector<unsigned int> &atom_ring_systs ,
                      const vector<OEAtomBase *> &bond_path ) {

  if( 5 != bond_path.size() ||
      ( OEElemNo::C == bond_path.back()->GetAtomicNum() &&
        OEElemNo::C == bond_path.front()->GetAtomicNum() ) ) {
    return false;
  }

  // For a change, atom_ring_systs indexes using atom.GetIdx() not the DACLIB
  // version. It should be safe for this, as the molecule won't be edited
  // or copied during bits of the program where this function might be called.
  // See if there are ring systems with more than 3 atoms in the bond_path.
  // They are deemed too inflexible to form a 6-membered pseudo-ring.
  for( size_t i = 0 , is = bond_path.size() ; i < is ; ++i ) {
    if( atom_ring_systs[bond_path[i]->GetIdx()] ) {
      int num_in_syst = 0;
      for( size_t j = 0 ; j < is ; ++j ) {
        if( atom_ring_systs[bond_path[j]->GetIdx()] ==
            atom_ring_systs[bond_path[i]->GetIdx()] ) {
          ++num_in_syst;
        }
      }
      if( num_in_syst > 2 ) {
        return false;
      }
    }
  }

  return true;

}

// ****************************************************************************
// rule 5 - a carbon had must have a bond path of length 2 to at least one
// heteroatom had.
// There's another possibility, though, shown in one of the tautomers of
// CHEMBL19253.  If we have a case CC(=X)C where X is hetero, or the hidden
// form CC(XH)=C, then a hydrogen can be passed from one carbon to the other
// via the intermediate enol form. This clearly matters in asymmetric systems.
// After wasting a lot of time trying to build rules to account for this in one
// go, it is now being dealt with by iterating round all the intermediate
// forms and expanding the t_skel as appropriate.
// Another possibility that Thalheim doesn't consider is an H shift where the
// carbon atom is 4 bonds from a hetero atom and they can do something akin to
// an internal hydrogen bond. An example is also in CHEMBL19253.  The test is
// if the C and non-C are in a bond path of length 4, and there aren't 3 atoms
// in the path in the same ring system, it's ok.
void apply_rule_5( const vector<unsigned int> &atom_ring_systs ,
                   vector<vector<OEAtomBase *> > &bond_paths ,
                   vector<OEAtomBase *> &hads , vector<int> &is_had ) {

  for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
    if( OEElemNo::C != hads[i]->GetAtomicNum() ) {
      continue;
    }
#ifdef NOTYET
    cout << "rule 5 for atom " << DACLIB::atom_index( *hads[i] ) + 1 << endl;
#endif
    bool ok_to_keep = false;
    for( size_t j = 0 , js = bond_paths.size() ; j < js ; ++j ) {
#ifdef NOTYET
      cout << "Path " << j << " :: ";
      for( size_t ii = 0 , iis = bond_paths[j].size() ; ii < iis ; ++ii ) {
        cout << " " << DACLIB::atom_index( *bond_paths[j][ii] ) + 1;
      }
      cout << endl;
#endif
      if( ( 3 != bond_paths[j].size() && 5 != bond_paths[j].size() ) ||
          ( hads[i] != bond_paths[j].front() &&
            hads[i] != bond_paths[j].back() ) ) {
        continue;
      }      
      if( 3 == bond_paths[j].size() &&
          ( OEElemNo::C != bond_paths[j].back()->GetAtomicNum() ||
          OEElemNo::C != bond_paths[j].front()->GetAtomicNum() ) ) {
#ifdef NOTYET
        cout << "passed normal, 2 bond, rule 5 with bond path :";
        for( size_t ii = 0 , iis = bond_paths[j].size() ; ii < iis ; ++ii ) {
          cout << " " << DACLIB::atom_index( *bond_paths[j][ii] ) + 1;
        }
        cout << endl;
#endif
        ok_to_keep = true;
        break;
      }
      ok_to_keep = rule_5_test_ext( atom_ring_systs , bond_paths[j] );
      if( ok_to_keep ) {
#ifdef NOTYET
        cout << "passed rule_5_test_ext with bond path :";
        for( size_t ii = 0 , iis = bond_paths[j].size() ; ii < iis ; ++ii ) {
          cout << " " << DACLIB::atom_index( *bond_paths[j][ii] ) + 1;
        }
        cout << endl;
#endif
        break;
      }
    }
    if( !ok_to_keep ) {
#ifdef NOTYET
      cout << "Removing HAD " << DACLIB::atom_index( *hads[i] ) + 1
           << " due to rule 5" << endl;
#endif
      is_had[DACLIB::atom_index( *hads[i] )] = 0;
      hads[i] = 0;
    } else {
#ifdef NOTYET
      cout << "HAD " << DACLIB::atom_index( *hads[i] ) + 1
           << " passes rule 5" << endl;
#endif
    }
  }

  hads.erase( remove( hads.begin() , hads.end() , static_cast<OEAtomBase *>( 0 ) ) ,
              hads.end() );

  // take out any paths that start or finish with a non-had, which can happen
  // if a had has been removed from the list
  remove_paths_of_non_hads( is_had , bond_paths );

}

// ****************************************************************************
// rule 6
// a heteroatom had must have a bond path of even length (all bond paths are
// even length the way I've built them) to at least one heteroatom had or
// a bond path of length 2 to at least 1 carbon had
void apply_rule_6( vector<vector<OEAtomBase *> > &bond_paths ,
                   vector<OEAtomBase *> &hads ,
                   vector<int> &is_had ) {

  for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
    if( OEElemNo::C == hads[i]->GetAtomicNum() ) {
      continue;
    }
#ifdef NOTYET
    cout << "Rule 6 for " << DACLIB::atom_index( *hads[i] ) + 1 << endl;
#endif
    bool ok_to_keep = false;
    for( size_t j = 0 , js = bond_paths.size() ; j < js ; ++j ) {
#ifdef NOTYET
      cout << "Path " << j << " :: ";
      for( size_t ii = 0 , iis = bond_paths[j].size() ; ii < iis ; ++ii ) {
        cout << " " << DACLIB::atom_index( *bond_paths[j][ii] ) + 1;
      }
      cout << endl;
#endif
      // do the carbon and length 3 case
      if( ( hads[i] == bond_paths[j].front() &&
            OEElemNo::C == bond_paths[j].back()->GetAtomicNum() &&
            3 == bond_paths[j].size() ) ||
          ( hads[i] == bond_paths[j].back() &&
            OEElemNo::C == bond_paths[j].front()->GetAtomicNum() &&
            3 == bond_paths[j].size() ) ) {
        ok_to_keep = true;
        break;
      }
      // do the hetero atom case
      if( ( hads[i] == bond_paths[j].front() &&
            OEElemNo::C != bond_paths[j].back()->GetAtomicNum() ) ||
          ( hads[i] == bond_paths[j].back() &&
            OEElemNo::C != bond_paths[j].front()->GetAtomicNum() ) ) {
        ok_to_keep = true;
        break;
      }
    }

    if( !ok_to_keep ) {
#ifdef NOTYET
      cout << "Removing HAD " << DACLIB::atom_index( *hads[i] ) + 1
           << " due to rule 6" << endl;
#endif
      is_had[DACLIB::atom_index( *hads[i] )] = 0;
      hads[i] = 0;
    }
  }

  hads.erase( remove( hads.begin() , hads.end() , static_cast<OEAtomBase *>( 0 ) ) ,
              hads.end() );

  // take out any paths that start or finish with a non-had, which can happen
  // if a had has been removed from the list
  remove_paths_of_non_hads( is_had , bond_paths );

}

// ****************************************************************************
// there are some extra rules to get rid of some of the hads already
// identified, that depend on the paths
void apply_had_pruning_rules( const vector<unsigned int> &atom_ring_systs ,
                              vector<vector<OEAtomBase *> > &bond_paths ,
                              vector<OEAtomBase *> &hads ,
                              vector<int> &is_had ) {

#ifdef NOTYET
  cout << "apply_had_pruning_rules" << endl;
#endif

  apply_rule_3( bond_paths , hads , is_had );
#ifdef NOTYET
  cout << "after rule 3 HAD list :";
  for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
    cout << " " << DACLIB::atom_index( *hads[i] ) + 1;
  }
  cout << endl;
#endif

  apply_rule_5( atom_ring_systs , bond_paths , hads , is_had );
#ifdef NOTYET
  cout << "after rule 5 HAD list :";
  for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
    cout << " " << DACLIB::atom_index( *hads[i] ) + 1;
  }
  cout << endl;
#endif

  apply_rule_6( bond_paths , hads , is_had );
#ifdef NOTYET
  cout << "after rule 6 HAD list :";
  for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
    cout << " " << DACLIB::atom_index( *hads[i] ) + 1;
  }
  cout << endl;
#endif

}

// ****************************************************************************
// actually, it's not really a single-double rule, it's a no consecutive single
// bonds rule. Consecutive multiple bonds are fine. If it's aromatic, to avoid
// problems with different kekule forms, use 5 for all aromatic bond orders.
void apply_single_double_rule( OEMolBase &mol ,
                               vector<vector<OEAtomBase *> > &bond_paths ,
                               vector<int> &keep_bond_paths ) {

  for( size_t i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    // only bother with ones still flagged as bad
    if( keep_bond_paths[i] ) {
      continue;
    }
#ifdef NOTYET
    cout << "Single-double rule for " << i << " ::";
    for( size_t ii = 0 , iis = bond_paths[i].size() ; ii < iis ; ++ii ) {
      cout << " " << DACLIB::atom_index( *bond_paths[i][ii] ) + 1;
    }
    cout << endl;
#endif
    OEBondBase *bond = mol.GetBond( bond_paths[i][0] , bond_paths[i][1] );
    unsigned int bo = bond->IsAromatic() ? 5 : bond->GetOrder();
    bool last_bond_single = (1 == bo);
    bool ok_to_keep = true;
    for( size_t j = 1 , js = bond_paths[i].size() - 1 ; j < js ; ++j ) {
      bond = mol.GetBond( bond_paths[i][j] , bond_paths[i][j+1] );
      bo = bond->IsAromatic() ? 5 : bond->GetOrder();
#ifdef NOTYET
      cout << "bond between " << DACLIB::atom_index( *bond_paths[i][j] ) + 1
           << " and " << DACLIB::atom_index( *bond_paths[i][j+1] ) + 1
          << " = " << bo << " and last_bond_single : " << last_bond_single
          << endl;
#endif
      if( last_bond_single && (1 == bo) ) {
#ifdef NOTYET
        cout << "Fails here" << endl;
#endif
        ok_to_keep = false;
        break;
      }
      last_bond_single = (1 == bo);
    }
    if( ok_to_keep ) {
#ifdef NOTYET
      cout << "Path " << i << " passes single-double rule :: ";
      for( size_t ii = 0 , iis = bond_paths[i].size() ; ii < iis ; ++ii ) {
        cout << " " << DACLIB::atom_index( *bond_paths[i][ii] ) + 1;
      }
      cout << endl;
#endif
      keep_bond_paths[i] = 1;
    } else {
#ifdef NOTYET
      cout << "Path " << i << " fails single-double rule :: ";
      for( size_t ii = 0 , iis = bond_paths[i].size() ; ii < iis ; ++ii ) {
        cout << " " << DACLIB::atom_index( *bond_paths[i][ii] ) + 1;
      }
      cout << endl;
#endif
    }
  }

}

// ****************************************************************************
void apply_conjugated_atom_rule( vector<vector<OEAtomBase *> > &bond_paths ,
                                 vector<int> &keep_bond_paths ) {

  // with the exception of the donor, which is the first or last atom, all atoms
  // have a multiple bond
  for( size_t i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    // only bother with ones still flagged as bad
    if( keep_bond_paths[i] ) {
      continue;
    }
    bool ok_to_keep = true;
    if( bond_paths[i].back()->GetTotalHCount() ) {
      for( size_t j = 0 , js = bond_paths[i].size() - 1 ; j < js ; ++j ) {
        if( !atom_has_multiple_bond( bond_paths[i][j] ) ) {
          ok_to_keep = false;
          break;
        }
      }
    } else if( bond_paths[i].front()->GetTotalHCount() ) {
      for( size_t j = 1 , js = bond_paths[i].size() ; j < js ; ++j ) {
        if( !atom_has_multiple_bond( bond_paths[i][j] ) ) {
          ok_to_keep = false;
          break;
        }
      }
    }
    if( ok_to_keep ) {
#ifdef NOTYET
      cout << "Path " << i << " passes conjugated atom rule :: ";
      for( size_t ii = 0 , iis = bond_paths[i].size() ; ii < iis ; ++ii ) {
        cout << " " << DACLIB::atom_index( *bond_paths[i][ii] ) + 1;
      }
      cout << endl;
#endif
      keep_bond_paths[i] = 1;
    }
  }

}

// ****************************************************************************
// rule 3 is that a HAD without an H-atom (acceptor) must have a bond-path
// to another HAD (donor) that can release a HAD. But the path (as opposed to
// the HAD) is only valid if the double bond to the acceptor, which makes it
// an acceptor, is in the path.  In the case C=CN=C, such as we see in
// the tautomer of CHEMBL19253 CC(C)(C)C(=C)N=C1C(C(=O)C1=NC2=CC=C(CC2O)C#N)O,
// the path C=CN isn't valid as the N can't accept the H without further
// knock-on effects.  It might be that the N is an acceptor from the other
// end, via a different path, but not this path.
void apply_rule_3_to_bond_paths( OEMolBase &mol ,
                                 const vector<vector<OEAtomBase *> > &bond_paths ,
                                 vector<int> &keep_bond_paths ) {

  for( size_t i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    if( !keep_bond_paths[i] ) {
      continue;
    }
    OEBondBase *bond = 0;
    if( !bond_paths[i].front()->GetTotalHCount() ) {
      bond = mol.GetBond( bond_paths[i][0] , bond_paths[i][1] );
    } else if( !bond_paths[i].back()->GetTotalHCount() ) {
      size_t last_ind = bond_paths[i].size() - 1;
      bond = mol.GetBond( bond_paths[i][last_ind] , bond_paths[i][last_ind-1] );
    }
    if( bond && 1 == bond->GetOrder() ) {
      keep_bond_paths[i] = 0;
    }
  }

}

// ****************************************************************************
void apply_rule_5_to_bond_paths( const vector<unsigned int> &atom_ring_systs ,
                                 const vector<vector<OEAtomBase *> > &bond_paths ,
                                 vector<int> &keep_bond_paths ) {

  for( size_t i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    if( !keep_bond_paths[i] ||
        ( OEElemNo::C != bond_paths[i].front()->GetAtomicNum() &&
          OEElemNo::C != bond_paths[i].back()->GetAtomicNum() ) ) {
      continue;
    }
    // we might reject it, but we'll see. Default to reject.
    keep_bond_paths[i] = 0;
    if( rule_5_test( bond_paths[i] ) ||
        rule_5_test_ext( atom_ring_systs , bond_paths[i] ) ) {
      keep_bond_paths[i] = 1;
    }
    // CHEMBL19253 showed that we need to knock potential cyclic allenes out
    // early.  So, if the carbon atom has an unsaturated bond to an atom
    // not in the path, reject the path.
    if( keep_bond_paths[i] ) {
      // if we've rejected it, it's failed a more general test
      OEAtomBase *c_atom = 0 , *next_atom = 0;
      if( OEElemNo::C == bond_paths[i].front()->GetAtomicNum() ) {
        c_atom = bond_paths[i][0];
        next_atom = bond_paths[i][1];
      } else {
        c_atom = bond_paths[i][bond_paths[i].size()-1];
        next_atom = bond_paths[i][bond_paths[i].size()-2];
      }
      OEAtomBase *unsat_nbour = atom_has_multiple_bond( c_atom );
      if( unsat_nbour && unsat_nbour != next_atom && c_atom->IsInRing() &&
          next_atom->IsInRing() && unsat_nbour->IsInRing() ) {
        keep_bond_paths[i] = 0;
      }
    }

#ifdef NOTYET
    if( !keep_bond_paths[i] ) {
      cout << "Path " << i << " :: " << keep_bond_paths[i] << " :: ";
      for( size_t ii = 0 , iis = bond_paths[i].size() ; ii < iis ; ++ii ) {
        cout << " " << DACLIB::atom_index( *bond_paths[i][ii] ) + 1;
      }
      cout << " fails rule 5" << endl;
    } else {
      cout << "Path " << i << " :: " << keep_bond_paths[i] << " :: ";
      for( size_t ii = 0 , iis = bond_paths[i].size() ; ii < iis ; ++ii ) {
        cout << " " << DACLIB::atom_index( *bond_paths[i][ii] ) + 1;
      }
      cout << " passes rule 5" << endl;
    }
#endif
  }

}

// ****************************************************************************
// as well as calculating the abes, put together a list of atoms that have
// H atoms on them that need to come off in the t_skel.
void calc_abes( OEMolBase &mol , const vector<OEAtomBase *> &hads ,
                const vector<vector<OEAtomBase *> > &bond_paths ,
                vector<int> &abes , int &num_act_h ,
                vector<int> &mobile_h ) {

  // get flags for all hads in paths
  vector<int> in_bond_path( abes.size() , 0 );
  for( size_t i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    for( size_t j = 0 , js = bond_paths[i].size() ; j < js ; ++j ) {
      in_bond_path[DACLIB::atom_index( (*bond_paths[i][j] ) )] = 1;
    }
  }

  // abes for mobile hydrogens
  for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
    if( !hads[i]->GetTotalHCount() ) {
      // it must have an H to be a donor!
      continue;
    }
    // if it's a carbon atom and it has a non-single bond, it's an H atom
    // acceptor but not a donor no matter how many H atoms it holds, unless
    // the other end of the double bond isn't a had or tsa i.e. in a bond path
    if( OEElemNo::C == hads[i]->GetAtomicNum() ) {
      OEAtomBase *nb = atom_has_multiple_bond( hads[i] );
      if( nb && in_bond_path[DACLIB::atom_index( *nb )] ) {
        continue;
      }
    }
    // also, if it's a hetero atom and it has a non-single bond, it's only
    // a donor if it has a neighbour that has another non-single bond
    bool is_donor = true;
    if( OEElemNo::C != hads[i]->GetAtomicNum() &&
        atom_has_multiple_bond( hads[i] ) ) {
      is_donor = false;
      for( OEIter<OEBondBase> bond = hads[i]->GetBonds( OENot<OEBondBase>( OEHasOrder( 1 ) ) ) ;
           bond ; ++bond ) {
        OEAtomBase *nb = bond->GetBgn() == hads[i] ? bond->GetEnd() : bond->GetBgn();
        for( OEIter<OEBondBase> nb_bond = nb->GetBonds( OENot<OEBondBase>( OEHasOrder( 1 ) ) ) ;
             nb_bond ; ++nb_bond ) {
          OEAtomBase *nb_nb = nb_bond->GetBgn() == nb ? nb_bond->GetEnd() : nb_bond->GetBgn();

          if( nb_nb != hads[i] ) {
            is_donor = true;
            break;
          }
        }
      }
    }
    if( !is_donor ) {
      continue;
    }

    ++abes[DACLIB::atom_index( *hads[i] )];
    mobile_h[DACLIB::atom_index( *hads[i] )] = 1;

  }
  num_act_h = static_cast<int>( count_if( abes.begin() , abes.end() ,
                                          bind( greater<int>() , _1 , 0 ) ) );

  // abes due to unsaturated bonds
  vector<int> done_bond( mol.GetMaxBondIdx() , 0 );
  for( size_t i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    for( size_t j = 0 , js = bond_paths[i].size() - 1 ; j < js ; ++j ) {
      OEBondBase *bond = mol.GetBond( bond_paths[i][j] , bond_paths[i][j+1] );
      if( bond && !done_bond[bond->GetIdx()] ) {
        int bo = bond->GetOrder();
        if( bo > 1 ) {
          abes[DACLIB::atom_index( *bond->GetBgn() )] += bo - 1;
          abes[DACLIB::atom_index( *bond->GetEnd() )] += bo - 1;
          done_bond[bond->GetIdx()] = bo;
        }
      }
    }
  }

#ifdef NOTYET
  cout << "GLOBAL ABES :";
  for( size_t ii = 0 , iis = abes.size() ; ii < iis ; ++ii ) {
    cout << " " << abes[ii];
  }
  cout << endl;
  cout << "Mobile H :";
  for( size_t ii = 0 , iis = mobile_h.size() ; ii < iis ; ++ii ) {
    cout << " " << mobile_h[ii];
  }
  cout << endl;
  cout << "num_act_h = " << num_act_h << endl;
#endif

}

// ****************************************************************************
// check if this is an amide carbon, the one in the middle of NC(=O)C. Also
// count imidic acids : N=C(O)C.
bool is_it_amide_c( OEMolBase &mol , OEAtomBase *atom ) {

#ifdef NOTYET
  cout << "is_it_amide_c for " << DACLIB::atom_index( *atom ) + 1 << endl;
#endif
  // it's not an amide if it's aromatic - it's a pyridone
  if( atom->IsAromatic() ) {
    return false;
  }

  OEIter<OEAtomBase> c_atom = atom->GetAtoms( OEHasAtomicNum( OEElemNo::C ) );
  if( !c_atom ) {
    return false;
  }
  OEIter<OEAtomBase> n_atom = atom->GetAtoms( OEHasAtomicNum( OEElemNo::N ) );
  if( !n_atom ) {
    return false;
  }
  OEIter<OEAtomBase> o_atom = atom->GetAtoms( OEHasAtomicNum( OEElemNo::O ) );
  // O atoms that are connected to 2 heavy atoms are ok.
  if( !o_atom || o_atom->GetHvyDegree() > 1 ) {
    return false;
  }

  OEBondBase *bond1 = mol.GetBond( o_atom , atom );
  OEBondBase *bond2 = mol.GetBond( n_atom , atom );

  if( ( bond1->GetOrder() > 1 && bond2->GetOrder() == 1 ) ||
      ( bond1->GetOrder() == 1 && bond2->GetOrder() > 1 ) ) {
#ifdef NOTYET
    cout << "returns true" << endl;
#endif
    return true;
  }

#ifdef NOTYET
  cout << "returns false def" << endl;
#endif
  return false;

}

// ****************************************************************************
// bond_path should be 3 atoms long, starting with an O or N atom. See if it's
// an amide.  It might be an imidic acid, as well, so treat that similarly.
// We also don't want to move the C=N bond in an imidic acid: N=C(O)C to
// NC(O)=C.
bool is_amide_path( OEMolBase &mol , vector<OEAtomBase *> &bond_path ) {

#ifdef NOTYET
  cout << "is_amide_path for : " << DACLIB::atom_index( *bond_path[0] ) + 1
      << " (" << bond_path[0]->GetAtomicNum() << ")"
      << " " << DACLIB::atom_index( *bond_path[1] ) + 1
      << " (" << bond_path[1]->GetAtomicNum() << ")"
      << " " << DACLIB::atom_index( *bond_path[2] ) + 1
      << " (" << bond_path[2]->GetAtomicNum() << ")"
      << endl;
#endif
  // it's an amide if the bond between 0 and 1 is a double, 1 is a carbon
  // and attached to 2 a carbon and an N not in bond_path
  if( OEElemNo:: C != bond_path[1]->GetAtomicNum() ||
      OEElemNo:: C != bond_path[2]->GetAtomicNum() ) {
    return false;
  }

  return is_it_amide_c( mol , bond_path[1] );

}

// ****************************************************************************
// check if this is an acid carbon, the one in the middle of OC(=O)C. As
// warfarin (CHEMBL1464) shows, we should ignore aromatic C and allow cases
// where the non-carbonyl C is unsaturated (which in warfarin amounts to the
// same rule).
bool is_it_acid_c( OEMolBase &mol , OEAtomBase *atom ) {

#ifdef NOTYET
  cout << "is_it_acid_c for " << DACLIB::atom_index( *atom ) + 1 << endl;
#endif

  // it's not if it's aromatic
  if( atom->IsAromatic() ) {
    return false;
  }
  OEIter<OEAtomBase> c_atom = atom->GetAtoms( OEHasAtomicNum( OEElemNo::C ) );
  if( !c_atom || atom_has_multiple_bond( c_atom ) ) {
    return false;
  }
  OEIter<OEAtomBase> o_atoms = atom->GetAtoms( OEHasAtomicNum( OEElemNo::O ) );
  if( !o_atoms ) {
    return false;
  }

  OEAtomBase *o_atom1 = o_atoms;
  ++o_atoms;
  if( !o_atoms ) {
    return false;
  }
  OEAtomBase *o_atom2 = o_atoms;

  OEBondBase *bond1 = mol.GetBond( o_atom1 , atom );
  OEBondBase *bond2 = mol.GetBond( o_atom2 , atom );

  if( ( bond1->GetOrder() > 1 && bond2->GetOrder() == 1 ) ||
      ( bond1->GetOrder() == 1 && bond2->GetOrder() > 1 ) ) {
#ifdef NOTYET
    cout << "returns true" << endl;
#endif
    return true;
  }

#ifdef NOTYET
  cout << "returns false" << endl;
#endif
  return false;

}

// ****************************************************************************
// bond_path should be 3 atoms long, starting with an O or N (because it comes
// from the same code as is_amide_path) atom. See if it's an
// acid. It might also be O=CO which is clearly pointless.
bool is_acid_path( OEMolBase &mol , vector<OEAtomBase *> &bond_path ) {

#ifdef NOTYET
  cout << "is_acid_path for : " << DACLIB::atom_index( *bond_path[0] ) + 1
      << " (" << bond_path[0]->GetAtomicNum() << ")"
      << " " << DACLIB::atom_index( *bond_path[1] ) + 1
      << " (" << bond_path[1]->GetAtomicNum() << ")"
      << " " << DACLIB::atom_index( *bond_path[2] ) + 1
      << " (" << bond_path[2]->GetAtomicNum() << ")"
      << endl;
#endif
  // OCO might be a poly-ether
  if( OEElemNo:: O == bond_path[0]->GetAtomicNum() &&
      OEElemNo:: C == bond_path[1]->GetAtomicNum() &&
      OEElemNo:: O == bond_path[2]->GetAtomicNum() ) {
    return is_it_acid_c( mol , bond_path[1] );
  }

  // it's an acid if the bond between 0 and 1 is a double, 1 is a carbon
  // and attached to 2 a carbon and an O not in bond_path
  if( OEElemNo:: O != bond_path[0]->GetAtomicNum() ||
      OEElemNo:: C != bond_path[1]->GetAtomicNum() ||
      OEElemNo:: C != bond_path[2]->GetAtomicNum() ) {
    return false;
  }

  return is_it_acid_c( mol , bond_path[1] );

}

// ****************************************************************************
// this one is optional, but removes any paths that are amides where they are
// not in another path and isolated by carbons (that latter relaxation is due
// to CHEMBL155287). It's so peptides don't take forever.
void apply_ignore_amides_rule( OEMolBase &mol ,
                               vector<vector<OEAtomBase *> > &bond_paths ,
                               vector<int> &keep_bond_paths ) {

  for( size_t i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    if( !keep_bond_paths[i] ) {
      continue; // it's already set to be hosed
    }
    if( 3 == bond_paths[i].size() ) {
      if( is_it_amide_c( mol , bond_paths[i][1] ) ) {
        // if either end atom is attached to a Het, we'll let it go. The
        // O is unlikely to be, but NUTPOACTSYO.
        if( atom_has_het_nbour( bond_paths[i][0] ) ||
            atom_has_het_nbour( bond_paths[i][2] ) ) {
          continue;
        }
        // see if this atom pops up anywhere other than the middle of another
        // 3-atom path (which is most likely the reverse of this path)
        bool found_it( false );
        for( size_t j = 0 , js = bond_paths.size() ; j < js ; ++j ) {
          if( !keep_bond_paths[j] ) {
            continue; // it's already set to be hosed
          }
          if( 3 == bond_paths[j].size() &&
              bond_paths[j][1] == bond_paths[i][1] ) {
            continue;
          }
          for( size_t k = 0 , ks = bond_paths[j].size() ; k < ks ; ++k ) {
            if( bond_paths[i][1] == bond_paths[j][k] ) {
              found_it = true;
              break;
            }
          }
          if( found_it ) {
            break;
          }
        }
        if( !found_it ) {
#ifdef NOTYET
          cout << "amide path :";
          for( size_t jj = 0 , jjs = bond_paths[i].size() ; jj < jjs ; ++jj ) {
            cout << " " << DACLIB::atom_index( *bond_paths[i][jj] ) + 1;
          }
          cout << endl;
#endif
          keep_bond_paths[i] = 0;
        }
      }
    }
  }

}

// ****************************************************************************
// Acids and amides obviously have a CC=O group which will appear as a
// possible tautomer position at the end of a path (the O should be the last
// atom) because it looks like an aldehyde if you don't look beyond the chain
// itself.  These should not be included as C=C(N)O is a bit of an abomination.
// The energy difference between it and CC(N)=O is about 30KCal/Mol according
// to a paper Martin Packer found.  C=C(O)(O) is even sillier
void apply_acid_amide_rule( OEMolBase &mol ,
                            vector<vector<OEAtomBase *> > &bond_paths ,
                            vector<int> &keep_bond_paths ) {

  for( size_t i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    if( !keep_bond_paths[i] ) {
      continue; // it's already set to be hosed
    }
#ifdef NOTYET
    cout << "apply_acid_amide_rule for path " << i << " :";
    for( size_t jj = 0 , jjs = bond_paths[i].size() ; jj < jjs ; ++jj ) {
      cout << " " << DACLIB::atom_index( *bond_paths[i][jj] ) + 1;
    }
    cout << endl;
#endif
    vector<OEAtomBase *> poss_bad_path;
    if( OEElemNo::O == bond_paths[i].front()->GetAtomicNum() ||
        OEElemNo::N == bond_paths[i].front()->GetAtomicNum() ) {
      poss_bad_path.insert( poss_bad_path.begin() , bond_paths[i].begin() ,
                              bond_paths[i].begin() + 3 );
      if( is_amide_path( mol , poss_bad_path ) ||
          is_acid_path( mol , poss_bad_path ) ) {
        keep_bond_paths[i] = 0;
      }
    }
    if( !keep_bond_paths[i] ) {
      continue; // already hosed it
    }
    if( OEElemNo::O == bond_paths[i].back()->GetAtomicNum() ||
        OEElemNo::N == bond_paths[i].back()->GetAtomicNum() ) {
      poss_bad_path.clear();
      poss_bad_path.insert( poss_bad_path.begin() , bond_paths[i].rbegin() ,
                              bond_paths[i].rbegin() + 3 );
      if( is_amide_path( mol , poss_bad_path ) ||
          is_acid_path( mol , poss_bad_path ) ) {
        keep_bond_paths[i] = 0;
      }
    }
#ifdef NOTYET
    cout << "result for path " << i << " : " << keep_bond_paths[i] << endl;
#endif
  }

}

// ****************************************************************************
// we've got all the possible bond paths between acceptors and donors, but some
// won't be allowed by the rules, so take those out. They're an eclectic set of
// rules.
// According to TT,:
// A bond path is an even number of bonds, alternating between single and
// double bonds. Or, with the exception of the donor, a conjugation of
// sp2-hybridised atoms. Or there can be cumulated double bonds.
// And in all of this, a double bond can be a triple bond as well.
// It is possible that this is a bit simplistic as it doesn't account for
// intermediate tautomers.  After a lot of mucking about, I've decided to deal
// with this by generating the t_skel in an iterative manner via the
// intermediates rather than trying to build rules that take account of them
// in one go.
// We've already ensured that it's an even number of bonds, at least 2.
void prune_bond_paths( OEMolBase &mol , bool ignore_amides ,
                       const vector<unsigned int> &atom_ring_systs ,
                       vector<vector<OEAtomBase *> > &bond_paths ) {

  // only want to hose a bond path if there is no reason to keep it.
  vector<int> keep_bond_paths( bond_paths.size() , 0 );
  apply_single_double_rule( mol , bond_paths , keep_bond_paths );
#ifdef NOTYET
  cout << "Results of apply single bond rule" << endl;
  for( size_t jj = 0 , jjs = bond_paths.size() ; jj < jjs ; ++jj ) {
    cout << "Path " << jj << " :: " << keep_bond_paths[jj] << " :: ";
    for( size_t ii = 0 , iis = bond_paths[jj].size() ; ii < iis ; ++ii ) {
      cout << " " << DACLIB::atom_index( *bond_paths[jj][ii] ) + 1;
    }
    cout << endl;
  }
#endif

  apply_conjugated_atom_rule( bond_paths , keep_bond_paths );
#ifdef NOTYET
  cout << "Results of apply conjugated atom rule" << endl;
  for( size_t jj = 0 , jjs = bond_paths.size() ; jj < jjs ; ++jj ) {
    cout << "Path " << jj << " :: " << keep_bond_paths[jj] << " :: ";
    for( size_t ii = 0 , iis = bond_paths[jj].size() ; ii < iis ; ++ii ) {
      cout << " " << DACLIB::atom_index( *bond_paths[jj][ii] ) + 1;
    }
    cout << endl;
  }
#endif

  apply_rule_3_to_bond_paths( mol , bond_paths , keep_bond_paths );
#ifdef NOTYET
  cout << "Results of rule 3 to bond paths" << endl;
  for( size_t jj = 0 , jjs = bond_paths.size() ; jj < jjs ; ++jj ) {
    cout << "Path " << jj << " :: " << keep_bond_paths[jj] << " :: ";
    for( size_t ii = 0 , iis = bond_paths[jj].size() ; ii < iis ; ++ii ) {
      cout << " " << DACLIB::atom_index( *bond_paths[jj][ii] ) + 1;
    }
    cout << endl;
  }
#endif

  apply_rule_5_to_bond_paths( atom_ring_systs , bond_paths , keep_bond_paths );
#ifdef NOTYET
  cout << "Results of rule 5 to bond paths" << endl;
  for( size_t jj = 0 , jjs = bond_paths.size() ; jj < jjs ; ++jj ) {
    cout << "Path " << jj << " :: " << keep_bond_paths[jj] << " :: ";
    for( size_t ii = 0 , iis = bond_paths[jj].size() ; ii < iis ; ++ii ) {
      cout << " " << DACLIB::atom_index( *bond_paths[jj][ii] ) + 1;
    }
    cout << endl;
  }
#endif

  if( ignore_amides ) {
    apply_ignore_amides_rule( mol , bond_paths , keep_bond_paths );
#ifdef NOTYET
    cout << "Results of ignore amide rule to bond paths" << endl;
    for( size_t jj = 0 , jjs = bond_paths.size() ; jj < jjs ; ++jj ) {
      cout << "Path " << jj << " :: " << keep_bond_paths[jj] << " :: ";
      for( size_t ii = 0 , iis = bond_paths[jj].size() ; ii < iis ; ++ii ) {
        cout << " " << DACLIB::atom_index( *bond_paths[jj][ii] ) + 1;
      }
      cout << endl;
    }
#endif
  }

  apply_acid_amide_rule( mol , bond_paths , keep_bond_paths );
#ifdef NOTYET
  cout << "Results of acid/amide rule to bond paths" << endl;
  for( size_t jj = 0 , jjs = bond_paths.size() ; jj < jjs ; ++jj ) {
    cout << "Path " << jj << " :: " << keep_bond_paths[jj] << " :: ";
    for( size_t ii = 0 , iis = bond_paths[jj].size() ; ii < iis ; ++ii ) {
      cout << " " << DACLIB::atom_index( *bond_paths[jj][ii] ) + 1;
    }
    cout << endl;
  }
#endif

  // do the prune
  for( size_t i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    if( !keep_bond_paths[i] ) {
      bond_paths[i].clear();
    }
  }
  bond_paths.erase( remove_if( bond_paths.begin() , bond_paths.end() ,
                               bind( &vector<OEAtomBase *>::empty , _1 ) ) ,
                    bond_paths.end() );

#ifdef NOTYET
  cout << "Pruned bond paths" << endl;
  for( size_t jj = 0 , jjs = bond_paths.size() ; jj < jjs ; ++jj ) {
    cout << "Path " << jj << " :: " << keep_bond_paths[jj] << " :: ";
    for( size_t ii = 0 , iis = bond_paths[jj].size() ; ii < iis ; ++ii ) {
      cout << " " << DACLIB::atom_index( *bond_paths[jj][ii] ) + 1;
    }
    cout << endl;
  }
#endif

}

// ****************************************************************************
void remove_hads_not_at_path_end( const vector<vector<OEAtomBase *> > &bond_paths ,
                                  vector<OEAtomBase *> &hads ,
                                  vector<unsigned int> &had_idxs ,
                                  vector<int> &is_had ) {

  for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
    bool had_found = false;
    for( size_t j = 0 , js = bond_paths.size() ; j < js ; ++j ) {
      if( bond_paths[j].front() == hads[i] || bond_paths[j].back() == hads[i] ) {
        had_found = true;
        break;
      }
    }
    if( !had_found ) {
#ifdef NOTYET
      cout << "Removing HAD " << DACLIB::atom_index( *hads[i] ) + 1
           << " as not at start or finish of any path." << endl;
#endif
      is_had[DACLIB::atom_index( *hads[i] )] = 0;
      hads[i] = 0;
      had_idxs[i] = numeric_limits<unsigned int>::max();
    }
  }

  hads.erase( remove( hads.begin() , hads.end() ,
                      static_cast<OEAtomBase *>( 0 ) ) , hads.end() );
  had_idxs.erase( remove( had_idxs.begin() , had_idxs.end() ,
                          numeric_limits<unsigned int>::max() ) , had_idxs.end() );

}

// ****************************************************************************
void remove_hads_not_in_paths( const vector<vector<OEAtomBase *> > &bond_paths ,
                               vector<OEAtomBase *> &hads ,
                               vector<unsigned int> &had_idxs ,
                               vector<int> &is_had ) {

  for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
    bool had_found = false;
    for( size_t j = 0 , js = bond_paths.size() ; j < js ; ++j ) {
      if( bond_paths[j].end() != find( bond_paths[j].begin() , bond_paths[j].end() , hads[i] ) ) {
        had_found = true;
        break;
      }
    }
    if( !had_found ) {
#ifdef NOTYET
      cout << "Removing HAD " << DACLIB::atom_index( *hads[i] ) + 1
           << " as not at start or finish of any path." << endl;
#endif
      is_had[DACLIB::atom_index( *hads[i] )] = 0;
      hads[i] = 0;
      had_idxs[i] = numeric_limits<unsigned int>::max();
    }
  }

  hads.erase( remove( hads.begin() , hads.end() ,
                      static_cast<OEAtomBase *>( 0 ) ) , hads.end() );
  had_idxs.erase( remove( had_idxs.begin() , had_idxs.end() ,
                          numeric_limits<unsigned int>::max() ) , had_idxs.end() );

}

// ****************************************************************************
void sort_hads_by_idx( vector<OEAtomBase *> &hads ,
                       vector<unsigned int> &had_idxs ) {

  vector<pair<OEAtomBase *,unsigned int> > had_pairs;
  had_pairs.reserve( hads.size() );
  for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
    had_pairs.push_back( make_pair( hads[i] , had_idxs[i] ) );
  }

  sort( had_pairs.begin() , had_pairs.end() ,
        bind( less<unsigned int>() ,
              bind( &pair<OEAtomBase *,unsigned int>::second , _1 ) ,
              bind( &pair<OEAtomBase *,unsigned int>::second , _2 ) ) );
  hads.clear();
  had_idxs.clear();
  for( size_t i = 0 , is = had_pairs.size() ; i < is ; ++i ) {
    hads.push_back( had_pairs[i].first );
    had_idxs.push_back( had_pairs[i].second );
  }

}

// ****************************************************************************
// HAD is an H-atom acceptor atom or H-atom donor atom. Not to be confused with
// hydrogen-bond donors and acceptors. See Thalheim et al p4.
void find_hads( OEMolBase &mol , bool ignore_amides ,
                vector<vector<OEAtomBase *> > &bond_paths ,
                vector<OEAtomBase *> &hads , vector<unsigned int> &had_idxs ,
                vector<int> &is_had ) {

  build_initial_had_list( mol , hads , is_had );
#ifdef NOTYET
  cout << "initial HAD list :";
  for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
    cout << " " << DACLIB::atom_index( *hads[i] ) + 1;
  }
  cout << endl;
#endif

  // get the number of the ring system each atom is in, which we'll need for
  // the extension of rule 5 to non-cyclic 1,5 shifts.
  vector<unsigned int> atom_ring_systs( mol.GetMaxAtomIdx() , 0 );
  OEDetermineRingSystems( mol , &atom_ring_systs[0] );
#ifdef NOTYET
  int num_ring_systs = OEDetermineRingSystems( mol , &atom_ring_systs[0] );
  cout << "num ring systems : " << num_ring_systs << endl;
#endif

  // A bond path is an even number of bonds, alternating between single and
  // double bonds. Or, with the exception of the donor, a conjugation of
  // sp2-hybridised atoms. Or there can be cumulated double bonds.
  // And in all of this, a double bond can be a triple bond as well.
  // first off, just get all the paths between hads, then
  // prune afterwards
  build_bond_paths( hads , is_had , bond_paths );

  // for the later rules, it's best if the bond paths are in descending order
  // of size
  stable_sort( bond_paths.begin() , bond_paths.end() ,
               bind( greater<size_t>() ,
                     bind( &vector<OEAtomBase *>::size , _1 ) ,
                     bind( &vector<OEAtomBase *>::size , _2 ) ) );

#ifdef NOTYET
  cout << "All initial paths" << endl;
  for( size_t jj = 0 , jjs = bond_paths.size() ; jj < jjs ; ++jj ) {
    cout << "Path " << jj << " :: ";
    for( size_t ii = 0 , iis = bond_paths[jj].size() ; ii < iis ; ++ii ) {
      cout << " " << DACLIB::atom_index( *bond_paths[jj][ii] ) + 1;
    }
    cout << endl;
  }
#endif

  // iterate round the pruning until we've reached a stable state. Rule 3,
  // for example, might pass atoms because there's a path to an atom that
  // is taken out by Rule 5, so Rule 3 comes into play next time round.
  // CHEMBL3306810 pulled this one up.
  had_idxs.clear();
  for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
    had_idxs.push_back( DACLIB::atom_index( *hads[i] ) );
  }
  sort_hads_by_idx( hads , had_idxs );

  while( true ) {
#ifdef NOTYET
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "Next Round" << endl;
#endif
    vector<OEAtomBase *> curr_hads = hads;
    prune_bond_paths( mol , ignore_amides , atom_ring_systs , bond_paths );

#ifdef NOTYET
    cout << "after prune_bond_paths HAD list :";
    for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
      cout << " " << DACLIB::atom_index( *hads[i] ) + 1;
    }
    cout << endl;
    cout << "intermediate bond paths" << endl;
    for( size_t jj = 0 , jjs = bond_paths.size() ; jj < jjs ; ++jj ) {
      cout << "Path " << jj << " :: ";
      for( size_t ii = 0 , iis = bond_paths[jj].size() ; ii < iis ; ++ii ) {
        cout << " " << DACLIB::atom_index( *bond_paths[jj][ii] ) + 1;
      }
      cout << endl;
    }
#endif

    // there are some extra rules to get rid of some of the hads already
    // identified, that depend on the paths
    apply_had_pruning_rules( atom_ring_systs , bond_paths , hads , is_had );

#ifdef NOTYET
    cout << "after had pruning rules HAD list :";
    for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
      cout << " " << DACLIB::atom_index( *hads[i] ) + 1;
    }
    cout << endl;
    cout << "after pruning rules bond paths" << endl;
    for( size_t jj = 0 , jjs = bond_paths.size() ; jj < jjs ; ++jj ) {
      cout << "Path " << jj << " :: ";
      for( size_t ii = 0 , iis = bond_paths[jj].size() ; ii < iis ; ++ii ) {
        cout << " " << DACLIB::atom_index( *bond_paths[jj][ii] ) + 1;
      }
      cout << endl;
    }
#endif

    // and if, at the end of all this, we have HADS that aren't in a bond path
    // they need to go. Create had_idxs before this, as remove_hads_not_in_paths
    // used elsewhere as well.
    had_idxs.clear();
    for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
      had_idxs.push_back( DACLIB::atom_index( *hads[i] ) );
    }
    remove_hads_not_at_path_end( bond_paths , hads , had_idxs , is_had );

    sort_hads_by_idx( hads , had_idxs );

#ifdef NOTYET
    cout << "final HAD list :";
    for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
      cout << " " << DACLIB::atom_index( *hads[i] ) + 1;
    }
    cout << endl;
    cout << "Final bond paths" << endl;
    for( size_t jj = 0 , jjs = bond_paths.size() ; jj < jjs ; ++jj ) {
      cout << "Path " << jj << " :: ";
      for( size_t ii = 0 , iis = bond_paths[jj].size() ; ii < iis ; ++ii ) {
        cout << " " << DACLIB::atom_index( *bond_paths[jj][ii] ) + 1;
      }
      cout << endl;
    }
#endif
    // If nothing changed, we're done. Otherwise, go round again.
    if( curr_hads == hads ) {
      break;
    }
  }

}

// ****************************************************************************
// see if there is an atom with a bonding electron (abes is non-zero)
// that doesn't have neighbour with a bonding electron
bool are_there_isolated_atoms( const vector<int> &abes ,
                               const vector<vector<unsigned int> > &nb_idxs ) {

  for( size_t i = 0 , is = abes.size() ; i < is ; ++i ) {
    // if this atom has no bonding electrons, it doesn't matter
    if( !abes[i] ) {
      continue;
    }
#ifdef NOTYET
    cout << "looking at atom " << i + 1 << " abe = " << abes[i] << endl;
#endif
    bool abe_nb = false;
    for( size_t j = 0 , js = nb_idxs[i].size() ; j < js ; ++j ) {
      if( abes[nb_idxs[i][j]] ) {
        abe_nb = true;
        break;
      }
    }
    if( !abe_nb ) {
#ifdef NOTYET
      cout << "isolated : " << i + 1 << endl;
#endif
      return true;
    }
  }

  return false;

}

// ****************************************************************************
// see if there is an atom with a bonding electron (abes is non-zero)
// that doesn't have neighbour with a bonding electron
bool are_there_isolated_atoms( const vector<int> &abes ,
                               const vector<vector<unsigned int> > &nb_idxs ,
                               unsigned int &isolated_idx ) {

  for( unsigned int i = 0 , is = static_cast<unsigned int>( abes.size() ) ; i < is ; ++i ) {
    // if this atom has no bonding electrons, it doesn't matter
    if( !abes[i] ) {
      continue;
    }
#ifdef NOTYET
    cout << "looking at atom " << i + 1 << " abe = " << abes[i] << endl;
#endif
    bool abe_nb = false;
    for( size_t j = 0 , js = nb_idxs[i].size() ; j < js ; ++j ) {
      if( abes[nb_idxs[i][j]] ) {
        abe_nb = true;
        break;
      }
    }
    if( !abe_nb ) {
#ifdef NOTYET
      cout << "isolated : " << i + 1 << endl;
#endif
      isolated_idx = i; // i doesn't have a neighbour that's has an abe

      return true;
    }
  }

  return false;

}

// ****************************************************************************
// find an atom which has an abe that has the minimum number of neighbours
// with abes
OEAtomBase *find_unsat_start_atom( const vector<int> &abes ,
                                   OEMolBase &t_skel_mol ) {

  OEAtomBase *ret_val = static_cast<OEAtomBase *>( 0 );
  unsigned int min_count = numeric_limits<unsigned int>::max();

  for( unsigned int i = 0 , is = static_cast<unsigned int>( abes.size() ) ; i < is ; ++i ) {
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
void build_connect_sets( const vector<int> &abes , OEMolBase &t_skel_mol ,
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
  for( size_t i = 0 , is = connect_sets.size() ; i < is ; ++i ) {
    cout << i << " :";
    for( size_t j = 0 , js = connect_sets[i].size() ; j < js ; ++j ) {
      cout << " " << DACLIB::atom_index( *connect_sets[i][j] ) + 1;
    }
    cout << endl;
  }
  cout << "abe_sets" << endl;
  for( size_t i = 0 , is = abe_sets.size() ; i < is ; ++i ) {
    cout << i << " : ";
    copy( abe_sets[i].begin() , abe_sets[i].end() , intOut );
    cout << endl;
  }
#endif

}

// ****************************************************************************
// get all the bond pairs except bad_bond_atom_pairs and sort into descending
// order of number of n'bours.
void build_bond_atom_pairs( OEMolBase &mol ,
                            const vector<vector<unsigned int> > &nb_idxs ,
                            const vector<pair<unsigned int,unsigned int> > &bad_bond_atom_pairs ,
                            vector<pair<OEAtomBase *,OEAtomBase *> > &bond_atom_pairs ,
                            vector<pair<unsigned int,unsigned int> > &bond_atom_pairs_idxs ) {

  for( unsigned int i = 0 , is = static_cast<unsigned int>( nb_idxs.size() ) ; i < is ; ++i ) {
    OEAtomBase *atom1 = mol.GetAtom( DACLIB::HasAtomIndex( i ) );
    for( size_t j = 0 , js = nb_idxs[i].size() ; j < js ; ++j ) {
      OEAtomBase *atom2 = mol.GetAtom( DACLIB::HasAtomIndex( nb_idxs[i][j] ) );
      // only want each bond in once
      if( i < nb_idxs[i][j] ) {
        pair<unsigned int,unsigned int> this_pair( make_pair( i , nb_idxs[i][j] ) );
        if( bad_bond_atom_pairs.end() != find( bad_bond_atom_pairs.begin() ,
                                               bad_bond_atom_pairs.end() ,
                                               this_pair ) ) {
#ifdef NOTYET
	  cout << "bad bond pair : " << i + 1 << " to " << nb_idxs[i][j] + 1 << endl;
#endif
          continue;
        }
        bond_atom_pairs_idxs.push_back( this_pair );
        bond_atom_pairs.push_back( make_pair( atom1 , atom2 ) );
      }
    }
  }

  vector<pair<int,unsigned int> > tmp;
  for( size_t i = 0 , is = bond_atom_pairs_idxs.size() ; i < is ; ++i ) {
    unsigned int max_conn = static_cast<unsigned int>( max( nb_idxs[bond_atom_pairs_idxs[i].first].size() ,
                                                       nb_idxs[bond_atom_pairs_idxs[i].second].size() ) );
    tmp.push_back( make_pair( i , max_conn ) );
  }

  stable_sort( tmp.begin() , tmp.end() ,
               bind( greater<unsigned int>() ,
                     bind( &pair<int,unsigned int>::second , _1 ) ,
                     bind( &pair<int,unsigned int>::second , _2 ) ) );

  vector<pair<OEAtomBase *,OEAtomBase *> > tmp_bond_atom_pairs;
  vector<pair<unsigned int,unsigned int> > tmp_bond_atom_pairs_idxs;
  tmp_bond_atom_pairs.reserve( bond_atom_pairs.size() );
  tmp_bond_atom_pairs_idxs.reserve( bond_atom_pairs_idxs.size() );
  for( size_t i = 0 , is = tmp.size() ; i < is ; ++i ) {
    tmp_bond_atom_pairs.push_back( bond_atom_pairs[tmp[i].first] );
    tmp_bond_atom_pairs_idxs.push_back( bond_atom_pairs_idxs[tmp[i].first] );
  }

  bond_atom_pairs = tmp_bond_atom_pairs;
  bond_atom_pairs_idxs = tmp_bond_atom_pairs_idxs;

#ifdef NOTYET
  cout << "bond_atom_pairs" << endl;
  for( size_t ii = 0 , iis = bond_atom_pairs_idxs.size() ; ii < iis ; ++ii ) {
    cout << bond_atom_pairs_idxs[ii].first + 1 << " to " << bond_atom_pairs_idxs[ii].second + 1 << endl;
  }
#endif

}

// ****************************************************************************
void extend_h_atom_combo( unsigned int conn_set_num_h ,
                          unsigned int num_poss_h_idxs ,
                          const vector<unsigned int> &poss_h_idxs ,
                          const vector<vector<unsigned int> > &nb_idxs ,
                          const vector<int> &abe_set ,
                          const vector<pair<unsigned int,unsigned int> > &bad_h_atom_pairs ,
                          vector<int> &curr_comb ,
                          vector<vector<unsigned int> > &h_atom_idxs ) {

  if( curr_comb.size() == conn_set_num_h ) {
    return;
  }

  unsigned int start_i = curr_comb.empty() ? 0 : curr_comb.back() + 1;
  for( unsigned int i = start_i ; i < num_poss_h_idxs ; ++i ) {

    unsigned int next_at = poss_h_idxs[i];

    // check whether this is a pair we've decided not to add H atoms to
    // both of
    bool duff_pair( false );
    for( size_t j = 0 , js = bad_h_atom_pairs.size() ; j < js ; ++j ) {
      unsigned int other_end = numeric_limits<unsigned int>::max();
      if( next_at == bad_h_atom_pairs[j].first ) {
        other_end = bad_h_atom_pairs[j].second;
      } else if( next_at == bad_h_atom_pairs[j].second ) {
        other_end = bad_h_atom_pairs[j].first;
      }
      for( size_t k = 0 , ks = curr_comb.size() ; k < ks ; ++k ) {
        if( poss_h_idxs[curr_comb[k]] == other_end ) {
#ifdef NOTYET
          cout << "rejecting " << next_at + 1 << " and " << poss_h_idxs[curr_comb[k]] + 1 << " as duff pair" << endl;
#endif
          duff_pair = true;
          break;
        }
      }
      if( duff_pair ) {
        break;
      }
    }
    if( duff_pair ) {
      continue;
    }

    curr_comb.push_back( i );
#ifdef NOTYET
    cout << i << " :";
    for( size_t ii = 0 ; ii < curr_comb.size() ; ++ii ) {
      cout << " " << curr_comb[ii];
    }
    cout << " :: ";
    for( size_t ii = 0 ; ii < curr_comb.size() ; ++ii ) {
      cout << " " << poss_h_idxs[curr_comb[ii]] + 1;
    }
    cout << endl;
#endif
    if( curr_comb.size() == conn_set_num_h ) {
      // final check for isolated atoms and accept if not
      vector<int> this_abes = abe_set;
      for( unsigned int j = 0 ; j < conn_set_num_h ; ++j ) {
        --this_abes[poss_h_idxs[curr_comb[j]]];
      }
      if( !are_there_isolated_atoms( this_abes , nb_idxs ) ) {
#ifdef NOTYET
        cout << "Final H comb : " << h_atom_idxs.size() << " :";
        for( unsigned int ii = 0 ; ii < conn_set_num_h ; ++ii ) {
          cout << " " << curr_comb[ii];
        }
        cout << endl;
  #endif
        h_atom_idxs.push_back( vector<unsigned int>( curr_comb.size() , 0 ) );
        for( size_t j = 0 , js = curr_comb.size() ; j < js ; ++j ) {
          h_atom_idxs.back()[j] = poss_h_idxs[curr_comb[j]];
        }
      }
    }
    // check to see if this comb has produced isolated tsas. If so, this is
    // not one that can work
    vector<int> this_abes = abe_set;
    for( size_t j = 0 , js = curr_comb.size() ; j < js ; ++j ) {
      --this_abes[poss_h_idxs[curr_comb[j]]];
    }
    // because nb_idxs are sorted, isolated_idx will always be the lowest
    // one for the atom - putting an H on this atom may isolate more than
    // 1 atom
    unsigned int isolated_idx;
    if( are_there_isolated_atoms( this_abes , nb_idxs , isolated_idx ) &&
        isolated_idx < poss_h_idxs[curr_comb.back()] ) {
#ifdef NOTYET
      cout << "isolated atom - not continuing with BAD COMBO : " << isolated_idx + 1 << " vs "
           << poss_h_idxs[curr_comb.back()] + 1 << endl;
#endif
    } else {
      extend_h_atom_combo( conn_set_num_h , num_poss_h_idxs ,
                           poss_h_idxs , nb_idxs , abe_set ,
                           bad_h_atom_pairs , curr_comb , h_atom_idxs );
    }
    curr_comb.pop_back();
  }

}

// ****************************************************************************
// build the neighbour indices for each atom for speed of access.
// sort them into order for branch pruning in find_h_atom_combos.
void build_nbor_idxs( OEMolBase &mol ,
                      vector<vector<unsigned int> > &nb_idxs ) {

  for( OEIter<OEAtomBase> atom = mol.GetAtoms() ; atom ; ++atom ) {
    unsigned int at_idx = DACLIB::atom_index( *atom );
    for( OEIter<OEAtomBase> nb = atom->GetAtoms() ; nb ; ++nb ) {
      nb_idxs[at_idx].push_back( DACLIB::atom_index( *nb ) );
    }
    sort( nb_idxs[at_idx].begin() , nb_idxs[at_idx].end() );
  }

}

// ****************************************************************************
// check for bad C-S bonds
void check_bad_C_S_bonds( OEMolBase &mol , OEAtomBase *had ,
                          unsigned int had_idx , vector<int> &bad_atoms ,
                          vector<pair<unsigned int,unsigned int> > &bad_h_atom_pairs ,
                          vector<pair<unsigned int,unsigned int> > &bad_bond_atom_pairs ) {

  if( OEElemNo::S == had->GetAtomicNum() && had->GetHvyDegree() > 1 ) {
    bad_atoms[had_idx] = 1;
    for( OEIter<OEAtomBase> nbor = had->GetAtoms( OEHasAtomicNum( OEElemNo::C ) ) ; nbor ; ++nbor ) {
      if( had_idx < DACLIB::atom_index( *nbor ) ) {
        bad_h_atom_pairs.push_back( make_pair( had_idx , DACLIB::atom_index( *nbor ) ) );
      } else {
        bad_h_atom_pairs.push_back( make_pair( DACLIB::atom_index( *nbor ) , had_idx ) );
      }
      bad_bond_atom_pairs.push_back( bad_h_atom_pairs.back() );
    }
  }

  if( OEElemNo::C == had->GetAtomicNum() ) {
    for( OEIter<OEAtomBase> nbour = had->GetAtoms( OEHasAtomicNum( OEElemNo::S ) ) ; nbour ; ++nbour ) {
      if( nbour->GetHvyDegree() > 1 ) {
        unsigned int nbour_idx = DACLIB::atom_index( *nbour );
        bad_atoms[nbour_idx] = 1;
        if( had_idx < nbour_idx ) {
          bad_h_atom_pairs.push_back( make_pair( had_idx , nbour_idx ) );
        } else {
          bad_h_atom_pairs.push_back( make_pair( nbour_idx , had_idx ) );
        }
        // if the bond was unsaturated on input, clearly it's ok we just want
        // to not make silly ones ourselves.
        OEBondBase *bond = mol.GetBond( had , nbour );
        if( !bad_h_atom_pairs.empty() && bond->GetOrder() == 1 &&
            !bond->IsAromatic() ) {
          bad_bond_atom_pairs.push_back( bad_h_atom_pairs.back() );
        }
      }
    }
  }

}

// ****************************************************************************
// don't add H atoms to connected Cs. CHEMBL1082532 provided the impetus.
void check_bad_C_C_nbours( OEAtomBase *had , unsigned int had_idx ,
                           vector<pair<unsigned int,unsigned int> > &bad_h_atom_pairs ) {

  if( OEElemNo::C == had->GetAtomicNum() ) {
    for( OEIter<OEAtomBase> nbour = had->GetAtoms( OEHasAtomicNum( OEElemNo::C ) ) ; nbour ; ++nbour ) {
      if( had_idx < DACLIB::atom_index( *nbour ) ) {
        bad_h_atom_pairs.push_back( make_pair( had_idx , DACLIB::atom_index( *nbour ) ) );
      } else {
        bad_h_atom_pairs.push_back( make_pair( DACLIB::atom_index( *nbour ) , had_idx ) );
      }
    }
  }

}

// ****************************************************************************
// Don't allow a bald NH-OH, but hydroxamic acids (O=CNO) as in CHEMBL155287
// are ok.
void check_bad_N_O_nbours( OEMolBase &mol , OEAtomBase *had ,
                           unsigned int had_idx ,
                           vector<pair<unsigned int,unsigned int> > &bad_h_atom_pairs ) {

#ifdef NOTYET
  cout << "checking Bad N_O nbours for " << had_idx + 1 << " in " << DACLIB::create_noncansmi( mol ) << endl;
#endif

  if( OEElemNo::N == had->GetAtomicNum() ) {
    bool c_hits( false );
    OEAtomBase *o_nbour = 0;
    for( OEIter<OEAtomBase> nbor = had->GetAtoms() ; nbor ; ++nbor ) {
      if( OEElemNo::O == nbor->GetAtomicNum() ) {
        o_nbour = nbor;
      }
      if( OEElemNo::C == nbor->GetAtomicNum() && is_it_amide_c( mol , nbor ) ) {
        c_hits = true;
      }
    }
    if( !c_hits && o_nbour ) {
      if( had_idx < DACLIB::atom_index( *o_nbour ) ) {
        bad_h_atom_pairs.push_back( make_pair( had_idx , DACLIB::atom_index( *o_nbour ) ) );
      } else {
        bad_h_atom_pairs.push_back( make_pair( DACLIB::atom_index( *o_nbour ) , had_idx ) );
      }
    }
  }

  if( OEElemNo::O == had->GetAtomicNum() ) {
    for( OEIter<OEAtomBase> nbor = had->GetAtoms( OEHasAtomicNum( OEElemNo::N ) ) ; nbor ; ++nbor ) {
      bool amide_c( false );
      for( OEIter<OEAtomBase> n_nbour = nbor->GetAtoms( OEHasAtomicNum( OEElemNo::C ) ) ; n_nbour ; ++n_nbour ) {
        if( is_it_amide_c( mol , n_nbour ) ) {
          amide_c = true;
        }
      }
      if( !amide_c ) {
        if( had_idx < DACLIB::atom_index( *nbor ) ) {
          bad_h_atom_pairs.push_back( make_pair( had_idx , DACLIB::atom_index( *nbor ) ) );
        } else {
          bad_h_atom_pairs.push_back( make_pair( DACLIB::atom_index( *nbor ) , had_idx ) );
        }
      }
    }
  }

}

// ****************************************************************************
void check_abes_dead_end( OEAtomBase *had , unsigned int had_idx ,
                          const vector<int> &abes ,
                          vector<pair<unsigned int,unsigned int> > &bad_h_atom_pairs ) {

#ifdef NOTYET
  cout << "check_abes_dead_end" << endl;
#endif
  if( 1 == abes[had_idx] && 1 == had->GetHvyDegree() ) {
    for( OEIter<OEAtomBase> nbor = had->GetAtoms() ; nbor ; ++nbor ) {
      if( 2 == nbor->GetHvyDegree() && 1 == abes[DACLIB::atom_index( *nbor )] &&
          ( OEElemNo::C == had->GetAtomicNum() || OEElemNo::C == nbor->GetAtomicNum() ) ) {
#ifdef NOTYET
        cout << had_idx + 1 << " and " << DACLIB::atom_index( *nbor ) + 1 << " abes dead end" << endl;
#endif
        if( had_idx < DACLIB::atom_index( *nbor ) ) {
          bad_h_atom_pairs.push_back( make_pair( had_idx , DACLIB::atom_index( *nbor ) ) );
        } else {
          bad_h_atom_pairs.push_back( make_pair( DACLIB::atom_index( *nbor ) , had_idx ) );
        }
      }
    }
  }

  if( 2 == abes[had_idx] && 2 == had->GetHvyDegree() ) {
    for( OEIter<OEAtomBase> nbor = had->GetAtoms() ; nbor ; ++nbor ) {
      if( 1 == nbor->GetHvyDegree() && 1 == abes[DACLIB::atom_index( *nbor )] &&
          ( OEElemNo::C == had->GetAtomicNum() || OEElemNo::C == nbor->GetAtomicNum() ) ) {
#ifdef NOTYET
        cout << had_idx + 1 << " and " << DACLIB::atom_index( *nbor ) + 1 << " abes dead end" << endl;
#endif
        if( had_idx < DACLIB::atom_index( *nbor ) ) {
          bad_h_atom_pairs.push_back( make_pair( had_idx , DACLIB::atom_index( *nbor ) ) );
        } else {
          bad_h_atom_pairs.push_back( make_pair( DACLIB::atom_index( *nbor ) , had_idx ) );
        }
      }
    }
  }

}

// ****************************************************************************
// don't want to make C=C(O)N groups. Aromatic ones are ok, they're pyridones.
// Ones that could be pyridones are trickier.
void check_C_C_to_amide( OEMolBase &mol , unsigned int atom_idx ,
                         const vector<int> &abes ,
                         vector<pair<unsigned int,unsigned int> > &bad_bond_atom_pairs ) {

  OEAtomBase *atom = mol.GetAtom( DACLIB::HasAtomIndex( atom_idx ) );

  if( OEElemNo:: C != atom->GetAtomicNum() || atom->GetHvyDegree() < 3 ||
      atom->IsAromatic() ) {
    return;
  }

  if( abes[atom_idx] ) {
    OEIter<OEAtomBase> c_nbour = atom->GetAtoms( OEHasAtomicNum( OEElemNo::C ) );
    if( !c_nbour || !abes[DACLIB::atom_index( *c_nbour )] ) {
      return;
    }
    OEIter<OEAtomBase> n_nbour = atom->GetAtoms( OEHasAtomicNum( OEElemNo::N ) );
    if( !n_nbour ) {
      return;
    }
    OEIter<OEAtomBase> o_nbour = atom->GetAtoms( OEHasAtomicNum( OEElemNo::O ) );
    // if O is connected to more than 1 heavy atom, it's ok.
    if( !o_nbour || o_nbour->GetHvyDegree() > 1 ) {
      return;
    }
    // if the bond was unsaturated on input, clearly it's ok we just want
    // to not make silly ones ourselves.
    OEBondBase *bond = mol.GetBond( atom , c_nbour );
    if( bond->GetOrder() == 1 && !bond->IsAromatic() ) {
      if( atom_idx < DACLIB::atom_index( *c_nbour ) ) {
        bad_bond_atom_pairs.push_back( make_pair( atom_idx , DACLIB::atom_index( *c_nbour ) ) );
      } else {
        bad_bond_atom_pairs.push_back( make_pair( DACLIB::atom_index( *c_nbour ) , atom_idx ) );
      }
    }
  }

}

// ****************************************************************************
// check for adding H atoms to things that shouldn't have them. At the
// moment, that's S atoms that are more than 1-connected, so we don't get
// credibility-damaging things like H-S(O)(O)=C, as happened in our old
// friend CHEMBL31034.
// Also, don't add Hs to both an N and and O if they're bonded to each other.
// That's also from CHEMBL31034
// From CHEMBL500474_small: Don't add an H to a 1-connected atom and a 2-
// that are bonded if they both have abes of 1 - that will create 2 single
// bonds in a line that is a dead-end tautomer-wise so not an acceptable
// tautomer. But CHEMBL155287 showed that this shouldn't be applied if the
// 2 atoms are both non-C.  I think that should work as a general rule -
// the extra het atom gives another chance for it to escape from this dead end.
// And flag bonds between C and S as unsavoury if the S has been deemed bad
// for H atom addition.
// Also, don't add H atoms to connected Cs - this almost certainly breaks the
// rule about non-terminal atoms changing hybridisation.
// Flag C=C bonds that would give C=C(O)N which we don't allow when generating
// bond paths so will create a tautomer that won't necessarily generate the
// same t_skel as that which created it. Seen in C1C=C(N=C(N1)N)C(C(N)O)N=CO
// from CHEMBL503551.
void flag_bad_atoms_for_adding_h_or_bond( OEMolBase &mol ,
                                          const vector<OEAtomBase *> &hads ,
                                          const vector<unsigned int> &had_idxs ,
                                          const vector<int> &abes ,
                                          vector<int> &bad_atoms ,
                                          vector<pair<unsigned int,unsigned int> > &bad_h_atom_pairs ,
                                          vector<pair<unsigned int,unsigned int> > &bad_bond_atom_pairs ) {

  for( size_t i = 0 , is = hads.size() ; i < is ; ++i ) {
    check_bad_C_S_bonds( mol , hads[i] , had_idxs[i] , bad_atoms ,
                         bad_h_atom_pairs , bad_bond_atom_pairs );
    check_bad_C_C_nbours( hads[i] , had_idxs[i] , bad_h_atom_pairs );
    check_bad_N_O_nbours( mol , hads[i] , had_idxs[i] , bad_h_atom_pairs );
    check_abes_dead_end( hads[i] , had_idxs[i] , abes , bad_h_atom_pairs );
  }

  // C=C bonds probably aren't hads, so go on abes
  for( unsigned int i = 0 , is = static_cast<unsigned int>( abes.size() ) ; i < is ; ++i ) {
    if( abes[i] ) {
      check_C_C_to_amide( mol , i , abes , bad_bond_atom_pairs );
    }
  }

  sort( bad_h_atom_pairs.begin() , bad_h_atom_pairs.end() );
  bad_h_atom_pairs.erase( unique( bad_h_atom_pairs.begin() , bad_h_atom_pairs.end() ) ,
                          bad_h_atom_pairs.end() );

  sort( bad_bond_atom_pairs.begin() , bad_bond_atom_pairs.end() );
  bad_bond_atom_pairs.erase( unique( bad_bond_atom_pairs.begin() , bad_bond_atom_pairs.end() ) ,
                             bad_bond_atom_pairs.end() );

#ifdef NOTYET
  cout << "bad atoms :";
  for( size_t i = 0 , is = bad_atoms.size() ; i < is ; ++i ) {
    cout << " " << bad_atoms[i];
  }
  cout << endl;
  cout << "bad h atom pairs :";
  for( size_t ii = 0 , iis = bad_h_atom_pairs.size() ; ii < iis ; ++ii ) {
    cout << " " << bad_h_atom_pairs[ii].first + 1 << "," << bad_h_atom_pairs[ii].second + 1;
  }
  cout << endl;
  cout << "bad bond atom pairs :";
  for( size_t ii = 0 , iis = bad_bond_atom_pairs.size() ; ii < iis ; ++ii ) {
    cout << " " << bad_bond_atom_pairs[ii].first + 1 << "->" << bad_bond_atom_pairs[ii].second + 1;
  }
  cout << endl;
#endif

}

// ****************************************************************************
// Find all combinations of atoms on which an H atom might be placed. They
// won't all be valid, as for some it won't be possible to generate a valid
// bonding arrangement.  At this point, though, we just want all possibilities.
// It's done by a recursive search via extend_h_atom_combo.
void find_h_atom_combos( const vector<OEAtomBase *> &conn_set ,
                         const vector<int> &abe_set ,
                         const vector<int> &mobile_h ,
                         const vector<int> &atom_at_end_of_bond_paths ,
                         const vector<vector<unsigned int> > &nb_idxs ,
                         const vector<int> &bad_atoms ,
                         const vector<pair<unsigned int,unsigned int> > &bad_h_atom_pairs ,
                         vector<vector<unsigned int> > &h_atom_idxs ) {

  unsigned int conn_set_num_h = 0;
  vector<unsigned int> poss_h_idxs;
  for( size_t i = 0 , is = conn_set.size() ; i < is ; ++i ) {
    unsigned int csi = DACLIB::atom_index( *conn_set[i] );
    if( bad_atoms[csi] ) {
      continue;
    }
    if( mobile_h[csi] ) {
      ++conn_set_num_h;
    }
    if( abe_set[csi] && atom_at_end_of_bond_paths[csi] ) {
      poss_h_idxs.push_back( csi );
    }
  }

#ifdef NOTYET
  cout << "Num conn set mobile h = " << conn_set_num_h << " on :";
  sort( poss_h_idxs.begin() , poss_h_idxs.end() );
  for( size_t ii = 0 , iis = poss_h_idxs.size() ; ii < iis ; ++ii ) {
    cout << " " << poss_h_idxs[ii] + 1;
  }
  cout << endl;
#endif

  vector<int> curr_comb;
  extend_h_atom_combo( conn_set_num_h , static_cast<unsigned int>( poss_h_idxs.size() ) ,
                       poss_h_idxs , nb_idxs , abe_set ,
                       bad_h_atom_pairs , curr_comb , h_atom_idxs );

#ifdef NOTYET
  cout << "valid h combs" << endl;
  for( size_t ii = 0 , iis = h_atom_idxs.size() ; ii < iis ; ++ii ) {
    cout << ii << " : ";
    for( size_t jj = 0 , jjs = h_atom_idxs[ii].size() ; jj < jjs ; ++jj ) {
      cout << " " << h_atom_idxs[ii][jj] + 1;
    }
    cout << endl;
  }
#endif

}

// ****************************************************************************
// find the bonds that should be set to 1 for this tautomer
void find_bonds_to_set_to_1( OEMolBase &mol ,
                             const vector<int> &abes ,
                             const vector<int> &ap_idxs ,
                             const vector<pair<unsigned int,unsigned int> > &bond_atom_pairs_idxs ,
                             vector<int> &bonds_to_1 ) {

  for( size_t i = 0 , is = ap_idxs.size() ; i < is ; ++i ) {
    unsigned int at1_idx = bond_atom_pairs_idxs[ap_idxs[i]].first;
    unsigned int at2_idx = bond_atom_pairs_idxs[ap_idxs[i]].second;
#ifdef NOTYET
    cout <<"at pair " << i << " : " << at1_idx + 1 << " -> " << at2_idx + 1
        << " :: " << abes[at1_idx] << " , " << abes[at2_idx] << endl;
#endif
    if( abes[at1_idx] && abes[at2_idx] ) {
      OEAtomBase *at1 = mol.GetAtom( DACLIB::HasAtomIndex( at1_idx ) );
      OEAtomBase *at2 = mol.GetAtom( DACLIB::HasAtomIndex( at2_idx ) );
      OEBondBase *bond = mol.GetBond( at1 , at2 );
      if( bond && bond->GetOrder() > 1 ) {
        // use bond->GetOrder() - 1 in case it's a triple bond that needs to
        // be returned to double.
        bonds_to_1[DACLIB::bond_index( *bond )] = bond->GetOrder() - 1;
      }
    }
  }

#ifdef NOTYET
  cout << "leaving Bonds to set to 1 : ";
  for( size_t ii = 0 , iis = bonds_to_1.size() ; ii < iis ; ++ii ) {
    if( bonds_to_1[ii] ) {
      OEBondBase *bond = mol.GetBond( DACLIB::HasBondIndex( ii ) );
      cout << DACLIB::atom_index( *bond->GetBgn() ) + 1 << "->"
           << DACLIB::atom_index( *bond->GetEnd() ) + 1
	   << "(" << bonds_to_1[ii] << ") ";
    }
  }
  cout << endl;
#endif

}

// ****************************************************************************
void calc_abes_nconns( const vector<OEAtomBase *> &conn_set ,
                       const vector<int> &abe_set ,
                       const vector<vector<unsigned int> > &nb_idxs ,
                       vector<int> &abes_nconns ) {

  for( size_t i = 0 , is = conn_set.size() ; i < is ; ++i ) {
    unsigned int at_idx = DACLIB::atom_index( *conn_set[i] );
    for( size_t j = 0 , js = nb_idxs[at_idx].size() ; j < js ; ++j ) {
      if( abe_set[nb_idxs[at_idx][j]] ) {
        ++abes_nconns[at_idx];
      }
    }
#ifdef NOTYET
    cout << "at_idx " << at_idx + 1 << " abes_nconn = " << abes_nconns[at_idx] << endl;
#endif
  }

}

// ****************************************************************************
void set_bond_unsatd( OEMolBase &mol , unsigned int at1_idx ,
                      unsigned int at2_idx ,
                      const vector<vector<unsigned int> > &nb_idxs ,
                      vector<unsigned int> &unsat_bond_idxs ,
                      vector<int> &abe_set ,
                      vector<int> &abes_nconns , int &num_abes ) {

#ifdef NOTYET
  cout << "setting " << at1_idx + 1 << "->" << at2_idx + 1 << " to unsatd" << endl;
#endif

  OEAtomBase *at1 = mol.GetAtom( DACLIB::HasAtomIndex( at1_idx ) );
  OEAtomBase *at2 = mol.GetAtom( DACLIB::HasAtomIndex( at2_idx ) );
  OEBondBase *bond = mol.GetBond( at1 , at2 );
  unsat_bond_idxs.push_back( DACLIB::bond_index( *bond ) );
  --abe_set[at1_idx];
  // if at1 now has no abes, then its neighbours need to have
  // abes_nconns reduced by 1
  if( !abe_set[at1_idx] ) {
    for( size_t k = 0 , ks = nb_idxs[at1_idx].size() ; k < ks ; ++k ) {
      if( abes_nconns[nb_idxs[at1_idx][k]] ) {
        --abes_nconns[nb_idxs[at1_idx][k]];
      }
    }
  }
  --abe_set[at2_idx];
  // likewise with at2
  if( !abe_set[at2_idx] ) {
    for( size_t k = 0 , ks = nb_idxs[at2_idx].size() ; k < ks ; ++k ) {
      if( abes_nconns[nb_idxs[at2_idx][k]] ) {
        --abes_nconns[nb_idxs[at2_idx][k]];
      }
    }
  }
  num_abes -= 2;

}

// ****************************************************************************
// assign double bonds from a terminal atom down the chain until there's a
// tri-node
void assign_term_double_bonds( OEMolBase &mol ,
                                const vector<vector<unsigned int> > &nb_idxs ,
                                vector<int> &abe_set ,int &num_abes ,
                                vector<int> &abes_nconns ,
                                vector<unsigned int> &unsat_bond_idxs ) {

#ifdef NOTYET
  cout << "entering assign_term_double_bonds" << endl;
  cout << "abes_nconns : ";
  copy( abes_nconns.begin() , abes_nconns.end() , intOut );
  cout << endl;
  cout << "abe_set : ";
  copy( abe_set.begin() , abe_set.end() , intOut );
  cout << endl;
  cout << "num_abes = " << num_abes << endl;
#endif

  // find any singly-connected abes atoms.  These must necessarily be at the
  // start of a chain and take a double bond.
  bool did_something( true );
  while( did_something && num_abes ) {
    did_something = false;
    for( unsigned int i = 0 , is = static_cast<unsigned int>( abes_nconns.size() ) ; i < is ; ++i ) {
      if( 1 == abes_nconns[i] && abe_set[i] ) {
        for( size_t j = 0 , js = nb_idxs[i].size() ; j < js ; ++j ) {
          unsigned int other_at = nb_idxs[i][j];
          if( abe_set[other_at] ) {
            set_bond_unsatd( mol , i , other_at , nb_idxs , unsat_bond_idxs ,
                             abe_set , abes_nconns , num_abes );
            did_something = true;
            break;
          }
        }
      }
    }
  }

  // at this point, there may be atoms with abes_nconns of 1, but no abes. An
  // example would be the C of Cc1ccccc1 which would have to be left over from
  // a larger molecule fragment. Take these out to make things tidier.
  for( size_t i = 0 , is = abes_nconns.size() ; i < is ; ++i ) {
    if( abes_nconns[i] && !abe_set[i] ) {
      abes_nconns[i] = 0;
    }
  }

#ifdef NOTYET
  cout << "leaving assign_term_double_bonds, num_abes = " << num_abes << endl;
  cout << "abes_nconns : ";
  copy( abes_nconns.begin() , abes_nconns.end() , intOut );
  cout << endl;
  cout << "abe_set : ";
  copy( abe_set.begin() , abe_set.end() , intOut );
  cout << endl;
  cout << "num_abes = " << num_abes << endl;
#endif

}

// ****************************************************************************
void assign_acyclic_double_bonds( OEMolBase &mol ,
                                  const vector<vector<unsigned int> > &nb_idxs ,
                                  vector<int> &abe_set ,int &num_abes ,
                                  vector<int> &abes_nconns ,
                                  vector<unsigned int> &unsat_bond_idxs ) {

#ifdef NOTYET
  cout << "entering assign_acyclic_double_bonds" << endl;
  cout << "abes_nconns : ";
  copy( abes_nconns.begin() , abes_nconns.end() , intOut );
  cout << endl;
  cout << "abe_set : ";
  copy( abe_set.begin() , abe_set.end() , intOut );
  cout << endl;
  cout << "num_abes = " << num_abes << endl;
#endif

  // Kearsley says to do binode to trinode before trinode to trinode
  // he also says to do variations on trinode to trinode in a specific order,
  // which isn't implemented here.
  for( int isum = 5 ; isum < 7 ; ++isum ) {
    for( size_t i = 0 , is = nb_idxs.size() ; i < is ; ++i ) {
      for( size_t j = 0 , js = nb_idxs[i].size() ; j < js ; ++j ) {
        unsigned int other_at = nb_idxs[i][j];
#ifdef NOTYET
        if( abes_nconns[i] || abes_nconns[other_at] ) {
          cout << "looking at " << i + 1 << "->" << other_at + 1 << " to go double : "
               << abes_nconns[i] << " + " << abes_nconns[other_at] << endl;
        }
#endif
        if( isum == abes_nconns[i] + abes_nconns[other_at] &&
            abe_set[i] && abe_set[other_at] ) {
          // see if removing this bond fragments the molecule. If it doesn't,
          // it's cyclic
          unsigned int old_i = abe_set[i];
          unsigned int old_nb_i = abe_set[other_at];
          abe_set[i] = 0;
          abe_set[other_at] = 0;
          vector<vector<OEAtomBase *> > conn_sets;
          vector<vector<int> > new_abe_sets;
          build_connect_sets( abe_set , mol , conn_sets , new_abe_sets );
          abe_set[i] = old_i;
          abe_set[other_at] = old_nb_i;
          if( 2 == conn_sets.size() ) {
            int num_abes1 = accumulate( new_abe_sets.front().begin() , new_abe_sets.front().end() , 0 );
            int num_abes2 = accumulate( new_abe_sets.back().begin() , new_abe_sets.back().end() , 0 );
            if( !(num_abes1 % 2 ) && !(num_abes2 % 2 ) ) {
              set_bond_unsatd( mol , static_cast<unsigned int>( i ) , other_at , nb_idxs , unsat_bond_idxs ,
                               abe_set , abes_nconns , num_abes );
              vector<int> curr_abes_nconns( DACLIB::max_atom_index( mol ) , 0 );
              for( size_t k = 0 , ks = abes_nconns.size() ; k < ks ; ++k ) {
                if( new_abe_sets.front()[k] && abes_nconns[k] ) {
                  curr_abes_nconns[k] = abes_nconns[k];
                }
              }
              assign_term_double_bonds( mol , nb_idxs , new_abe_sets.front() ,
                                        num_abes , curr_abes_nconns , unsat_bond_idxs );
              // update the global abe_set and abes_nconns
              for( size_t k = 0 , ks = conn_sets.front().size() ; k < ks ; ++k ) {
                unsigned int at_idx = DACLIB::atom_index( *conn_sets.front()[k] );
                abe_set[at_idx] = new_abe_sets.front()[at_idx];
                abes_nconns[at_idx] = curr_abes_nconns[at_idx];
              }
              curr_abes_nconns = vector<int>( DACLIB::max_atom_index( mol ) , 0 );
              for( size_t k = 0 , ks = abes_nconns.size() ; k < ks ; ++k ) {
                if( new_abe_sets.back()[k] && abes_nconns[k] ) {
                  curr_abes_nconns[k] = abes_nconns[k];
                }
              }
              assign_term_double_bonds( mol , nb_idxs , new_abe_sets.back() ,
                                        num_abes , curr_abes_nconns , unsat_bond_idxs );
              // update the global abe_set and abes_nconns
              for( size_t k = 0 , ks = conn_sets.back().size() ; k < ks ; ++k ) {
                unsigned int at_idx = DACLIB::atom_index( *conn_sets.back()[k] );
                abe_set[at_idx] = new_abe_sets.back()[at_idx];
                abes_nconns[at_idx] = curr_abes_nconns[at_idx];
              }
            }
          }
        }
      }
    }
  }

#ifdef NOTYET
  cout << "leaving assign_acyclic_double_bonds" << endl;
  cout << "abes_nconns : ";
  copy( abes_nconns.begin() , abes_nconns.end() , intOut );
  cout << endl;
  cout << "abe_set : ";
  copy( abe_set.begin() , abe_set.end() , intOut );
  cout << endl;
  cout << "num_abes = " << num_abes << endl;
#endif

}

// ****************************************************************************
void assign_cyclic_double_bonds( OEMolBase &mol ,
                                 const vector<vector<unsigned int> > &nb_idxs ,
                                 vector<int> &abe_set , int &num_abes ,
                                 vector<int> &abes_nconns ,
                                 vector<unsigned int> &unsat_bond_idxs ) {

#ifdef NOTYET
  cout << "entering sub assign_cyclic_double_bonds" << endl;
  cout << "abes_nconns : ";
  copy( abes_nconns.begin() , abes_nconns.end() , intOut );
  cout << endl;
  cout << "abe_set : ";
  copy( abe_set.begin() , abe_set.end() , intOut );
  cout << endl;
  cout << "num_abes = " << num_abes << endl;
#endif

  // at this point, we should just have rings of alternating single/double
  // bonds. In the vast majority of cases we can just pick one at random, and
  // generate a valid Kekule form.  However, good old CHEMBL19253 with its
  // 4-membered ring demonstrated that for a non-aromatic ring, the order
  // matters - the difference between OC1=CC=C1O and OC1C=CC(O)=1. So pad
  // the end of unsat_bond_idxs with all possibilities separated by
  // MAX_UINT and unpack them later
  bool did_something = true;
  while( did_something && num_abes ) {
    did_something = false;
    int final_num_abes = num_abes;
    vector<int> final_abe_set;
    for( unsigned int i = 0 , is = static_cast<unsigned int>( abes_nconns.size() ) ; i < is ; ++i ) {
      if( abes_nconns[i] && abe_set[i] ) {
        int number_done = 0;
        for( size_t j = 0 , js = nb_idxs[i].size() ; j < js ; ++j ) {
          unsigned int other_at = nb_idxs[i][j];
          vector<int> tmp_abe_set( abe_set );
          if( abe_set[other_at] ) {
#ifdef NOTYET
            cout << "arbitrary unsat : " << i + 1 << "->" << other_at + 1 << endl;
#endif
            vector<int> tmp_abes_nconns( abes_nconns );
            int tmp_num_abes = num_abes;
            vector<unsigned int> tmp_unsat_bond_idxs;
            set_bond_unsatd( mol , i , other_at , nb_idxs , tmp_unsat_bond_idxs ,
                             tmp_abe_set , tmp_abes_nconns , tmp_num_abes );
            assign_term_double_bonds( mol , nb_idxs , tmp_abe_set ,
                                      tmp_num_abes , tmp_abes_nconns ,
                                      tmp_unsat_bond_idxs );
#ifdef NOTYET
            cout << "Makes result : " << tmp_num_abes << endl;
#endif
            if( 0 == tmp_num_abes ) {
              ++number_done;
              // the sort-of porphyrin CHEMBL443682 shows that this isn't always
              // successful due to rings within rings
              unsat_bond_idxs.push_back( numeric_limits<unsigned int>::max() );
              unsat_bond_idxs.insert( unsat_bond_idxs.end() ,
                                      tmp_unsat_bond_idxs.begin() ,
                                      tmp_unsat_bond_idxs.end() );
              final_num_abes = tmp_num_abes;
              final_abe_set = tmp_abe_set;
            }
          }
          // we only want to try each way round the ring from 1 starting
          // point
          if( number_done ) {
            did_something = true;
            num_abes = final_num_abes;
            abe_set = final_abe_set;
          }
        }
      }
      if( did_something ) {
        break;
      }
    }
  }

#ifdef NOTYET
  cout << "leaving sub assign_cyclic_double_bonds" << endl;
  cout << "abes_nconns : ";
  copy( abes_nconns.begin() , abes_nconns.end() , intOut );
  cout << endl;
  cout << "abe_set : ";
  copy( abe_set.begin() , abe_set.end() , intOut );
  cout << endl;
  cout << "num_abes = " << num_abes << endl;
#endif

}

// ****************************************************************************
void update_abes( const vector<vector<int> > &abe_sets ,
                  vector<int> &abes ) {

  if( abe_sets.empty() ) {
    return;
  }
  abes = vector<int>( abe_sets.front().size() , 0 );
  for( size_t i = 0 , is = abe_sets.size() ; i < is ; ++i ) {
    for( size_t j = 0 , js = abe_sets[i].size() ; j < js ; ++j ) {
      if( abe_sets[i][j] ) {
        abes[j] = abe_sets[i][j];
      }
    }
  }

}

// ****************************************************************************
bool check_conn_sets_for_valid_soln( vector<vector<int> > &abe_sets ) {

  for( size_t i = 0 , is = abe_sets.size() ; i < is ; ++i ) {
    int num_abes = accumulate( abe_sets[i].begin() , abe_sets[i].end() , 0 );
    if( num_abes % 2 ) {
#ifdef NOTYET
      cout << "Failure num_abes = " << num_abes << " :" << endl;
#endif
      return false; // an odd number of abes can't be satisfied by double bonds.
    }
  }

  return true;
}

// ****************************************************************************
// take the already-fragment abes and assign unsaturation from any terminal
// chains.
bool assign_terminal_double_bonds( OEMolBase &mol ,
                                   const vector<vector<unsigned int> > &nb_idxs ,
                                   vector<int> &abes ,
                                   vector<unsigned int> &term_unsat_bond_idxs ) {

#ifdef NOTYET
  cout << "entering assign_terminal_double_bonds : " << endl;
  cout << "abes : ";
  copy( abes.begin() , abes.end() , intOut );
  cout << endl;
#endif

  // split the abes into different connect sets. At this stage, this will
  // be due to adding H atoms, as we do terminal chains first
  vector<vector<OEAtomBase *> > conn_sets;
  vector<vector<int> > abe_sets;
  build_connect_sets( abes , mol , conn_sets , abe_sets );

  // check that a solution is possible for all pieces
  if( !check_conn_sets_for_valid_soln( abe_sets ) ) {
    return false;
  }

  // do each piece in turn
  for( size_t i = 0 , is = conn_sets.size() ; i < is ; ++i ) {
#ifdef NOTYET
    cout << "conn set " << i << endl;
#endif
    int num_abes = accumulate( abe_sets[i].begin() , abe_sets[i].end() , 0 );
    vector<int> abes_nconns( DACLIB::max_atom_index( mol ) , 0 );
    // find the number of abes connected to each atom in the conn set
    calc_abes_nconns( conn_sets[i] , abe_sets[i] , nb_idxs , abes_nconns );
    assign_term_double_bonds( mol , nb_idxs , abe_sets[i] ,
                              num_abes , abes_nconns , term_unsat_bond_idxs );
    if( are_there_isolated_atoms( abe_sets[i] , nb_idxs ) ) {
      term_unsat_bond_idxs.clear();
#ifdef NOTYET
      cout << "Failure, isolated atoms : " << endl;
#endif
      return false;
    }
  }

  // update the global abes to reflect changes just made
  update_abes( abe_sets , abes );

#ifdef NOTYET
  cout << "leaving assign_terminal_double_bonds : " << endl;
  for( size_t ii = 0 , iis = term_unsat_bond_idxs.size() ; ii < iis ; ++ii ) {
    OEBondBase *bond = mol.GetBond( DACLIB::HasBondIndex( term_unsat_bond_idxs[ii] ) );
    cout << DACLIB::atom_index( *bond->GetBgn() ) + 1 << "->"
         << DACLIB::atom_index( *bond->GetEnd() ) + 1 << " ";
  }
  cout << endl;
  cout << "abes : ";
  copy( abes.begin() , abes.end() , intOut );
  cout << endl;
#endif

  return true;

}

// ****************************************************************************
// now do the acyclic double bonds, that aren't terminal chains
bool assign_acyclic_double_bonds( OEMolBase &mol ,
                                  const vector<vector<unsigned int> > &nb_idxs ,
                                  vector<int> &abes ,
                                  vector<unsigned int> &acyc_unsat_bond_idxs ) {

  // split the abes into different connect sets
  vector<vector<OEAtomBase *> > conn_sets;
  vector<vector<int> > abe_sets;
  build_connect_sets( abes , mol , conn_sets , abe_sets );

  // check that a solution is possible for all pieces
  if( !check_conn_sets_for_valid_soln( abe_sets ) ) {
    return false;
  }

  // do each piece in turn
  for( size_t i = 0 , is = conn_sets.size() ; i < is ; ++i ) {
#ifdef NOTYET
    cout << "conn set " << i << endl;
#endif
    int num_abes = accumulate( abe_sets[i].begin() , abe_sets[i].end() , 0 );
    vector<int> abes_nconns( DACLIB::max_atom_index( mol ) , 0 );
    // find the number of abes connected to each atom in the conn set
    calc_abes_nconns( conn_sets[i] , abe_sets[i] , nb_idxs , abes_nconns );
    assign_acyclic_double_bonds( mol , nb_idxs , abe_sets[i] ,
                                 num_abes , abes_nconns , acyc_unsat_bond_idxs );
    if( are_there_isolated_atoms( abe_sets[i] , nb_idxs ) ) {
      acyc_unsat_bond_idxs.clear();
#ifdef NOTYET
      cout << "Failure, isolated atoms : " << endl;
#endif
      return false;
    }
  }

  // update the global abes to reflect changes just made
  update_abes( abe_sets , abes );

#ifdef NOTYET
  cout << "leaving assign_acyclic_double_bonds : " << endl;
  for( size_t ii = 0 , iis = acyc_unsat_bond_idxs.size() ; ii < iis ; ++ii ) {
    OEBondBase *bond = mol.GetBond( DACLIB::HasBondIndex( acyc_unsat_bond_idxs[ii] ) );
    cout << DACLIB::atom_index( *bond->GetBgn() ) + 1 << "->"
         << DACLIB::atom_index( *bond->GetEnd() ) + 1 << " ";
  }
  cout << endl;
  cout << "abes : ";
  copy( abes.begin() , abes.end() , intOut );
  cout << endl;
#endif
  return true;

}

// ****************************************************************************
// now do the cyclic double bonds, which are a bit more complicated due to the
// need to allow, in rare cases, for both valid assignments of atoms.
// We return a vector of unsat_bond_idxs for each conn_set that's remaining,
// each of which will have 2 sets of alternative unsaturations in it,
// separated by a MAX_UINT.
bool assign_cyclic_double_bonds( OEMolBase &mol ,
                                 const vector<vector<unsigned int> > &nb_idxs ,
                                 vector<int> &abes ,
                                 vector<vector<unsigned int> > &cyc_unsat_bond_idxs ) {

#ifdef NOTYET
  cout << "entering assign_cyclic_double_bonds" << endl;
#endif

  // split the abes into different connect sets
  vector<vector<OEAtomBase *> > conn_sets;
  vector<vector<int> > abe_sets;
  build_connect_sets( abes , mol , conn_sets , abe_sets );

  // check that a solution is possible for all pieces
  if( !check_conn_sets_for_valid_soln( abe_sets ) ) {
    return false;
  }

  cyc_unsat_bond_idxs = vector<vector<unsigned int> >( conn_sets.size() , vector<unsigned int>() );

  // do each piece in turn
  for( size_t i = 0 , is = conn_sets.size() ; i < is ; ++i ) {
#ifdef NOTYET
    cout << "conn set " << i << " of " << is << endl;
#endif
    int num_abes = accumulate( abe_sets[i].begin() , abe_sets[i].end() , 0 );
    vector<int> abes_nconns( DACLIB::max_atom_index( mol ) , 0 );
    // find the number of abes connected to each atom in the conn set
    calc_abes_nconns( conn_sets[i] , abe_sets[i] , nb_idxs , abes_nconns );
    assign_cyclic_double_bonds( mol , nb_idxs , abe_sets[i] ,
                                num_abes , abes_nconns , cyc_unsat_bond_idxs[i] );
    if( are_there_isolated_atoms( abe_sets[i] , nb_idxs ) ) {
      cyc_unsat_bond_idxs.clear();
#ifdef NOTYET
      cout << "Failure, isolated atoms : " << endl;
#endif
      return false;
    }
  }

  // update the global abes to reflect changes just made
  update_abes( abe_sets , abes );

#ifdef NOTYET
  cout << "leaving assign_cyclic_double_bonds : " << endl;
  for( size_t ii = 0 , iis = cyc_unsat_bond_idxs.size() ; ii < iis ; ++ii ) {
    for( size_t jj = 0 , jjs = cyc_unsat_bond_idxs[ii].size() ; jj < jjs ; ++jj ) {
      if( numeric_limits<unsigned int>::max() == cyc_unsat_bond_idxs[ii][jj] ) {
        cout << " : ";
      } else {
        OEBondBase *bond = mol.GetBond( DACLIB::HasBondIndex( cyc_unsat_bond_idxs[ii][jj] ) );
        cout << DACLIB::atom_index( *bond->GetBgn() ) + 1 << "->"
             << DACLIB::atom_index( *bond->GetEnd() ) + 1 << " ";
      }
    }
    cout << endl;
  }
  cout << "abes : ";
  copy( abes.begin() , abes.end() , intOut );
  cout << endl;
#endif

  return true;

}

// ****************************************************************************
// these_unsat_bond_idxs have multiple possibilities in each vector. Create
// all combinations of these, into exp_unsat_bond_idxs.
void expand_multi_unsat_bonds( vector<vector<unsigned int> > &these_unsat_bond_idxs ,
                               vector<vector<unsigned int> > &exp_unsat_bond_idxs ) {

  for( size_t i = 0 , is = these_unsat_bond_idxs.size() ; i < is ; ++i ) {
    vector<unsigned int> &this_ubi = these_unsat_bond_idxs[i];
    vector<unsigned int>::iterator p = find( this_ubi.begin() , this_ubi.end() ,
                                             numeric_limits<unsigned int>::max() );
    vector<vector<unsigned int> > these_exp_bits;
    if( p == this_ubi.end() ) {
      these_exp_bits.push_back( this_ubi );
    } else {
      vector<unsigned int> stem_ubi( this_ubi.begin() , p );
      ++p;

      vector<vector<unsigned int> > next_bits;
      while( true ) {
        vector<unsigned int>::iterator q = find( p + 1 , this_ubi.end() ,
                                                 numeric_limits<unsigned int>::max() );
        vector<unsigned int> next_bit( p , q );
        sort( next_bit.begin() , next_bit.end() );
        if( next_bits.end() == find( next_bits.begin() , next_bits.end() ,
                                     next_bit ) ) {
          next_bits.push_back( next_bit );
        }
        if( q == this_ubi.end() ) {
          break;
        }
        p = q + 1;
      }
      for( size_t j = 0 , js = next_bits.size() ; j < js ; ++j ) {
        these_exp_bits.push_back( stem_ubi );
        these_exp_bits.back().insert( these_exp_bits.back().end() ,
                                      next_bits[j].begin() , next_bits[j].end() );
      }
    }
    if( exp_unsat_bond_idxs.empty() ) {
      exp_unsat_bond_idxs = these_exp_bits;
    } else {
      // form all combinations of each of these_exp_bits with what's already
      // there
      vector<vector<unsigned int> > last_exp_unsat_bond_idxs = exp_unsat_bond_idxs;
      exp_unsat_bond_idxs.clear();
      exp_unsat_bond_idxs.reserve( last_exp_unsat_bond_idxs.size() * these_exp_bits.size() );
      for( size_t j = 0 , js = these_exp_bits.size() ; j < js ; ++j ) {
        size_t curr_top = exp_unsat_bond_idxs.size();
        exp_unsat_bond_idxs.insert( exp_unsat_bond_idxs.end() , last_exp_unsat_bond_idxs.begin() , last_exp_unsat_bond_idxs.end() );
        for( size_t k = 0 , ks = last_exp_unsat_bond_idxs.size() ; k < ks ; ++k ) {
          vector<unsigned int> &t = exp_unsat_bond_idxs[curr_top+k];
          t.insert( t.end() , these_exp_bits[j].begin() , these_exp_bits[j].end() );
        }
      }
    }
  }

#ifdef NOTYET
  cout << "leaving expand_multi_unsat_bonds" << endl;
  for( size_t ii = 0 , iis = exp_unsat_bond_idxs.size() ; ii < iis ; ++ii ) {
    cout << ii << " : ";
    copy( exp_unsat_bond_idxs[ii].begin() , exp_unsat_bond_idxs[ii].end() , uintOut );
    cout << endl;
  }
#endif

}

// ****************************************************************************
// take this abes set, which is assumed to have had the mobile h atoms added
// already, and find the bonds that need to be made unsaturated to fill the
// valences. This is using ideas culled from Kearsley, Computers. Chem.,
// 17, 1-10 (1993).
void build_unsatd_bond_combos( OEMolBase &mol , const vector<int> &abes ,
                               const vector<vector<unsigned int> > &nb_idxs ,
                               const vector<int> &all_ap_idxs ,
                               vector<pair<unsigned int,unsigned int> > &bond_atom_pairs_idxs ,
                               vector<unsigned int> &h_atom_idxs ,
                               vector<vector<unsigned int> > &unsat_bond_idxs ,
                               vector<int> &bonds_to_1 ) {

#ifdef NOTYET
  cout << "build_unsatd_bond_combos" << endl;
  cout << "ABES : ";
  copy( abes.begin() , abes.end() , intOut );
  cout << endl;
  cout << DACLIB::create_cansmi( mol ) << endl;
  sort( h_atom_idxs.begin() , h_atom_idxs.end() );
  cout << "Atoms for Hs :";
  for( size_t ii = 0 , iis = h_atom_idxs.size() ; ii < iis ; ++ii ) {
    cout << " " << h_atom_idxs[ii] + 1;
  }
  cout << endl;
#endif

  vector<int> this_abes( abes );
  find_bonds_to_set_to_1( mol , this_abes , all_ap_idxs ,
                          bond_atom_pairs_idxs , bonds_to_1 );
  // amend available bonding electrons for this H atom combo
  for( size_t j = 0 , js = h_atom_idxs.size() ; j < js ; ++j ) {
    --this_abes[h_atom_idxs[j]];
  }

  unsat_bond_idxs.clear();

  vector<unsigned int> term_unsat_bond_idxs;
  if( !assign_terminal_double_bonds( mol , nb_idxs , this_abes , term_unsat_bond_idxs ) ) {
    return;
  }

  vector<unsigned int> acyc_unsat_bond_idxs;
  if( !assign_acyclic_double_bonds( mol , nb_idxs , this_abes , acyc_unsat_bond_idxs ) ) {
    return;
  }

  vector<vector<unsigned int> > cyc_unsat_bond_idxs;
  if( !assign_cyclic_double_bonds( mol , nb_idxs , this_abes , cyc_unsat_bond_idxs ) ) {
    return;
  }

  // if there are still abes to do, it's a fail
  if( accumulate( this_abes.begin() , this_abes.end() , 0 ) ) {
#ifdef NOTYET
    cout << "Fail on final num_abes : "
         << accumulate( this_abes.begin() , this_abes.end() , 0 ) << endl;
#endif
    return;
  }

  // now put all the pieces together
  if( cyc_unsat_bond_idxs.empty() ) {
    unsat_bond_idxs.push_back( term_unsat_bond_idxs );
    unsat_bond_idxs.front().insert( unsat_bond_idxs.front().end() ,
                                    acyc_unsat_bond_idxs.begin(),
                                    acyc_unsat_bond_idxs.end() );
  } else {
    expand_multi_unsat_bonds( cyc_unsat_bond_idxs , unsat_bond_idxs );
    // tack the terminal and acyclic bits on
    for( size_t i = 0 , is = unsat_bond_idxs.size() ; i < is ; ++i ) {
      unsat_bond_idxs[i].insert( unsat_bond_idxs[i].end() ,
                                 term_unsat_bond_idxs.begin(),
                                 term_unsat_bond_idxs.end() );
      unsat_bond_idxs[i].insert( unsat_bond_idxs[i].end() ,
                                 acyc_unsat_bond_idxs.begin(),
                                 acyc_unsat_bond_idxs.end() );
    }
  }

#ifdef NOTYET
  cout << "leaving build_unsatd_bond_combos : ";
  cout << "Success :";
  for( size_t ii = 0 , iis = h_atom_idxs.size() ; ii < iis ; ++ii ) {
    cout << " " << h_atom_idxs[ii] + 1;
  }
  cout << endl;
  for( size_t jj = 0 , jjs = unsat_bond_idxs.size() ; jj < jjs ; ++jj ) {
    copy( unsat_bond_idxs[jj].begin() , unsat_bond_idxs[jj].end() , uintOut );
    cout << endl;
    for( size_t ii = 0 , iis = unsat_bond_idxs[jj].size() ; ii < iis ; ++ii ) {
      if( numeric_limits<unsigned int>::max() == unsat_bond_idxs[jj][ii] ) {
        continue;
      }
      OEBondBase *bond = mol.GetBond( DACLIB::HasBondIndex( unsat_bond_idxs[jj][ii] ) );
      cout << DACLIB::atom_index( *bond->GetBgn() ) + 1 << "->"
           << DACLIB::atom_index( *bond->GetEnd() ) + 1 << " ";
    }
    cout << endl;
  }
  cout << "bonds_to_1 : ";
  copy( bonds_to_1.begin() , bonds_to_1.end() , intOut );
  cout << endl;
  for( size_t ii = 0 , iis = bonds_to_1.size() ; ii < iis ; ++ii ) {
    if( bonds_to_1[ii] ) {
      OEBondBase *bond = mol.GetBond( DACLIB::HasBondIndex( ii ) );
      cout << DACLIB::atom_index( *bond->GetBgn() ) + 1 << "->"
           << DACLIB::atom_index( *bond->GetEnd() ) + 1 << " ";
    }
  }
  cout << endl;
#endif

}

// ****************************************************************************
void build_bond_atom_pair_indices( const OEMolBase &mol ,
                                   const vector<pair<OEAtomBase *,OEAtomBase *> > &bond_atom_pairs ,
                                   vector<int> &bond_ap_idxs ) {

  vector<pair<unsigned int,int> > tmp_idxs;
  tmp_idxs.reserve( bond_atom_pairs.size() );
  for( size_t i = 0 , is = bond_atom_pairs.size() ; i < is ; ++i ) {
    OEBondBase *bond = mol.GetBond( bond_atom_pairs[i].first ,
                                    bond_atom_pairs[i].second );
    tmp_idxs.push_back( make_pair( i , bond->GetOrder() ) );
  }

  stable_sort( tmp_idxs.begin() , tmp_idxs.end() ,
               bind( greater<int>() ,
                     bind( &pair<unsigned int,int>::second , _1 ) ,
                     bind( &pair<unsigned int,int>::second , _2 ) ) );

  for( size_t i = 0 , is = tmp_idxs.size() ; i < is ; ++i ) {
    bond_ap_idxs.push_back( tmp_idxs[i].first );
  }

}

// ****************************************************************************
// There's possible multiple entries in these_unsat_bond_idxs for
// the ith bonds_to_1 and h_atom_idxs.
void unpack_multi_unsat_bonds( vector<vector<unsigned int> > &these_unsat_bond_idxs ,
                               unsigned int i ,
                               vector<vector<unsigned int> > &unsat_bond_idxs ,
                               vector<vector<unsigned int> > &h_atom_idxs ,
                               vector<vector<int> > &bonds_to_1 ) {

  // if there's only 1 entry, don't make life difficult
  if( 1 == these_unsat_bond_idxs.size() ) {
    unsat_bond_idxs[i] = these_unsat_bond_idxs.front();
    return;
  }

  // Take a copy of the ith entries, the ones of interest, then clear them.
  // find_unsatd_bonds removes the empties all in one go at the end.
  unsat_bond_idxs[i].clear();
  vector<unsigned int> this_hai( h_atom_idxs[i] );
  h_atom_idxs[i].clear();
  vector<int> this_b_to_1( bonds_to_1[i] );
  bonds_to_1[i].clear();

  // now put the appropriate number of copies of h_atom_idxs and bonds_to_1
  // for the new exp_unsat_bond_idxs.
  for( size_t j = 0 , js = these_unsat_bond_idxs.size() ; j < js ; ++j ) {
    unsat_bond_idxs.push_back( these_unsat_bond_idxs[j] );
    h_atom_idxs.push_back( this_hai );
    bonds_to_1.push_back( this_b_to_1 );
  }

}

// ****************************************************************************
// take the given set of atoms, which form a tautomer network, and the sets of
// atoms to add Hs to, and find the bonds to add to make up the valences in
// abes.  If a particular H combination can't be reconciled with a set of double
// bonds, remove it from the input set.
void find_unsaturated_bonds( OEMolBase &mol ,
                             const vector<int> &abes ,
                             const vector<vector<unsigned int> > &nb_idxs ,
                             vector<pair<OEAtomBase *,OEAtomBase *> > &bond_atom_pairs ,
                             vector<pair<unsigned int,unsigned int> > &bond_atom_pairs_idxs ,
                             vector<vector<unsigned int> > &h_atom_idxs ,
                             vector<vector<unsigned int> > &unsat_bond_idxs ,
                             vector<vector<int> > &bonds_to_1 ) {

#ifdef NOTYET
  cout << "find_unsaturated_bonds" << endl;
  cout << "Number of bond pairs in this connect set : " << bond_atom_pairs.size() << endl;
#endif

  unsat_bond_idxs = vector<vector<unsigned int> >( h_atom_idxs.size() , vector<unsigned int>() );
  bonds_to_1 = vector<vector<int> >( h_atom_idxs.size() , vector<int>( DACLIB::max_bond_index( mol ) , 0 ) );

  for( unsigned int i = 0 , is = static_cast<unsigned int>( h_atom_idxs.size() ) ; i < is ; ++i ) {
#ifdef NOTYET
    cout << "doing h atom combo : " << i << " :";
    for( size_t ii = 0 , iis = h_atom_idxs[i].size() ; ii < iis ; ++ii ) {
      cout << " " << h_atom_idxs[i][ii] + 1;
    }
    cout << endl;
#endif

    // get indices for the bond atom pairs, sorted in descending order of
    // original bond order. This is to get round the kekule issue seen in
    // CHEMBL3306810, where the re-saturation of the aromatic system
    // gave the other kekule form so it didn't see that it was always
    // setting the same bonds to double, and should have been taken out of
    // the t_skel.
    vector<int> all_ap_idxs;
    build_bond_atom_pair_indices( mol , bond_atom_pairs , all_ap_idxs );
    vector<vector<unsigned int> > these_unsat_bond_idxs;
    build_unsatd_bond_combos( mol , abes , nb_idxs , all_ap_idxs ,
                              bond_atom_pairs_idxs , h_atom_idxs[i] ,
                              these_unsat_bond_idxs , bonds_to_1[i] );

    if( these_unsat_bond_idxs.empty() ) {
      // this combo is a bust
      h_atom_idxs[i].clear();
      bonds_to_1[i].clear();
    } else {
      unpack_multi_unsat_bonds( these_unsat_bond_idxs , i , unsat_bond_idxs ,
                                h_atom_idxs , bonds_to_1 );
    }
  }

  // take out the empty ones
  h_atom_idxs.erase( remove_if( h_atom_idxs.begin() , h_atom_idxs.end() ,
                                bind( &vector<unsigned int>::empty , _1 ) ) ,
                     h_atom_idxs.end() );
  unsat_bond_idxs.erase( remove_if( unsat_bond_idxs.begin() , unsat_bond_idxs.end() ,
                                    bind( &vector<unsigned int>::empty , _1 ) ) ,
                         unsat_bond_idxs.end() );
  bonds_to_1.erase( remove_if( bonds_to_1.begin() , bonds_to_1.end() ,
                               bind( &vector<int>::empty , _1 ) ) ,
                    bonds_to_1.end() );

}

// ****************************************************************************
// atoms_for_hs indexes into hads and gives the atoms that should have an H
// atom added as part of the tautomer generation. At the end of this,
// is_had, hads and had_idxs might be changed.
void find_atoms_for_hs_and_unsat_bonds( OEMolBase &inmol ,
                                        const vector<int> &abes ,
                                        const vector<int> &atom_at_end_of_bond_path ,
                                        const vector<OEAtomBase *> connect_set ,
                                        vector<int> &abe_set ,
                                        vector<int> &mobile_h ,
                                        vector<OEAtomBase *> &hads ,
                                        vector<unsigned int> &had_idxs ,
                                        vector<vector<unsigned int> > &atoms_for_hs ,
                                        vector<vector<unsigned int> > &unsat_bond_idxs ,
                                        vector<vector<int> > &bonds_to_1 ) {

#ifdef NOTYET
  for( OEIter<OEBondBase> bond = inmol.GetBonds() ; bond ; ++bond ) {
    cout << DACLIB::bond_index( *bond ) << " : "
         << DACLIB::atom_index( *bond->GetBgn() ) + 1 << "->"
         << DACLIB::atom_index( *bond->GetEnd() ) + 1 << " : "
         << bond->GetOrder() << " : " << bond->IsAromatic() << endl;
  }
#endif

#ifdef NOTYET
  cout << "find_atoms_for_hs_and_unsat_bonds" << endl;
  cout << "mobile_h : ";
  copy( mobile_h.begin() , mobile_h.end() , intOut );
  cout << endl;
#endif

  // build the neighbour indices for each atom for speed of access later.
  vector<vector<unsigned int> > nb_idxs( DACLIB::max_atom_index( inmol ) , vector<unsigned int>() );
  build_nbor_idxs( inmol , nb_idxs );

  // find atoms that aren't good for adding H atoms to, either individually
  // or as pairs.
  vector<int> bad_atoms( DACLIB::max_atom_index( inmol ) , 0 );
  vector<pair<unsigned int,unsigned int> > bad_h_atom_pairs , bad_bond_atom_pairs;
  flag_bad_atoms_for_adding_h_or_bond( inmol , hads , had_idxs , abes , bad_atoms ,
                                       bad_h_atom_pairs , bad_bond_atom_pairs );

  // get all the bonds except bad_bond_atom_pairs, sorted in descending order
  // of number of n'bours.
  vector<pair<OEAtomBase *,OEAtomBase *> > bond_atom_pairs;
  vector<pair<unsigned int,unsigned int> > bond_atom_pairs_idxs;
  build_bond_atom_pairs( inmol , nb_idxs , bad_bond_atom_pairs ,
                         bond_atom_pairs , bond_atom_pairs_idxs );

  // get a set of possible atoms to add H's to.
  find_h_atom_combos( connect_set , abe_set , mobile_h ,
                      atom_at_end_of_bond_path , nb_idxs ,
                      bad_atoms , bad_h_atom_pairs , atoms_for_hs );
  // find_unsaturated_bonds may take things out of atoms_for_hs
  // if it's not possible to satisfy the valences with unsaturation.
  // This is one of the most time-consuming parts of the code, normally.
  find_unsaturated_bonds( inmol , abe_set , nb_idxs , bond_atom_pairs ,
                          bond_atom_pairs_idxs , atoms_for_hs ,
                          unsat_bond_idxs , bonds_to_1 );

}

// ****************************************************************************
void find_atoms_at_end_of_bond_paths( OEMolBase &inmol ,
                                      const vector<vector<OEAtomBase *> > &bond_paths ,
                                      vector<int> &atom_at_end_of_bond_path ) {

  atom_at_end_of_bond_path = vector<int>( DACLIB::max_atom_index( inmol ) , 0 );
  for( size_t i = 0 , is = bond_paths.size() ; i < is ; ++i ) {
    atom_at_end_of_bond_path[DACLIB::atom_index( *bond_paths[i].front() )] = 1;
    atom_at_end_of_bond_path[DACLIB::atom_index( *bond_paths[i].back() )] = 1;
  }

#ifdef NOTYET
  cout << "atoms at ends of bond paths :";
  for( size_t i = 0 , is = atom_at_end_of_bond_path.size() ; i < is ; ++i ) {
    if( atom_at_end_of_bond_path[i] ) {
      cout << " " << i + 1;
    }
  }
  cout << endl;
#endif

}

// ****************************************************************************
// take the combinations of atoms for hs and bonds for unsaturation and remove
// from hads and bond_paths any atoms that don't appear in any of them, which
// will possibly then feed back to atoms_for_hs and unsat_bond_idxs.
void update_hads_and_bond_paths( OEMolBase &mol ,
                                 vector<vector<unsigned int> > &atoms_for_hs ,
                                 vector<vector<int> > &bonds_to_1 ,
                                 vector<vector<unsigned int> > &unsat_bond_idxs ,
                                 vector<OEAtomBase *> hads ,
                                 vector<unsigned int> &had_idxs ,
                                 vector<int> &is_had ,
                                 vector<vector<OEAtomBase *> > &bond_paths ) {

#ifdef NOTYET
  cout << "entering update_hads_and_bond_paths" << endl;
  cout << "Atoms for Hs" << endl;
  for( unsigned ii = 0 , iis = atoms_for_hs.size() ; ii < iis ; ++ii ) {
    cout << ii << " : ";
    for( unsigned jj = 0 , jjs = atoms_for_hs[ii].size() ; jj < jjs ; ++jj ) {
      cout << atoms_for_hs[ii][jj] + 1 << " ";
    }
    cout << endl;
  }
  cout << "bonds_to_1 :" << endl;
  for( size_t ii = 0 , iis = bonds_to_1.size() ; ii < iis ; ++ii ) {
    for( size_t jj = 0 , jjs = bonds_to_1[ii].size() ; jj < jjs ; ++jj ) {
      if( bonds_to_1[ii][jj] ) {
        OEBondBase *b = mol.GetBond( DACLIB::HasBondIndex( jj ) );
        cout << " (" << jj << ") " << DACLIB::atom_index( *b->GetBgn() ) + 1 << "->"
             << DACLIB::atom_index( *b->GetEnd() ) + 1;
      }
    }
    cout << endl;
  }
  cout << endl;

  cout << "Unsat bonds" << endl;
  for( size_t ii = 0 , iis = unsat_bond_idxs.size() ; ii < iis ; ++ii ) {
    cout << ii << " : ";
    for( unsigned jj = 0 , jjs = unsat_bond_idxs[ii].size() ; jj < jjs ; ++jj ) {
      OEBondBase *bond = mol.GetBond( DACLIB::HasBondIndex( unsat_bond_idxs[ii][jj] ) );
      cout << DACLIB::atom_index( *bond->GetBgn() ) + 1 << " -> "
           << DACLIB::atom_index( *bond->GetEnd() ) + 1 << " : ";
    }
    cout << endl;
  }
  cout << "HADS :";
  for( size_t ii = 0 , iis = had_idxs.size() ; ii < iis ; ++ii ) {
    cout << " " << had_idxs[ii] + 1;
  }
  cout << endl;
  cout << "is_had : ";
  copy( is_had.begin() , is_had.end() , intOut );
  cout << endl;
  cout << "Bond paths" << endl;
  for( size_t jj = 0 , jjs = bond_paths.size() ; jj < jjs ; ++jj ) {
    cout << "Path " << jj << " ::";
    for( size_t ii = 0 , iis = bond_paths[jj].size() ; ii < iis ; ++ii ) {
      cout << " " << DACLIB::atom_index( *bond_paths[jj][ii] ) + 1;
    }
    cout << endl;
  }
#endif

  vector<int> used_atom( DACLIB::max_atom_index( mol ) , 0 );
  for( size_t i = 0 , is = atoms_for_hs.size() ; i < is ; ++i ) {
    for( size_t j = 0 , js = atoms_for_hs[i].size() ; j < js ; ++j ) {
      used_atom[atoms_for_hs[i][j]] = 1;
    }
  }

  for( size_t i = 0 , is = unsat_bond_idxs.size() ; i < is ; ++i ) {
    for( size_t j = 0 , js = unsat_bond_idxs[i].size() ; j < js ; ++j ) {
      OEBondBase *bond = mol.GetBond( DACLIB::HasBondIndex( unsat_bond_idxs[i][j] ) );
      used_atom[DACLIB::atom_index( *bond->GetBgn() )] = 1;
      used_atom[DACLIB::atom_index( *bond->GetEnd() )] = 1;
    }
  }

  for( size_t i = 0 , is = used_atom.size() ; i < is ; ++i ) {
    if( !used_atom[i] ) {
#ifdef NOTYET
      if( is_had[i] ) {
        cout << "had " << i + 1 << " not used" << endl;
      }
#endif
      is_had[i] = 0;
    }
  }

  for( size_t i = 0 , is = had_idxs.size() ; i < is ; ++i ) {
    if( !is_had[had_idxs[i]] ) {
      had_idxs[i] = numeric_limits<unsigned int>::max();
      hads[i] = 0;
    }
  }
  had_idxs.erase( remove( had_idxs.begin() , had_idxs.end() ,
                          numeric_limits<unsigned int>::max() ) ,
                  had_idxs.end() );
  hads.erase( remove( hads.begin() , hads.end() ,
                      static_cast<OEAtomBase *>( 0 ) ) , hads.end() );

  remove_paths_with_non_hads_both_ends( is_had , bond_paths );
  remove_hads_not_in_paths( bond_paths , hads , had_idxs , is_had );

  // now remove any non-hads from atoms_for_hs. If this leaves an empty atoms_for_hs
  // vector, the corresponding unsat_bond_idxs must be cleared as well.
  vector<int> atom_removed( DACLIB::max_atom_index( mol ) , 0 );
  for( size_t i = 0 , is = atoms_for_hs.size() ; i < is ; ++i ) {
    for( size_t j = 0 , js = atoms_for_hs[i].size() ; j < js ; ++j ) {
      if( !is_had[atoms_for_hs[i][j]] ) {
        atom_removed[atoms_for_hs[i][j]] = 1;
        atoms_for_hs[i][j] = numeric_limits<unsigned int>::max();
      }
    }
    atoms_for_hs[i].erase( remove( atoms_for_hs[i].begin() , atoms_for_hs[i].end() ,
                                   numeric_limits<unsigned int>::max() ) ,
                           atoms_for_hs[i].end() );
    if( atoms_for_hs[i].empty() ) {
      unsat_bond_idxs[i].clear();
    }
  }

  // and remove from unsat_bond_idxs any bond that includes one of the atoms
  // removed from atoms_for_hs.
  // Empty the corresponding atoms_for_hs if that gives an empty vector.
  for( size_t i = 0 , is = unsat_bond_idxs.size() ; i < is ; ++i ) {
    for( size_t j = 0 , js = unsat_bond_idxs[i].size() ; j < js ; ++j ) {
      OEBondBase *bond = mol.GetBond( DACLIB::HasBondIndex( unsat_bond_idxs[i][j] ) );
      unsigned int bb = DACLIB::atom_index( *bond->GetBgn() );
      unsigned int be = DACLIB::atom_index( *bond->GetEnd() );
      if( atom_removed[bb] || atom_removed[be] ) {
        cout << "bond " << unsat_bond_idxs[i][j] << " : " << be + 1 << "->" << bb + 1 << " removed from set " << i
             << " by update_hads_and_bond_paths" << endl;
        cout << "bonds_to_1 for this bond is " << bonds_to_1[i][unsat_bond_idxs[i][j]] << endl;
        if( bonds_to_1[i][unsat_bond_idxs[i][j]] ) {
          // if we're not going to unsaturate this bond, don't set it to 1
          // first if it's flagged for it.
          --bonds_to_1[i][unsat_bond_idxs[i][j]];
          cout << "bonds_to_1 for this bond is " << bonds_to_1[i][unsat_bond_idxs[i][j]] << endl;
        }
        unsat_bond_idxs[i][j] = numeric_limits<unsigned int>::max();
      }
    }
    unsat_bond_idxs[i].erase( remove( unsat_bond_idxs[i].begin() , unsat_bond_idxs[i].end() ,
                                      numeric_limits<unsigned int>::max() ) ,
                              unsat_bond_idxs[i].end() );
    if( unsat_bond_idxs[i].empty() ) {
      bonds_to_1[i].clear();
      atoms_for_hs[i].clear();
    }

  }

  // now take out any empty atoms_for_hs and unsat_bond_idxs - they should be
  // in step.
  atoms_for_hs.erase( remove_if( atoms_for_hs.begin() , atoms_for_hs.end() ,
                                 bind( &vector<unsigned int>::empty , _1 ) ) ,
                      atoms_for_hs.end() );
  bonds_to_1.erase( remove_if( bonds_to_1.begin() , bonds_to_1.end() ,
                               bind( &vector<int>::empty , _1 ) ) ,
                    bonds_to_1.end() );
  unsat_bond_idxs.erase( remove_if( unsat_bond_idxs.begin() , unsat_bond_idxs.end() ,
                                    bind( &vector<unsigned int>::empty , _1 ) ) ,
                         unsat_bond_idxs.end() );

#ifdef NOTYET
  cout << "leaving update_hads_and_bond_paths" << endl;
  cout << "Atoms for Hs" << endl;
  for( unsigned ii = 0 , iis = atoms_for_hs.size() ; ii < iis ; ++ii ) {
    cout << ii << " : ";
    for( unsigned jj = 0 , jjs = atoms_for_hs[ii].size() ; jj < jjs ; ++jj ) {
      cout << atoms_for_hs[ii][jj] + 1 << " ";
    }
    cout << endl;
  }
  cout << "bonds_to_1 :" << endl;
  for( size_t ii = 0 , iis = bonds_to_1.size() ; ii < iis ; ++ii ) {
    for( size_t jj = 0 , jjs = bonds_to_1[ii].size() ; jj < jjs ; ++jj ) {
      if( bonds_to_1[ii][jj] ) {
        OEBondBase *b = mol.GetBond( DACLIB::HasBondIndex( jj ) );
        cout << " (" << jj << ") " << DACLIB::atom_index( *b->GetBgn() ) + 1 << "->"
             << DACLIB::atom_index( *b->GetEnd() ) + 1;
      }
    }
    cout << endl;
  }
  cout << endl;

  cout << "Unsat bonds" << endl;
  for( size_t ii = 0 , iis = unsat_bond_idxs.size() ; ii < iis ; ++ii ) {
    cout << ii << " : ";
    for( unsigned jj = 0 , jjs = unsat_bond_idxs[ii].size() ; jj < jjs ; ++jj ) {
      OEBondBase *bond = mol.GetBond( DACLIB::HasBondIndex( unsat_bond_idxs[ii][jj] ) );
      cout << DACLIB::atom_index( *bond->GetBgn() ) + 1 << " -> "
           << DACLIB::atom_index( *bond->GetEnd() ) + 1 << " : ";
    }
    cout << endl;
  }
  cout << "HADS :";
  for( size_t ii = 0 , iis = had_idxs.size() ; ii < iis ; ++ii ) {
    cout << " " << had_idxs[ii] + 1;
  }
  cout << endl;
  cout << "Bond paths" << endl;
  for( size_t jj = 0 , jjs = bond_paths.size() ; jj < jjs ; ++jj ) {
    cout << "Path " << jj << " ::";
    for( size_t ii = 0 , iis = bond_paths[jj].size() ; ii < iis ; ++ii ) {
      cout << " " << DACLIB::atom_index( *bond_paths[jj][ii] ) + 1;
    }
    cout << endl;
  }
#endif

}

// ****************************************************************************
// take the set of bonds to be set to 1 for all tautomers and merge into a
// a global set which will be used to generate the t_skel
void merge_bonds_to_1( const vector<vector<int> > &all_bonds_to_1 ,
                       vector<int> &bonds_to_1 ) {

  if( all_bonds_to_1.empty() ) {
    return;
  }

  for( size_t i = 0 , is = all_bonds_to_1.size() ; i < is ; ++i ) {
    for( size_t j = 0 , js = all_bonds_to_1[i].size() ; j < js ; ++j ) {
      if( all_bonds_to_1[i][j] ) {
        bonds_to_1[j] = 1;
      }
    }
  }

}

// ****************************************************************************
// if a mobile_h has an h added to it for all atoms_for_hs, it's a bit
// pointless and adds an unnecessary atom to the global t_skel. Take them out
// if they occur.
void check_mobile_h_counts( vector<int> &mobile_h ,
                            vector<vector<unsigned int> > &atoms_for_hs ) {

  vector<unsigned int> atoms_for_hs_cnts( mobile_h.size() , 0 );
  for( size_t j = 0 , js = atoms_for_hs.size() ; j < js ; ++j ) {
    for( size_t k = 0 , ks = atoms_for_hs[j].size() ; k < ks ; ++k ) {
      ++atoms_for_hs_cnts[atoms_for_hs[j][k]];
    }
  }

  for( size_t i = 0 , is = atoms_for_hs_cnts.size() ; i < is ; ++i ) {
    if( mobile_h[i] && atoms_for_hs.size() == atoms_for_hs_cnts[i] ) {
      mobile_h[i] = 0;
      for( size_t j = 0 , js = atoms_for_hs.size() ; j < js ; ++j ) {
        atoms_for_hs[j].erase( remove( atoms_for_hs[j].begin() ,
                                       atoms_for_hs[j]. end() , i ) ,
                               atoms_for_hs[j].end() );
      }
    }
  }

}

// ****************************************************************************
// find the bits that need to be changed for this input molecule - the atoms
// with mobile hs and the bonds that need to be set to 1 in the t_skel, and
// the combinations of atoms and bonds to be altered in the t_skel to generate
// all the tautomers.
void find_tautomers_details( OEMolBase &mol , bool ignore_amides ,
                             vector<vector<int> > &all_mobile_h ,
                             vector<vector<int> > &all_t_skel_bonds_to_1 ,
                             vector<vector<vector<int> > > &all_bonds_to_1 ,
                             vector<vector<vector<unsigned int > > > &all_atoms_for_hs ,
                             vector<vector<vector<unsigned int> > > &all_unsat_bond_idxs ) {

  // HAD is an H-atom acceptor atom or H-atom donor atom.
  vector<OEAtomBase *> hads;
  vector<unsigned int> had_idxs;
  vector<int> is_had;

  vector<vector<OEAtomBase *> > bond_paths;
  find_hads( mol , ignore_amides , bond_paths , hads , had_idxs , is_had );

  if( hads.empty() ) {
    return;
  }

  // ABEs are Atom Bonding Electrons and are the free valences when the mobile H
  // atoms are removed and the bondpath orders are set to 1.
  vector<int> abes( DACLIB::max_atom_index( mol ) , 0 );
  vector<int> mobile_h( DACLIB::max_atom_index( mol ) , 0 );
  int num_act_h = 0;
  calc_abes( mol , hads , bond_paths , abes , num_act_h , mobile_h );

  // only atoms at an end of a bond path can receive an H atom.
  vector<int> atom_at_end_of_bond_path;
  find_atoms_at_end_of_bond_paths( mol , bond_paths ,
                                   atom_at_end_of_bond_path );

  // Get the connect sets together - these are contiguous bits
  // of abe atoms.  We want to treat these independently, as we don't want to
  // move H atoms from one tautomer system to another.  This happened in
  // the peptide (peptide.smi)
  // NC(=O)C(CC(=O)O)NC(=O)C(CCCCN)NC(=O)C(CO)NC(=O)C2CCCN2C(=O)C3CCC(=O)N3
  // And also, to prevent a combinatarial explosion when doing big molecules
  // such as peptides.  Do the tautomer expansion for each connect set
  // separately. If you want all the tautomers at the end, that will be
  // a separate step.
  vector<vector<OEAtomBase *> > connect_sets;
  vector<vector<int> > abe_sets;
  build_connect_sets( abes , mol , connect_sets , abe_sets );

  // get the combinations of atoms to which h atoms should be added and
  // unsaturated bonds made for each tautomer.  Ignoring symmetry, there'll
  // be one tautomer for each atoms_for_hs/unsat_bond_idxs entry.
  // Do it separately for each connect set. It would be more efficient if each
  // connect set was done entirely independently, including generating
  // tautomers and iterating round to to see if the t_skel changes with the
  // new tautomers. However, that would take a significant re-factoring of the
  // code, so I haven't bothered for now.  Doing it this way means that if a
  // particular connect set generates a new set of tautomers, then all the
  // tautomers of the other connect sets in the molecule will be generated
  // even if that is just repeating itself and not adding anything new. The
  // alternative is more complicated, and would need to deal with a connect set
  // splitting into 2 which would be a book-keeping complication.  This new
  // method is already a massive performance score over the previous method of
  // generating all combinations of all changes across all the connect sets,
  // so let's not bother.

  for( size_t i = 0 , is = connect_sets.size() ; i < is ; ++i ) {
#ifdef NOTYET
    cout << "Connect set " << i << " : ";
    for( size_t j = 0 , js = connect_sets[i].size() ; j < js ; ++j ) {
      cout << " " << DACLIB::atom_index( *connect_sets[i][j] ) + 1;
    }
    cout << endl;
    cout << "had_idxs : ";
    for( size_t jj = 0 , jjs = had_idxs.size() ; jj < jjs ; ++jj ) {
      cout << " " << had_idxs[jj] + 1;
    }
    cout << endl;
    cout << "abes : ";
    copy( abes.begin() , abes.end() , intOut );
    cout << endl;
#endif
    // the hads and had_idxs and all the other stuff can be changed in the
    // loop, so take a working copy
    vector<OEAtomBase *> these_hads( hads );
    vector<unsigned int> these_had_idxs( had_idxs );
    vector<int> these_is_had( is_had );
    vector<vector<OEAtomBase *> > these_bond_paths( bond_paths );

    vector<vector<unsigned int> > these_atoms_for_hs;
    vector<vector<unsigned int> > these_unsat_bond_idxs;
    vector<vector<int> > these_bonds_to_1;
    vector<int> this_mobile_h( DACLIB::max_atom_index( mol ) , 0 );
    for( size_t j = 0 , js = connect_sets[i].size() ; j < js ; ++j ) {
      unsigned int at_idx = DACLIB::atom_index( *connect_sets[i][j] );
      this_mobile_h[at_idx] = mobile_h[at_idx];
    }

    // find the atoms that need H atoms adding to them and bonds that need
    // their orders increasing to make all the different tautomers. This can
    // change hads, had_idxs and is_had because if an h is added to an atom in all
    // successful tautomers, it's clearly not part of the tautomer skeleton
    // and so is removed.
    find_atoms_for_hs_and_unsat_bonds( mol , abes , atom_at_end_of_bond_path ,
                                       connect_sets[i] , abe_sets[i] ,
                                       this_mobile_h , these_hads ,
                                       these_had_idxs , these_atoms_for_hs ,
                                       these_unsat_bond_idxs ,
                                       these_bonds_to_1 );

    // see if we remove and add an H to the same atom in all cases and
    // update this_mobile_h and atoms_for_hs if so.
    check_mobile_h_counts( this_mobile_h , these_atoms_for_hs );

    // now we know the atoms that are to take hs and bonds that need to be
    // unsaturated, adjust the hads and bond paths. This is because we may have
    // taken out hads that appear in all combinations of atoms_for_hs or
    // unsat_bond_idxs, and this can feed back into atoms_for_hs and
    // unsat_bond_idxs
    update_hads_and_bond_paths( mol , these_atoms_for_hs , these_bonds_to_1 ,
                                these_unsat_bond_idxs , these_hads ,
                                these_had_idxs , these_is_had ,
                                these_bond_paths );

    all_mobile_h.push_back( this_mobile_h );
    all_bonds_to_1.push_back( these_bonds_to_1 );
    all_atoms_for_hs.push_back( these_atoms_for_hs );
    all_unsat_bond_idxs.push_back( these_unsat_bond_idxs );
    vector<int> these_t_skel_bonds_to_1( DACLIB::max_bond_index( mol ) , 0 );
    merge_bonds_to_1( these_bonds_to_1 , these_t_skel_bonds_to_1 );
    all_t_skel_bonds_to_1.push_back( these_t_skel_bonds_to_1 );

#ifdef NOTYET
    cout << "leaving find_tautomers_details" << endl;
    cout << "atoms for hs :" << endl;
    for( size_t ii = 0 , iis = these_atoms_for_hs.size() ; ii < iis ; ++ii ) {
      cout << ii << " :";
      for( size_t jj = 0 , jjs = these_atoms_for_hs[ii].size() ; jj < jjs ; ++jj ) {
        cout << " " << these_atoms_for_hs[ii][jj] + 1;
      }
      cout << " ::";
      for( size_t jj = 0 , jjs = these_bonds_to_1[ii].size() ; jj < jjs ; ++jj ) {
        if( these_bonds_to_1[ii][jj] ) {
          OEBondBase *bond = mol.GetBond( DACLIB::HasBondIndex( static_cast<unsigned int>( jj ) ) );
          cout << " " << DACLIB::atom_index( *bond->GetBgn() ) + 1 << "->"
               << DACLIB::atom_index( *bond->GetEnd() ) + 1;
        }
      }
      cout << " ::";
      for( size_t jj = 0 , jjs = these_unsat_bond_idxs[ii].size() ; jj < jjs ; ++jj ) {
        OEBondBase *bond = mol.GetBond( DACLIB::HasBondIndex( these_unsat_bond_idxs[ii][jj] ) );
        cout << " " << DACLIB::atom_index( *bond->GetBgn() ) + 1 << "->"
             << DACLIB::atom_index( *bond->GetEnd() ) + 1;
      }
      cout << endl;
    }
#endif
  }

}

// ****************************************************************************
int update_global_t_skel_details( const vector<int> &mobile_h ,
                                  const vector<int> &bonds_to_1 ,
                                  vector<int> &global_mobile_h ,
                                  vector<int> &global_bonds_to_1 ) {

  int num_changes = 0;
  for( size_t j = 0 , js = mobile_h.size() ; j < js ; ++j ) {
    if( mobile_h[j] ) {
      if( !global_mobile_h[j] ) {
        ++num_changes;
      }
      global_mobile_h[j] = max( global_mobile_h[j] , mobile_h[j] );
    }
  }
  for( size_t j = 0 , js = bonds_to_1.size() ; j < js ; ++j ) {
    if( bonds_to_1[j] ) {
      if( !global_bonds_to_1[j] ) {
        ++num_changes;
      }
      global_bonds_to_1[j] = 1;
    }
  }

#ifdef NOTYET
  cout << "update_global_t_skel_details, num_changes = " << num_changes << endl;
#endif
  return num_changes;

}

// ****************************************************************************
// because of the way mobile_h and bond_paths_to_1 are built up iteratively,
// we may have mobile_h greater than the number of Hs on the atom in mol, so
// adjust down.
void check_global_mobile_h( OEMolBase &mol ,
                            vector<int> &mobile_h ) {

  for( OEIter<OEAtomBase> atom = mol.GetAtoms() ; atom ; ++atom ) {
    unsigned int at_idx = DACLIB::atom_index( *atom );
    if( static_cast<int>( atom->GetTotalHCount() ) < mobile_h[at_idx] ) {
      mobile_h[at_idx] = atom->GetTotalHCount();
    }
  }

}

// ****************************************************************************
// because of the way mobile_h and bond_paths_to_1 are built up iteratively,
// we may have arrived at a case, as seen in CHEMBL501944_small, where the
// rule has been broken about carbon atom hads that are unsaturated never
// being H atom donors. So check that and amend mobile H accordingly.
void check_mobile_h_on_c( OEMolBase &mol , const vector<int> &bonds_to_1 ,
                          vector<int> &mobile_h ) {

  for( unsigned int i = 0 , is = static_cast<unsigned int>( bonds_to_1.size() ) ; i < is ; ++i ) {
    if( bonds_to_1[i] ) {
      OEBondBase *bond = mol.GetBond( DACLIB::HasBondIndex( i ) );
      if( bond && bond->GetOrder() > 1 ) {
        if( OEElemNo::C == bond->GetBgn()->GetAtomicNum() ) {
          mobile_h[DACLIB::atom_index( *bond->GetBgn() )] = 0;
        }
        if( OEElemNo::C == bond->GetEnd()->GetAtomicNum() ) {
          mobile_h[DACLIB::atom_index( *bond->GetEnd() )] = 0;
        }
      }
    }
  }

}

// ****************************************************************************
// terminal atoms that are unsaturated and due to have both an H removed and
// the bond set to 1 shouldn't lose their H atoms. As seen in
// CHEMBL153534 : CC1=CC(=CN1C)C2=CSC(=NC(=N)N)N2.
// But our old friend CHEMBL19253 says that HN=C=C shouldn't have this
// applied. The new check_redundant_bonds_to_1 should have made this
// 2nd check unnecessary
void check_mobile_h_on_term_atoms( OEMolBase &mol ,
                                   const vector<int> &bonds_to_1 ,
                                   vector<int> &mobile_h ) {

  for( OEIter<OEAtomBase> atom = mol.GetAtoms( OEHasHvyDegree( 1 ) ) ; atom ; ++atom ) {
    unsigned int at_ind = DACLIB::atom_index( *atom );
    if( mobile_h[at_ind] ) {
      OEIter<OEBondBase> bond = atom->GetBonds();
      if( bond->GetOrder() > 1 && bonds_to_1[DACLIB::bond_index( *bond )] ) {
        OEAtomBase *nbr = atom == bond->GetBgn() ? bond->GetEnd() : bond->GetBgn();
        int num_unsat = 0;
        for( OEIter<OEBondBase> nbr_bond = nbr->GetBonds() ; nbr_bond ; ++nbr_bond ) {
          if( nbr_bond->GetOrder() > 1 || nbr_bond->IsAromatic() ) {
            ++num_unsat;
          }
        }
        if( 1 == num_unsat ) {
#ifdef NOTYET
          cout << "Mobile H for " << at_ind + 1 << " rescinded : "
               << bond->GetOrder() << " and "
               << bonds_to_1[DACLIB::bond_index( *bond )] << " :: "
               << DACLIB::atom_index( *bond->GetBgn() ) + 1
               << " -> " << DACLIB::atom_index( *bond->GetEnd() ) + 1
               << endl;
#endif
          mobile_h[at_ind] = 0;
        }
      }
    }
  }

}

// ****************************************************************************
// if removing mobile Hs from mol reduces the number of Hs on an atom in mol
// to less than all the tautomers have, then it's not right.
// It needs to be all non-canonical tautomers, not all canonical ones, see
// CHEMBL153534.
void check_minimum_h_on_tauts( OEMolBase &mol ,
                               const vector<pTautGen> &all_taut_gens ,
                               vector<int> &mobile_h ) {

  vector<unsigned int> min_h_count( DACLIB::max_atom_index( mol ) , 100 );
  for( size_t i = 0 , is = all_taut_gens.size() ; i < is ; ++i ) {
    vector<pOEMolBase> all_tauts = all_taut_gens[i]->generate_all_tautomers();
    for( size_t j = 0 , js = all_tauts.size() ; j < js ; ++j ) {
      for( OEIter<OEAtomBase> atom = all_tauts[j]->GetAtoms() ; atom ; ++atom ) {
        unsigned int at_idx = DACLIB::atom_index( *atom );
        min_h_count[at_idx] = min( min_h_count[at_idx] , atom->GetTotalHCount() );
      }
    }
  }

  for( OEIter<OEAtomBase> atom = mol.GetAtoms() ; atom ; ++atom ) {
    unsigned int at_idx = DACLIB::atom_index( *atom );
#ifdef NOTYET
    cout << at_idx + 1 << " : " << atom->GetTotalHCount() << " : "
         << mobile_h[at_idx] << " : " << min_h_count[at_idx] << endl;
#endif
    // because atom->GetTotalHCount() returns an unsigned int we need to be
    // careful with the maths because, whilst I don't think it can happen,
    // if mobile_h[at_idx] is greater than atom->GetTotalHCount() it could
    // all go horribly wrong with unsigned ints. NUTPOACTSYO.
    int final_h_count = atom->GetTotalHCount() - mobile_h[at_idx];
    if( 100 != min_h_count[at_idx] &&
        final_h_count < static_cast<int>( min_h_count[at_idx] ) ) {
#ifdef NOTYET
      cout << "BAD mobile_h for " << at_idx + 1 << " : " << final_h_count
           << " vs " << min_h_count[at_idx] << endl;
#endif
      mobile_h[at_idx] = static_cast<int>( min_h_count[at_idx] ) -
          static_cast<int>( atom->GetTotalHCount() );
    }
  }

}

// ****************************************************************************
void check_mobile_h( OEMolBase &mol ,
                     const vector<pTautGen> &all_taut_gens ,
                     const vector<int> &bonds_to_1 ,
                     vector<int> &mobile_h ) {

  // to start with, we can't remove more Hs than there are. Some tautomers
  // might have more Hs than mol, but since we're making a global t_skel,
  // from mol, that's not relevant.
  check_global_mobile_h( mol , mobile_h );
  check_mobile_h_on_c( mol , bonds_to_1 , mobile_h );
  check_mobile_h_on_term_atoms( mol , bonds_to_1 , mobile_h );
  check_minimum_h_on_tauts( mol , all_taut_gens , mobile_h );

}

// ****************************************************************************
// see if the global_bonds_to_1 includes things that are always going to be
// unsaturated. These need to be taken out of the global t_skel, be left in
// the individual sets for creating tautomers.
// CHEMBL19253, the gift that keeps on giving - the NC=C=C tautomers.
// But there's an issue with kekule bonds, so don't do it for aromatic
// systems - CHEMBL31359
void check_redundant_bonds_to_1( OEMolBase &mol ,
                                 vector<pOEMolBase > &all_t_skel_master_mols ,
                                 const vector<vector<vector<unsigned int> > > &all_unsat_bond_idxs ,
                                 const vector<int> &is_aromatic ,
                                 const vector<vector<vector<int> > > &all_bonds_to_1 ,
                                 vector<int> &global_bonds_to_1 ) {

  vector<unsigned int> min_bond_orders( DACLIB::max_bond_index( mol ) ,
                                        numeric_limits<unsigned int>::max() );
  for( size_t i = 0 , is = all_unsat_bond_idxs.size() ; i < is ; ++i ) {
    vector<unsigned int> orig_bond_orders( DACLIB::max_bond_index( mol ) , 0 );
    for( OEIter<OEBondBase> bond = all_t_skel_master_mols[i]->GetBonds() ; bond ; ++bond ) {
      orig_bond_orders[DACLIB::bond_index( *bond )] = bond->GetOrder();
    }
    for( size_t j = 0 , js = all_unsat_bond_idxs[i].size() ; j < js ; ++j ) {
      vector<unsigned int> these_bond_orders( orig_bond_orders );
      for( size_t k = 0 , ks = these_bond_orders.size() ; k < ks ; ++k ) {
        these_bond_orders[k] -= all_bonds_to_1[i][j][k];
      }
      for( size_t k = 0 , ks = all_unsat_bond_idxs[i][j].size() ; k < ks ; ++k ) {
        ++these_bond_orders[all_unsat_bond_idxs[i][j][k]];
      }
      for( size_t k = 0 , ks = these_bond_orders.size() ; k < ks ; ++k ) {
        min_bond_orders[k] = min( min_bond_orders[k] , these_bond_orders[k] );
      }
    }
  }

  for( OEIter<OEBondBase> bond = mol.GetBonds() ; bond ; ++bond ) {
    unsigned int bi = DACLIB::bond_index( *bond );
    if( !is_aromatic[bi] && min_bond_orders[bi] == bond->GetOrder() ) {
#ifdef NOTYET
      cout << "REDUNDANT bond in global_t_skel : " << bi
           << " : " << DACLIB::atom_index( *bond->GetBgn() ) + 1
           << " -> " << DACLIB::atom_index( *bond->GetEnd() ) << endl;
#endif
      global_bonds_to_1[bi] = 0;
    }
  }

}

// ****************************************************************************
void flag_aromatic_bonds( OEMolBase &mol ,
                         vector<int> &is_aromatic ) {

  for( OEIter<OEBondBase> bond = mol.GetBonds( OEIsAromaticBond() ) ; bond ; ++bond ) {
    is_aromatic[DACLIB::bond_index( *bond )] = 1;
  }

}

// ****************************************************************************
// Save this generator and make the connect set tautomers for it, and add them
// to the ones already found.
void store_new_taut_gen( pTautGen &new_taut_gen ,
                         vector<pTautGen> &all_taut_gens ,
                         vector<pOEMolBase> &new_tauts ) {
#ifdef NOTYET
  cout << "store new taut_gen for " << new_taut_gen->global_t_skel_smi() << endl;
#endif

  all_taut_gens.push_back( new_taut_gen );

  vector<pOEMolBase> these_tauts = new_taut_gen->generate_conn_set_tauts();
  new_tauts.insert( new_tauts.end() , these_tauts.begin() , these_tauts.end() );

}

// ****************************************************************************
// new_taut_gen and old_taut_gen have the same global t_skel, but they may
// have arrived there from different starting points and therefore may have
// different tautomers, different mobile_h etc.  If the new one creates
// tautomers not seen in the old, then keep it.
void store_existing_taut_gen( pTautGen &new_taut_gen ,
                              pTautGen &old_taut_gen ,
                              vector<pTautGen> &all_taut_gens ,
                              vector<pOEMolBase> &new_tauts ) {

#ifdef NOTYET
  cout << "storing existing taut_gen for " << new_taut_gen->global_t_skel_smi() << endl;
#endif

  new_taut_gen->prune( *old_taut_gen );
  if( !new_taut_gen->global_t_skel_smi().empty() ) {
    // empty global_t_skel_smi means it pruned down to nothing
    store_new_taut_gen( new_taut_gen , all_taut_gens, new_tauts );
  }

}

// ****************************************************************************
void update_tauts_list( vector<pOEMolBase> &new_tauts ,
                        vector<string> &all_taut_smis ,
                        vector<pOEMolBase> &all_tauts ) {

  int next_new = 0;
  for( size_t j = 0 , js = new_tauts.size() ; j < js ; ++j ) {
    string ts = DACLIB::create_cansmi( *new_tauts[j] );
    if( all_taut_smis.end() == find( all_taut_smis.begin() , all_taut_smis.end() , ts ) ) {
#ifdef NOTYET
      cout << "New tautomer : " << ts << endl;
#endif
      new_tauts[next_new++] = new_tauts[j];
      all_taut_smis.push_back( ts );
      all_tauts.push_back( new_tauts[j] );
    } else {
      new_tauts[j].reset();
    }
  }
  new_tauts.erase( new_tauts.begin() + next_new , new_tauts.end() );

}

#ifdef NOTYET
// This isn't being used - it seemed like a good idea, but wasn't in practice.
// It seemed a shame to ditch it, though. One day we might want to know how
// to find the symmetric atoms and bonds in a molecule.
// ****************************************************************************
void find_symmetric_atoms( OEMolBase &mol ,
                           vector<pair<unsigned int,unsigned int> > &equiv_atoms ) {

  OEPerceiveSymmetry( mol );

  for( OEIter<OEAtomBase> at1 = mol.GetAtoms() ; at1 ; ++at1 ) {
    for( OEIter<OEAtomBase> at2 = mol.GetAtoms() ; at2 ; ++at2 ) {
      if( at1 == at2 || at1->GetSymmetryClass() != at2->GetSymmetryClass() ) {
        continue;
      }
      unsigned int at1_idx = DACLIB::atom_index( *at1 );
      unsigned int at2_idx = DACLIB::atom_index( *at2 );
      if( at1_idx > at2_idx ) {
        swap( at1_idx , at2_idx );
      }
      pair<unsigned int,unsigned int> at_pair( at1_idx , at2_idx );
      if( equiv_atoms.end() == find( equiv_atoms.begin() , equiv_atoms.end() ,
                                     at_pair ) ) {
        equiv_atoms.push_back( at_pair );
      }
    }
  }

  cout << "equiv_atoms : " << endl;
  for( size_t i = 0 , is = equiv_atoms.size() ; i < is ; ++i ) {
    cout << equiv_atoms[i].first + 1 << " and " << equiv_atoms[i].second + 1 << endl;
  }

}

// ****************************************************************************
void store_bond_pair_equiv( OEMolBase &mol ,
                            OEAtomBase *bond1_1 , OEAtomBase *bond1_2 ,
                            OEAtomBase *bond2_1 , OEAtomBase *bond2_2 ,
                            vector<pair<unsigned int,unsigned int> > &equiv_bonds ) {

  OEBondBase *b1 = mol.GetBond( bond1_1 , bond1_2 );
  unsigned int b1_idx = DACLIB::bond_index( *b1 );
  OEBondBase *b2 = mol.GetBond( bond2_1 , bond2_2 );
  unsigned int b2_idx = DACLIB::bond_index( *b2 );
  if( b1_idx > b2_idx ) {
    swap( b1_idx , b2_idx );
  }
  pair<unsigned int,unsigned int> bd_pair( b1_idx , b2_idx );
  if( equiv_bonds.end() == find( equiv_bonds.begin() , equiv_bonds.end() ,
                                 bd_pair ) ) {
    equiv_bonds.push_back( bd_pair );
  }

}

// ****************************************************************************
void find_symmetric_bonds( OEMolBase &mol ,
                           vector<pair<unsigned int,unsigned int> > &equiv_atoms ,
                           vector<pair<unsigned int,unsigned int> > &equiv_bonds ) {

  // there's no equivalent of OEAtomBase::GetSymmetryClass for bonds, more's
  // the pity.

  for( size_t i = 0 , is = equiv_atoms.size() ; i < is ; ++i ) {
    OEAtomBase *at1 = mol.GetAtom( DACLIB::HasAtomIndex( equiv_atoms[i].first ) );
    OEAtomBase *at2 = mol.GetAtom( DACLIB::HasAtomIndex( equiv_atoms[i].second ) );
    for( OEIter<OEAtomBase> at1_nb = at1->GetAtoms() ; at1_nb ; ++at1_nb ) {
      for( OEIter<OEAtomBase> at2_nb = at2->GetAtoms() ; at2_nb ; ++at2_nb ) {
        if( at1_nb == at2_nb ) {
          // at1 and at2 have a common neighbour, so the bonds must be equivalent
          store_bond_pair_equiv( mol , at1 , at1_nb , at2 , at2_nb ,
                                 equiv_bonds );
        }
        // if at1_nb and at2_nb are also an equiv_atoms pair, then the bond
        // is equivalent
        unsigned int at1_nb_idx = DACLIB::atom_index( *at1_nb );
        unsigned int at2_nb_idx = DACLIB::atom_index( *at2_nb );
        if( at1_nb_idx > at2_nb_idx ) {
          swap( at1_nb_idx , at2_nb_idx );
        }
        if( equiv_atoms.end() != find( equiv_atoms.begin() , equiv_atoms.end() ,
                                       make_pair( at1_nb_idx , at2_nb_idx ) ) ) {
          store_bond_pair_equiv( mol , at1 , at1_nb , at2 , at2_nb ,
                                 equiv_bonds );
        }
      }
    }
  }

  cout << "equiv_bonds : " << endl;
  for( size_t i = 0 , is = equiv_bonds.size() ; i < is ; ++i ) {
    cout << equiv_bonds[i].first + 1 << " and " << equiv_bonds[i].second + 1 << endl;
  }

}

// ****************************************************************************
void apply_symmetry_to_t_skel( OEMolBase &t_skel_mol ,
                               vector<int> &mobile_h ,
                               vector<int> &bonds_to_1 ) {

  vector<pair<unsigned int,unsigned int> > equiv_atoms;
  find_symmetric_atoms( t_skel_mol , equiv_atoms );

  for( size_t i = 0 , is = equiv_atoms.size() ; i < is ; ++i ) {
    cout << "mobile_h equivs : " << mobile_h[equiv_atoms[i].first]
         << " and " << mobile_h[equiv_atoms[i].second]
         << " for " << equiv_atoms[i].first + 1 << " and "
         << equiv_atoms[i].second + 1 << endl;
    if( mobile_h[equiv_atoms[i].first] ) {
      mobile_h[equiv_atoms[i].second] = mobile_h[equiv_atoms[i].first];
    }
    if( mobile_h[equiv_atoms[i].second] ) {
      mobile_h[equiv_atoms[i].first] = mobile_h[equiv_atoms[i].second];
    }
  }

  cout << bonds_to_1.size() << endl;
  vector<pair<unsigned int,unsigned int> > equiv_bonds;
  find_symmetric_bonds( t_skel_mol , equiv_atoms , equiv_bonds );

  for( size_t i = 0 , is = equiv_bonds.size() ; i < is ; ++i ) {
    OEBondBase *b1 = t_skel_mol.GetBond( DACLIB::HasBondIndex( equiv_bonds[i].first ) );
    OEBondBase *b2 = t_skel_mol.GetBond( DACLIB::HasBondIndex( equiv_bonds[i].second ) );
    cout << "bonds_to_1 equivs : " << bonds_to_1[equiv_bonds[i].first]
         << " and " << bonds_to_1[equiv_bonds[i].second]
         << " for " << DACLIB::atom_index( *b1->GetBgn() ) + 1 << "->"
         << DACLIB::atom_index( *b1->GetEnd() ) + 1 << " and "
         << DACLIB::atom_index( *b2->GetBgn() ) + 1 << "->"
         << DACLIB::atom_index( *b2->GetEnd() ) + 1 << endl;
    if( bonds_to_1[equiv_bonds[i].first] ) {
      bonds_to_1[equiv_bonds[i].second] = bonds_to_1[equiv_bonds[i].first];
    }
    if( bonds_to_1[equiv_bonds[i].second] ) {
      bonds_to_1[equiv_bonds[i].first] = bonds_to_1[equiv_bonds[i].second];
    }
  }

}
#endif

// ****************************************************************************
// merge all the t_skels found into 1 global one for the whole molecule.
void create_global_t_skel( pOEMolBase &master_mol ,
                           const vector<int> &global_is_aromatic ,
                           vector<pTautGen> &all_taut_gens ,
                           bool skip_h_check ,
                           OEMolBase &final_t_skel_mol ) {

#ifdef NOTYET
  cout << "create_global_t_skel" << endl;
#endif

  vector<int> global_mobile_h( DACLIB::max_atom_index( *master_mol ) , 0 );
  vector<int> global_bonds_to_1( DACLIB::max_bond_index( *master_mol ) , 0 );
  vector<pOEMolBase> all_t_skel_master_mols;
  vector<vector<vector<int> > > all_bonds_to_1;
  vector<vector<vector<unsigned int> > > all_unsat_bond_idxs;

  for( size_t i = 0 , is = all_taut_gens.size() ; i < is ; ++i ) {
    all_t_skel_master_mols.push_back( all_taut_gens[i]->mol() );
    all_bonds_to_1.push_back( all_taut_gens[i]->get_all_bonds_to_1() );
    all_unsat_bond_idxs.push_back( all_taut_gens[i]->get_all_unsat_bond_idxs() );
    const vector<vector<int> > &mh = all_taut_gens[i]->mobile_h();
    for( size_t j = 0 , js = mh.size() ; j < js ; ++j ) {
      for( size_t k = 0 , ks = mh[j].size() ; k < ks ; ++k ) {
        global_mobile_h[k] = max( global_mobile_h[k] , mh[j][k] );
      }
    }
    const vector<vector<int> > &b_to_1 = all_taut_gens[i]->t_skel_bonds_to_1();
    for( size_t j = 0 , js = b_to_1.size() ; j < js ; ++j ) {
      for( size_t k = 0 , ks = b_to_1[j].size() ; k < ks ; ++k ) {
        if( b_to_1[j][k] ) {
          global_bonds_to_1[k] = 1;
        }
      }
    }
  }

  // Check for redundant global_bonds_to_1 in final set.
  check_redundant_bonds_to_1( *master_mol , all_t_skel_master_mols ,
                              all_unsat_bond_idxs , global_is_aromatic ,
                              all_bonds_to_1 , global_bonds_to_1 );

  // Because of the way mobile_h and bond_paths_to_1 are built up iteratively,
  // we may have arrived at a case, as seen in CHEMBL501944_small, where the
  // rule has been broken about carbon atom hads that are unsaturated never
  // being H atom donors. So check that and amend mobile H accordingly.
  // Also, if it's something like C=N, don't take the H off the N if it's
  // terminal and the double bond is also going.
  // However, doing this involved enumerating all tautomers.  If the search has
  // timed out, that's not sensible because it's pretty much guaranteed that
  // the number of tautomers is what caused the timeout.  By timing out,
  // we've already agreed that we're not getting the complete picture.
  if( !skip_h_check ) {
    check_mobile_h( *master_mol , all_taut_gens , global_bonds_to_1 ,
                    global_mobile_h );
  }

  final_t_skel_mol = *master_mol;
  build_t_skel_mol( global_mobile_h , global_bonds_to_1 , final_t_skel_mol );

}

// ****************************************************************************
// Takes an input SMILES string and returns an OEMolBase which is the tautomer
// skeleton (t_skel) for it.  Also returns all the information needed to
// generate all the tautomers that lead to the t_skel.  If all you care about
// is the t_skel itself, you only need the first and last arguments, the rest
// can be thrown away by the calling function.
// ignore_amides says whether or not to ignore isolated amides (such as amide
// bonds in peptides).  This can have a dramatic effect on performance, and
// the amide/imidic acid tautomer is a bit controversial.
void generate_t_skel( const string &in_smi , const string &mol_name ,
                      bool ignore_amides , float max_time ,
                      vector<pTautGen> &all_taut_gens ,
                      OEMolBase &final_t_skel_mol ) {

  // for the t_skel, we need a complete set of all the atoms that need an H
  // to be removed (the mobile_h vector) and bonds that need to be set to 1
  // (the bonds_to_1 vector).  To generate a tautomer, we need the mobile_h
  // and bonds_to_1 vectors for the individual input structure, and the
  // atoms_for_hs and unsat_bond_idxs for each individual tautomer.
  // It's essential that each new tautomer molecule is derived from the master
  // molecule created from in_smi, because otherwise the atom and bond indices
  // won't be consistent.

  pOEMolBase full_master_mol( OENewMolBase( OEMolBaseType::OEDefault ) );
  if( !OESmilesToMol( *full_master_mol , in_smi ) ) {
    cerr << "Error parsing SMILES " << in_smi << "." << endl;
    return;
  }
  full_master_mol->SetTitle( mol_name.c_str() );

  // Deal with mixtures by taking each part separately and putting the
  // t_skels together at the end. For historical reasons, each component
  // is called master_mol, and pulled out of full_master_mol. You can
  // probably work out why it ended up like this!

  // create atom indices for the full master mol so that mixtures
  // still have each atom with a unique index.  This may not be important,
  // but I'm a pill 'n' condom sort of guy.
  DACLIB::create_atom_indices( *full_master_mol );
  DACLIB::create_bond_indices( *full_master_mol );

  // timeouts to the whole mixture rather than individual components.
  OEStopwatch watch;
  watch.Start();

  std::vector<unsigned int> parts( full_master_mol->GetMaxAtomIdx());
  unsigned int pcount = OEDetermineComponents( *full_master_mol , &parts[0] );
  OEPartPred pred( &parts[0] , full_master_mol->GetMaxAtomIdx() );
  for( unsigned int mol_part = 1; mol_part <= pcount ; ++mol_part ) {

    pred.SelectPart( mol_part );
    pOEMolBase master_mol( OENewMolBase( OEMolBaseType::OEDefault ) );
    OESubsetMol( *master_mol , *full_master_mol , pred );

    // If there are fewer than 3 atoms (metal ions being a case in point)
    // there's no possibility of a tautomer, and some of the code asssumes
    // at least one bond with messy results if there isn't one.
    if( master_mol->NumAtoms() < 3 ) {
      final_t_skel_mol += *master_mol;
      continue;
    }

#ifdef NOTYET
    cout << "Component mol_part : " << DACLIB::create_noncansmi( *master_mol ) << endl;
#endif
    vector<OEAtomBase *> rad_atoms;
    DACLIB::radical_atoms( *master_mol , rad_atoms );
    if( !rad_atoms.empty() ) {
      cout << full_master_mol->GetTitle() << " : cowardly refusal to tackle radical species : "
           << DACLIB::create_noncansmi( *full_master_mol ) << endl;
      final_t_skel_mol += *master_mol;
      continue;
    }

    DACLIB::apply_daylight_aromatic_model( *master_mol );
    vector<int> global_is_aromatic( DACLIB::max_bond_index( *master_mol ) , 0 );

    vector<pOEMolBase > new_tauts( 1 , master_mol );
    vector<string> all_taut_smis( 1 , DACLIB::create_cansmi( *master_mol ) );
    vector<pOEMolBase > all_tauts( 1 , master_mol );

    while( 1 ) {

#ifdef NOTYET
      cout << "NEXT LOOP ROUND, number of tautomers to do : " << new_tauts.size() << endl;
      cout << "Number of taut_gens : " << all_taut_gens.size()
           << " number of tauts to consider : " << new_tauts.size() << endl;
#endif

      vector<pOEMolBase > next_new_tauts;

      for( size_t i = 0 , is = new_tauts.size() ; i < is ; ++i ) {
#ifdef NOTYET
        cout << endl << "find all tautomers for " << i + 1 << " : " << DACLIB::create_noncansmi( *new_tauts[i] )
             << " : " << DACLIB::create_cansmi( *new_tauts[i] )
             << " (of " << new_tauts.size() << ")" << endl;
#endif

        if( watch.Elapsed() > max_time ) {
          cout << "Timeout after " << watch.Elapsed() << "s" << endl;
          break;
        }
        flag_aromatic_bonds( *new_tauts[i] , global_is_aromatic );

        // find the tautomer details for each tautomer system in this tautomer
        vector<vector<int> > these_all_mobile_h;
        vector<vector<int> > these_t_skel_bonds_to_1;
        vector<vector<vector<int> > > these_all_bonds_to_1;
        vector<vector<vector<unsigned int > > > these_all_atoms_for_hs;
        vector<vector<vector<unsigned int> > > these_all_unsat_bond_idxs;

        find_tautomers_details( *new_tauts[i] , ignore_amides , these_all_mobile_h ,
                                these_t_skel_bonds_to_1 , these_all_bonds_to_1 ,
                                these_all_atoms_for_hs , these_all_unsat_bond_idxs );

        // Make a tautomer generator for this set of results. The t_skels for
        // each one will be made at the same time, as well as the global
        // t_skel
        pTautGen this_taut_gen( new TautomerGenerator( new_tauts[i] ,
                                                       these_all_mobile_h ,
                                                       these_t_skel_bonds_to_1 ,
                                                       these_all_bonds_to_1 ,
                                                       these_all_atoms_for_hs ,
                                                       these_all_unsat_bond_idxs ) );
        unsigned int old_tg = 0;
        for( ; old_tg < all_taut_gens.size() ; ++old_tg ) {
          if( *this_taut_gen == *all_taut_gens[old_tg] ) {
            break;
          }
        }
        if( old_tg == all_taut_gens.size() ) {
          // if this t_skel hasn't been seen yet, then it represents a possible new
          // set of tautomers, so store the details and generate the tautomers so
          // we can go round again.
          store_new_taut_gen( this_taut_gen , all_taut_gens , next_new_tauts );
        } else {
          // we've already got this t_skel, but it has obviously come to us
          // from a different starting tautomer.  It might, therefore, generate
          // different tautomers.  If it does, we'll treat it as a new t_skel,
          // if it doesn't then it's not interesting so we'll discard it.
          store_existing_taut_gen( this_taut_gen , all_taut_gens[old_tg] ,
                                   all_taut_gens , next_new_tauts );
        }

      }
      if( next_new_tauts.empty() ) {
        break;
      }
      new_tauts = next_new_tauts;
      update_tauts_list( new_tauts , all_taut_smis , all_tauts );

#ifdef NOTYET
      cout << "Unique new tauts : " << new_tauts.size() << endl;
      for( size_t ii = 0 , iis = all_taut_smis.size() ; ii < iis ; ++ii ) {
        cout << all_taut_smis[ii] << endl;
      }
#endif
      if( new_tauts.empty() ) {
        break;
      }
      // EARLY BREAK
      // cout << "early break" << endl;
      // break;
      if( watch.Elapsed() > max_time ) {
        // cout << "timeout break" << endl;
        break;
      }
    }

#ifdef NOTYET
    cout << "Number of t_skels found : " << all_taut_gens.size() << endl;
    for( size_t ii = 0 , iis = all_taut_gens.size() ; ii < iis ; ++ii ) {
      cout << ii << " : " << all_taut_gens[ii]->global_t_skel_smi() << endl;
    }
#endif

    // Build the t_skel mol that encompasses all the tautomers found.
    // in create_global_t_skel, there's an option to see whether we're being too
    // enthusiastic with Hydrogen atoms.
    bool skip_h_check( watch.Elapsed() > max_time );
    pOEMolBase part_t_skel_mol( OENewMolBase( OEMolBaseType::OEDefault ) );
    create_global_t_skel( master_mol , global_is_aromatic , all_taut_gens ,
                          skip_h_check , *part_t_skel_mol );
    final_t_skel_mol += *part_t_skel_mol;
  }

#ifdef NOTYET
  cout << "Global t_skel : " << DACLIB::create_noncansmi( final_t_skel_mol ) << endl;
#endif

}

// ****************************************************************************
// returns the total number of tautomers, which includes duplicates
unsigned int make_taut_skeleton( const string &in_smi , const string &mol_name ,
                                 float max_time , string &t_skel_smi ) {

#ifdef NOTYET
  cout << "in_smi : " << in_smi << endl;
#endif

  vector<pTautGen> taut_gens;
  OEGraphMol final_t_skel_mol;

  generate_t_skel( in_smi , mol_name , true , max_time , taut_gens ,
                   final_t_skel_mol );

  unsigned int num_tauts = 0;
  for( size_t i = 0 , is = taut_gens.size() ; i < is ; ++i ) {
    num_tauts += taut_gens[i]->num_all_tautomers();
  }
  t_skel_smi = DACLIB::create_cansmi( final_t_skel_mol );

  return num_tauts;

#ifdef NOTYET
  cout << "t_skel_smi : " << t_skel_smi << endl;
#endif

}

// ****************************************************************************
// take the individual t_skel_mols and their associated atoms_for_hs and
// unsat_bond_idxs and generate all tautomer SMILES strings.
void enumerate_all_tautomers( vector<pTautGen> &taut_gens ,
                              vector<string> &taut_smis ) {

#ifdef NOTYET
  cout << BOOST_CURRENT_FUNCTION << endl;
#endif

  taut_smis.clear();

  for( size_t i = 0 , is = taut_gens.size() ; i < is ; ++i ) {
    vector<string> these_taut_smis = taut_gens[i]->generate_all_tautomer_smiles();
    taut_smis.insert( taut_smis.end() , these_taut_smis.begin() ,
                      these_taut_smis.end() );
  }

  sort( taut_smis.begin() , taut_smis.end() );
  taut_smis.erase( unique( taut_smis.begin() , taut_smis.end() ) ,
                   taut_smis.end() );

}

// ****************************************************************************
void make_taut_skeleton_and_tauts( const string &in_smi , const string &mol_name ,
                                   string &t_skel_smi ,
                                   vector<string> &taut_smis ) {

  vector<pTautGen> taut_gens;
  OEGraphMol final_t_skel_mol;

#ifdef NOTYET
  cout << BOOST_CURRENT_FUNCTION << endl;
#endif
  generate_t_skel( in_smi , mol_name , true , numeric_limits<float>::max() ,
                   taut_gens , final_t_skel_mol );
  t_skel_smi = DACLIB::create_cansmi( final_t_skel_mol );
  enumerate_all_tautomers( taut_gens , taut_smis );

#ifdef NOTYET
  cout << taut_smis.size() << " tautomers : " << endl;
  for( size_t i = 0 , is = taut_smis.size() ; i < is ; ++i ) {
    cout << "Taut_Smi : " << taut_smis[i] << endl;
  }
#endif

}
