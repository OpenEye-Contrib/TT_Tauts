//
// File TautomerGenerator.H
// David Cosgrove
// AstraZeneca
// 16th December 2016
//
// This class holds the means to generate all the tautmers from an input
// structure, and its tautomer skeleton.
// The analysis to arrive at its data is done elsewhere, this is for the most
// part a convenient container to avoid what was turning into a complicated
// and impenetrable set of vectors of vectors of vectors of vectors.

#ifndef TAUTOMERGENERATOR_H
#define TAUTOMERGENERATOR_H

#include <string>
#include <vector>

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>

namespace OEChem {
class OEMolBase;
}

typedef boost::shared_ptr<OEChem::OEMolBase> pOEMolBase;

// ****************************************************************************
class TautomerGenerator {

public :

  // for when there's no tautomers, we just want a copy of the input molecule
  TautomerGenerator( pOEMolBase &input_mol );

  TautomerGenerator( pOEMolBase &input_mol ,
                     const std::vector<std::vector<int> > &ts_bnds_to_1 ,
                     const std::vector<std::vector<std::vector<int> > > &bnds_to_1 ,
                     const std::vector<std::vector<std::pair<unsigned int,unsigned int> > > &h_mov ,
                     const std::vector<std::vector<std::vector<unsigned int> > > &us_bonds );

  pOEMolBase mol() { return mol_; }
  std::vector<std::string> t_skel_smis() const { return t_skel_smis_; }
  std::string global_t_skel_smi() const { return global_t_skel_smi_; }

  const std::vector<std::vector<std::pair<unsigned int,unsigned int> > > &h_moves() const { return h_moves_; }
  const std::vector<std::vector<int> > &t_skel_bonds_to_1() const { return t_skel_bonds_to_1_; }

  // Generate a set of tautomers, applying each h_moves_ etc. in turn.
  // This will give a set of molecules where all tautomers of connect set 1 are
  // followed by all tautomers of connect set 2 etc.
  // This will not generate all possible tautomers, which needs all combinations
  // of h_moves_ etc to be applied and is a separate function.
  std::vector<pOEMolBase> generate_conn_set_tauts() const;
  std::vector<std::string> generate_conn_set_taut_smiles() const;
  unsigned int num_conn_set_tauts() const; // the maximum, including duplicates

  // this one does to all combinations of connect sets to return all
  // tautomers possible from the input molecule.
  std::vector<pOEMolBase> generate_all_tautomers();
  std::vector<std::string> generate_all_tautomer_smiles();
  unsigned int num_all_tautomers() const; // the maximum, including duplicates

  // return accumulation of everything in bonds_to_1_
  std::vector<std::vector<int> > get_all_bonds_to_1() const;
  // return accumulation of everything in unsat_bond_idxs_
  std::vector<std::vector<unsigned int> > get_all_unsat_bond_idxs() const;

  bool operator==( const TautomerGenerator &rhs );
  bool operator!=( const TautomerGenerator &rhs );

  // take the TautomerGenerator old_one passed in, and take out any
  // transformations in this one that produces a tautomer that's already
  // in old_one
  void prune( const TautomerGenerator &old_one );

private :

  pOEMolBase mol_;
  std::string can_smi_;

  // each connect set/tautomer system has its own outer vector. If there is
  // only 1 tautomer system there'll only be 1 vector in t_skel_bonds_to_1,
  // for example. Each connect set is dealt with separately, to generate
  // its own t_skel.
  std::vector<std::vector<int> > t_skel_bonds_to_1_;

  std::vector<std::vector<std::vector<int> > > bonds_to_1_;
  std::vector<std::vector<std::pair<unsigned int,unsigned int> > > h_moves_;
  std::vector<std::vector<std::vector<unsigned int> > > unsat_bond_idxs_;

  std::vector<pOEMolBase> t_skel_mols_;
  std::vector<std::string> t_skel_smis_;

  mutable std::vector<std::string> conn_set_taut_smis_;

  // these are the global t_skel details, formed by combining the individual
  // ones above
  pOEMolBase global_t_skel_mol_;
  std::string global_t_skel_smi_;

  // remove any instances where we set a bond to zero and then set it back
  // to what it was before.
  void normalise_bond_changes();

  // make the individual tautomer skeletons and the global one
  void generate_t_skels();
  void make_global_t_skel();
  // prune() can change the h_moves_ etc, after which we need to rebuild the
  // t_skel stuff.
  void rebuild_t_skels();

};

void build_t_skel_mol( unsigned int h_from ,
                       const std::vector<int> &bonds_to_1 ,
                       OEChem::OEMolBase &t_skel_mol );
void build_t_skel_mol( const std::vector<int> &mobile_h ,
                       const std::vector<int> &bonds_to_1 ,
                       OEChem::OEMolBase &t_skel_mol );
void build_t_skel_mol( const std::vector<std::pair<unsigned int,unsigned int> > &h_moves ,
                       const std::vector<int> &bonds_to_1 ,
                       OEChem::OEMolBase &t_skel_mol );
void remove_h_from_t_skel( unsigned int h_from ,
                           OEChem::OEMolBase &t_skel_mol );
void set_bonds_to_1( const std::vector<int> &bonds_to_1 ,
                     OEChem::OEMolBase &mol );
void generate_tautomers( const OEChem::OEMolBase &master_mol ,
                         const std::vector<std::pair<unsigned int, unsigned int> > &h_moves ,
                         const std::vector<std::vector<unsigned int> > &unsat_bond_idxs ,
                         const std::vector<std::vector<int> > &bonds_to_1 ,
                         std::vector<pOEMolBase > &tauts );

#endif // TAUTOMERGENERATOR_H
