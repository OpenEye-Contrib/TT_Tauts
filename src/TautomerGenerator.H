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

  TautomerGenerator( pOEMolBase &input_mol ,
                     std::vector<std::vector<int> > &mob_h ,
                     std::vector<std::vector<int> > &ts_bnds_to_1 ,
                     std::vector<std::vector<std::vector<int> > > &bnds_to_1 ,
                     std::vector<std::vector<std::vector<unsigned int> > > &ats_for_hs ,
                     std::vector<std::vector<std::vector<unsigned int> > > &us_bonds );

  pOEMolBase mol() { return mol_; }
  std::vector<std::string> t_skel_smis() const { return t_skel_smis_; }
  std::string global_t_skel_smi() const { return global_t_skel_smi_; }

  const std::vector<std::vector<int> > &mobile_h() const { return mobile_h_; }
  const std::vector<std::vector<int> > &t_skel_bonds_to_1() const { return t_skel_bonds_to_1_; }

  // Generate a set of tautomers, applying each mobile_h_ etc. in turn.
  // This will give a set of molecules where all tautomers of connect set 1 are
  // followed by all tautomerse of connect set 2 etc.
  // This will not generate all possible tautomers, which needs all combinations
  // of mobile_h_ etc to be applied and is a separate function.
  std::vector<pOEMolBase> generate_conn_set_tauts() const;
  std::vector<std::string> generate_conn_set_taut_smiles() const;
  unsigned int num_conn_set_tauts() const; // the maximum, including duplicates

  // this one does to all combinations of connect sets to return all
  // tautomers possible from the input molecule.
  std::vector<pOEMolBase> generate_all_tautomers() const;
  std::vector<std::string> generate_all_tautomer_smiles() const;
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
  // only 1 tautomer system there'll only be 1 vector in mobile_h_, for example.
  // Each connect set is dealt with separately, to generate its own t_skel.
  std::vector<std::vector<int> > mobile_h_;
  std::vector<std::vector<int> > t_skel_bonds_to_1_;

  std::vector<std::vector<std::vector<int> > > bonds_to_1_;
  std::vector<std::vector<std::vector<unsigned int > > > atoms_for_hs_;
  std::vector<std::vector<std::vector<unsigned int> > > unsat_bond_idxs_;

  std::vector<pOEMolBase> t_skel_mols_;
  std::vector<std::string> t_skel_smis_;

  mutable std::vector<std::string> conn_set_taut_smis_;

  // these are the global t_skel details, formed by combining the individual
  // ones above
  pOEMolBase global_t_skel_mol_;
  std::string global_t_skel_smi_;

  // make the individual tautomer skeletons and the global one
  void generate_t_skels();
  void make_global_t_skel();

};

void build_t_skel_mol( const std::vector<int> &mobile_h ,
                       const std::vector<int> &bonds_to_1 ,
                       OEChem::OEMolBase &t_skel_mol );
void remove_h_from_t_skel( const std::vector<int> &mobile_h ,
                           OEChem::OEMolBase &t_skel_mol );
void set_bonds_to_1( const std::vector<int> &bonds_to_1 ,
                     OEChem::OEMolBase &mol );
void generate_tautomers( const OEChem::OEMolBase &master_mol ,
                         const std::vector<int> &mobile_h ,
                         const std::vector<std::vector<unsigned int> > &atoms_for_hs ,
                         const std::vector<std::vector<unsigned int> > &unsat_bond_idxs ,
                         const std::vector<std::vector<int> > &bonds_to_1 ,
                         std::vector<pOEMolBase > &tauts );

#endif // TAUTOMERGENERATOR_H
