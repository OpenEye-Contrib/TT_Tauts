//
// file TTTautsSettings.H
// David Cosgrove
// AstraZeneca
// 8th June 2015
//
// Arguments parser for program tt_tauts.

#ifndef TTTAUTSSETTINGS_H
#define TTTAUTSSETTINGS_H

#include <iosfwd>
#include <string>
#include <boost/program_options/options_description.hpp>

// **************************************************************************

class TTTautsSettings {

public :

  TTTautsSettings( int argc , char **argv );

  void print_usage( std::ostream &os ) const;

  const std::string in_mol_file() const { return in_mol_file_; }
  const std::string out_mol_file() const { return out_mol_file_; }
  bool standardise_mols() const {return !dont_standardise_mols_;}
  float max_time() const {return max_time_;}
  int max_tauts() const {return max_tauts_;}
  const std::string &usage_text() { return usage_text_; }

private :

  std::string in_mol_file_;
  std::string out_mol_file_;
  bool dont_standardise_mols_;
  float max_time_;
  int max_tauts_;

  std::string usage_text_;

  void build_program_options( boost::program_options::options_description &desc );

};

#endif // TTTAUTSSETTINGS_H

