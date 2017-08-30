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

  TTTautsSettings();
  virtual ~TTTautsSettings() {}

  void parse_options(int argc, char **argv);
  void print_usage( std::ostream &os ) const;

  const std::string in_mol_file() const { return in_mol_file_; }
  const std::string out_mol_file() const { return out_mol_file_; }
  bool standardise_mols() const {return !dont_standardise_mols_;}
  float max_time() const {return max_time_;}
  int max_tauts() const {return max_tauts_;}
  const std::string &usage_text() { return usage_text_; }

protected:

  boost::program_options::options_description desc_;

private :

  std::string in_mol_file_;
  std::string out_mol_file_;
  bool dont_standardise_mols_;
  float max_time_;
  int max_tauts_;

  mutable std::string usage_text_;

  virtual void build_program_options();

};

#endif // TTTAUTSSETTINGS_H

