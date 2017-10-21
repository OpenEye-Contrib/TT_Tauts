//
// file TTTautsBatchSettings.h
// David Cosgrove
// CozChemIx Limited
//
// Arguments parser for tt_tauts_batch, derived from TTTautsSettings.H

#ifndef TTTAUTSBATCHSETTINGS_H
#define TTTAUTSBATCHSETTINGS_H

#include "TTTautsSettings.H"

class TTTautsBatchSettings : public TTTautsSettings {

public :

    TTTautsBatchSettings();
    virtual ~TTTautsBatchSettings() {}

    int start_mol() const { return start_mol_; }
    int stop_mol() const { return stop_mol_; }
    int max_atoms() const { return max_atoms_; }

private :

    // start and stop numbers for molecules in file. Molecules start_mol_ to
    // stop_mol_ inclusive will be processed.
    int start_mol_;
    int stop_mol_;
    int max_atoms_; // maximum number of heavy atoms

    void build_program_options();

};

#endif // TTTAUTSBATCHSETTINGS_H
