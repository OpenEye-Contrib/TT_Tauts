//
// file TTTauts.H
// David Cosgrove
// AstraZeneca
// 8th June 2015
//
// This is the interface for the class TTTauts, the main class in program tt_tauts,
// which, amongst other things, defines the GUI.

#ifndef TTTAUTS_H
#define TTTAUTS_H

#include <string>
#include <vector>

#include <QMainWindow>

class QScrollArea;
class QSlider;

namespace DACLIB {
class QTMolDisplay2D;
class MolGridDisplay;
}

namespace OEChem {
class OEAtomBase;
}

// ****************************************************************************

class TTTauts : public QMainWindow {

  Q_OBJECT

public :

  TTTauts();
  void parse_args( int argc , char **argv );

private slots :

  void slot_quit();
  void slot_read_molecules();
  void slot_slider_changed();

private :

  std::vector<std::string> mol_names_;
  std::vector<std::string> in_smiles_;
  std::vector<std::string> t_skel_smiles_;
  float max_time_;
  int max_tauts_;
  bool standardise_mols_;

  DACLIB::QTMolDisplay2D *mol_disp_ , *t_skel_disp_;
  DACLIB::MolGridDisplay *tauts_disp_;
  QSlider *mol_slider_;
  QScrollArea *tauts_disp_area_;
  QAction *file_quit_ , *file_read_act_;

  void build_actions();
  void build_menubar();
  void build_widget();

  void read_molecule_file( const std::string &filename );

  void make_t_skeleton( unsigned int mol_num );

};

#endif // TTTAUTS_H

