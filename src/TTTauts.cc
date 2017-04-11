//
// file TTTauts.cc
// David Cosgrove
// AstraZeneca
// 8th June 2015
//

#include "TTTauts.H"
#include "TTTautsSettings.H"

#include "FileExceptions.H"
#include "MolGridDisplay.H"
#include "QTMolDisplay2D.H"

#include <QAction>
#include <QLayout>
#include <QMenu>
#include <QMenuBar>
#include <QScrollArea>
#include <QSlider>
#include <QSplitter>
#include <QStatusBar>

#include <oechem.h>

#include <boost/lexical_cast.hpp>

using namespace std;
using namespace OEChem;
using namespace OESystem;

// in make_taut_skeleton.cpp
void make_taut_skeleton_and_tauts( const string &in_smi , const string &mol_name ,
                                   string &t_skel_smi ,
                                   vector<string> &taut_smis,
                                   bool &timed_out ,
                                   float max_time = std::numeric_limits<float>::max() );

// ****************************************************************************
TTTauts::TTTauts() : QMainWindow() {

  build_actions();
  build_menubar();
  build_widget();

}

// ****************************************************************************
void TTTauts::parse_args( int argc , char **argv ) {

  TTTautsSettings ttts( argc , argv );
  if( !ttts.in_mol_file().empty() ) {
    read_molecule_file( ttts.in_mol_file() );
  }

}

// ****************************************************************************
void TTTauts::slot_quit() {

  exit( 0 );

}

// ****************************************************************************
void TTTauts::slot_read_molecules() {

}

// ****************************************************************************
void TTTauts::slot_slider_changed() {

  int slider_val = mol_slider_->value();

  if( mol_names_.empty() ) {
    return;
  }

  OEMol mol , skel;
  OEParseSmiles( mol , in_smiles_[slider_val] );
  mol.SetTitle( mol_names_[slider_val] );
  mol_disp_->set_display_molecule( &mol.SCMol() );

  make_t_skeleton( slider_val );
  OEParseSmiles( skel , t_skel_smiles_[slider_val] );
  skel.SetTitle( mol_names_[slider_val] );
  t_skel_disp_->set_display_molecule( &skel.SCMol() );

}

// ****************************************************************************
void TTTauts::build_actions() {

  file_quit_ = new QAction( "Quit" , this );
  file_quit_->setShortcut( QString( "Ctrl+Q" ) );
  file_quit_->setStatusTip( "Exit the application" );
  connect( file_quit_ , &QAction::triggered , this , &TTTauts::slot_quit );

  file_read_act_ = new QAction( "Read Molecules" , this );
  file_read_act_->setStatusTip( "Read molecule file for input structures." );
  connect( file_read_act_ , &QAction::triggered ,
           this , &TTTauts::slot_read_molecules );

}

// ****************************************************************************
void TTTauts::build_menubar() {

  QMenu *file_menu = menuBar()->addMenu( "File" );
  file_menu->addAction( file_read_act_ );
  file_menu->addSeparator();
  file_menu->addAction( file_quit_ );

}

// ****************************************************************************
void TTTauts::build_widget() {

  mol_slider_ = new QSlider;
  mol_slider_->setInvertedAppearance( true );
  mol_slider_->setInvertedControls( true );
  mol_slider_->setSingleStep( 1 );
  mol_slider_->setPageStep( 1 );
  mol_slider_->setEnabled( false );
  connect( mol_slider_ , &QSlider::valueChanged ,
           this , &TTTauts::slot_slider_changed );

  mol_disp_ = new DACLIB::QTMolDisplay2D;
  mol_disp_->number_atoms_by_default( true );
  t_skel_disp_ = new DACLIB::QTMolDisplay2D;
  t_skel_disp_->number_atoms_by_default( true );

  tauts_disp_ = new DACLIB::MolGridDisplay;

  tauts_disp_area_ = new QScrollArea;
  tauts_disp_area_->setWidget( tauts_disp_ );

  QHBoxLayout *hbox = new QHBoxLayout;
  hbox->addWidget( mol_slider_ );
  QSplitter *split = new QSplitter;
  split->addWidget( mol_disp_ );
  split->addWidget( t_skel_disp_ );
  split->addWidget( tauts_disp_area_ );
  split->setStretchFactor( 0 , 1 );
  split->setStretchFactor( 1 , 1 );
  split->setStretchFactor( 2 , 2 );
  hbox->addWidget( split );

#ifdef NOTYET
  hbox->addWidget( mol_disp_ );
  hbox->addWidget( t_skel_disp_ );
  hbox->addWidget( tauts_disp_area_ );
#endif

  QWidget *main_wid = new QWidget;
  main_wid->setLayout( hbox );

  setCentralWidget( main_wid );

}

// ****************************************************************************
void TTTauts::read_molecule_file( const string &filename ) {

  oemolistream ims;
  if( !ims.open( filename ) ) {
    throw DACLIB::FileReadOpenError( filename.c_str() );
  }

  OEMol mol;
  while( ims >> mol ) {
    string smi;
    OECreateSmiString( smi , mol , OESMILESFlag::AtomStereo | OESMILESFlag::BondStereo | OESMILESFlag::Canonical );
    in_smiles_.push_back( smi );
    mol_names_.push_back( mol.GetTitle() );
    if( mol_names_.back().empty() ) {
      mol_names_.back() = string( "Str_" ) + boost::lexical_cast<string>( mol_names_.size() );
    }
    t_skel_smiles_.push_back( smi );
  }

  if( mol_names_.size() > 1 ) {
    mol_slider_->setRange( 0 , int( mol_names_.size() ) - 1 );
    mol_slider_->setEnabled( true );
  }

  slot_slider_changed();
  statusBar()->showMessage( QString( "Number of molecules now : %1." ).arg( mol_names_.size() ) );

}

// ****************************************************************************
void TTTauts::make_t_skeleton( unsigned int mol_num ) {

#ifdef NOTYET
  cout << "make_t_skeleton for " << mol_names_[mol_num] << endl;
#endif

  string tss;
  vector<string> taut_smis;
  bool timed_out = false;
  make_taut_skeleton_and_tauts( in_smiles_[mol_num] , mol_names_[mol_num] ,
                                tss , taut_smis , timed_out );
  t_skel_smiles_[mol_num] = tss;

  cout << "input SMILES : " << in_smiles_[mol_num] << " t_skel_smi : " << tss << " for " << mol_names_[mol_num];
  if( timed_out ) {
    cout << " but timed out.";
  }
  cout << endl;
  if( taut_smis.size() < 100 ) {
    tauts_disp_->set_smiles( taut_smis );

    QSize td_size_hint = tauts_disp_->sizeHint();
    QSize new_size( td_size_hint.width() > tauts_disp_area_->width() ? td_size_hint.width() : tauts_disp_area_->width() ,
                    td_size_hint.height() > tauts_disp_area_->height() ? td_size_hint.height() : tauts_disp_area_->height());

    tauts_disp_->resize( new_size );
  } else {
    cout << "Found " << taut_smis.size() << " tautomers - too many to show sensibly." << endl;
    taut_smis.clear();
    tauts_disp_->set_smiles( taut_smis );
  }

}
