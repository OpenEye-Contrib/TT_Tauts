//
// file MolGridDisplay.cc
// David Cosgrove
// AstraZeneca
// 26th June 2015
//

#include "MolGridDisplay.H"
#include "QTMolDisplay2D.H"

#include <QGridLayout>

#include <iostream>

#include <oechem.h>

using namespace std;
using namespace OEChem;

namespace DACLIB {

// ****************************************************************************
MolGridDisplay::MolGridDisplay( QWidget *par , Qt::WindowFlags f ) :
  QWidget( par , f ) {

  build_widget();

}

// ****************************************************************************
QSize MolGridDisplay::sizeHint() const {

  if( disps_.empty() ) {
    return QSize( 200 , 200 );
  }

  QSize msh = disps_.front()->sizeHint();
  int w = grid_->columnCount() * msh.width();
  int h = grid_->rowCount() * msh.height();

  return QSize( w , h );

}

// ****************************************************************************
void MolGridDisplay::set_smiles( const vector<string> &new_smis ) {

  smiles_ = new_smis;
  cout << smiles_.size() << " mols to show" << endl;
  unsigned int n_col = static_cast<unsigned int>( ( sqrt( double( smiles_.size() ) ) ) );
  if( n_col * n_col < smiles_.size() ) {
    ++n_col; // prefer wider grid over longer
  }
#ifdef NOTYET
  cout << "num cols : " << n_col << endl;
#endif

  for( size_t i = 0 , is = disps_.size() ; i < is ; ++i ) {
    disps_[i]->hide();
    grid_->removeWidget( disps_[i] );
  }

  for( unsigned int i = 0 , is = static_cast<unsigned int>( smiles_.size() ) ; i < is ; ++i ) {
#ifdef NOTYET
    cout << i << " : " << smiles_[i] << endl;
#endif
    if( i == disps_.size() ) {
      disps_.push_back( new DACLIB::QTMolDisplay2D );
    }
    OEMolBase *mol = OENewMolBase( OEMolBaseType::OEDefault );
    OEParseSmiles( *mol , smiles_[i] );
    disps_[i]->set_display_molecule( mol );
    delete mol; // QTMolDisplay2D takes a copy

    disps_[i]->show();
    unsigned int r = i / n_col;
    unsigned int c = i % n_col;
    grid_->addWidget( disps_[i] , r , c );
#ifdef NOTYET
    cout << "done " << i << endl;
#endif
  }

}

// ****************************************************************************
void MolGridDisplay::build_widget() {

  grid_ = new QGridLayout;
  setLayout( grid_ );

}

} // EO namespace DACLIB
