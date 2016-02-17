//
// file MolGridDisplay.H
// David Cosgrove
// AstraZeneca
// 26th June 2015
//
// This class takes a vector of SMILES strings and displays them in a grid of
// QTMolDisplay2D objects.

#ifndef MOLGRIDDISPLAY_H
#define MOLGRIDDISPLAY_H

#include <string>
#include <vector>

#include <QWidget>

// ****************************************************************************

class QGridLayout;
namespace DACLIB {
class QTMolDisplay2D;
}

// ****************************************************************************

namespace DACLIB {

class MolGridDisplay : public QWidget {

public :

  MolGridDisplay( QWidget *par = 0 , Qt::WindowFlags f = 0 );

  QSize sizeHint() const;

  void set_smiles( const std::vector<std::string> &new_smis );

private :

  std::vector<std::string> smiles_;

  QGridLayout *grid_;
  std::vector<QTMolDisplay2D *> disps_;

  void build_widget();

};

} // EO namespace DACLIB

#endif // MOLGRIDDISPLAY_H

