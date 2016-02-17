//
// File tt_tauts.cc
// David Cosgrove
// AstraZeneca
// 8th June 2015
//
// This is the main line for the program tt_tauts.cc. 
// The program implements the tautomer enumeration algoirthm of Thalheim et al.,
// 'A Branch-and-Bound Approach for Tautomer Enumeration',
// Molecular Informatics, 2015 - At present I have the ASAP version, without a full
// reference.  DOI: 10.1002/minf.201400128.
// It's behind a paywall.
//
// The purpose of this program is for testing.  It reads a molecule file, generates
// the tautomer skeleton and all tautomers, and displays them for inspection.

#include <iostream>

#include "TTTauts.H"

#include <QApplication>
#include <QDesktopWidget>
#include <QScreen>

#include <oechem.h>

using namespace std;

extern string BUILD_TIME; // in build_time.cc

// *****************************************************************************
int main( int argc , char **argv ) {

  cout << "tt_tauts" << endl
        << "Built " << BUILD_TIME << " using OEToolkits version "
        << OEChem::OEChemGetRelease() << "." << endl << endl
        << "Built with Qt version " << QT_VERSION_STR << endl
        << "Running with Qt version " << qVersion() << endl << endl;

  QApplication a( argc , argv );
  TTTauts *ttt = new TTTauts;
  ttt->resize( QApplication::desktop()->screenGeometry().width() * 3 / 5 ,
               QApplication::desktop()->screenGeometry().height() * 3 / 5 );
  ttt->move( QApplication::desktop()->screenGeometry().width() , 0 );

  ttt->show();
  ttt->parse_args( argc , argv );

  return a.exec();

}
