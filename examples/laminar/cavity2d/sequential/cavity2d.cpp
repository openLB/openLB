/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006 - 2012 Mathias J. Krause, Jonas Fietz,
 *  Jonas Latt, Jonas Kratzke
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

/* cavity2d.cpp:
 * This example illustrates a flow in a cuboid, lid-driven cavity.
 * It also shows how to use the XML parameter files and has an
 * example description file for OpenGPI. This version is for sequential
 * use. A version for parallel use is also available.
 */


#include "olb2D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb2D.hh"   // include full template code
#endif
#include <cmath>
#include <iostream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9<>

void prepareLattice( UnitConverter<T,DESCRIPTOR> const& converter,
                     BlockLatticeStructure2D<T,DESCRIPTOR>& lattice,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& bc ) {

  const int nx = lattice.getNx();
  const int ny = lattice.getNy();
  const T omega = converter.getLatticeRelaxationFrequency();

  // link lattice with dynamics for collision step
  lattice.defineDynamics( 0,nx-1, 0,ny-1, &bulkDynamics );

  // left boundary
  bc.addVelocityBoundary0N(   0,   0,   1,ny-2, omega );
  // right boundary
  bc.addVelocityBoundary0P( nx-1,nx-1,   1,ny-2, omega );
  // bottom boundary
  bc.addVelocityBoundary1N(   1,nx-2,   0,   0, omega );
  // top boundary
  bc.addVelocityBoundary1P(   1,nx-2,ny-1,ny-1, omega );

  // corners
  bc.addExternalVelocityCornerNN(   0,   0, omega );
  bc.addExternalVelocityCornerNP(   0,ny-1, omega );
  bc.addExternalVelocityCornerPN( nx-1,   0, omega );
  bc.addExternalVelocityCornerPP( nx-1,ny-1, omega );
}

void setBoundaryValues( UnitConverter<T,DESCRIPTOR> const& converter,
                        BlockLatticeStructure2D<T,DESCRIPTOR>& lattice, int iT ) {

  if ( iT==0 ) {

    const int nx = lattice.getNx();
    const int ny = lattice.getNy();

    // set initial values: v = [0,0]
    for ( int iX=0; iX<nx; ++iX ) {
      for ( int iY=0; iY<ny; ++iY ) {
        T vel[] = { T(), T()};
        lattice.get( iX,iY ).defineRhoU( ( T )1, vel );
        lattice.get( iX,iY ).iniEquilibrium( ( T )1, vel );
      }
    }

    // set non-zero velocity for upper boundary cells
    for ( int iX=1; iX<nx-1; ++iX ) {
      T u = converter.getCharLatticeVelocity();
      T vel[] = { u, T() };
      lattice.get( iX,ny-1 ).defineRhoU( ( T )1, vel );
      lattice.get( iX,ny-1 ).iniEquilibrium( ( T )1, vel );
    }

    // Make the lattice ready for simulation
    lattice.initialize();
  }
}

void getResults( BlockLatticeStructure2D<T,DESCRIPTOR>& lattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT, Timer<T>* timer,
                 const T logT, const T imSave, const T vtkSave,
                 std::string filenameGif, std::string filenameVtk,
                 const int timerPrintMode, const int timerTimeSteps, bool converged ) {

  // Get statistics
  if ( iT%converter.getLatticeTime( logT )==0 || converged ) {
    lattice.getStatistics().print( iT, converter.getPhysTime( iT ) );
  }

//  if ( iT%timerTimeSteps==0 || converged ) {
  if ( iT%timerTimeSteps==0 ) {
    timer->print( iT,timerPrintMode );
  }

  BlockVTKwriter2D<T> vtkWriter( filenameVtk );
  BlockLatticePhysVelocity2D<T,DESCRIPTOR> velocity( lattice,converter );
  BlockLatticePhysPressure2D<T,DESCRIPTOR> pressure( lattice,converter );
  vtkWriter.addFunctor( velocity );
  vtkWriter.addFunctor( pressure );

  // Writes the Gif files
  if ( ( iT%converter.getLatticeTime( imSave )==0 && iT>0 ) || converged ) {
    BlockEuklidNorm2D<T,DESCRIPTOR> normVel( velocity );
    BlockGifWriter<T> gifWriter;
    gifWriter.write( normVel, 0, 3, iT, filenameVtk );
//    gifWriter.write(normVel, iT, "vel");
  }

  // Writes the VTK files
  if ( ( iT%converter.getLatticeTime( vtkSave )==0 && iT>0 ) || converged ) {
    vtkWriter.write( iT );
  }
}


int main( int argc, char* argv[] ) {

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  OstreamManager clout( std::cout,"main" );

  string fName( "cavity2d.xml" );
  XMLreader config( fName );

  std::string olbdir = "../../";  //config["Application"]["OlbDir"].get<std::string>();
  std::string outputdir = "./tmp/"; //config["Output"]["OutputDir"].get<std::string>();
  singleton::directories().setOlbDir( olbdir );
  singleton::directories().setOutputDir( outputdir );

  // call creator functions using xml data
  UnitConverter<T,DESCRIPTOR>* converter = createUnitConverter<T,DESCRIPTOR>( config );
  // Prints the converter log as console output
  converter->print();
  // Writes the converter log in a file
  converter->write("cavity2d");

  int N = converter->getLatticeLength(1) + 1; // number of voxels in x,y,z direction
  Timer<T>* timer = createTimer<T>( config, *converter, N*N );

  // === 3rd Step: Prepare Lattice ===
  T logT = 0.1;   //config["Output"]["Log"]["SaveTime"].get<T>();
  T imSave = 1;   //config["Output"]["VisualizationImages"]["SaveTime"].get<T>();
  T vtkSave = 1;  //config["Output"]["VisualizationVTK"]["SaveTime"].get<T>();
  T maxPhysT = 100; //config["Application"]["PhysParam"]["MaxTime"].get<T>();
  int timerSkipType = 0;  //config["Output"]["Timer"]["SkipType"].get<T>();
  int timerPrintMode = 0; //config["Output"]["Timer"]["PrintMode"].get<int>();
  int timerTimeSteps = 1;

  if ( timerSkipType == 0 ) {
    timerTimeSteps = converter->getLatticeTime( .1 );
  }
//  else {
//    config["Output"]["Timer"]["TimeSteps"].read( timerTimeSteps );
//  }


  std::string filenameGif = "cavity2dimage";  //config["Output"]["VisualizationImages"]["Filename"].get<std::string>();
  std::string filenameVtk = "cavity2dvtk";    //config["Output"]["VisualizationVTK"]["Filename"].get<std::string>();

  BlockLattice2D<T, DESCRIPTOR> lattice( N, N );

  ConstRhoBGKdynamics<T, DESCRIPTOR> bulkDynamics (
    converter->getLatticeRelaxationFrequency(),
    instances::getBulkMomenta<T,DESCRIPTOR>()
  );

  OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
  boundaryCondition = createInterpBoundaryCondition2D<T,DESCRIPTOR,ConstRhoBGKdynamics<T,DESCRIPTOR> >( lattice );

  prepareLattice( *converter, lattice, bulkDynamics, *boundaryCondition );

  // === 4th Step: Main Loop with Timer ===

  int interval = converter->getLatticeTime( 1 /*config["Application"]["ConvergenceCheck"]["interval"].get<T>()*/ );
  T epsilon = 1e-3; //config["Application"]["ConvergenceCheck"]["residuum"].get<T>();
  util::ValueTracer<T> converge( interval, epsilon );

  timer->start();
  for ( int iT=0; iT <= converter->getLatticeTime( maxPhysT ); ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << endl;
      getResults( lattice, *converter, iT, timer, logT, imSave, vtkSave, filenameGif, filenameVtk, timerPrintMode, timerTimeSteps, converge.hasConverged() );

      break;
    }

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( *converter, lattice, iT );
    // === 6th Step: Collide and Stream Execution ===
    lattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults( lattice, *converter, iT, timer, logT, imSave, vtkSave, filenameGif, filenameVtk, timerPrintMode, timerTimeSteps, converge.hasConverged() );
    converge.takeValue( lattice.getStatistics().getAverageEnergy(), true );
  }

  timer->stop();
  timer->printSummary();
  delete converter;
  delete timer;
  delete boundaryCondition;
}
