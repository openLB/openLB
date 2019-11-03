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
 * example description file for OpenGPI. This version is for parallel
 * use. A version for sequential use is also available.
 */


#include "olb2D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb2D.hh"   // include full template code
#endif
#include <vector>
#include <cmath>
#include <iostream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9<>

void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry2D<T>& superGeometry ) {
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;


  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,1,1 );
  superGeometry.clean();

  T eps = converter.getConversionFactorLength();
  Vector<T,2> extend( T( 1 ) + 2*eps, 2*eps );
  Vector<T,2> origin( T() - eps, T( 1 ) - eps );
  IndicatorCuboid2D<T> lid( extend, origin );
  // Set material number for lid
  superGeometry.rename( 2,3,1,lid );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice2D<T, DESCRIPTOR>& sLattice,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& sBoundaryCondition,
                     SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // link lattice with dynamics for collision step

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 1, &bulkDynamics );

  // Material=2,3 -->bulk dynamics, velocity boundary
  sLattice.defineDynamics( superGeometry, 2, &bulkDynamics );
  sLattice.defineDynamics( superGeometry, 3, &bulkDynamics );
  sBoundaryCondition.addVelocityBoundary( superGeometry, 2, omega );
  sBoundaryCondition.addVelocityBoundary( superGeometry, 3, omega );

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues( UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperLattice2D<T, DESCRIPTOR>& sLattice,
                        int iT, SuperGeometry2D<T>& superGeometry ) {

  if ( iT==0 ) {
    // set initial values: v = [0,0]
    AnalyticalConst2D<T,T> rhoF( 1 );
    std::vector<T> velocity( 2,T() );
    AnalyticalConst2D<T,T> uF( velocity );

    auto bulkIndicator = superGeometry.getMaterialIndicator({1, 2, 3});
    sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );
    sLattice.defineRhoU( bulkIndicator, rhoF, uF );

    // set non-zero velocity for upper boundary cells
    velocity[0] = converter.getCharLatticeVelocity();
    AnalyticalConst2D<T,T> u( velocity );
    sLattice.defineU( superGeometry, 3, u );

    // Make the lattice ready for simulation
    sLattice.initialize();
  }
}

void getResults( SuperLattice2D<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT, Timer<T>* timer,
                 const T logT, const T maxPhysT, const T imSave, const T vtkSave,
                 std::string filenameGif, std::string filenameVtk,
                 const int timerPrintMode,
                 const int timerTimeSteps, SuperGeometry2D<T>& superGeometry, bool converged ) {

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter2D<T> vtmWriter( filenameVtk );


  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Get statistics
  if ( iT%converter.getLatticeTime( logT )==0 || converged ) {
    sLattice.getStatistics().print( iT, converter.getPhysTime( iT ) );
  }

  if ( iT%timerTimeSteps==0 || converged ) {
    timer->print( iT,timerPrintMode );
  }

  // Writes the VTK files
  if ( ( iT%converter.getLatticeTime( vtkSave )==0 && iT>0 ) || converged ) {
    SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice, converter );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.write( iT );
  }

  // Writes the Gif files
  if ( ( iT%converter.getLatticeTime( imSave )==0 && iT>0 ) || converged ) {
    SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity( sLattice, converter );
    SuperEuklidNorm2D<T,DESCRIPTOR> normVel( velocity );
    BlockReduction2D2D<T> planeReduction( normVel, 600, BlockDataSyncMode::ReduceOnly );
    // write output of velocity as JPEG
    heatmap::write(planeReduction, iT);
  }

  // Output for x-velocity along y-position at the last time step
  if ( iT == converter.getLatticeTime( maxPhysT ) || converged ) {
    // Gives access to velocity information on lattice
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocityField( sLattice, converter );
    // Interpolation functor with velocityField information
    AnalyticalFfromSuperF2D<T> interpolation( velocityField, true, 1 );

    Vector<int,17> y_coord( {128, 125, 124, 123, 122, 109, 94, 79, 64, 58, 36, 22, 13, 9, 8, 7, 0} );
    // Ghia, Ghia and Shin, 1982: "High-Re Solutions for Incompressible Flow Using the Navier-Stokes Equations and a Multigrid Method";  Table 1
    Vector<T,17> vel_ghia_RE1000( { 1.0,     0.65928, 0.57492, 0.51117, 0.46604,
                                    0.33304, 0.18719, 0.05702,-0.06080,-0.10648,
                                    -0.27805,-0.38289,-0.29730,-0.22220,-0.20196,
                                    -0.18109, 0.0
                                  } );
    Vector<T,17> vel_ghia_RE100( {1.0,     0.84123, 0.78871, 0.73722, 0.68717,
                                  0.23151, 0.00332,-0.13641,-0.20581,-0.21090,
                                  -0.15662,-0.10150,-0.06434,-0.04775,-0.04192,
                                  -0.03717, 0.0
                                 } );
    Vector<T,17> vel_simulation;

    // Gnuplot interface to create plots
    static Gnuplot<T> gplot( "centerVelocityX" );
    // Define comparison values
    Vector<T,17> comparison = vel_ghia_RE1000;

    for ( int nY = 0; nY < 17; ++nY ) {
      // 17 data points evenly distributed between 0 and 1 (height)
      T position[2] = {0.5, y_coord[nY]/128.0};
      T velocity[2] = {T(), T()};
      // Interpolate velocityField at "position" and save it in "velocity"
      interpolation( velocity, position );
      // Save value of velocity (in x-direction) in "vel_simulation" for every position "nY"
      vel_simulation[nY] = velocity[0];
      // Set data for plot output
      gplot.setData( position[1], {vel_simulation[nY],comparison[nY]}, {"simulated","Ghia"} );
    }
    // Create PNG file
    gplot.writePNG();
    // Console output with results
    clout << "absoluteErrorL2(line)=" << ( vel_simulation - comparison ).norm() / 17. << "; relativeErrorL2(line)=" << ( vel_simulation - comparison ).norm() / comparison.norm() << std::endl;
  }
}



int main( int argc, char* argv[] ) {

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  OstreamManager clout( std::cout,"main" );

  string fName( "cavity2d.xml" );
  XMLreader config( fName );

  std::string olbdir, outputdir;
  config["Application"]["OlbDir"].read( olbdir );
  config["Output"]["OutputDir"].read( outputdir );
  singleton::directories().setOlbDir( olbdir );
  singleton::directories().setOutputDir( outputdir );

  UnitConverter<T,DESCRIPTOR>* converter = createUnitConverter<T,DESCRIPTOR>( config );
  // Prints the converter log as console output
  converter->print();
  // Writes the converter log in a file
  converter->write("cavity2d");

  int N = converter->getLatticeLength(1) + 1; // number of voxels in x,y,z direction
  Timer<T>* timer = createTimer<T>( config, *converter, N*N );


  T logT = config["Output"]["Log"]["SaveTime"].get<T>();
  T imSave = config["Output"]["VisualizationImages"]["SaveTime"].get<T>();
  T vtkSave = config["Output"]["VisualizationVTK"]["SaveTime"].get<T>();
  T maxPhysT = config["Application"]["PhysParameters"]["PhysMaxTime"].get<T>();
  int timerSkipType = config["Output"]["Timer"]["SkipType"].get<T>();
  int timerPrintMode = config["Output"]["Timer"]["PrintMode"].get<int>();
  int timerTimeSteps = 1;

  if ( timerSkipType == 0 ) {
    timerTimeSteps = converter->getLatticeTime( 1. /*config["Output"]["Timer"]["PhysTime"].get<T>()*/ );
  } else {
//    config["Output"]["Timer"]["TimeSteps"].read( timerTimeSteps );
  }

  std::string filenameGif = config["Output"]["VisualizationImages"]["Filename"].get<std::string>();
  std::string filenameVtk = config["Output"]["VisualizationVTK"]["Filename"].get<std::string>();

  // === 2rd Step: Prepare Geometry ===
  Vector<T,2> extend( 1,1 );
  Vector<T,2> origin( 0,0 );
  IndicatorCuboid2D<T> cuboid( extend, origin );

#ifdef PARALLEL_MODE_MPI
  CuboidGeometry2D<T> cuboidGeometry( cuboid, converter->getConversionFactorLength(), singleton::mpi().getSize() );
#else
  CuboidGeometry2D<T> cuboidGeometry( cuboid, converter->getConversionFactorLength(), 7 );
#endif

  cuboidGeometry.print();

  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );
  SuperGeometry2D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );
  prepareGeometry( *converter, superGeometry );

  // === 3rd Step: Prepare Lattice ===

  SuperLattice2D<T, DESCRIPTOR> sLattice( superGeometry );

  ConstRhoBGKdynamics<T, DESCRIPTOR> bulkDynamics (
    converter->getLatticeRelaxationFrequency(),
    instances::getBulkMomenta<T,DESCRIPTOR>()
  );

  sOnLatticeBoundaryCondition2D<T,DESCRIPTOR> sBoundaryCondition( sLattice );
  createInterpBoundaryCondition2D<T,DESCRIPTOR,ConstRhoBGKdynamics<T,DESCRIPTOR> > ( sBoundaryCondition );

  prepareLattice( *converter, sLattice, bulkDynamics, sBoundaryCondition, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  int interval = converter->getLatticeTime( 1 /*config["Application"]["ConvergenceCheck"]["interval"].get<T>()*/ );
  T epsilon = 1e-3; //config["Application"]["ConvergenceCheck"]["residuum"].get<T>();
  util::ValueTracer<T> converge( interval, epsilon );

  timer->start();
  for ( int iT=0; iT <= converter->getLatticeTime( maxPhysT ); ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << endl;
      getResults( sLattice, *converter, iT, timer, logT, maxPhysT, imSave, vtkSave, filenameGif, filenameVtk,
                  timerPrintMode, timerTimeSteps, superGeometry, converge.hasConverged() );
      break;
    }
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( *converter, sLattice, iT, superGeometry );
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, *converter, iT, timer, logT, maxPhysT, imSave, vtkSave, filenameGif, filenameVtk,
                timerPrintMode, timerTimeSteps, superGeometry, converge.hasConverged() );
    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
  }
  timer->stop();
  timer->printSummary();
  delete converter;
  delete timer;
}
