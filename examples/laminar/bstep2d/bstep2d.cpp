/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006, 2007, 2012 Jonas Latt, Mathias J. Krause
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

/* bstep2d.cpp:
 * The implementation of a backward facing step. It is furthermore
 * shown how to use checkpointing to save the state of the
 * simulation regularly.
 */


#include "olb2D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb2D.hh"   // include full template code
#endif

#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;


typedef double T;
#define DESCRIPTOR D2Q9<>


// Parameters for the simulation setup
const T lx1   = 5.0;    // length of step in meter
const T ly1   = 0.75;   // height of step in meter
const T lx0   = 20.0;   // length of channel in meter
const T ly0   = 1.5;    // height of channel in meter
const int N = 60;        // resolution of the model
const int M = 50;        // resolution of the model
const T maxPhysT = 40.; // max. simulation time in s, SI unit


// Stores geometry information in form of material numbers
SuperGeometry2D<T> prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // set number of cuboids/blocks
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 6;
#endif

  // setup channel
  Vector<T,2> extendChannel( lx0,ly0 );
  Vector<T,2> originChannel( 0, 0 );
  std::shared_ptr<IndicatorF2D<T>> channel = std::make_shared<IndicatorCuboid2D<T>>( extendChannel, originChannel );

  // setup step
  Vector<T,2> extendStep( lx1,ly1 );
  Vector<T,2> originStep( 0, 0 );
  std::shared_ptr<IndicatorF2D<T>> step = std::make_shared<IndicatorCuboid2D<T>>( extendStep, originStep );

  CuboidGeometry2D<T>* cuboidGeometry = new CuboidGeometry2D<T>( *(channel-step), converter.getConversionFactorLength(), noOfCuboids );
  HeuristicLoadBalancer<T>* loadBalancer = new HeuristicLoadBalancer<T>( *cuboidGeometry );
  // Instantiation of a superGeometry
  SuperGeometry2D<T> superGeometry( *cuboidGeometry, *loadBalancer, 2 );

  // material numbers from zero to 2 inside geometry defined by indicator
  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,1,1 );

  Vector<T,2> extendBC( 0,ly0 );
  Vector<T,2> originBC;
  IndicatorCuboid2D<T> inflow( extendBC, originBC );
  // Set material number for inflow
  superGeometry.rename( 2,3,1,inflow );
  originBC[0] = lx0;
  IndicatorCuboid2D<T> outflow( extendBC, originBC );
  // Set material number for outflow
  superGeometry.rename( 2,4,1,outflow );
  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
  return superGeometry;
}

// Set up the geometry of the simulation
void prepareLattice( UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice2D<T,DESCRIPTOR>& sLattice,
                     Dynamics<T,DESCRIPTOR>& bulkDynamics,
                     sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& bc,
                     SuperGeometry2D<T>& superGeometry ) {
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  auto bulkIndicator = superGeometry.getMaterialIndicator({1, 3, 4});

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T,DESCRIPTOR>() );
  // Material=1 -->bulk dynamics
  // Material=3 -->bulk dynamics (inflow)
  // Material=4 -->bulk dynamics (outflow)
  sLattice.defineDynamics( bulkIndicator, &bulkDynamics );
  // Material=2 -->bounce back
  sLattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T,DESCRIPTOR>() );

  // Setting of the boundary conditions
  bc.addVelocityBoundary( superGeometry, 3, converter.getLatticeRelaxationFrequency() );
  bc.addPressureBoundary( superGeometry, 4, converter.getLatticeRelaxationFrequency() );

  // Initial conditions
  AnalyticalConst2D<T,T> ux( 0. );
  AnalyticalConst2D<T,T> uy( 0. );
  AnalyticalConst2D<T,T> rho( 1. );
  AnalyticalComposed2D<T,T> u( ux,uy );

  //Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rho, u );
  sLattice.iniEquilibrium( bulkIndicator, rho, u );

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperLattice2D<T,DESCRIPTOR>& sLattice, int iT,
                        SuperGeometry2D<T>& superGeometry ) {
  OstreamManager clout( std::cout,"setBoundaryValues" );

  // time for smooth start-up
  int iTmaxStart = converter.getLatticeTime( maxPhysT*0.2 );
  int iTupdate = 5;

  if ( iT%iTupdate == 0 && iT<= iTmaxStart ) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));
    // Smooth start curve, polynomial
    PolynomialStartScale<T,int> StartScale( iTmaxStart, T( 1 ) );
    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1] = {iT};
    T frac[1] = {};
    StartScale( frac,iTvec );
    T maxVelocity = converter.getCharLatticeVelocity()*3./2.*frac[0];
    T distance2Wall = converter.getConversionFactorLength()/2.;
    Poiseuille2D<T> poiseuilleU( superGeometry, 3, maxVelocity, distance2Wall );
    // define lattice speed on inflow
    sLattice.defineU( superGeometry, 3, poiseuilleU );
  }
}

// write data to termimal and file system
void getResults( SuperLattice2D<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry2D<T>& superGeometry, Timer<T>& timer,
                 SuperPlaneIntegralFluxVelocity2D<T>& velocityFlux,
                 SuperPlaneIntegralFluxPressure2D<T>& pressureFlux ) {
  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter2D<T> vtmWriter( "bstep2d" );

  if ( iT==0 ) {
    // Writes geometry, cuboid no. and rank no. to file system
    SuperLatticeGeometry2D<T,DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid2D<T,DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank2D<T,DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Writes every 0.2  seconds
  if ( iT%converter.getLatticeTime( 0.2 )==0 ) {
    SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice, converter );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    // write vtk to file system
    vtmWriter.write( iT );
    sLattice.communicate();
    SuperEuklidNorm2D<T,DESCRIPTOR> normVel( velocity );
    BlockReduction2D2D<T> planeReduction( normVel, 600, BlockDataSyncMode::ReduceOnly );
    // write output as JPEG
    heatmap::write(planeReduction, iT);
  }

  // Writes every 0.1 simulated
  if ( iT%converter.getLatticeTime( 0.1 )==0 ) {
    velocityFlux.print();
    pressureFlux.print();

    // write to terminal
    timer.update( iT );
    timer.printStep();
    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

  // Saves lattice data
  if ( iT%converter.getLatticeTime( 1 )==0 && iT>0 ) {
    clout << "Checkpointing the system at t=" << iT << std::endl;
    sLattice.save( "bstep2d.checkpoint" );
    // The data can be reloaded using
    //     sLattice.load("bstep2d.checkpoint");
  }
}


int main( int argc, char* argv[] ) {
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );  // set output directory
  OstreamManager clout( std::cout, "main" );

  UnitConverter<T,DESCRIPTOR> converter(
    (T)   1./N,     // physDeltaX: spacing between two lattice cells in __m__
    (T)   1./(M*N), // physDeltaT: time step in __s__
    (T)   1.,       // charPhysLength: reference length of simulation geometry
    (T)   1.,       // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1./500.,  // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.        // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("bstep2d");

  // === 2nd Step: Prepare Geometry ===
  // Instantiation of a superGeometry
  SuperGeometry2D<T> superGeometry( prepareGeometry(converter) );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice2D<T,DESCRIPTOR> sLattice( superGeometry );

  BGKdynamics<T,DESCRIPTOR> bulkDynamics (
    converter.getLatticeRelaxationFrequency(),
    instances::getBulkMomenta<T,DESCRIPTOR>()
  );

  // choose between local and non-local boundary condition
  sOnLatticeBoundaryCondition2D<T,DESCRIPTOR> sBoundaryCondition( sLattice );
  // createInterpBoundaryCondition2D<T,DESCRIPTOR>(sBoundaryCondition);
  createLocalBoundaryCondition2D<T,DESCRIPTOR>( sBoundaryCondition );

  prepareLattice( converter, sLattice, bulkDynamics, sBoundaryCondition, superGeometry );

  // instantiate reusable functors
  SuperPlaneIntegralFluxVelocity2D<T> velocityFlux( sLattice, converter, superGeometry, {19.,1.}, {0.,1.} );
  SuperPlaneIntegralFluxPressure2D<T> pressureFlux( sLattice, converter, superGeometry, {19.,1.}, {0.,1.} );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( int iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, sLattice, iT, superGeometry );
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer, velocityFlux, pressureFlux );
  }

  timer.stop();
  timer.printSummary();
}
