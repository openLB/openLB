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

/* bstep3d.cpp:
 * The implementation of a backward facing step. It is furthermore
 * shown how to use checkpointing to save the state of the
 * simulation regularly.
 */


#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb3D.hh"   // include full template code
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;
#define DESCRIPTOR D3Q19<>


// Parameters for the simulation setup
const T lx1   = 5.0;     // length of step
const T ly1   = 0.75;    // height of step
const T lx0   = 18.0;    // length of channel
const T ly0   = 1.5;     // height of channel
const T lz0   = 1.5;     // width of channel
const int N = 20;         // resolution of the model
const int M = 25;         // resolution of the model
const T maxPhysT = 40.;  // max. simulation time in s, SI unit


// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2 );

  superGeometry.rename( 2,1,1,1,1 );

  Vector<T,3> extend( lx1, ly1, lz0 );
  Vector<T,3> origin;
  IndicatorCuboid3D<T> cuboid2( extend, origin );

  superGeometry.rename( 1,2,cuboid2 );

  // Set material number for inflow
  extend = {2*converter.getConversionFactorLength(), ly0, lz0};
  origin[0] -= converter.getConversionFactorLength()/2.;
  IndicatorCuboid3D<T> inflow( extend, origin );
  superGeometry.rename( 2,3,1,inflow );

  // Set material number for outflow
  origin[0] = lx0 - converter.getConversionFactorLength()*1.5;
  IndicatorCuboid3D<T> outflow( extend, origin );
  superGeometry.rename( 2,4,1,outflow );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice3D<T,DESCRIPTOR>& sLattice,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& bc,
                     SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  // Material=3 -->bulk dynamics (inflow)
  // Material=4 -->bulk dynamics (outflow)
  auto bulkIndicator = superGeometry.getMaterialIndicator({1, 3, 4});
  sLattice.defineDynamics( bulkIndicator, &bulkDynamics );

  // Material=2 -->bounce back
  sLattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>() );

  // Setting of the boundary conditions
  bc.addVelocityBoundary( superGeometry, 3, omega );
  bc.addPressureBoundary( superGeometry, 4, omega );

  // Initial conditions
  AnalyticalConst3D<T,T> ux( 0. );
  AnalyticalConst3D<T,T> uy( 0. );
  AnalyticalConst3D<T,T> uz( 0. );
  AnalyticalConst3D<T,T> rho( 1. );
  AnalyticalComposed3D<T,T> u( ux,uy,uz );

  //Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rho, u );
  sLattice.iniEquilibrium( bulkIndicator, rho, u );

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperLattice3D<T,DESCRIPTOR>& sLattice, int iT,
                        SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime( maxPhysT*0.2 );
  int iTupdate = 5;

  if ( iT%iTupdate==0 && iT<= iTmaxStart ) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));

    // Smooth start curve, polynomial
    PolynomialStartScale<T,int> startScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1]= {iT};
    T frac[1]= {};
    startScale( frac,iTvec );
    std::vector<T> maxVelocity( 3,0 );
    maxVelocity[0] = 2.25*frac[0]*converter.getCharLatticeVelocity();

    T distance2Wall = converter.getConversionFactorLength()/2.;
    RectanglePoiseuille3D<T> poiseuilleU( superGeometry, 3, maxVelocity, distance2Wall, distance2Wall, distance2Wall );
    sLattice.defineU( superGeometry, 3, poiseuilleU );

    clout << "step=" << iT << "; maxVel=" << maxVelocity[0] << std::endl;
  }
}

// Output to console and files
void getResults( SuperLattice3D<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter,
                 BlockReduction3D2D<T>& planeReduction,
                 int iT,
                 SuperGeometry3D<T>& superGeometry, Timer<T>& timer)
{
  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "bstep3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  const int  vtkIter  = converter.getLatticeTime( 0.2 );
  const int  statIter = converter.getLatticeTime( 0.1 );
  const int  saveIter = converter.getLatticeTime( 1. );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Writes the ppm files
  if ( iT%vtkIter==0 ) {
    vtmWriter.write( iT );
    planeReduction.update();
    // write output as JPEG
    heatmap::write(planeReduction, iT);
  }

  // Writes output on the console
  if ( iT%statIter==0 && iT>=0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

  // Saves lattice data
  if ( iT%( saveIter/2 )==0 && iT>0 ) {
    clout << "Checkpointing the system at t=" << iT << endl;
    sLattice.save( "bstep3d.checkpoint" );
    // The data can be reloaded using
    //     sLattice.load("bstep3d.checkpoint");
  }
}

int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  UnitConverter<T,DESCRIPTOR> converter(
    (T)   1./N,     // physDeltaX: spacing between two lattice cells in __m__
    (T)   1./(M*N), // physDeltaT: time step in __s__
    (T)   1.,       // charPhysLength: reference length of simulation geometry
    (T)   1.,       // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1./100.,  // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.        // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("bstep3d");

  // === 2nd Step: Prepare Geometry ===
  Vector<T,3> extend( lx0, ly0, lz0 );
  Vector<T,3> origin;
  IndicatorCuboid3D<T> cuboid( extend, origin );

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 7;
#endif
  CuboidGeometry3D<T> cuboidGeometry( cuboid, converter.getConversionFactorLength(), noOfCuboids );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry3D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( converter, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

  BGKdynamics<T, DESCRIPTOR> bulkDynamics (
    converter.getLatticeRelaxationFrequency(),
    instances::getBulkMomenta<T,DESCRIPTOR>()
  );

  // choose between local and non-local boundary condition
  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sBoundaryCondition( sLattice );
  // createInterpBoundaryCondition3D<T,DESCRIPTOR>(sBoundaryCondition);
  createLocalBoundaryCondition3D<T,DESCRIPTOR>( sBoundaryCondition );

  prepareLattice( converter, sLattice, bulkDynamics, sBoundaryCondition, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << endl;
  Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  // Set up persistent measuring functors for result extraction
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperEuklidNorm3D<T, DESCRIPTOR> normVel( velocity );

  BlockReduction3D2D<T> planeReduction(
    normVel,
    Hyperplane3D<T>().centeredIn(cuboidGeometry.getMotherCuboid()).normalTo({0,0,1}),
    600,
    BlockDataSyncMode::ReduceOnly);

  for ( int iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, sLattice, iT, superGeometry );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, planeReduction, iT, superGeometry, timer );
  }

  timer.stop();
  timer.printSummary();
}
