/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Peter Weisbrod
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

/* phaseSeparation3d.cpp:
 * In this example the simulation is initialized with a given
 * density plus a small random number all over the domain. This
 * condition is unstable and leads to liquid-vapor phase separation.
 * Boundaries are assumed to be periodic. This example shows the
 * usage of multiphase flow.
 */


#include "olb3D.h"
#include "olb3D.hh"   // use only generic version!
#include <cstdlib>
#include <iostream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;
#define DESCRIPTOR ShanChenDynOmegaForcedD3Q19Descriptor


// Parameters for the simulation setup
const int maxIter  = 2000;
const int nx   = 76;
const int ny   = 76;
const int nz   = 76;


// Stores geometry information in form of material numbers
void prepareGeometry( SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // Sets material number for fluid
  superGeometry.rename( 0,1 );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics1,
                     SuperGeometry3D<T>& superGeometry ) {

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 1, &bulkDynamics1 );

  // Initial conditions
  AnalyticalConst3D<T,T> noise( .01 );
  std::vector<T> v( 3,T() );
  AnalyticalConst3D<T,T> zeroVelocity( v );
  AnalyticalConst3D<T,T> oldRho( .125 );
  AnalyticalRandom3D<T,T> random;
  AnalyticalIdentity3D<T,T> newRho( random*noise+oldRho );

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( superGeometry, 1, newRho, zeroVelocity );
  sLattice.iniEquilibrium( superGeometry, 1, newRho, zeroVelocity );

  // Make the lattice ready for simulation
  sLattice.initialize();
}

// Output to console and files
void getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice, int iT,
                 SuperGeometry3D<T>& superGeometry, Timer<T>& timer ) {

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "phaseSeparation3d" );
  SuperLatticeVelocity3D<T, DESCRIPTOR> velocity( sLattice );
  SuperLatticeDensity3D<T, DESCRIPTOR> density( sLattice );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( density );

  const int vtkIter  = 20;
  const int statIter = 20;

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

  // Writes the vtk files
  if ( iT%vtkIter==0 ) {
    clout << "Writing VTK and JPEG..." << std::endl;
    vtmWriter.write( iT );

    BlockReduction3D2D<T> planeReduction( density, {0, 0, 1} );
    // write output as JPEG
    heatmap::write(planeReduction, iT);
  }

  // Writes output on the console
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,iT );
  }
}

int main( int argc, char *argv[] ) {

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  const T omega1 = 1.0;
  const T G      = -1.;

  // === 2rd Step: Prepare Geometry ===

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidGeometry3D<T> cuboidGeometry( 0, 0, 0, 1, nx, ny, nz, noOfCuboids );

  // Periodic boundaries in x- and y- and z-direction
  cuboidGeometry.setPeriodicity( true, true, true );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry3D<T> superGeometry( cuboidGeometry,loadBalancer,2 );

  prepareGeometry( superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

  ForcedShanChenBGKdynamics<T, DESCRIPTOR> bulkDynamics1 (
    omega1, instances::getExternalVelocityMomenta<T,DESCRIPTOR>() );

  std::vector<T> rho0;
  rho0.push_back( 1 );
  rho0.push_back( 1 );
  CarnahanStarling<T,T> interactionPotential( G );
  ShanChenForcedSingleComponentGenerator3D<T,DESCRIPTOR> coupling( G,rho0,interactionPotential );

  sLattice.addLatticeCoupling( coupling, sLattice );

  prepareLattice( sLattice, bulkDynamics1, superGeometry );

  // === 4th Step: Main Loop ===
  int iT = 0;
  clout << "starting simulation..." << endl;
  Timer<T> timer( maxIter, superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( iT = 0; iT < maxIter; ++iT ) {

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    // in this application no boundary conditions have to be adjusted

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    sLattice.communicate();
    sLattice.executeCoupling();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, iT, superGeometry, timer );
  }

  timer.stop();
  timer.printSummary();
}
