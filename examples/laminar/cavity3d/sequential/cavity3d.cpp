/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2014 Mathias J. Krause
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

/* cavity3d.cpp:
 * This example illustrates a flow in a cuboid, lid-driven cavity.
 * This version is for sequential use. A version for parallel use
 * is also available.
 */


#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb3D.hh"   // include full template code
#endif

#include <cmath>
#include <iostream>
#include <fstream>


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
#define DESCRIPTOR D3Q19<>

const int N = 30; // resolution of the model
//const int M = 1; // time discretization refinement
const T maxT = (T) 100.; // max. simulation time in s, SI unit

const T interval = 1.0; // Time intervall in seconds for convergence check
const T epsilon = 1e-3; // Residuum for convergence check

void prepareGeometry( UnitConverter<T, DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator, BlockGeometry3D<T>& blockGeometry ) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // Sets material number for fluid and boundary
  blockGeometry.rename( 0,2,indicator );
  blockGeometry.rename( 2,1,1,1,1 );

  Vector<T,3> origin( T(), converter.getCharPhysLength(), T() );
  Vector<T,3> extend( converter.getCharPhysLength(), converter.getConversionFactorLength(), converter.getCharPhysLength() );
  IndicatorCuboid3D<T> lid( extend,origin );

  blockGeometry.rename( 2,3,1,lid );

  // Removes all not needed boundary voxels outside the surface
  blockGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  blockGeometry.innerClean();
  blockGeometry.checkForErrors();

  clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice( UnitConverter<T, DESCRIPTOR> const& converter,
                     BlockLatticeStructure3D<T,DESCRIPTOR>& lattice,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& bc,
                     BlockGeometry3D<T>& blockGeometry) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  lattice.defineDynamics( blockGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  lattice.defineDynamics( blockGeometry, 1, &bulkDynamics );

  // Material=2 -->bounce back
  //lattice.defineDynamics(superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>());

  // Material=2,3 -->bulk dynamics, velocity boundary
  lattice.defineDynamics( blockGeometry, 2, &bulkDynamics );
  lattice.defineDynamics( blockGeometry, 3, &bulkDynamics );
  bc.addVelocityBoundary( blockGeometry, 2, omega );
  bc.addVelocityBoundary( blockGeometry, 3, omega );

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues( UnitConverter<T, DESCRIPTOR> const& converter,
                        BlockLatticeStructure3D<T,DESCRIPTOR>& lattice, BlockGeometry3D<T>& blockGeometry, int iT ) {

  OstreamManager clout( std::cout,"setBoundaryValues" );

  if ( iT==0 ) {

    AnalyticalConst3D<T,T> rhoF( 1 );
    std::vector<T> velocity( 3,T() );
    AnalyticalConst3D<T,T> uF( velocity );

    lattice.iniEquilibrium( blockGeometry, 1, rhoF, uF );
    lattice.iniEquilibrium( blockGeometry, 2, rhoF, uF );
    lattice.iniEquilibrium( blockGeometry, 3, rhoF, uF );

    lattice.defineRhoU( blockGeometry, 1, rhoF, uF );
    lattice.defineRhoU( blockGeometry, 2, rhoF, uF );
    lattice.defineRhoU( blockGeometry, 3, rhoF, uF );

    velocity[0]=converter.getCharLatticeVelocity();
    AnalyticalConst3D<T,T> u( velocity );

    lattice.defineU( blockGeometry,3,u );

    // Make the lattice ready for simulation
    lattice.initialize();
  }
}

void getResults( BlockLatticeStructure3D<T,DESCRIPTOR>& lattice,
                 UnitConverter<T, DESCRIPTOR> const& converter, BlockGeometry3D<T>& blockGeometry, int iT, Timer<T>& timer, bool converged ) {

  OstreamManager clout( std::cout,"getResults" );
  BlockVTKwriter3D<T> vtkWriter( "cavity3d" );

  const T logT     = ( T )1.;
  const T vtkSave  = ( T )1.;

  if ( iT==0 ) {
    BlockLatticeGeometry3D<T, DESCRIPTOR> geometry( lattice, blockGeometry );
    vtkWriter.write( geometry );
  }

  // Get statistics
  if ( (iT%converter.getLatticeTime( logT )==0 && iT>0) || converged ) {
    timer.update( iT );
    timer.printStep( 2 );
    lattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

  // Writes the VTK files
  if ( (iT%converter.getLatticeTime( vtkSave )==0 && iT>0) || converged ) {

    BlockLatticePhysVelocity3D<T, DESCRIPTOR> velocity( lattice, 0, converter );
    BlockLatticePhysPressure3D<T, DESCRIPTOR> pressure( lattice, 0, converter );
    vtkWriter.addFunctor( velocity );
    vtkWriter.addFunctor( pressure );

    vtkWriter.write( iT );
  }
}



int main( int argc, char **argv ) {

  // === 1st Step: Initialization ===

  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},     // resolution: number of voxels per charPhysL
    (T)   0.509, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   1.0,   // charPhysLength: reference length of simulation geometry
    (T)   1.0,   // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.001, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0    // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("cavity3d");


  // === 2nd Step: Prepare Geometry ===

  // Instantiation of a unit cube by an indicator
  Vector<T,3> origin;
  Vector<T,3> extend( converter.getCharPhysLength() );
  IndicatorCuboid3D<T> cube( extend,origin );

  Cuboid3D<T> cuboid( cube, converter.getConversionFactorLength() );

  // Instantiation of a block geometry
  BlockGeometry3D<T> blockGeometry( cuboid );

  prepareGeometry( converter, cube, blockGeometry );


  // === 3rd Step: Prepare Lattice ===

  BlockLattice3D<T, DESCRIPTOR> lattice( blockGeometry.getNx(), blockGeometry.getNy(), blockGeometry.getNz(), blockGeometry );

  ConstRhoBGKdynamics<T, DESCRIPTOR> bulkDynamics (
    converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T,DESCRIPTOR>() );

  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>*
  boundaryCondition = createInterpBoundaryCondition3D<T,DESCRIPTOR,ConstRhoBGKdynamics<T,DESCRIPTOR> >( lattice );

  prepareLattice( converter, lattice, bulkDynamics, *boundaryCondition, blockGeometry );

  // === 4th Step: Main Loop with Timer ===
  util::ValueTracer<T> converge( converter.getLatticeTime(interval), epsilon );

  Timer<T> timer( converter.getLatticeTime( maxT ), std::pow<int>(converter.getResolution(),3) );
  timer.start();
  int iT;

  for ( iT=0; iT < converter.getLatticeTime( maxT ); ++iT ) {

    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << endl;
      getResults( lattice, converter, blockGeometry, iT, timer, converge.hasConverged() );
      break;
    }

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, lattice, blockGeometry, iT );

    // === 6th Step: Collide and Stream Execution ===
    lattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( lattice, converter, blockGeometry, iT, timer, converge.hasConverged() );
    converge.takeValue( lattice.getStatistics().getAverageEnergy(), true );
  }

  timer.stop();
  timer.printSummary();
}
