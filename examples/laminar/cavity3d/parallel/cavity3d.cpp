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
 * This version is for parallel use. A version for sequential use
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

void prepareGeometry( UnitConverter<T, DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator, SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // Sets material number for fluid and boundary
  superGeometry.rename( 0,2,indicator );
  superGeometry.rename( 2,1,1,1,1 );

  T eps = converter.getConversionFactorLength();
  Vector<T,3> origin( -eps, converter.getCharPhysLength() - eps, -eps );
  Vector<T,3> extend( converter.getCharPhysLength() + 2*eps, 2*eps, converter.getCharPhysLength() + 2*eps );
  IndicatorCuboid3D<T> lid( extend,origin );

  superGeometry.rename( 2,3,1,lid );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice( UnitConverter<T, DESCRIPTOR> const& converter,
                     SuperLattice3D<T,DESCRIPTOR>& lattice,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& bc,
                     SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  lattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  lattice.defineDynamics( superGeometry, 1, &bulkDynamics );

  // Material=2,3 -->bulk dynamics, velocity boundary
  lattice.defineDynamics( superGeometry, 2, &bulkDynamics );
  lattice.defineDynamics( superGeometry, 3, &bulkDynamics );
  bc.addVelocityBoundary( superGeometry, 2, omega );
  bc.addVelocityBoundary( superGeometry, 3, omega );

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues( UnitConverter<T, DESCRIPTOR> const& converter,
                        SuperLattice3D<T,DESCRIPTOR>& lattice, SuperGeometry3D<T>& superGeometry, int iT ) {

  OstreamManager clout( std::cout,"setBoundaryValues" );

  if ( iT==0 ) {

    AnalyticalConst3D<T,T> rhoF( T( 1 ) );
    AnalyticalConst3D<T,T> uF( T( 0 ), T( 0 ), T( 0 ) );

    auto bulkIndicator = superGeometry.getMaterialIndicator({1, 2, 3});
    lattice.iniEquilibrium( bulkIndicator, rhoF, uF );
    lattice.defineRhoU( bulkIndicator, rhoF, uF );

    clout << converter.getCharLatticeVelocity() << std::endl;
    AnalyticalConst3D<T,T> uTop( converter.getCharLatticeVelocity(), T( 0 ), T( 0 ) );
    lattice.defineU( superGeometry,3,uTop );

    // Make the lattice ready for simulation
    lattice.initialize();
  }
}

void getResults( SuperLattice3D<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T, DESCRIPTOR> const& converter, SuperGeometry3D<T>& superGeometry,
                 int iT, Timer<T>& timer, bool converged ) {

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "cavity3d" );

  const T logT     = ( T )1.;
  const T vtkSave  = ( T )1.;
  const T imSave   = ( T )5.;

  if ( iT==0 ) {
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Get statistics
  if ( (iT%converter.getLatticeTime( logT )==0 && iT>0) || converged ) {
    timer.update( iT );
    timer.printStep( 2 );
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

  // Writes the VTK
  if ( (iT%converter.getLatticeTime( vtkSave )==0 && iT>0) || converged ) {
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );

    vtmWriter.write( iT );
  }

  // Writes the JPEG files
  if ( (iT%converter.getLatticeTime( imSave )==0 && iT>0) || converged ) {
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
    // define vector which span the plane
    Vector<T,3> u( 1,0,0 );
    Vector<T,3> v( 0,1,0 );
    T tmp = T( converter.getCharPhysLength() / 2. );
    T origin[3] = {tmp,tmp,tmp};

    SuperEuklidNorm3D<T, DESCRIPTOR> normVel( velocity );
    BlockReduction3D2D<T> planeReduction( normVel, origin, u, v, 600 );

    // write a heatmap
    heatmap::plotParam<T> plotParam;
    plotParam.maxValue = 1.;
    plotParam.name = "velocity";
    heatmap::write(planeReduction, iT, plotParam);
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
  Vector<T,3> origin( T( 0 ) );
  Vector<T,3> extend( converter.getCharPhysLength() );
  IndicatorCuboid3D<T> cube( extend,origin );

  // Instantiation of a cuboid geometry with weights
  int noCuboids = singleton::mpi().getSize();
  if ( noCuboids < 7 ) {
    noCuboids = 7;
  }
  CuboidGeometry3D<T> cuboidGeometry( cube, converter.getConversionFactorLength(), noCuboids );

  // Instantiation of a load balancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a super geometry
  SuperGeometry3D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( converter, cube, superGeometry );


  // === 3rd Step: Prepare Lattice ===
  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

  ConstRhoBGKdynamics<T, DESCRIPTOR> bulkDynamics (
    converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T,DESCRIPTOR>() );

  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sBoundaryCondition( sLattice );
  createInterpBoundaryCondition3D<T,DESCRIPTOR,ConstRhoBGKdynamics<T,DESCRIPTOR> >( sBoundaryCondition );

  prepareLattice( converter, sLattice, bulkDynamics, sBoundaryCondition, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  util::ValueTracer<T> converge( converter.getLatticeTime(interval), epsilon );

  Timer<T> timer( converter.getLatticeTime( maxT ), std::pow<int>(converter.getResolution(),3) );
  timer.start();

  for ( int iT = 0; iT <= converter.getLatticeTime( maxT ); ++iT ) {

    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << endl;
      getResults( sLattice, converter, superGeometry, iT, timer, converge.hasConverged() );
      break;
    }

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, sLattice, superGeometry, iT );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, superGeometry, iT, timer, converge.hasConverged() );
    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
  }

  timer.stop();
  timer.printSummary();
}
