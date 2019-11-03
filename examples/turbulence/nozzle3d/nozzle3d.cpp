/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2015 Mathias J. Krause, Patrick Nathan
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

/* nozzle3d.cpp:
 * This example examines a turbulent flow in a nozzle injection tube. At the
 * main inlet, either a block profile or a power 1/7 profile is imposed as a
 * Dirchlet velocity boundary condition, whereas at the outlet a
 * Dirichlet pressure condition is set by p=0 (i.e. rho=1).
 *
 * The example shows the usage of turbulent models.
 *
 * The results of a simular simulation setup are publish in TODO
 */

#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used
#include "olb3D.hh"     // Include full template code
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;

// Choose your turbulent model of choice
//#define RLB
#define Smagorinsky //default
//#define ConsitentStrainSmagorinsky
//#define ShearSmagorinsky
//#define Krause

#ifdef ShearSmagorinsky
#define DESCRIPTOR ShearSmagorinskyD3Q19Descriptor
#else
#define DESCRIPTOR D3Q19<>
#endif

// Parameters for the simulation setup
const int N = 3;                 // resolution of the model, for RLB N>=5, others N>=2, but N>=5 recommended
const int M = 1;                 // time discretization refinement
const int inflowProfileMode = 0; // block profile (mode=0), power profile (mode=1)
const T maxPhysT = 200.;         // max. simulation time in s, SI unit

template <typename T, typename _DESCRIPTOR>
class TurbulentVelocity3D : public AnalyticalF3D<T,T> {

protected:
  // block profile (mode=0), power profile (mode=1)
  int _mode;
  T rho;
  T nu;
  T u0;
  T p0;
  T charL;
  T dx;

public:
  TurbulentVelocity3D( UnitConverter<T,_DESCRIPTOR> const& converter, int mode=0 ) : AnalyticalF3D<T,T>( 3 ) {
    _mode = mode;
    u0 = converter.getCharLatticeVelocity();
    rho = converter.getPhysDensity();
    nu = converter.getPhysViscosity();
    charL = converter.getCharPhysLength();
    p0 = converter.getCharPhysPressure();
    dx = converter.getConversionFactorLength();

    this->getName() = "turbulentVelocity3d";
  };

  bool operator()( T output[], const S input[] ) override {
    T y = input[1];
    T z = input[2];
    // block profile inititalization
    T u_calc = u0;
    // power profile inititalization
    if ( _mode==1 ) {
      T obst_y = 5.5+dx;
      T obst_z = 5.5+dx;
      T obst_r = 0.5;

      T B      = 5.5;
      T kappa  = 0.4;
      T ReTau  = 183.6;

      u_calc = u0/7.*( 2.*nu*ReTau/( charL*kappa )*log( fabs( 2.*ReTau/charL*( obst_r - sqrt( pow( y - obst_y, 2. )
                       + pow( z - obst_z, 2. ) ) )*1.5*( 1 + sqrt( pow( y - obst_y, 2. )
                           + pow( z - obst_z, 2. ) )/obst_r )/( 1 + 2.*pow( sqrt( pow( y - obst_y, 2. )
                               + pow( z - obst_z, 2. ) )/obst_r, 2. ) ) ) + B ) );
    }
    T a = -1., b = 1.;
    T nRandom = rand()/( T )RAND_MAX*( b-a ) + a;

    output[0] = u_calc+  0.15*u0*nRandom;
    output[1] = 0.15*u0*nRandom;
    output[2] = 0.15*u0*nRandom;
    return true;
  };
};


void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator, SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // Sets material number for fluid and boundary
  superGeometry.rename( 0,2,indicator );

  Vector<T,3> origin( T(),
                      5.5*converter.getCharPhysLength()+converter.getConversionFactorLength(),
                      5.5*converter.getCharPhysLength()+converter.getConversionFactorLength() );

  Vector<T,3> extend( 4.*converter.getCharPhysLength()+5*converter.getConversionFactorLength(),
                      5.5*converter.getCharPhysLength()+converter.getConversionFactorLength(),
                      5.5*converter.getCharPhysLength()+converter.getConversionFactorLength() );

  IndicatorCylinder3D<T> inletCylinder( extend, origin, converter.getCharPhysLength() );
  superGeometry.rename( 2,1,inletCylinder );


  origin[0]=4.*converter.getCharPhysLength();
  origin[1]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();
  origin[2]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();

  extend[0]=54.*converter.getCharPhysLength();
  extend[1]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();
  extend[2]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();

  IndicatorCylinder3D<T> injectionTube( extend, origin, 5.5*converter.getCharPhysLength() );
  superGeometry.rename( 2,1,injectionTube );

  origin[0]=converter.getConversionFactorLength();
  origin[1]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();
  origin[2]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();

  extend[0]=T();
  extend[1]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();
  extend[2]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();

  IndicatorCylinder3D<T> cylinderIN( extend, origin, converter.getCharPhysLength() );
  superGeometry.rename( 1,3,cylinderIN );


  origin[0]=54.*converter.getCharPhysLength()-converter.getConversionFactorLength();
  origin[1]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();
  origin[2]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();

  extend[0]=54.*converter.getCharPhysLength();
  extend[1]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();
  extend[2]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();

  IndicatorCylinder3D<T> cylinderOUT( extend, origin, 5.5*converter.getCharPhysLength() );
  superGeometry.rename( 1,4,cylinderOUT );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( SuperLattice3D<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& bc,
                     sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& offBc,
                     SuperGeometry3D<T>& superGeometry ) {



  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  // Material=3 -->bulk dynamics (inflow)
  // Material=4 -->bulk dynamics (outflow)
  sLattice.defineDynamics( superGeometry.getMaterialIndicator({1, 3, 4}), &bulkDynamics );

  // Material=2 -->bounce back
  sLattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>() );

  bc.addVelocityBoundary( superGeometry, 3, omega );
  bc.addPressureBoundary( superGeometry, 4, omega );

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues( UnitConverter<T,DESCRIPTOR> const&converter,
                        SuperLattice3D<T,DESCRIPTOR>& lattice, SuperGeometry3D<T>& superGeometry, int iT ) {

  OstreamManager clout( std::cout,"setBoundaryValues" );

  if ( iT==0 ) {
    AnalyticalConst3D<T,T> rhoF( 1 );
    std::vector<T> velocity( 3,T() );
    AnalyticalConst3D<T,T> uF( velocity );

    // Seeding of random fluctuations and definition of the velocity field
    srand( time( nullptr ) );
    TurbulentVelocity3D<T,DESCRIPTOR> uSol( converter, inflowProfileMode );

    lattice.iniEquilibrium( superGeometry.getMaterialIndicator({1, 2, 4}), rhoF, uF );
    lattice.iniEquilibrium( superGeometry, 3, rhoF, uSol );

    lattice.defineU( superGeometry, 3, uSol );
    lattice.defineRho( superGeometry, 4, rhoF );

    // Make the lattice ready for simulation
    lattice.initialize();
  }
}

void getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry3D<T>& superGeometry, Timer<T>& timer ) {

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "nozzle3d" );

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
  if ( iT%converter.getLatticeTime( maxPhysT/100. )==0 ) {
    // Create the data-reading functors...
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.write( iT );

    SuperEuklidNorm3D<T, DESCRIPTOR> normVel( velocity );
    BlockReduction3D2D<T> planeReduction( normVel, {0, 1, 0} );
    // write output as JPEG
    heatmap::write(planeReduction, iT);
  }

  // Writes output on the console
  if ( iT%converter.getLatticeTime( maxPhysT/200. )==0 ) {
    timer.update( iT );
    timer.printStep();
    sLattice.getStatistics().print( iT, converter.getPhysTime( iT ) );
  }
}



int main( int argc, char* argv[] ) {

  // === 1st Step: Initialization ===

  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  // clout.setMultiOutput(true);

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},        // resolution: number of voxels per charPhysL
    (T)   0.500018, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   1,        // charPhysLength: reference length of simulation geometry
    (T)   1,        // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.0002,   // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0       // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("nozzle3d");

  Vector<T,3> origin;
  Vector<T,3> extend( 54.*converter.getCharPhysLength(), 11.*converter.getCharPhysLength()+2.*converter.getConversionFactorLength(), 11.*converter.getCharPhysLength()+2.*converter.getConversionFactorLength() );

  IndicatorCuboid3D<T> cuboid( extend,origin );

  CuboidGeometry3D<T> cuboidGeometry( cuboid, converter.getConversionFactorLength(), singleton::mpi().getSize() );
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // === 2nd Step: Prepare Geometry ===

  SuperGeometry3D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );
  prepareGeometry( converter, cuboid, superGeometry );

  // === 3rd Step: Prepare Lattice ===

  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

  Dynamics<T, DESCRIPTOR>* bulkDynamics;
  const T omega = converter.getLatticeRelaxationFrequency();
#if defined(RLB)
  bulkDynamics = new RLBdynamics<T, DESCRIPTOR>( omega, instances::getBulkMomenta<T, DESCRIPTOR>() );
#elif defined(Smagorinsky)
  bulkDynamics = new SmagorinskyBGKdynamics<T, DESCRIPTOR>( omega, instances::getBulkMomenta<T, DESCRIPTOR>(),
      0.15);
#elif defined(ShearSmagorinsky)
  bulkDynamics = new ShearSmagorinskyBGKdynamics<T, DESCRIPTOR>( omega, instances::getBulkMomenta<T, DESCRIPTOR>(),
      0.15);
#elif defined(Krause)
  bulkDynamics = new KrauseBGKdynamics<T, DESCRIPTOR>( omega, instances::getBulkMomenta<T, DESCRIPTOR>(),
      0.15);
#else //ConsitentStrainSmagorinsky
  bulkDynamics = new ConStrainSmagorinskyBGKdynamics<T, DESCRIPTOR>( omega, instances::getBulkMomenta<T, DESCRIPTOR>(),
      0.05);
#endif

  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sBoundaryCondition( sLattice );
  createInterpBoundaryCondition3D<T, DESCRIPTOR> ( sBoundaryCondition );

  sOffLatticeBoundaryCondition3D<T, DESCRIPTOR> sOffBoundaryCondition( sLattice );
  createBouzidiBoundaryCondition3D<T, DESCRIPTOR> ( sOffBoundaryCondition );

  prepareLattice( sLattice, converter, *bulkDynamics, sBoundaryCondition, sOffBoundaryCondition, superGeometry );

  // === 4th Step: Main Loop with Timer ===

  Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( int iT = 0; iT <= converter.getLatticeTime( maxPhysT ); ++iT ) {
    // === 5ath Step: Apply filter
#ifdef ADM
    SuperLatticeADM3D<T, DESCRIPTOR> admF( sLattice, 0.01, 2 );
    admF.execute( superGeometry, 1 );
#endif
    // === 5bth Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, sLattice, superGeometry, iT );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer );
  }
  timer.stop();
  timer.printSummary();
  delete bulkDynamics;
}
