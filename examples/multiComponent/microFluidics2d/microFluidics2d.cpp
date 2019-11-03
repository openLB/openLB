/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2019 Sam Avis
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

/*  microfluidics2d.cpp:
 *  This example shows a microfluidic channel creating droplets of
 *  two fluid components. Poiseuille velocity profiles are imposed
 *  at the various channel inlets, while a constant density outlet
 *  is imposed at the end of the channel to allow the droplets to
 *  exit the simulation.
 *
 *  This example demonstrates the use of three fluid components
 *  with the free energy model. It also shows the use of open
 *  boundary conditions, specifically velocity inlet and density
 *  outlet boundaries.
 */

#include "olb2D.h"
#include "olb2D.hh"   // use only generic version!
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9<CHEM_POTENTIAL,FORCE>

// Parameters for the simulation setup
const int N  = 100;
const T nx   = 800.;
const T ny   = 100.;
const T dx = ny / N;

const T inSize = 175.;
const T xl1 = inSize * 2./7.;
const T yl1 = ny / 4.;
const T xl2 = inSize / 7.;
const T yl2 = ny;
const T xl3 = inSize * 3./7.;
const T yl3 = ny / 4.;
const T xl4 = inSize / 7.;
const T yl4 = ny;
const T xl5 = nx - inSize;
const T yl5 = ny / 2.;

const T inlet1Velocity = 0.00056; // [lattice units]
const T inlet2Velocity = 0.00055; // [lattice units]
const T inlet3Velocity = 0.0014;  // [lattice units]
const T outletDensity = 1.;       // [lattice units]
const T alpha = 1.;               // Interfacial width          [lattice units]
const T kappa1 = 0.0132;          // For surface tensions       [lattice units]
const T kappa2 = 0.0012;          // For surface tensions       [lattice units]
const T kappa3 = 0.0013;          // For surface tensions       [lattice units]
const T gama = 1.;                // For mobility of interfaces [lattice units]
const T h1 = 0.;                  // Contact angle 90 degrees   [lattice units]
const T h2 = 0.;                  // Contact angle 90 degrees   [lattice units]
const T h3 = 0.;                  // Contact angle 90 degrees   [lattice units]

const int maxIter  = 1000000;
const int vtkIter  = 1000;
const int statIter = 2000;


void prepareGeometry( SuperGeometry2D<T>& superGeometry ) {
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  std::shared_ptr<IndicatorF2D<T>> section1 = std::make_shared<IndicatorCuboid2D<T>>( xl1, yl1, std::vector<T>{xl1/2., ny/2.} );
  std::shared_ptr<IndicatorF2D<T>> section2 = std::make_shared<IndicatorCuboid2D<T>>( xl2, yl2, std::vector<T>{xl1 + xl2/2., ny/2.} );
  std::shared_ptr<IndicatorF2D<T>> section3 = std::make_shared<IndicatorCuboid2D<T>>( xl3, yl3, std::vector<T>{xl1 + xl2 + xl3/2., ny/2.} );
  std::shared_ptr<IndicatorF2D<T>> section4 = std::make_shared<IndicatorCuboid2D<T>>( xl4, yl4, std::vector<T>{xl1 + xl2 + xl3 + xl4/2., ny/2.} );
  std::shared_ptr<IndicatorF2D<T>> section5 = std::make_shared<IndicatorCuboid2D<T>>( xl5, yl5, std::vector<T>{xl1 + xl2 + xl3 + xl4 + xl5/2., ny/2.} );
  IndicatorIdentity2D<T> channel( section1 + section2 + section3 + section4 + section5 );

  superGeometry.rename( 0, 2, channel );
  superGeometry.rename( 2,1,1,1 );

  // Inlets and outlet
  IndicatorCuboid2D<T> inlet1 ( dx, yl1, {0., ny/2.} );
  IndicatorCuboid2D<T> inlet21( xl2 - dx, dx, {xl1 + xl2/2., 0.} );
  IndicatorCuboid2D<T> inlet22( xl2 - dx, dx, {xl1 + xl2/2., ny} );
  IndicatorCuboid2D<T> inlet31( xl4 - dx, dx, {xl1 + xl2 + xl3 + xl4/2., 0.} );
  IndicatorCuboid2D<T> inlet32( xl4 - dx, dx, {xl1 + xl2 + xl3 + xl4/2., ny} );
  IndicatorCuboid2D<T> outlet( dx, yl5, {nx, ny/2.} );
  superGeometry.rename( 2, 3, 1, inlet1 );
  superGeometry.rename( 2, 4, 1, inlet21 );
  superGeometry.rename( 2, 5, 1, inlet22 );
  superGeometry.rename( 2, 6, 1, inlet31 );
  superGeometry.rename( 2, 7, 1, inlet32 );
  superGeometry.rename( 2, 8, 1, outlet );

  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice( SuperLattice2D<T, DESCRIPTOR>& sLattice1,
                     SuperLattice2D<T, DESCRIPTOR>& sLattice2,
                     SuperLattice2D<T, DESCRIPTOR>& sLattice3,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics1,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics2,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics3,
                     sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& sOnBC1,
                     sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& sOnBC2,
                     sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& sOnBC3,
                     UnitConverter<T, DESCRIPTOR>& converter,
                     SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;
 
  // define lattice dynamics
  sLattice1.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice3.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  sLattice1.defineDynamics( superGeometry, 1, &bulkDynamics1 );
  sLattice2.defineDynamics( superGeometry, 1, &bulkDynamics2 );
  sLattice3.defineDynamics( superGeometry, 1, &bulkDynamics3 );

  sLattice1.defineDynamics( superGeometry, 2, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 2, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice3.defineDynamics( superGeometry, 2, &instances::getNoDynamics<T, DESCRIPTOR>() );

  sLattice1.defineDynamics( superGeometry, 3, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 3, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice3.defineDynamics( superGeometry, 3, &instances::getNoDynamics<T, DESCRIPTOR>() );

  sLattice1.defineDynamics( superGeometry, 4, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 4, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice3.defineDynamics( superGeometry, 4, &instances::getNoDynamics<T, DESCRIPTOR>() );

  sLattice1.defineDynamics( superGeometry, 5, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 5, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice3.defineDynamics( superGeometry, 5, &instances::getNoDynamics<T, DESCRIPTOR>() );

  sLattice1.defineDynamics( superGeometry, 6, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 6, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice3.defineDynamics( superGeometry, 6, &instances::getNoDynamics<T, DESCRIPTOR>() );

  sLattice1.defineDynamics( superGeometry, 7, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 7, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice3.defineDynamics( superGeometry, 7, &instances::getNoDynamics<T, DESCRIPTOR>() );

  sLattice1.defineDynamics( superGeometry, 8, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 8, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice3.defineDynamics( superGeometry, 8, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // add wall boundary
  sOnBC1.addFreeEnergyWallBoundary( superGeometry, 2, alpha, kappa1, kappa2, kappa3, h1, h2, h3, 1 );
  sOnBC2.addFreeEnergyWallBoundary( superGeometry, 2, alpha, kappa1, kappa2, kappa3, h1, h2, h3, 2 );
  sOnBC3.addFreeEnergyWallBoundary( superGeometry, 2, alpha, kappa1, kappa2, kappa3, h1, h2, h3, 3 );

  // add inlet boundaries
  T omega = converter.getLatticeRelaxationFrequency();
  auto inlet1Indicator = superGeometry.getMaterialIndicator(3);
  sOnBC1.addFreeEnergyInletBoundary( inlet1Indicator, omega, "velocity", 1 );
  sOnBC2.addFreeEnergyInletBoundary( inlet1Indicator, omega, "velocity", 2 );
  sOnBC3.addFreeEnergyInletBoundary( inlet1Indicator, omega, "velocity", 3 );

  auto inlet2Indicator = superGeometry.getMaterialIndicator({4, 5});
  sOnBC1.addFreeEnergyInletBoundary( inlet2Indicator, omega, "velocity", 1 );
  sOnBC2.addFreeEnergyInletBoundary( inlet2Indicator, omega, "velocity", 2 );
  sOnBC3.addFreeEnergyInletBoundary( inlet2Indicator, omega, "velocity", 3 );

  auto inlet3Indicator = superGeometry.getMaterialIndicator({6, 7});
  sOnBC1.addFreeEnergyInletBoundary( inlet3Indicator, omega, "velocity", 1 );
  sOnBC2.addFreeEnergyInletBoundary( inlet3Indicator, omega, "velocity", 2 );
  sOnBC3.addFreeEnergyInletBoundary( inlet3Indicator, omega, "velocity", 3 );

  // add outlet boundary
  auto outletIndicator = superGeometry.getMaterialIndicator(8);
  sOnBC1.addFreeEnergyOutletBoundary( outletIndicator, omega, "density", 1 );
  sOnBC2.addFreeEnergyOutletBoundary( outletIndicator, omega, "density", 2 );
  sOnBC3.addFreeEnergyOutletBoundary( outletIndicator, omega, "density", 3 );

  // bulk initial conditions
  std::vector<T> v( 2,T() );
  AnalyticalConst2D<T,T> zeroVelocity( v );

  AnalyticalConst2D<T,T> zero ( 0. );
  AnalyticalConst2D<T,T> one ( 1. );
  SmoothIndicatorCuboid2D<T,T> section1( {xl1/2., ny/2.}, xl1+dx, ny, 0. );
  SmoothIndicatorCuboid2D<T,T> section2( {xl1 + (xl2 + xl3)/2., ny/2.}, xl2 + xl3, ny, 0. );

  AnalyticalIdentity2D<T,T> c1( section1 );
  AnalyticalIdentity2D<T,T> c2( section2 );
  AnalyticalIdentity2D<T,T> rho( one );
  AnalyticalIdentity2D<T,T> phi( c1 - c2 );
  AnalyticalIdentity2D<T,T> psi( rho - c1 - c2 );

  auto allIndicator = superGeometry.getMaterialIndicator({1, 2, 3, 4, 5, 6});
  sLattice1.iniEquilibrium( allIndicator, rho, zeroVelocity );
  sLattice2.iniEquilibrium( allIndicator, phi, zeroVelocity );
  sLattice3.iniEquilibrium( allIndicator, psi, zeroVelocity );

  // inlet boundary conditions
  Poiseuille2D<T> inlet1U( superGeometry, 3, 1.5*inlet1Velocity, 0. );
  sLattice1.defineU( inlet1Indicator, inlet1U );
  sLattice2.defineRho( inlet1Indicator, phi );
  sLattice3.defineRho( inlet1Indicator, psi );

  Poiseuille2D<T> inlet21U( superGeometry, 4, 1.5*inlet2Velocity, 0. );
  Poiseuille2D<T> inlet22U( superGeometry, 5, 1.5*inlet2Velocity, 0. );
  sLattice1.defineU( superGeometry, 4, inlet21U );
  sLattice1.defineU( superGeometry, 5, inlet22U );
  sLattice2.defineRho( inlet2Indicator, phi );
  sLattice3.defineRho( inlet2Indicator, psi );

  Poiseuille2D<T> inlet31U( superGeometry, 6, 1.5*inlet3Velocity, 0. );
  Poiseuille2D<T> inlet32U( superGeometry, 7, 1.5*inlet3Velocity, 0. );
  sLattice1.defineU( superGeometry, 6, inlet31U );
  sLattice1.defineU( superGeometry, 7, inlet32U );
  sLattice2.defineRho( inlet3Indicator, phi );
  sLattice3.defineRho( inlet3Indicator, psi );

  // outlet initial / boundary conditions
  AnalyticalConst2D<T,T> rhoOutlet( outletDensity );
  AnalyticalIdentity2D<T,T> phiOutlet( zero );
  AnalyticalIdentity2D<T,T> psiOutlet( rhoOutlet );
  sLattice1.defineRho( outletIndicator, rhoOutlet );
  sLattice2.defineRho( outletIndicator, phiOutlet );
  sLattice3.defineRho( outletIndicator, psiOutlet );

  // initialise lattices
  sLattice1.initialize();
  sLattice2.initialize();
  sLattice3.initialize();

  sLattice1.communicate();
  sLattice2.communicate();
  sLattice3.communicate();

  clout << "Prepare Lattice ... OK" << std::endl;
}


void prepareCoupling(SuperLattice2D<T, DESCRIPTOR>& sLattice1,
                     SuperLattice2D<T, DESCRIPTOR>& sLattice2,
                     SuperLattice2D<T, DESCRIPTOR>& sLattice3,
                     SuperGeometry2D<T>& superGeometry) {
  OstreamManager clout( std::cout,"prepareCoupling" );
  clout << "Add lattice coupling" << endl;

  // Bulk couplings
  FreeEnergyChemicalPotentialGenerator2D<T,DESCRIPTOR> coupling2( alpha, kappa1, kappa2, kappa3 );
  FreeEnergyForceGenerator2D<T,DESCRIPTOR> coupling3;

  // Inlet / outlet couplings
  FreeEnergyDensityOutletGenerator2D<T,DESCRIPTOR> coupling1( outletDensity );
  FreeEnergyInletOutletGenerator2D<T,DESCRIPTOR> coupling4;

  // The DensityOutlet coupling must be added to the first lattice and come before the ChemicalPotential coupling
  // The InletOutlet couplings must be added to the second lattice and come after the Force coupling.
  sLattice1.addLatticeCoupling<DESCRIPTOR>( superGeometry, 8, coupling1, {&sLattice2, &sLattice3} );

  sLattice1.addLatticeCoupling<DESCRIPTOR>( superGeometry, 1, coupling2, {&sLattice2, &sLattice3} );
  sLattice2.addLatticeCoupling<DESCRIPTOR>( superGeometry, 1, coupling3, {&sLattice1, &sLattice3} );

  sLattice2.addLatticeCoupling<DESCRIPTOR>( superGeometry, 3, coupling4, {&sLattice1, &sLattice3} );
  sLattice2.addLatticeCoupling<DESCRIPTOR>( superGeometry, 4, coupling4, {&sLattice1, &sLattice3} );
  sLattice2.addLatticeCoupling<DESCRIPTOR>( superGeometry, 5, coupling4, {&sLattice1, &sLattice3} );
  sLattice2.addLatticeCoupling<DESCRIPTOR>( superGeometry, 6, coupling4, {&sLattice1, &sLattice3} );
  sLattice2.addLatticeCoupling<DESCRIPTOR>( superGeometry, 7, coupling4, {&sLattice1, &sLattice3} );
  sLattice2.addLatticeCoupling<DESCRIPTOR>( superGeometry, 8, coupling4, {&sLattice1, &sLattice3} );

  clout << "Add lattice coupling ... OK!" << endl;
}


void getResults( SuperLattice2D<T, DESCRIPTOR>& sLattice1,
                 SuperLattice2D<T, DESCRIPTOR>& sLattice2,
                 SuperLattice2D<T, DESCRIPTOR>& sLattice3, int iT,
                 SuperGeometry2D<T>& superGeometry, Timer<T>& timer,
                 UnitConverter<T, DESCRIPTOR> converter) {

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter2D<T> vtmWriter( "microFluidics2d" );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry( sLattice1, superGeometry );
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice1 );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice1 );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Get statistics
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();
    sLattice1.getStatistics().print( iT, converter.getPhysTime(iT) );
    sLattice2.getStatistics().print( iT, converter.getPhysTime(iT) );
    sLattice3.getStatistics().print( iT, converter.getPhysTime(iT) );
  }

  // Writes the VTK files
  if ( iT%vtkIter==0 ) {
    SuperLatticeVelocity2D<T, DESCRIPTOR> velocity( sLattice1 );
    SuperLatticeDensity2D<T, DESCRIPTOR> density1( sLattice1 );
    density1.getName() = "rho";
    SuperLatticeDensity2D<T, DESCRIPTOR> density2( sLattice2 );
    density2.getName() = "phi";
    SuperLatticeDensity2D<T, DESCRIPTOR> density3( sLattice3 );
    density3.getName() = "density-fluid-3";
    
    AnalyticalConst2D<T,T> half_( 0.5 );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> half(half_, sLattice1);

    SuperIdentity2D<T,T> c1 (half*(density1+density2-density3));
    c1.getName() = "density-fluid-1";
    SuperIdentity2D<T,T> c2 (half*(density1-density2-density3));
    c2.getName() = "density-fluid-2";

    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( density1 );
    vtmWriter.addFunctor( density2 );
    vtmWriter.addFunctor( density3 );
    vtmWriter.addFunctor( c1 );
    vtmWriter.addFunctor( c2 );
    vtmWriter.write( iT );
  }
}


int main( int argc, char *argv[] ) {

  // === 1st Step: Initialization ===

  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converter(
    (T)   N, // resolution
    (T)   1., // lattice relaxation time (tau)
    (T)   ny, // charPhysLength: reference length of simulation geometry
    (T)   1.e-6, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.1, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1. // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();

  // === 2nd Step: Prepare Geometry ===
  std::vector<T> extend = { nx, ny };
  std::vector<T> origin = { 0, 0 };
  IndicatorCuboid2D<T> cuboid(extend,origin);
#ifdef PARALLEL_MODE_MPI
  CuboidGeometry2D<T> cGeometry( cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize() );
#else
  CuboidGeometry2D<T> cGeometry( cuboid, converter.getPhysDeltaX() );
#endif

  // Instantiation of loadbalancer
  HeuristicLoadBalancer<T> loadBalancer( cGeometry );
  loadBalancer.print();

  // Instantiation of superGeometry
  SuperGeometry2D<T> superGeometry( cGeometry,loadBalancer );

  prepareGeometry( superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice2D<T, DESCRIPTOR> sLattice1( superGeometry );
  SuperLattice2D<T, DESCRIPTOR> sLattice2( superGeometry );
  SuperLattice2D<T, DESCRIPTOR> sLattice3( superGeometry );

  ForcedBGKdynamics<T, DESCRIPTOR> bulkDynamics1 (
    converter.getLatticeRelaxationFrequency(),
    instances::getBulkMomenta<T,DESCRIPTOR>() );

  FreeEnergyBGKdynamics<T, DESCRIPTOR> bulkDynamics23 (
    converter.getLatticeRelaxationFrequency(), gama,
    instances::getBulkMomenta<T,DESCRIPTOR>() );

  sOnLatticeBoundaryCondition2D<T, DESCRIPTOR> sOnBC1( sLattice1 );
  sOnLatticeBoundaryCondition2D<T, DESCRIPTOR> sOnBC2( sLattice2 );
  sOnLatticeBoundaryCondition2D<T, DESCRIPTOR> sOnBC3( sLattice3 );
  createLocalBoundaryCondition2D<T, DESCRIPTOR> (sOnBC1);
  createLocalBoundaryCondition2D<T, DESCRIPTOR> (sOnBC2);
  createLocalBoundaryCondition2D<T, DESCRIPTOR> (sOnBC3);

  prepareLattice( sLattice1, sLattice2, sLattice3, bulkDynamics1, bulkDynamics23,
                  bulkDynamics23, sOnBC1, sOnBC2, sOnBC3, converter, superGeometry );

  prepareCoupling( sLattice1, sLattice2, sLattice3, superGeometry );

  SuperExternal2D<T,DESCRIPTOR,CHEM_POTENTIAL> sExternal1 (superGeometry, sLattice1, sLattice1.getOverlap() );
  SuperExternal2D<T,DESCRIPTOR,CHEM_POTENTIAL> sExternal2 (superGeometry, sLattice2, sLattice2.getOverlap() );
  SuperExternal2D<T,DESCRIPTOR,CHEM_POTENTIAL> sExternal3 (superGeometry, sLattice3, sLattice3.getOverlap() );

  // === 4th Step: Main Loop with Timer ===
  int iT = 0;
  clout << "starting simulation..." << endl;
  Timer<T> timer( maxIter, superGeometry.getStatistics().getNvoxel() );
  timer.start();
  
  for ( iT=0; iT<maxIter; ++iT ) {
    // Computation and output of the results
    getResults( sLattice1, sLattice2, sLattice3, iT, superGeometry, timer, converter );

    // Collide and stream execution
    sLattice1.collideAndStream();
    sLattice2.collideAndStream();
    sLattice3.collideAndStream();

    // MPI communication for lattice data
    sLattice1.communicate();
    sLattice2.communicate();
    sLattice3.communicate();
 
    // Execute coupling between the two lattices
    sLattice1.executeCoupling();
    sExternal1.communicate();
    sExternal2.communicate();
    sExternal3.communicate();
    sLattice2.executeCoupling();
  }

  timer.stop();
  timer.printSummary();

}
