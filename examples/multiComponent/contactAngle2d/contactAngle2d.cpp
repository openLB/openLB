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

/* contactAngle2d.cpp
 * In this example a semi-circular droplet of fluid is initialised
 * within a different fluid at a solid boundary. The contact angle
 * is measured as the droplet comes to equilibrium. This is compared
 * with the analytical angle (100 degrees) predicted by the
 * parameters set for the boundary.
 *
 * This example demonstrates how to use the wetting solid boundaries
 * for the free-energy model with two fluid components.
 */

#include "olb2D.h"
#include "olb2D.hh"   // use only generic version!
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9<CHEM_POTENTIAL,FORCE>

// Parameters for the simulation setup
const int N = 75;
const T nx  = 75.;
const T ny  = 50.;
const T radius = 0.25 * nx;

const T alpha = 1.;      // Interfacial width         [lattice units]
const T kappa1 = 0.005;  // For surface tensions      [lattice units]
const T kappa2 = 0.005;  // For surface tensions      [lattice units]
const T gama = 10.;      // For mobility of interface [lattice units]
const T h1 =  0.0001448; // Contact angle 100 degrees [lattice units]
const T h2 = -0.0001448; // Contact angle 100 degrees [lattice units]

const int maxIter = 70000;
const int vtkIter  = 1000;
const int statIter = 1000;
const bool calcAngle = true;

T angle_prev = 90.;


void prepareGeometry( SuperGeometry2D<T>& superGeometry, 
                      UnitConverter<T, DESCRIPTOR>& converter) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2 );

  Vector<T,2> extend(nx+2., ny-1.*converter.getPhysDeltaX() );
  Vector<T,2> origin( -1., 0.5*converter.getPhysDeltaX() );
  IndicatorCuboid2D<T> inner ( extend, origin );
  superGeometry.rename( 2,1,inner );

  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice( SuperLattice2D<T, DESCRIPTOR>& sLattice1,
                     SuperLattice2D<T, DESCRIPTOR>& sLattice2,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics1,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics2,
                     UnitConverter<T, DESCRIPTOR>& converter,
                     SuperGeometry2D<T>& superGeometry,
                     sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& sOnBC1,
                     sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& sOnBC2 ) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;
 
  // Define lattice Dynamics
  sLattice1.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice1.defineDynamics( superGeometry, 1, &bulkDynamics1 );
  sLattice2.defineDynamics( superGeometry, 1, &bulkDynamics2 );
  sLattice1.defineDynamics( superGeometry, 2, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 2, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Add wall boundary
  sOnBC1.addFreeEnergyWallBoundary( superGeometry, 2, alpha, kappa1, kappa2, h1, h2, 1 );
  sOnBC2.addFreeEnergyWallBoundary( superGeometry, 2, alpha, kappa1, kappa2, h1, h2, 2 );

  // Bulk initial conditions
  // Define circular domain for fluid 2
  std::vector<T> v( 2,T() );
  AnalyticalConst2D<T,T> zeroVelocity( v );

  AnalyticalConst2D<T,T> one( 1.0 );
  SmoothIndicatorCircle2D<T,T> circle( {nx/2., 0.}, radius, 10.*alpha );

  AnalyticalIdentity2D<T,T> rho( one );
  AnalyticalIdentity2D<T,T> phi( one - circle - circle );
  
  sLattice1.iniEquilibrium( superGeometry, 1, rho, zeroVelocity );
  sLattice2.iniEquilibrium( superGeometry, 1, phi, zeroVelocity );

  sLattice1.iniEquilibrium( superGeometry, 2, rho, zeroVelocity );
  sLattice2.iniEquilibrium( superGeometry, 2, phi, zeroVelocity );

  sLattice1.initialize();
  sLattice2.initialize();

  sLattice1.communicate();
  sLattice2.communicate();

  clout << "Prepare Lattice ... OK" << std::endl;
}


void prepareCoupling( SuperLattice2D<T, DESCRIPTOR>& sLattice1,
                      SuperLattice2D<T, DESCRIPTOR>& sLattice2,
                      SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareCoupling" );
  clout << "Add lattice coupling" << endl;

  // Add the lattice couplings (not to the solid nodes)
  // The chemical potential coupling must come before the force coupling
  FreeEnergyChemicalPotentialGenerator2D<T, DESCRIPTOR> coupling1(
    alpha, kappa1, kappa2);
  FreeEnergyForceGenerator2D<T, DESCRIPTOR> coupling2;

  sLattice1.addLatticeCoupling( superGeometry, 1, coupling1, sLattice2 );
  sLattice2.addLatticeCoupling( superGeometry, 1, coupling2, sLattice1 );

  clout << "Add lattice coupling ... OK!" << endl;
}


void getResults( SuperLattice2D<T, DESCRIPTOR>& sLattice1,
                 SuperLattice2D<T, DESCRIPTOR>& sLattice2, int iT,
                 SuperGeometry2D<T>& superGeometry, Timer<T>& timer,
                 UnitConverter<T, DESCRIPTOR> converter ) {

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter2D<T> vtmWriter( "contactAngle2d" );

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
  }

  // Writes the VTK files
  if ( iT%vtkIter==0 ) {
    AnalyticalConst2D<T,T> half_( 0.5 );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> half(half_, sLattice1);
    
    SuperLatticeVelocity2D<T, DESCRIPTOR> velocity( sLattice1 );
    SuperLatticeDensity2D<T, DESCRIPTOR> rho( sLattice1 );
    rho.getName() = "rho";
    SuperLatticeDensity2D<T, DESCRIPTOR> phi( sLattice2 );
    phi.getName() = "phi";

    SuperIdentity2D<T,T> c1 (half*(rho+phi));
    c1.getName() = "density-fluid-1";
    SuperIdentity2D<T,T> c2 (half*(rho-phi));
    c2.getName() = "density-fluid-2";

    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( rho );
    vtmWriter.addFunctor( phi );
    vtmWriter.addFunctor( c1 );
    vtmWriter.addFunctor( c2 );
    vtmWriter.write( iT );

    // Evaluation of contact angle
    if (calcAngle) {
      int Ny = (int)( N * ny / nx );
      T dx = converter.getPhysDeltaX();
      AnalyticalFfromSuperF2D<T,T> interpolPhi( phi, true, 1 );

      T base1 = 0.;
      T base2 = 0.;
      T height1 = 0.;
      T height2 = 0.;

      double pos[2] = {0., dx};
      for(int ix=0; ix<N; ix++) {
        T phi1, phi2;
        pos[0] = ix * dx;
        interpolPhi( &phi1, pos );
        if (phi1 < 0.) {
          pos[0] = (ix-1) * dx;
          interpolPhi( &phi2, pos );
          base1 = 2. * ( 0.5*N - ix + phi1/(phi1-phi2) );
          break;
        }
      }

      pos[1] = 3.*dx;
      for(int ix=0; ix<N; ix++) {
        T phi1, phi2;
        pos[0] = ix * dx;
        interpolPhi( &phi1, pos );
        if (phi1 < 0.) {
          pos[0] = (ix-1) * dx;
          interpolPhi( &phi2, pos );
          base2 = 2. * ( 0.5*N - ix + phi1/(phi1-phi2) );
          break;
        }
      }

      pos[0] = nx / 2.;
      for(int iy=2; iy<Ny; iy++) {
        T phi1, phi2;
        pos[1] = iy * dx;
        interpolPhi( &phi1, pos );
        if (phi1 > 0.) {
          pos[1] = (iy-1) * dx;
          interpolPhi( &phi2, pos );
          height1 = iy - 1. - phi1/(phi1-phi2);
          height2 = iy - 3. - phi1/(phi1-phi2);
          break;
        }
      }

      // Calculate simulated contact angle
      T pi = 3.14159265;
      T height = height1 + 1.;
      T base = base1 + 2 * (radius - height1) / base1;
      T radius = (4.*height2*height2 + base2*base2) / ( 8.*height2 );
      T angle_rad = pi + atan( 0.5*base / (radius - height) );
      T angle = angle_rad * 180. / pi;
      if ( angle > 180. ) angle -= 180.;

      // Calculate theoretical contact angle
      T ak1 = alpha * kappa1;
      T ak2 = alpha * kappa2;
      T k12 = kappa1 + kappa2;
      T num1 = pow(ak1 + 4 * h1, 1.5) - pow(ak1 - 4 * h1, 1.5);
      T num2 = pow(ak2 + 4 * h2, 1.5) - pow(ak2 - 4 * h2, 1.5);
      T angle_an = 180 / pi * acos(num2 / (2 * k12 * sqrt(ak2)) - \
                                   num1 / (2 * k12 * sqrt(ak1)));

      clout << "----->>>>> Contact angle: " << angle << " ; ";
      clout << "Analytical contact angle: " << angle_an <<  std::endl;
      clout << "----->>>>> Difference to previous: " << angle-angle_prev << std::endl;
      angle_prev = angle;   
    }
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
    (T)   nx, // charPhysLength: reference length of simulation geometry
    (T)   0.0001, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1.002e-8, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1. // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();

  // === 2nd Step: Prepare Geometry ===
  std::vector<T> extend = { nx, ny };
  std::vector<T> origin = { 0., 0. };
  IndicatorCuboid2D<T> cuboid(extend,origin);
#ifdef PARALLEL_MODE_MPI
  CuboidGeometry2D<T> cGeometry( cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize() );
#else
  CuboidGeometry2D<T> cGeometry( cuboid, converter.getPhysDeltaX() );
#endif

  // Set periodic boundaries to the domain
  cGeometry.setPeriodicity( true, false );

  // Instantiation of loadbalancer
  HeuristicLoadBalancer<T> loadBalancer( cGeometry );
  loadBalancer.print();

  // Instantiation of superGeometry
  SuperGeometry2D<T> superGeometry( cGeometry,loadBalancer );

  prepareGeometry( superGeometry, converter );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice2D<T, DESCRIPTOR> sLattice1( superGeometry );
  SuperLattice2D<T, DESCRIPTOR> sLattice2( superGeometry );

  ForcedBGKdynamics<T, DESCRIPTOR> bulkDynamics1 (
    converter.getLatticeRelaxationFrequency(),
    instances::getBulkMomenta<T,DESCRIPTOR>() );

  FreeEnergyBGKdynamics<T, DESCRIPTOR> bulkDynamics2 (
    converter.getLatticeRelaxationFrequency(), gama,
    instances::getBulkMomenta<T,DESCRIPTOR>() );

  sOnLatticeBoundaryCondition2D<T,DESCRIPTOR> sOnBC1( sLattice1 );
  sOnLatticeBoundaryCondition2D<T,DESCRIPTOR> sOnBC2( sLattice2 );
  createLocalBoundaryCondition2D<T,DESCRIPTOR> (sOnBC1);
  createLocalBoundaryCondition2D<T,DESCRIPTOR> (sOnBC2);

  prepareLattice( sLattice1, sLattice2, bulkDynamics1, bulkDynamics2,
                  converter, superGeometry, sOnBC1, sOnBC2 );

  prepareCoupling( sLattice1, sLattice2, superGeometry);

  SuperExternal2D<T,DESCRIPTOR,CHEM_POTENTIAL> sExternal1 (superGeometry, sLattice1, sLattice1.getOverlap() );
  SuperExternal2D<T,DESCRIPTOR,CHEM_POTENTIAL> sExternal2 (superGeometry, sLattice2, sLattice2.getOverlap() );

  // === 4th Step: Main Loop with Timer ===
  int iT = 0;
  clout << "starting simulation..." << endl;
  Timer<T> timer( maxIter, superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( iT=0; iT<=maxIter; ++iT ) {
    // Computation and output of the results
    getResults( sLattice1, sLattice2, iT, superGeometry, timer, converter );

    // Collide and stream execution
    sLattice1.collideAndStream();
    sLattice2.collideAndStream();

    // MPI communication for lattice data
    sLattice1.communicate();
    sLattice2.communicate();

    // Execute coupling between the two lattices
    sLattice1.executeCoupling();
    sExternal1.communicate();
    sExternal2.communicate();
    sLattice2.executeCoupling();
  }

  timer.stop();
  timer.printSummary();

}

