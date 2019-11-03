/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2014 Mathias J. Krause, Thomas Henn,
 *  Cyril Masquelier
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

/* venturi3d.cpp:
 * This example examines a steady flow in a venturi tube. At the
 * main inlet, a Poiseuille profile is imposed as Dirichlet velocity
 * boundary condition, whereas at the outlet and the minor inlet
 * a Dirichlet pressure condition is set by p=0 (i.e. rho=1).
 *
 * The example shows the usage of the Indicator functors to
 * build up a geometry and explains how to set boundary conditions
 * automatically.
 */

#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used
#include "olb3D.hh"     // Include full template code
#endif

#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
#define DESCRIPTOR D3Q19<>

T maxPhysT = 200.0; // max. simulation time in s, SI unit

SuperGeometry3D<T> prepareGeometry( ) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  std::string fName("venturi3d.xml");
  XMLreader config(fName);

  std::shared_ptr<IndicatorF3D<T> > inflow = createIndicatorCylinder3D<T>(config["Geometry"]["Inflow"]["IndicatorCylinder3D"], false);
  std::shared_ptr<IndicatorF3D<T> > outflow0 = createIndicatorCylinder3D<T>(config["Geometry"]["Outflow0"]["IndicatorCylinder3D"], false);
  std::shared_ptr<IndicatorF3D<T> > outflow1 = createIndicatorCylinder3D<T>(config["Geometry"]["Outflow1"]["IndicatorCylinder3D"], false);

  std::shared_ptr<IndicatorF3D<T> > venturi = createIndicatorF3D<T>(config["Geometry"]["Venturi"], false);

  // Build CoboidGeometry from IndicatorF (weights are set, remove and shrink is done)
  int N = config["Application"]["Discretization"]["Resolution"].get<int>();
  CuboidGeometry3D<T>* cuboidGeometry = new CuboidGeometry3D<T>( *venturi, 1./N, 20*singleton::mpi().getSize() );

  // Build LoadBalancer from CuboidGeometry (weights are respected)
  HeuristicLoadBalancer<T>* loadBalancer = new HeuristicLoadBalancer<T>( *cuboidGeometry );

  // Default instantiation of superGeometry
  SuperGeometry3D<T> superGeometry( *cuboidGeometry, *loadBalancer, 2 );

  // Set boundary voxels by rename material numbers
  superGeometry.rename( 0,2, venturi );
  superGeometry.rename( 2,1,1,1,1 );
  superGeometry.rename( 2,3,1, inflow );
  superGeometry.rename( 2,4,1, outflow0 );
  superGeometry.rename( 2,5,1, outflow1 );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();


  superGeometry.getStatistics().print();
  superGeometry.communicate();

  clout << "Prepare Geometry ... OK" << std::endl;
  return superGeometry;
}


void prepareLattice( SuperLattice3D<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T, DESCRIPTOR> const& converter,
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
  sLattice.defineDynamics( superGeometry, 1, &bulkDynamics );

  // Material=2 -->bounce back
  sLattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>() );

  // Material=3 -->bulk dynamics (inflow)
  sLattice.defineDynamics( superGeometry, 3, &bulkDynamics );

  // Material=4 -->bulk dynamics (outflow)
  sLattice.defineDynamics( superGeometry, 4, &bulkDynamics );

  // Material=5 -->bulk dynamics (2nd outflow)
  sLattice.defineDynamics( superGeometry, 5, &bulkDynamics );

  // Setting of the boundary conditions
  bc.addVelocityBoundary( superGeometry, 3, omega );
  bc.addPressureBoundary( superGeometry, 4, omega );
  bc.addPressureBoundary( superGeometry, 5, omega );

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing sinuidal inflow for the first iTMax timesteps
void setBoundaryValues( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                        SuperGeometry3D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime( maxPhysT*0.8 );
  int iTperiod = 50;

  if ( iT==0 ) {
    // Make the lattice ready for simulation
    sLattice.initialize();
  }

  else if ( iT%iTperiod==0 && iT<= iTmaxStart ) {
    //clout << "Set Boundary Values ..." << std::endl;

    //SinusStartScale<T,int> startScale(iTmaxStart, (T) 1);
    PolynomialStartScale<T,int> startScale( iTmaxStart, T( 1 ) );
    int iTvec[1]= {iT};
    T frac = T();
    startScale( &frac,iTvec );

    // Creates and sets the Poiseuille inflow profile using functors
    CirclePoiseuille3D<T> poiseuilleU( superGeometry, 3, frac*converter.getCharLatticeVelocity(), converter.getConversionFactorLength() );
    sLattice.defineU( superGeometry, 3, poiseuilleU );

    //clout << "step=" << iT << "; scalingFactor=" << frac << std::endl;
  }
  //clout << "Set Boundary Values ... ok" << std::endl;
}

void getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T, DESCRIPTOR>& converter, int iT,
                 SuperGeometry3D<T>& superGeometry, Timer<T>& timer ) {

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "venturi3d" );

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

  // Writes the vtm files
  if ( iT%converter.getLatticeTime( 1. )==0 ) {
    // Create the data-reading functors...
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.write( iT );

    SuperEuklidNorm3D<T, DESCRIPTOR> normVel( velocity );
    BlockReduction3D2D<T> planeReduction( normVel, {0, 0, 1} );

    // write output as JPEG
    heatmap::write(planeReduction, iT);

    // write output as JPEG and changing properties
    heatmap::plotParam<T> jpeg_Param;
    jpeg_Param.name             = "outflow";
    jpeg_Param.contourlevel     = 5;
    jpeg_Param.colour           = "blackbody";
    jpeg_Param.zoomOrigin       = {0.6, 0.3};
    jpeg_Param.zoomExtend       = {0.4, 0.7};
    heatmap::write(planeReduction, iT, jpeg_Param);
  }

  // Writes output on the console
  if ( iT%converter.getLatticeTime( 1. )==0 ) {
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

  std::string fName("venturi3d.xml");
  XMLreader config(fName);

  UnitConverter<T, DESCRIPTOR>* converter = createUnitConverter<T, DESCRIPTOR>(config);

  // Prints the converter log as console output
  converter->print();
  // Writes the converter log in a file
  converter->write("venturi3d");

  // === 2nd Step: Prepare Geometry ===

  SuperGeometry3D<T> superGeometry( prepareGeometry() );

  // === 3rd Step: Prepare Lattice ===

  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

  RLBdynamics<T, DESCRIPTOR> bulkDynamics( converter->getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() );

  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sBoundaryCondition( sLattice );
  createInterpBoundaryCondition3D<T, DESCRIPTOR> ( sBoundaryCondition );

  sOffLatticeBoundaryCondition3D<T, DESCRIPTOR> sOffBoundaryCondition( sLattice );
  createBouzidiBoundaryCondition3D<T, DESCRIPTOR> ( sOffBoundaryCondition );

  prepareLattice( sLattice, *converter, bulkDynamics, sBoundaryCondition, sOffBoundaryCondition, superGeometry );

  Timer<T> timer( converter->getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();
  getResults( sLattice, *converter, 0, superGeometry, timer );

  // === 4th Step: Main Loop with Timer ===
  for ( int iT = 0; iT <= converter->getLatticeTime( maxPhysT ); ++iT ) {

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLattice, *converter, iT, superGeometry );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, *converter, iT, superGeometry, timer );
  }

  timer.stop();
  timer.printSummary();
}
