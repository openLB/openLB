/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2007, 2012 Jonas Latt, Mathias J. Krause
 *  Vojtech Cvrcek, Peter Weisbrod
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

/* poiseuille2d.cpp:
 * This example examines a 2D Poseuille flow
 * It illustrates the computation of error norms.
 */


#include "olb2D.h"
#include "olb2D.hh"   // use only generic version!
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

//#define MRT
#ifdef MRT
#define DESCRIPTOR ForcedMRTD2Q9Descriptor
#else
#define DESCRIPTOR  D2Q9<FORCE>
#endif

typedef enum {forced, nonForced} FlowType;

typedef enum {bounceBack, local, interpolated, freeSlip, partialSlip} BoundaryType;


// Parameters for the simulation setup
FlowType flowType = forced;
BoundaryType boundaryType = interpolated;
const T lx  = 2.;             // length of the channel
const T ly  = 1.;             // height of the channel
int N = 50;                   // resolution of the model
const T Re = 10.;             // Reynolds number
const T maxPhysT = 20.;       // max. simulation time in s, SI unit
const T physInterval = 0.25;  // interval for the convergence check in s
const T residuum = 1e-5;      // residuum for the convergence check
const T tuner = 0.99;         // for partialSlip only: 0->bounceBack, 1->freeSlip


// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2 );

  superGeometry.rename( 2,1,1,1 );

  if (flowType == nonForced) {
    Vector<T,2> extend;
    Vector<T,2> origin;
    T physSpacing = converter.getPhysDeltaX();

    // Set material number for inflow
    extend[1] = ly;
    extend[0] = physSpacing / 2;
    origin[0] -= physSpacing / 4;
    IndicatorCuboid2D<T> inflow( extend, origin );
    superGeometry.rename( 2,3,1,inflow );

    // Set material number for outflow
    origin[0] = lx - physSpacing / 4;
    IndicatorCuboid2D<T> outflow( extend, origin );
    superGeometry.rename( 2,4,1,outflow );
  }

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
                     SuperLattice2D<T, DESCRIPTOR>& sLattice,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& sBoundaryCondition,
                     SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 1, &bulkDynamics );

  if (boundaryType == bounceBack) {
    sLattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>() );
  } else if (boundaryType == freeSlip) {
    sLattice.defineDynamics(superGeometry, 2, &instances::getNoDynamics<T, DESCRIPTOR>());
    sBoundaryCondition.addSlipBoundary( superGeometry, 2 );
  } else if (boundaryType == partialSlip) {
    sLattice.defineDynamics(superGeometry, 2, &instances::getNoDynamics<T, DESCRIPTOR>());
    sBoundaryCondition.addPartialSlipBoundary(tuner, superGeometry, 2 );
  } else {
    sLattice.defineDynamics( superGeometry, 2, &bulkDynamics );
    sBoundaryCondition.addVelocityBoundary( superGeometry, 2, omega );
  }

  if (flowType == nonForced) {
    // Material=3 -->bulk dynamics
    sLattice.defineDynamics( superGeometry, 3, &bulkDynamics );

    // Material=4 -->bulk dynamics
    sLattice.defineDynamics( superGeometry, 4, &bulkDynamics );

    sBoundaryCondition.addVelocityBoundary( superGeometry, 3, omega );
    sBoundaryCondition.addPressureBoundary( superGeometry, 4, omega );
  }

  // Initial conditions
  T Lx = converter.getLatticeLength( lx );
  T Ly = converter.getLatticeLength( ly );

  if (flowType == forced) {
    std::vector<T> poiseuilleForce( 2,T() );
    poiseuilleForce[0] = 8.*converter.getLatticeViscosity()*converter.getCharLatticeVelocity() / ( Ly*Ly );
    AnalyticalConst2D<T,T> force( poiseuilleForce );

    // Initialize force
    sLattice.defineField<FORCE>(superGeometry, 1, force);
    sLattice.defineField<FORCE>(superGeometry, 2, force);
  } else {
    T p0 =8.*converter.getLatticeViscosity()*converter.getCharLatticeVelocity()*Lx/( Ly*Ly );
    AnalyticalLinear2D<T,T> rho( -p0/lx*invCs2<T,DESCRIPTOR>(), 0 , p0*invCs2<T,DESCRIPTOR>()+1 );

    T maxVelocity = converter.getCharLatticeVelocity();
    T distance2Wall = converter.getConversionFactorLength();
    Poiseuille2D<T> u( superGeometry, 3, maxVelocity, distance2Wall );

    // Initialize all values of distribution functions to their local equilibrium
    sLattice.defineRhoU( superGeometry, 1, rho, u );
    sLattice.iniEquilibrium( superGeometry, 1, rho, u );
    sLattice.defineRhoU( superGeometry, 2, rho, u );
    sLattice.iniEquilibrium( superGeometry, 2, rho, u );
    sLattice.defineRhoU( superGeometry, 3, rho, u );
    sLattice.iniEquilibrium( superGeometry, 3, rho, u );
    sLattice.defineRhoU( superGeometry, 4, rho, u );
    sLattice.iniEquilibrium( superGeometry, 4, rho, u );
  }

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Compute error norms
void error( SuperGeometry2D<T>& superGeometry,
            SuperLattice2D<T, DESCRIPTOR>& sLattice,
            UnitConverter<T,DESCRIPTOR> const& converter,
            Dynamics<T, DESCRIPTOR>& bulkDynamics ) {

  OstreamManager clout( std::cout,"error" );

  int tmp[]= { };
  T result[2]= { };

  // velocity error
  const T maxVelocity = converter.getCharPhysVelocity();
  const T radius = ly/2.;
  std::vector<T> axisPoint( 2,T() );
  axisPoint[0] = lx/2.;
  axisPoint[1] = ly/2.;
  std::vector<T> axisDirection( 2,T() );
  axisDirection[0] = 1;
  axisDirection[1] = 0;
  Poiseuille2D<T> uSol( axisPoint, axisDirection, maxVelocity, radius );
  SuperLatticePhysVelocity2D<T,DESCRIPTOR> u( sLattice,converter );
  auto indicatorF = superGeometry.getMaterialIndicator(1);

  SuperAbsoluteErrorL1Norm2D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm2D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm2D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm2D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;

  // strainRate error
  PoiseuilleStrainRate2D<T,T,DESCRIPTOR> sSol( converter, ly );
  SuperLatticePhysStrainRate2D<T,DESCRIPTOR> s( sLattice,converter );

  SuperAbsoluteErrorL1Norm2D<T> absStrainRateErrorNormL1(s, sSol, indicatorF);
  absStrainRateErrorNormL1(result, tmp);
  clout << "strainRate-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm2D<T> relStrainRateErrorNormL1(s, sSol, indicatorF);
  relStrainRateErrorNormL1(result, tmp);
  clout << "; strainRate-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm2D<T> absStrainRateErrorNormL2(s, sSol, indicatorF);
  absStrainRateErrorNormL2(result, tmp);
  clout << "strainRate-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relStrainRateErrorNormL2(s, sSol, indicatorF);
  relStrainRateErrorNormL2(result, tmp);
  clout << "; strainRate-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absStrainRateErrorNormLinf(s, sSol, indicatorF);
  absStrainRateErrorNormLinf(result, tmp);
  clout << "strainRate-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm2D<T> relStrainRateErrorNormLinf(s, sSol, indicatorF);
  relStrainRateErrorNormLinf(result, tmp);
  clout << "; strainRate-Linf-error(rel)=" << result[0] << std::endl;

  if (flowType == nonForced) {
    // pressure error
    int Lx = converter.getLatticeLength( lx );
    int Ly = converter.getLatticeLength( ly );
    T p0 = 8.*converter.getLatticeViscosity()*converter.getCharLatticeVelocity()*Lx/T( Ly*Ly );
    AnalyticalLinear2D<T,T> pressureSol( -converter.getPhysPressure( p0 )/lx , 0 , converter.getPhysPressure( p0 ) );
    SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice,converter );

    SuperAbsoluteErrorL1Norm2D<T> absPressureErrorNormL1(pressure, pressureSol, indicatorF);
    absPressureErrorNormL1(result, tmp);
    clout << "pressure-L1-error(abs)=" << result[0];
    SuperRelativeErrorL1Norm2D<T> relPressureErrorNormL1(pressure, pressureSol, indicatorF);
    relPressureErrorNormL1(result, tmp);
    clout << "; pressure-L1-error(rel)=" << result[0] << std::endl;

    SuperAbsoluteErrorL2Norm2D<T> absPressureErrorNormL2(pressure, pressureSol, indicatorF);
    absPressureErrorNormL2(result, tmp);
    clout << "pressure-L2-error(abs)=" << result[0];
    SuperRelativeErrorL2Norm2D<T> relPressureErrorNormL2(pressure, pressureSol, indicatorF);
    relPressureErrorNormL2(result, tmp);
    clout << "; pressure-L2-error(rel)=" << result[0] << std::endl;

    SuperAbsoluteErrorLinfNorm2D<T> absPressureErrorNormLinf(pressure, pressureSol, indicatorF);
    absPressureErrorNormLinf(result, tmp);
    clout << "pressure-Linf-error(abs)=" << result[0];
    SuperRelativeErrorLinfNorm2D<T> relPressureErrorNormLinf(pressure, pressureSol, indicatorF);
    relPressureErrorNormLinf(result, tmp);
    clout << "; pressure-Linf-error(rel)=" << result[0] << std::endl;
  }
}

// Output to console and files
void getResults( SuperLattice2D<T,DESCRIPTOR>& sLattice, Dynamics<T, DESCRIPTOR>& bulkDynamics,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry2D<T>& superGeometry, Timer<T>& timer, bool hasConverged ) {

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter2D<T> vtmWriter( "forcedPoiseuille2d" );
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  const int vtmIter  = converter.getLatticeTime( maxPhysT/20. );
  const int statIter = converter.getLatticeTime( maxPhysT/20. );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
    superGeometry.rename( 0,2 );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 || hasConverged ) {
    vtmWriter.write( iT );

    SuperEuklidNorm2D<T, DESCRIPTOR> normVel( velocity );
    BlockReduction2D2D<T> planeReduction( normVel, 600, BlockDataSyncMode::ReduceOnly );
    // write output as JPEG
    heatmap::write(planeReduction, iT);
  }

  if ( hasConverged ) {
    Gnuplot<T> gplot( "centerVelocity" );
    T Ly = converter.getLatticeLength( ly );
    for ( int iY=0; iY<=Ly; ++iY ) {
      T dx = 1. / T(converter.getResolution());
      const T maxVelocity = converter.getPhysVelocity( converter.getCharLatticeVelocity() );
      T point[2]= {T(),T()};
      point[0] = lx/2.;
      point[1] = ( T )iY/Ly;
      const T radius = ly/2.;
      std::vector<T> axisPoint( 2,T() );
      axisPoint[0] = lx/2.;
      axisPoint[1] = ly/2.;
      std::vector<T> axisDirection( 2,T() );
      axisDirection[0] = 1;
      axisDirection[1] = 0;
      Poiseuille2D<T> uSol( axisPoint, axisDirection, maxVelocity, radius );
      T analytical[2] = {T(),T()};
      uSol( analytical,point );
      SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( sLattice, converter );
      AnalyticalFfromSuperF2D<T> intpolateVelocity( velocity, true );
      T numerical[2] = {T(),T()};
      intpolateVelocity( numerical,point );
      gplot.setData( iY*dx, {analytical[0],numerical[0]}, {"analytical","numerical"} );
    }
    // Create PNG file
    gplot.writePNG();
  }

  // Writes output on the console
  if ( iT%statIter==0 || hasConverged ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Error norms
    error( superGeometry, sLattice, converter, bulkDynamics );
  }
}

int main( int argc, char* argv[] ) {

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  if (argc > 1) {
    if (argv[1][0]=='-'&&argv[1][1]=='h') {
      OstreamManager clout( std::cout,"help" );
      clout<<"Usage: program [Resolution] [FlowType] [BoundaryType]"<<std::endl;
      clout<<"FlowType: 0=forced, 1=nonForced"<<std::endl;
      clout<<"BoundaryType: 0=bounceBack, 1=local, 2=interpolated, 3=freeSlip, 4=partialSlip"<<std::endl;
      clout<<"Default: FlowType=forced, Resolution=50, BoundaryType=interpolated"<<std::endl;
      return 0;
    }
  }

  if (argc > 1) {
    N = atoi(argv[1]);
    if (N < 1) {
      std::cerr << "Fluid domain is too small" << std::endl;
      return 1;
    }
  }

  if (argc > 2) {
    int flowTypeNumber = atoi(argv[2]);
    if (flowTypeNumber < 0 || flowTypeNumber > (int)nonForced) {
      std::cerr << "Unknown fluid flow type" << std::endl;
      return 2;
    }
    flowType = (FlowType) flowTypeNumber;
  }

  if (argc > 3) {
    int boundaryTypeNumber = atoi(argv[3]);
    if (boundaryTypeNumber < 0 || boundaryTypeNumber > (int) partialSlip) {
      std::cerr << "Unknown boundary type" << std::endl;
      return 3;
    }
    boundaryType = (BoundaryType) boundaryTypeNumber;
  }

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},     // resolution: number of voxels per charPhysL
    (T)   0.8,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   1,     // charPhysLength: reference length of simulation geometry
    (T)   1,     // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1./Re, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0    // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("forcedPoiseuille2d");


  // === 2nd Step: Prepare Geometry ===
  Vector<T,2> extend( lx, ly );
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid( extend, origin );

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidGeometry2D<T> cuboidGeometry( cuboid, converter.getConversionFactorLength(), noOfCuboids );


  if (flowType == forced) {
    // Periodic boundaries in x-direction
    cuboidGeometry.setPeriodicity( true, false );
  }

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry2D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( converter, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice2D<T, DESCRIPTOR> sLattice( superGeometry );

  Dynamics<T, DESCRIPTOR>* bulkDynamics;

#if defined(MRT)
  if (flowType == forced) {
    bulkDynamics = new ForcedMRTdynamics<T, DESCRIPTOR>( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() );
  } else {
    bulkDynamics = new MRTdynamics<T, DESCRIPTOR>( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() );
  }
#else
  if (flowType == forced) {
    bulkDynamics = new ForcedBGKdynamics<T, DESCRIPTOR>( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() );
  } else {
    bulkDynamics = new BGKdynamics<T, DESCRIPTOR>( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() );
  }
#endif


  // choose between local and non-local boundary condition
  sOnLatticeBoundaryCondition2D<T, DESCRIPTOR> sBoundaryCondition( sLattice );
  if (boundaryType == local) {
    createLocalBoundaryCondition2D<T, DESCRIPTOR> (sBoundaryCondition);
  } else {
    createInterpBoundaryCondition2D<T, DESCRIPTOR> ( sBoundaryCondition );
  }

  prepareLattice( converter, sLattice, *bulkDynamics, sBoundaryCondition, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << endl;
  Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( converter.getLatticeTime( physInterval ), residuum );
  timer.start();

  for ( int iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << endl;
      getResults( sLattice, *bulkDynamics, converter, iT, superGeometry, timer, converge.hasConverged() );

      break;
    }

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    // in this application no boundary conditions have to be adjusted

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, *bulkDynamics, converter, iT, superGeometry, timer, converge.hasConverged()  );
    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
  }

  timer.stop();
  timer.printSummary();
}
