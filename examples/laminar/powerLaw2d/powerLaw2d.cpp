/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2016 Vojtech Cvrcek, Mathias J. Krause
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

/* powerLaw2d.cpp:
 * This example examines a steady flow of a non-newtonian fluid in a channel.
 * At the inlet, a profile for non-newtonian fluid is imposed on the velocity,
 * where as the outlet implements an outflow condition grad_x p = 0.
 * The power law model is for n=1 and m=charNu in fact the classical poiseuille flow.
 * One can validate the error with using functors in void error.
 *
 *
 */

#include "olb2D.h"
#include "olb2D.hh"   // include full template code

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
#define DESCRIPTOR DynOmegaD2Q9Descriptor

// Parameters for the simulation setup
int N = 40;            // resolution of the model
T Re = 10.;      // Reynolds number
T tau = 0.8;
T lx = 2.;             // channel lenght
T ly = 1.;             // channel width
T maxU = 1;     // Max velocity
T Tmax = 20;      // max. phys. time in s
T Tprint = 1;      // Phys time at which the status of the system is print
// set the changes for n and m in powerLawBGKdynamics.h
T n = .2;              // parameter in power law model (n=1 Newtonian fluid)
bool bcTypePeriodic = false; //true works only with one core

const T residuum = 1e-5;      // residuum for the convergence check

void prepareGeometry( PowerLawUnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry2D<T>& superGeometry ) {
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T,2> extend( lx, ly );
  Vector<T,2> origin;

  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,1,1 );

  // Set material number for inflow
  extend[0] = 1.2*converter.getConversionFactorLength();
  origin[0] = -converter.getConversionFactorLength();
  IndicatorCuboid2D<T> inflow( extend, origin );
  if (bcTypePeriodic)
    superGeometry.rename( 1,3,inflow );
  else
    superGeometry.rename( 2,3,1,inflow );
  // Set material number for outflow
  origin[0] = lx-.5*converter.getConversionFactorLength();
  IndicatorCuboid2D<T> outflow( extend, origin );
  if (bcTypePeriodic)
    superGeometry.rename( 1,4,outflow );
  else
    superGeometry.rename( 2,4,1,outflow );
  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();

  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice2D<T,DESCRIPTOR>& sLattice,
                     PowerLawUnitConverter<T,DESCRIPTOR> const& converter,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics,
                     Dynamics<T, DESCRIPTOR>& inDynamics,
                     Dynamics<T, DESCRIPTOR>& outDynamics,
                     sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& sBoundaryCondition,
                     SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry.getMaterialIndicator(0), &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics( superGeometry.getMaterialIndicator(1), &bulkDynamics );

  // Material=2 -->bounce back
  sLattice.defineDynamics( superGeometry.getMaterialIndicator(2), &instances::getBounceBack<T, DESCRIPTOR>() );

  // Material=3 -->bulk dynamics (inflow)
  if (bcTypePeriodic)
    sLattice.defineDynamics( superGeometry.getMaterialIndicator(3), &inDynamics );
  else {
    sLattice.defineDynamics( superGeometry.getMaterialIndicator(3), &bulkDynamics );
    // Setting of the boundary conditions
    sBoundaryCondition.addVelocityBoundary( superGeometry, 3, omega );
  }

  // Material=4 -->bulk dynamics (outflow)
  if (bcTypePeriodic)
    sLattice.defineDynamics( superGeometry.getMaterialIndicator(4), &outDynamics );
  else {
    sLattice.defineDynamics( superGeometry.getMaterialIndicator(4), &bulkDynamics );
    // Setting of the boundary conditions
    sBoundaryCondition.addPressureBoundary( superGeometry, 4, omega );
  }
  clout << "Prepare Lattice ... OK" << std::endl;
}


void setBoundaryValues( SuperLattice2D<T, DESCRIPTOR>& sLattice,
                        PowerLawUnitConverter<T,DESCRIPTOR> const& converter,
                        int iT, SuperGeometry2D<T>& superGeometry ) {

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // Set initial and steady boundary conditions
  if ( iT==0 ) {

    // Define the analytical solutions for pressure and velocity
    T maxVelocity = converter.getCharLatticeVelocity();
    T distance2Wall = converter.getConversionFactorLength()/2.;

    T p0 = converter.getPhysConsistencyCoeff()*pow( converter.getCharPhysVelocity(),n )*pow( ( n + 1. )/n,n )*pow( 2./( ly-distance2Wall*2 ),n + 1. );

    AnalyticalLinear2D<T,T> rho( converter.getLatticeDensityFromPhysPressure( -p0 ) - 1., 0, converter.getLatticeDensityFromPhysPressure(  p0*(lx + distance2Wall*2.)/2. ) );

    PowerLaw2D<T> u( superGeometry, 3, maxVelocity, distance2Wall, ( n + 1. )/n );

    // Set the analytical solutions for pressure and velocity
    AnalyticalConst2D<T,T> omega0( converter.getLatticeRelaxationFrequency() );
    sLattice.defineField<descriptors::OMEGA>( superGeometry, 1, omega0 );
    sLattice.defineField<descriptors::OMEGA>( superGeometry, 3, omega0 );
    sLattice.defineField<descriptors::OMEGA>( superGeometry, 4, omega0 );

    // Set the analytical solutions for pressure and velocity
    // Initialize all values of distribution functions to their local equilibrium

    sLattice.defineRhoU( superGeometry, 1, rho, u );
    sLattice.iniEquilibrium( superGeometry, 1, rho, u );

    sLattice.iniEquilibrium( superGeometry, 3, rho, u );
    sLattice.defineRhoU( superGeometry, 3, rho, u );

    sLattice.iniEquilibrium( superGeometry, 4, rho, u );
    sLattice.defineRhoU( superGeometry, 4, rho, u );

    // Make the lattice ready for simulation
    sLattice.initialize();
  }
}

// Compute error norms
void error( SuperGeometry2D<T>& superGeometry,
            SuperLattice2D<T, DESCRIPTOR>& sLattice,
            PowerLawUnitConverter<T,DESCRIPTOR> const& converter,
            Dynamics<T, DESCRIPTOR>& bulkDynamics ) {
  OstreamManager clout( std::cout,"error" );

  int input[1] = { };
  T result[1]  = { };

  T distance2Wall = converter.getConversionFactorLength()/2.;

  PowerLaw2D<T> uSol( superGeometry,3,converter.getCharPhysVelocity(),distance2Wall,( n + 1. )/n );
  SuperLatticePhysVelocity2D<T,DESCRIPTOR> u( sLattice,converter );

  T p0 = converter.getPhysConsistencyCoeff()*pow( converter.getCharPhysVelocity(),n )*pow( ( n + 1. )/n,n )*pow( 2./( ly-distance2Wall*2 ),n + 1. );
  AnalyticalLinear2D<T,T> pressureSol( -p0, 0, p0*(lx + distance2Wall*2.)/2. );
  SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice,converter );

  auto indicatorF = superGeometry.getMaterialIndicator(1);

  // velocity error
  SuperAbsoluteErrorL1Norm2D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, input);
  clout << "velocity-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm2D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, input);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm2D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, input);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, input);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, input);
  clout << "velocity-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm2D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, input);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;

  // pressure error
  SuperAbsoluteErrorL1Norm2D<T> absPressureErrorNormL1(pressure, pressureSol, indicatorF);
  absPressureErrorNormL1(result, input);
  clout << "pressure-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm2D<T> relPressureErrorNormL1(pressure, pressureSol, indicatorF);
  relPressureErrorNormL1(result, input);
  clout << "; pressure-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm2D<T> absPressureErrorNormL2(pressure, pressureSol, indicatorF);
  absPressureErrorNormL2(result, input);
  clout << "pressure-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relPressureErrorNormL2(pressure, pressureSol, indicatorF);
  relPressureErrorNormL2(result, input);
  clout << "; pressure-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absPressureErrorNormLinf(pressure, pressureSol, indicatorF);
  absPressureErrorNormLinf(result, input);
  clout << "pressure-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm2D<T> relPressureErrorNormLinf(pressure, pressureSol, indicatorF);
  relPressureErrorNormLinf(result, input);
  clout << "; pressure-Linf-error(rel)=" << result[0] << std::endl;
}

// Output to console and files
void getResults( SuperLattice2D<T, DESCRIPTOR>& sLattice,
                 Dynamics<T, DESCRIPTOR>& bulkDynamics,
                 PowerLawUnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry2D<T>& superGeometry, Timer<double>& timer ) {
  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter2D<T> vtmWriter( "powerLaw2d" );
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure( sLattice, converter );

  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  if ( iT==0 ) {
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry( sLattice,superGeometry );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  if ( iT%converter.getLatticeTime( Tprint )==0 ) {
    vtmWriter.write( iT );

    SuperEuklidNorm2D<T, DESCRIPTOR> normVel( velocity );
    BlockReduction2D2D<T> planeReduction( normVel, 600, BlockDataSyncMode::ReduceOnly );
    // write output of velocity as JPEG
    heatmap::write(planeReduction, iT);
  }

  // Writes output on the console
  if ( iT%converter.getLatticeTime( Tprint )==0 ) {
    timer.update( iT );
    timer.printStep();
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
    error( superGeometry, sLattice, converter, bulkDynamics );
  }
  return;
}


int main( int argc, char* argv[] ) {

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  /*
  if ( argc > 1 ) {
    N = atoi( argv[1] );
  }
  if ( argc > 2 ) {
    n = atof( argv[2] );
  }
  */

  singleton::directories().setOutputDir( "./tmp/" );

  XMLreader config("input.xml");
  config["geometry"]["N"].read(N);
  config["geometry"]["lx"].read(lx);
  config["geometry"]["ly"].read(ly);
  config["dynamics"]["maxU"].read(maxU);
  config["dynamics"]["Re"].read(Re);
  config["dynamics"]["tau"].read(tau);
  config["dynamics"]["n"].read(n);
  config["time"]["Tmax"].read(Tmax);
  config["time"]["Tprint"].read(Tprint);

  /*
  T physDeltaX = (ly/N);
  T m = ly*maxU*pow( maxU/(2*ly),1-n )/Re;
  T physViscosity = m*pow( maxU/(2*ly),n-1 );
  T physDeltaT = (tau - 0.5) / DESCRIPTOR::invCs2 * pow(physDeltaX,2) / physViscosity;
  PowerLawUnitConverter<T, DESCRIPTOR> const converter(
    (T)   physDeltaX,
    (T)   physDeltaT,
   (T)   ly,
   (T)   maxU,
   (T)   m,
   (T)   n,     // power-law index
    (T)   1.0    // physDensity: physical density in __kg / m^3__
  );
  */
  PowerLawUnitConverterFrom_Resolution_RelaxationTime_Reynolds_PLindex<T, DESCRIPTOR> const converter(
  int {N},     // resolution: number of voxels per charPhysL
  (T)   tau,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
  (T)   ly,     // charPhysLength: reference length of simulation geometry
  (T)   maxU,     // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
  (T)   Re,        // Reynolds number
  (T)   n,     // power-law index
  (T)   1.0    // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("powerLaw2d");


  // === 2rd Step: Prepare Geometry ===
  // Instantiation of a cuboidGeometry with weights

  Vector<T,2> extend( lx, ly );
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid( extend, origin );

#ifdef PARALLEL_MODE_MPI
  CuboidGeometry2D<T> cuboidGeometry( cuboid, converter.getConversionFactorLength(), singleton::mpi().getSize() );
#else
  CuboidGeometry2D<T> cuboidGeometry( cuboid, converter.getConversionFactorLength(), 7 );
#endif

  // Periodic boundaries in x-direction
  if (bcTypePeriodic)
    cuboidGeometry.setPeriodicity( true, false );

  //cuboidGeometry.printExtended();

  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );
  SuperGeometry2D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );
  prepareGeometry( converter, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice2D<T, DESCRIPTOR> sLattice( superGeometry );

  T distance2Wall = converter.getConversionFactorLength()/2.;
  T p0 = converter.getPhysConsistencyCoeff()*pow( converter.getCharPhysVelocity(),n )*pow( ( n + 1. )/n,n )*pow( 2./( ly-distance2Wall*2 ),n + 1. );

  clout << "Dimensionalized version-1." << std::endl;
  PowerLawBGKdynamics<T, DESCRIPTOR> bulkDynamics( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>(), converter.getLatticeConsistencyCoeff(), n );

  PeriodicPressureDynamics<T, DESCRIPTOR, PowerLawBGKdynamics<T,DESCRIPTOR>> outDynamics( bulkDynamics,converter.getLatticeDensityFromPhysPressure( p0*(lx + distance2Wall*2.))-1,1,0);
  PeriodicPressureDynamics<T, DESCRIPTOR, PowerLawBGKdynamics<T,DESCRIPTOR>> inDynamics( bulkDynamics,-converter.getLatticeDensityFromPhysPressure( p0*(lx + distance2Wall*2. ))+1,-1,0);
  std::cout << -converter.getLatticeDensityFromPhysPressure( p0 )+1 << std::endl;

  sOnLatticeBoundaryCondition2D<T, DESCRIPTOR> sBoundaryCondition( sLattice );
  createLocalBoundaryCondition2D<T, DESCRIPTOR, PowerLawBGKdynamics<T,DESCRIPTOR> > ( sBoundaryCondition );

  prepareLattice( sLattice, converter, bulkDynamics, inDynamics, outDynamics, sBoundaryCondition, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  Timer<double> timer( converter.getLatticeTime( Tmax ), superGeometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( converter.getLatticeTime( 0.1*Tprint ), residuum );
  timer.start();

  for ( int iT=0; iT<converter.getLatticeTime( Tmax ); ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << endl;
      getResults( sLattice, bulkDynamics, converter, iT, superGeometry, timer );

      break;
    }

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLattice, converter, iT, superGeometry );
    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, bulkDynamics, converter, iT, superGeometry, timer );
    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    if (bcTypePeriodic)
      sLattice.stripeOffDensityOffset ( sLattice.getStatistics().getAverageRho()-(T)1 );
  }
  timer.stop();
  timer.printSummary();

  // === 7th Step: Gnuplot ===
  Gnuplot<T> gplot( "centerVelocity" );
  T Ly = converter.getLatticeLength( ly );
  for ( int iY=0; iY<=Ly; ++iY ) {
    T dx = 1. / T(converter.getResolution());
    // const T maxVelocity = converter.getPhysVelocity( converter.getCharLatticeVelocity() );
    T point[2]= {T(),T()};
    point[0] = .9*lx;
    point[1] = ( T )iY/Ly;
    // const T radius = ly/2.;
    std::vector<T> axisPoint( 2,T() );
    axisPoint[0] = lx/2.;
    axisPoint[1] = ly/2.;
    std::vector<T> axisDirection( 2,T() );
    axisDirection[0] = 1;
    axisDirection[1] = 0;
    T distance2Wall = converter.getConversionFactorLength()/2.;
    PowerLaw2D<T> uSol( superGeometry,3,converter.getCharPhysVelocity(),distance2Wall,( n + 1. )/n );
    //Poiseuille2D<T> uSol( axisPoint, axisDirection, maxVelocity, radius ); // <---
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
