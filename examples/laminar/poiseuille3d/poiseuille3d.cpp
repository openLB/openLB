/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2018 Marc Haußmann, Mathias J. Krause
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

/* poiseuille3d.cpp:
 * This example examines a 3D Poseuille flow
 * It illustrates the computation of error norms.
 */


#include "olb3D.h"
#include "olb3D.hh"

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
#define DESCRIPTOR ForcedMRTD3Q19Descriptor
#else
#define DESCRIPTOR  D3Q19<FORCE>
#endif

typedef enum {forced, nonForced} FlowType;

typedef enum {bounceBack, local, interpolated, bouzidi, freeSlip, partialSlip} BoundaryType;


// Parameters for the simulation setup
FlowType flowType = forced;
BoundaryType boundaryType = bouzidi;
const T length  = 2.;         // length of the pie
const T diameter  = 1.;       // diameter of the pipe
int N = 21;                   // resolution of the model
const T physU = 1.;           // physical velocity
const T Re = 10.;             // Reynolds number
const T physRho = 1.;         // physical density
const T tau = 0.8;            // lattice relaxation time
const T maxPhysT = 20.;       // max. simulation time in s, SI unit
const T residuum = 1e-5;      // residuum for the convergence check
const T tuner = 0.97;         // for partialSlip only: 0->bounceBack, 1->freeSlip

// Scaled Parameters
const T radius  = diameter/2.;            // radius of the pipe
const T physInterval = 0.0125*maxPhysT;   // interval for the convergence check in s


// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout(std::cout, "prepareGeometry");

  clout << "Prepare Geometry ..." << std::endl;

  Vector<T, 3> center0(-converter.getPhysDeltaX() * 0.2, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  if (flowType == forced) {
    center0[0] -= 3.*converter.getPhysDeltaX();
    center1[0] += 3.*converter.getPhysDeltaX();
  }
  IndicatorCylinder3D<T> pipe(center0, center1, radius);

  superGeometry.rename(0, 2);

  superGeometry.rename(2, 1, pipe);

  if (flowType == nonForced) {
    Vector<T, 3> origin(0, radius, radius);
    Vector<T, 3> extend = origin;

    // Set material number for inflow
    origin[0] = -converter.getPhysDeltaX() * 2;
    extend[0] = converter.getPhysDeltaX() * 2;
    IndicatorCylinder3D<T> inflow(origin, extend, radius);
    superGeometry.rename(2, 3, 1, inflow);

    // Set material number for outflow
    origin[0] = length - 2 * converter.getPhysDeltaX();
    extend[0] = length + 2 * converter.getPhysDeltaX();
    IndicatorCylinder3D<T> outflow(extend, origin, radius);
    superGeometry.rename(2, 4, 1, outflow);
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
void prepareLattice(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                    UnitConverter<T, DESCRIPTOR>const& converter,
                    Dynamics<T, DESCRIPTOR>& bulkDynamics,
                    sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& onBc,
                    sOffLatticeBoundaryCondition3D<T, DESCRIPTOR>& offBc,
                    SuperGeometry3D<T>& superGeometry)
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 1, &bulkDynamics );

  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length, radius, radius);

  std::vector<T> origin = { length, radius, radius};
  std::vector<T> axis = { 1, 0, 0 };

  CirclePoiseuille3D<T> poiseuilleU(origin, axis, converter.getCharLatticeVelocity(), radius);

  if (boundaryType == bounceBack) {
    sLattice.defineDynamics( superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>() );
  }
  else if (boundaryType == freeSlip) {
    sLattice.defineDynamics(superGeometry, 2, &instances::getNoDynamics<T, DESCRIPTOR>());
    onBc.addSlipBoundary( superGeometry, 2 );
  }
  else if (boundaryType == partialSlip) {
    sLattice.defineDynamics(superGeometry, 2, &instances::getNoDynamics<T, DESCRIPTOR>());
    onBc.addPartialSlipBoundary(tuner, superGeometry, 2 );
  }
  else if (boundaryType == bouzidi) {
    sLattice.defineDynamics(superGeometry, 2, &instances::getNoDynamics<T, DESCRIPTOR>() );

    center0[0] -= 0.5*converter.getPhysDeltaX();
    center1[0] += 0.5*converter.getPhysDeltaX();
    if (flowType == forced) {
      center0[0] -= 3.*converter.getPhysDeltaX();
      center1[0] += 3.*converter.getPhysDeltaX();
    }
    IndicatorCylinder3D<T> pipe(center0, center1, radius);
    offBc.addZeroVelocityBoundary(superGeometry, 2, pipe);
  }
  else {
    sLattice.defineDynamics( superGeometry, 2, &bulkDynamics );
    onBc.addVelocityBoundary( superGeometry, 2, omega );
  }

  if (flowType == nonForced) {
    if (boundaryType == bouzidi) {
      sLattice.defineDynamics(superGeometry, 3, &instances::getNoDynamics<T, DESCRIPTOR>() );
      IndicatorCylinder3D<T> pipe(center0, center1, radius);
      offBc.addVelocityBoundary(superGeometry, 3, pipe);
      offBc.defineU(superGeometry,3,poiseuilleU);
    }
    else {
      // Material=3 -->bulk dynamics
      sLattice.defineDynamics( superGeometry, 3, &bulkDynamics );
      onBc.addVelocityBoundary( superGeometry, 3, omega );
    }
    // Material=4 -->bulk dynamics
    sLattice.defineDynamics( superGeometry, 4, &bulkDynamics );
    onBc.addPressureBoundary( superGeometry, 4, omega );
  }

  if (flowType == forced) {
    // Initial conditions
    T D = converter.getLatticeLength(diameter);

    std::vector<T> poiseuilleForce(3, T());
    poiseuilleForce[0] = 4. * converter.getLatticeViscosity() * converter.getCharLatticeVelocity() / (D * D / 4. );
    AnalyticalConst3D<T,T> force( poiseuilleForce );

    // Initialize force
    sLattice.defineField<FORCE>(superGeometry, 1, force);
    sLattice.defineField<FORCE>(superGeometry, 2, force );


    AnalyticalConst3D<T, T> rhoF(1);

    sLattice.defineRhoU(superGeometry, 1, rhoF, poiseuilleU);
    sLattice.iniEquilibrium(superGeometry, 1, rhoF, poiseuilleU);
    sLattice.defineRhoU(superGeometry, 2, rhoF, poiseuilleU);
    sLattice.iniEquilibrium(superGeometry, 2, rhoF, poiseuilleU);
  }
  else {
    // Initial conditions
    T p0 = 4. * converter.getPhysViscosity() * converter.getCharPhysVelocity() * length / (radius * radius);

    p0 = converter.getLatticePressure(p0);
    AnalyticalLinear3D<T, T> rho(-p0 / length * descriptors::invCs2<T,DESCRIPTOR>(), 0, 0, p0 * descriptors::invCs2<T,DESCRIPTOR>() + 1);

    std::vector<T> velocity(3, T());
    AnalyticalConst3D<T, T> uF(velocity);

    // Initialize all values of distribution functions to their local equilibrium
    sLattice.defineRhoU(superGeometry, 0, rho, uF);
    sLattice.iniEquilibrium(superGeometry, 0, rho, uF);
    sLattice.defineRhoU(superGeometry, 1, rho, poiseuilleU);
    sLattice.iniEquilibrium(superGeometry, 1, rho, poiseuilleU);
    sLattice.defineRhoU(superGeometry, 2, rho, poiseuilleU);
    sLattice.iniEquilibrium(superGeometry, 2, rho, poiseuilleU);
    sLattice.defineRhoU(superGeometry, 3, rho, poiseuilleU);
    sLattice.iniEquilibrium(superGeometry, 3, rho, poiseuilleU);
    sLattice.defineRhoU(superGeometry, 4, rho, poiseuilleU);
    sLattice.iniEquilibrium(superGeometry, 4, rho, poiseuilleU);
  }

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Compute error norms
void error( SuperGeometry3D<T>& superGeometry,
            SuperLattice3D<T, DESCRIPTOR>& sLattice,
            UnitConverter<T,DESCRIPTOR> const& converter,
            Dynamics<T, DESCRIPTOR>& bulkDynamics,
            SuperLatticePhysWallShearStress3D<T,DESCRIPTOR>& wss)
{
  OstreamManager clout( std::cout,"error" );

  int tmp[]= { };
  T result[2]= { };

  // velocity error
  const T maxVelocity = converter.getCharPhysVelocity();
  std::vector<T> axisPoint = {length, radius, radius};
  std::vector<T> axisDirection = { 1, 0, 0 };
  CirclePoiseuille3D<T> uSol(axisPoint, axisDirection, maxVelocity, radius);
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> u( sLattice,converter );
  auto indicatorF = superGeometry.getMaterialIndicator(1);

  SuperAbsoluteErrorL1Norm3D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm3D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm3D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;

  // strainRate error
  CirclePoiseuilleStrainRate3D<T, DESCRIPTOR> sSol( converter, radius );
  SuperLatticePhysStrainRate3D<T,DESCRIPTOR> s( sLattice,converter );

  SuperAbsoluteErrorL1Norm3D<T> absStrainRateErrorNormL1(s, sSol, indicatorF);
  absStrainRateErrorNormL1(result, tmp);
  clout << "strainRate-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm3D<T> relStrainRateErrorNormL1(s, sSol, indicatorF);
  relStrainRateErrorNormL1(result, tmp);
  clout << "; strainRate-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absStrainRateErrorNormL2(s, sSol, indicatorF);
  absStrainRateErrorNormL2(result, tmp);
  clout << "strainRate-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relStrainRateErrorNormL2(s, sSol, indicatorF);
  relStrainRateErrorNormL2(result, tmp);
  clout << "; strainRate-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absStrainRateErrorNormLinf(s, sSol, indicatorF);
  absStrainRateErrorNormLinf(result, tmp);
  clout << "strainRate-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm3D<T> relStrainRateErrorNormLinf(s, sSol, indicatorF);
  relStrainRateErrorNormLinf(result, tmp);
  clout << "; strainRate-Linf-error(rel)=" << result[0] << std::endl;

  // wallShearStress error
  AnalyticalConst3D<T,T> wssSol(4. * converter.getPhysViscosity() * converter.getPhysDensity() * maxVelocity / diameter);
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> wssSolLattice (wssSol, sLattice);

  auto indicatorB = superGeometry.getMaterialIndicator(2);

  SuperAbsoluteErrorL1Norm3D<T> absWallShearStressErrorNormL1(wss, wssSol, indicatorB);
  absWallShearStressErrorNormL1(result, tmp);
  clout << "wss-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm3D<T> relWallShearStressErrorNormL1(wss, wssSol, indicatorB);
  relWallShearStressErrorNormL1(result, tmp);
  clout << "; wss-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absWallShearStressErrorNormL2(wss, wssSol, indicatorB);
  absWallShearStressErrorNormL2(result, tmp);
  clout << "wss-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relWallShearStressErrorNormL2(wss, wssSol, indicatorB);
  relWallShearStressErrorNormL2(result, tmp);
  clout << "; wss-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absWallShearStressErrorNormLinf(wss, wssSol, indicatorB);
  absWallShearStressErrorNormLinf(result, tmp);
  clout << "wss-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm3D<T> relWallShearStressErrorNormLinf(wss, wssSol, indicatorB);
  relWallShearStressErrorNormLinf(result, tmp);
  clout << "; wss-Linf-error(rel)=" << result[0] << std::endl;

  if (flowType == nonForced) {
    // pressure error
    T p0 = 4. * converter.getPhysViscosity() * maxVelocity * length / (radius * radius);
    AnalyticalLinear3D<T, T> pressureSol(-p0 / length, 0, 0, p0);
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);

    SuperAbsoluteErrorL1Norm3D<T> absPressureErrorNormL1(pressure, pressureSol, indicatorF);
    absPressureErrorNormL1(result, tmp);
    clout << "pressure-L1-error(abs)=" << result[0];
    SuperRelativeErrorL1Norm3D<T> relPressureErrorNormL1(pressure, pressureSol, indicatorF);
    relPressureErrorNormL1(result, tmp);
    clout << "; pressure-L1-error(rel)=" << result[0] << std::endl;

    SuperAbsoluteErrorL2Norm3D<T> absPressureErrorNormL2(pressure, pressureSol, indicatorF);
    absPressureErrorNormL2(result, tmp);
    clout << "pressure-L2-error(abs)=" << result[0];
    SuperRelativeErrorL2Norm3D<T> relPressureErrorNormL2(pressure, pressureSol, indicatorF);
    relPressureErrorNormL2(result, tmp);
    clout << "; pressure-L2-error(rel)=" << result[0] << std::endl;

    SuperAbsoluteErrorLinfNorm3D<T> absPressureErrorNormLinf(pressure, pressureSol, indicatorF);
    absPressureErrorNormLinf(result, tmp);
    clout << "pressure-Linf-error(abs)=" << result[0];
    SuperRelativeErrorLinfNorm3D<T> relPressureErrorNormLinf(pressure, pressureSol, indicatorF);
    relPressureErrorNormLinf(result, tmp);
    clout << "; pressure-Linf-error(rel)=" << result[0] << std::endl;
  }
}

// Output to console and files
void getResults( SuperLattice3D<T,DESCRIPTOR>& sLattice, Dynamics<T, DESCRIPTOR>& bulkDynamics,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry3D<T>& superGeometry, Timer<T>& timer, bool hasConverged,
                 SuperLatticePhysWallShearStress3D<T,DESCRIPTOR>& wss)
{

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "poiseuille3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );
  vtmWriter.addFunctor( wss );

  const int vtmIter  = converter.getLatticeTime( maxPhysT/20. );
  const int statIter = converter.getLatticeTime( maxPhysT/20. );

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

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 || hasConverged ) {
    vtmWriter.write( iT );

    SuperEuklidNorm3D<T, DESCRIPTOR> normVel( velocity );
    BlockReduction3D2D<T> planeReduction( normVel, {0,0,1}, 600, BlockDataSyncMode::ReduceOnly );
    // write output as JPEG
    heatmap::write(planeReduction, iT);

  }

  if ( hasConverged ) {
    Gnuplot<T> gplot( "centerVelocity" );
    T D = converter.getLatticeLength( diameter );
    for ( int iY=0; iY<=D; ++iY ) {
      T dx = 1. / T(converter.getResolution());
      T point[3]= {T(),T(),T()};
      point[0] = length/2.;
      point[1] = ( T )converter.getPhysLength(iY);
      point[2] = ( T )radius;
      const T maxVelocity = converter.getCharPhysVelocity();
      std::vector<T> axisPoint = {length, radius, radius};
      std::vector<T> axisDirection = { 1, 0, 0 };
      CirclePoiseuille3D<T> uSol(axisPoint, axisDirection, maxVelocity, radius);
      T analytical[3] = {T(),T(),T()};
      uSol( analytical,point );
      SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
      AnalyticalFfromSuperF3D<T> intpolateVelocity( velocity, true, 1 );
      T numerical[3] = {T(),T(),T()};
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
    error( superGeometry, sLattice, converter, bulkDynamics, wss );
  }
}

int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  if (argc > 1) {
    if (argv[1][0]=='-'&&argv[1][1]=='h') {
      OstreamManager clout( std::cout,"help" );
      clout<<"Usage: program [Resolution] [FlowType] [BoundaryType]"<<std::endl;
      clout<<"FlowType: 0=forced, 1=nonForced"<<std::endl;
      clout<<"BoundaryType: 0=bounceBack, 1=local, 2=interpolated, 3=bouzidi, 4=freeSlip, 5=partialSlip"<<std::endl;
      clout<<"Default: FlowType=forced, Resolution=21, BoundaryType=bouzidi"<<std::endl;
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
    int {N},                  // resolution: number of voxels per charPhysL
    (T)   tau,                // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   diameter,           // charPhysLength: reference length of simulation geometry
    (T)   physU,              // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   diameter*physU/Re,  // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   physRho             // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("poiseuille3d");


  // === 2nd Step: Prepare Geometry ===

  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  IndicatorLayer3D<T> extendedDomain(pipe, converter.getPhysDeltaX());

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = 2*singleton::mpi().getSize();
#else // ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = 6;
#endif // ifdef PARALLEL_MODE_MPI
  CuboidGeometry3D<T> cuboidGeometry(extendedDomain, converter.getPhysDeltaX(), noOfCuboids);
  if (flowType == forced) {
    // Periodic boundaries in x-direction
    cuboidGeometry.setPeriodicity( true, false, false );
  }

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  // Instantiation of a superGeometry
  SuperGeometry3D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(converter, superGeometry);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

  std::unique_ptr<Dynamics<T, DESCRIPTOR>> bulkDynamics;

#if defined(MRT)
  if (flowType == forced) {
    bulkDynamics.reset(new ForcedMRTdynamics<T, DESCRIPTOR>( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() ));
  }
  else {
    bulkDynamics.reset(new MRTdynamics<T, DESCRIPTOR>( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() ));
  }
#else
  if (flowType == forced) {
    bulkDynamics.reset(new ForcedBGKdynamics<T, DESCRIPTOR>( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() ));
  }
  else {
    bulkDynamics.reset(new BGKdynamics<T, DESCRIPTOR>( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() ));
  }
#endif


  // choose between local and non-local boundary condition
  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sOnBoundaryCondition( sLattice );
  sOffLatticeBoundaryCondition3D<T, DESCRIPTOR> sOffBoundaryCondition(sLattice);
  createBouzidiBoundaryCondition3D<T, DESCRIPTOR>(sOffBoundaryCondition);

  if (boundaryType == local) {
    createLocalBoundaryCondition3D<T, DESCRIPTOR> (sOnBoundaryCondition);
  }
  else {
    createInterpBoundaryCondition3D<T, DESCRIPTOR> ( sOnBoundaryCondition );
  }

  prepareLattice(sLattice, converter, *bulkDynamics, sOnBoundaryCondition, sOffBoundaryCondition, superGeometry);

  // set up size-increased indicator and instantiate wall shear stress functor (wss)
  Vector<T, 3> center0Extended(-converter.getPhysDeltaX() * 0.2, radius, radius);
  Vector<T, 3> center1Extended(length, radius, radius);
  if (flowType == forced) {
    center0Extended[0] -= 4.*converter.getPhysDeltaX();
    center1Extended[0] += 4.*converter.getPhysDeltaX();
  }
  IndicatorCylinder3D<T> pipeExtended(center0Extended, center1Extended, radius);
  IndicatorLayer3D<T> indicatorExtended (pipeExtended, 0.9*converter.getConversionFactorLength()*N/11.);
  SuperLatticePhysWallShearStress3D<T,DESCRIPTOR> wss(sLattice, superGeometry, 2, converter, indicatorExtended);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << endl;
  Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( converter.getLatticeTime( physInterval ), residuum );
  timer.start();

  for ( int iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << endl;
      getResults( sLattice, *bulkDynamics, converter, iT, superGeometry, timer, converge.hasConverged(), wss );

      break;
    }

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    // in this application no boundary conditions have to be adjusted

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, *bulkDynamics, converter, iT, superGeometry, timer, converge.hasConverged(), wss  );
    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
  }

  timer.stop();
  timer.printSummary();
}
