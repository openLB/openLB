/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2019 Fabian Klemens, Marc Haußmann
 *  Mathias J. Krause, Vojtech Cvrcek, Peter Weisbrod
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

/* porousPoiseuille3d.cpp:
 * This example examines a 3D Poseuille flow with porous media.
 * Two porous media LB methods can be used here:
 * Spaid and Phelan (doi:10.1063/1.869392), or
 * Guo and Zhao (doi:10.1103/PhysRevE.66.036304)
 */


#include "olb3D.h"
#include "olb3D.hh"

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace olb;

typedef double T;

#define SPAID_PHELAN
//#define GUO_ZHAO

T epsilon = 1.;
T K = 1e-3;

#ifdef SPAID_PHELAN
#define DESCRIPTOR descriptors::PorousD3Q19Descriptor
#elif defined GUO_ZHAO
#define DESCRIPTOR descriptors::GuoZhaoD3Q19Descriptor
#endif


// Parameters for the simulation setup
const T length  = 2.;         // length of the pie
const T diameter  = 1.;       // diameter of the pipe
int N = 21;                   // resolution of the model
const T physU = 1.;           // physical velocity
const T Re = 1.;             // Reynolds number
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
  IndicatorCylinder3D<T> pipe(center0, center1, radius);

  superGeometry.rename(0, 2);

  superGeometry.rename(2, 1, pipe);

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

  AnalyticalConst3D<T,T> zero(0.);
  AnalyticalConst3D<T,T> one(1.);

  T nu = (tau-0.5)/3.;
  T h = converter.getPhysDeltaX();
#ifdef SPAID_PHELAN
  T d = 1. - (h*h*nu*tau/K);
  clout << "Lattice Porosity: " << d << std::endl;
  clout << "Kmin: " << h*h*nu*tau << std::endl;
  if (K < h*h*nu*tau) {
    clout << "WARNING: Chosen K is too small!" << std::endl;
    exit(1);
  }
#endif

#ifdef SPAID_PHELAN
  AnalyticalConst3D<T,T> porosity(d);
  sLattice.defineField<descriptors::POROSITY>(superGeometry, 1, porosity);
#elif defined GUO_ZHAO
  AnalyticalConst3D<T,T> Nu(nu);
  AnalyticalConst3D<T,T> k(K/(h*h));
  AnalyticalConst3D<T,T> eps(epsilon);

  sLattice.defineField<descriptors::EPSILON>(superGeometry, 1,  eps);
  sLattice.defineField<descriptors::NU>(superGeometry, 1, 1, Nu);
  sLattice.defineField<descriptors::K>(superGeometry, 1, k);
#endif

  sLattice.defineDynamics(superGeometry, 2, &instances::getNoDynamics<T, DESCRIPTOR>() );

  center0[0] -= 0.5*converter.getPhysDeltaX();
  center1[0] += 0.5*converter.getPhysDeltaX();
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  offBc.addZeroVelocityBoundary(superGeometry, 2, pipe);

  sLattice.defineDynamics( superGeometry, 2, &bulkDynamics );
  onBc.addVelocityBoundary( superGeometry, 2, omega );

  sLattice.defineDynamics(superGeometry, 3, &instances::getNoDynamics<T, DESCRIPTOR>() );
  offBc.addVelocityBoundary(superGeometry, 3, pipe);
  offBc.defineU(superGeometry,3,poiseuilleU);

  // Material=4 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 4, &bulkDynamics );
  onBc.addPressureBoundary( superGeometry, 4, omega );

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

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

/// Compute error norms
void error( SuperGeometry3D<T>& superGeometry,
            SuperLattice3D<T, DESCRIPTOR>& sLattice,
            UnitConverter<T,DESCRIPTOR> const& converter,
            AnalyticalF3D<T,T>& porVel) {

  OstreamManager clout( std::cout,"error" );


  int tmp[]= { };
  T result[2]= { };
  auto indicatorF = superGeometry.getMaterialIndicator(1);
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> u( sLattice, converter );

  SuperAbsoluteErrorL1Norm3D<T> absVelocityErrorNormL1(u, porVel, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm3D<T> relVelocityErrorNormL1(u, porVel, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absVelocityErrorNormL2(u, porVel, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relVelocityErrorNormL2(u, porVel, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absVelocityErrorNormLinf(u, porVel, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm3D<T> relVelocityErrorNormLinf(u, porVel, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;
}


// Output to console and files
void getResults( SuperLattice3D<T,DESCRIPTOR>& sLattice, Dynamics<T, DESCRIPTOR>& bulkDynamics,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry3D<T>& superGeometry, Timer<T>& timer, bool hasConverged )
{

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "porousPoiseuille3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

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
  }


  // Writes output on the console
  if ( iT%statIter==0 || hasConverged ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Error norms
    AnalyticalFfromSuperF3D<T> intpolatePressure( pressure, true );

    T point1[3] = {0, radius, radius};
    T point2[3] = {0, radius, radius};

    point1[0] = length*0.5 - length*0.01;
    point2[0] = length*0.5 + length*0.01;

    T p1, p2;
    intpolatePressure( &p1,point1 );
    intpolatePressure( &p2,point2 );

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1-p2;
    clout << "; pressureDrop=" << pressureDrop;

    AnalyticalFfromSuperF3D<T> intpolateVelocity( velocity, true );
    T midVel[3];
    T mid[3] = {length*0.5, radius, radius};
    intpolateVelocity(midVel, mid);
    T mu = converter.getPhysViscosity()*converter.getPhysDensity();
    T l = point2[0] - point1[0];
    T vel = midVel[0];
    T pressureGradient = pressureDrop/l;

    AnalyticalPorousVelocity3D<T> porVel(superGeometry, 3, K, mu, pressureGradient, radius, epsilon);

    /// Darcy law
    T expectedPressureDrop = vel*mu*l/K;
    T darcyFlux =  K*pressureGradient/mu;
#ifdef GUO_ZHAO
    expectedPressureDrop *= epsilon;
    darcyFlux /= epsilon;
#endif

    clout << "; expected(darcy)=" << expectedPressureDrop << std::endl;
    clout << "peakVelocity=" <<  midVel[0];
    clout << "; expected(darcy)=" << darcyFlux << std::endl;
    clout << "peakVelocity(analytical)=" << porVel.getPeakVelocity();
    clout << "; peakVelocity-error(rel)=" << abs(porVel.getPeakVelocity()-vel)/porVel.getPeakVelocity() << std::endl;

    error(superGeometry, sLattice, converter, porVel);

    // Gnuplot
    Gnuplot<T> gplot( "velocityProfile" );
    T uAnalytical[3] = {};
    T uNumerical[3] = {};
    for (int i=0; i<101; i++) {
      T yInput = diameter*i/100.;
      T input[3] = {length*0.5, yInput, radius};
      porVel(uAnalytical, input);
      intpolateVelocity(uNumerical, input);
      gplot.setData( yInput, {uAnalytical[0], uNumerical[0]}, {"analytical","numerical"} );
    }

    // Create PNG file
    gplot.writePNG();
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
      clout<<"Usage: program [Resolution] [Permeability]" <<std::endl;
      clout<<"Default: Resolution=21, Permeability=1e-3" <<std::endl;
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
    K = atof(argv[2]);
    if (K < 0) {
      std::cerr << "Permeabilty must be greater than 0" << std::endl;
      return 2;
    }
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
  converter.write("porousPoiseuille3d");


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

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  // Instantiation of a superGeometry
  SuperGeometry3D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(converter, superGeometry);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

  std::unique_ptr<Dynamics<T, DESCRIPTOR>> bulkDynamics;

#ifdef SPAID_PHELAN
  bulkDynamics.reset(new PorousBGKdynamics<T, DESCRIPTOR>( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() ));
#elif defined GUO_ZHAO
  bulkDynamics.reset(new GuoZhaoBGKdynamics<T, DESCRIPTOR>( converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>() ));
#endif


  // choose between local and non-local boundary condition
  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sOnBoundaryCondition( sLattice );
  sOffLatticeBoundaryCondition3D<T, DESCRIPTOR> sOffBoundaryCondition(sLattice);
  createBouzidiBoundaryCondition3D<T, DESCRIPTOR>(sOffBoundaryCondition);

  createInterpBoundaryCondition3D<T, DESCRIPTOR> ( sOnBoundaryCondition );

  prepareLattice(sLattice, converter, *bulkDynamics, sOnBoundaryCondition, sOffBoundaryCondition, superGeometry);

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
