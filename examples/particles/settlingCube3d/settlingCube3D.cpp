/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2015 Fabian Klemens, Jonas Latt, Mathias J. Krause
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

/* settlingCube3d.cpp:
 * The case examines the settling of a cubical silica particle
 * under the influence of gravity.
 * The object is surrounded by water in a rectangular domain
 * limited by no-slip boundary conditions.
 * For the calculation of forces an DNS approach is chosen
 * which also leads to a back-coupling of the particle on the fluid,
 * inducing a flow.
 *
 * The simulation is based on the homogenised lattice Boltzmann approach
 * (HLBM) introduced in "Particle flow simulations with homogenised
 * lattice Boltzmann methods" by Krause et al.
 * and extended in "Towards the simulation of arbitrarily shaped 3D particles
 * using a homogenised lattice Boltzmann method" by Trunk et al.
 * for the simulation of 3D particles.
 *
 * This example demonstrates the usage of HLBM in the OpenLB framework.
 */

#include "olb3D.h"
#include "olb3D.hh"     // use generic version only!

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
#define DESCRIPTOR D3Q19<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR>

#define WriteVTK

// Discretization Settings
int res = 30;
T const charLatticeVelocity = 0.01;

// Time Settings
T const maxPhysT = 0.5;       // max. simulation time in s
T const iTwrite = 0.02;       // write out intervall in s

// Domain Settings
T const lengthX = 0.01;     
T const lengthY = 0.01;     
T const lengthZ = 0.05;

// Fluid Settings 
T const physDensity = 1000;
T const physViscosity = 1E-5;

//Particle Settings
T centerX = lengthX*.5;
T centerY = lengthY*.5;
T centerZ = lengthZ*.9;
T const cubeDensity = 2500;              
T const cubeEdgeLength = 0.0025; 
Vector<T,3> cubeCenter = {centerX,centerY,centerZ};
Vector<T,3> cubeOrientation = {0.,15.,0.}; 
Vector<T,3> cubeVelocity = {0.,0.,0.};
Vector<T,3> externalAcceleration = {.0, .0, -9.81 * (1. - physDensity / cubeDensity)};

// Characteristic Quantities
T const charPhysLength = lengthX;
T const charPhysVelocity = 0.15;    // Assumed maximal velocity


// Prepare geometry
void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry3D<T>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 2);
  superGeometry.rename(2, 1, 1, 1, 1);

  superGeometry.clean();
  superGeometry.innerClean();

  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}


// Set up the geometry of the simulation
void prepareLattice(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, UnitConverter<T,DESCRIPTOR> const& converter,
  Dynamics<T, DESCRIPTOR>& designDynamics,
  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& sBoundaryCondition,
  SuperGeometry3D<T>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  clout << "setting Velocity Boundaries ..." << std::endl;

  /// Material=0 -->do nothing
  sLattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>());
  sLattice.defineDynamics(superGeometry, 1, &designDynamics);
  sLattice.defineDynamics(superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>());

  clout << "Prepare Lattice ... OK" << std::endl;
}


//Set Boundary Values
void setBoundaryValues(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                       UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                       SuperGeometry3D<T>& superGeometry)
{
  OstreamManager clout(std::cout, "setBoundaryValues");

  if (iT == 0) {
    AnalyticalConst3D<T, T> zero(0.);
    AnalyticalConst3D<T, T> one(1.);
    sLattice.defineField<descriptors::POROSITY>(superGeometry.getMaterialIndicator({0,1,2}), one);    
    // Set initial condition
    AnalyticalConst3D<T, T> ux(0.);
    AnalyticalConst3D<T, T> uy(0.);
    AnalyticalConst3D<T, T> uz(0.);
    AnalyticalConst3D<T, T> rho(1.);
    AnalyticalComposed3D<T, T> u(ux, uy, uz);

    //Initialize all values of distribution functions to their local equilibrium
    sLattice.defineRhoU(superGeometry, 1, rho, u);
    sLattice.iniEquilibrium(superGeometry, 1, rho, u);

    // Make the lattice ready for simulation
    sLattice.initialize();
  }
}


/// Computes the pressure drop between the voxels before and after the cylinder
void getResults(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                SuperGeometry3D<T>& superGeometry, Timer<double>& timer,
                ParticleDynamics3D<T, DESCRIPTOR> particleDynamics)
{
  OstreamManager clout(std::cout, "getResults");

#ifdef WriteVTK
  SuperVTMwriter3D<T> vtkWriter("sedimentation");
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
  SuperLatticePhysExternalPorosity3D<T, DESCRIPTOR> externalPor(sLattice, converter);
  vtkWriter.addFunctor(velocity);
  vtkWriter.addFunctor(pressure);
  vtkWriter.addFunctor(externalPor);

  if (iT == 0) {
    /// Writes the converter log file
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry(sLattice, superGeometry);
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank3D<T, DESCRIPTOR> rank(sLattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
  }

  if (iT % converter.getLatticeTime(iTwrite) == 0) {
    vtkWriter.write(iT);
  }
#endif

  /// Writes output on the console
  if (iT % converter.getLatticeTime(iTwrite) == 0) {
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    particleDynamics.print(); 
  }
}

int main(int argc, char* argv[])
{
  /// === 1st Step: Initialization ===
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  UnitConverterFromResolutionAndLatticeVelocity<T,DESCRIPTOR> converter(
    (int)   res,                  //resolution
    ( T )   charLatticeVelocity,  //charLatticeVelocity
    ( T )   charPhysLength,       //charPhysLength
    ( T )   charPhysVelocity,     //charPhysVelocity
    ( T )   physViscosity,        //physViscosity
    ( T )   physDensity           //physDensity
  );
  converter.print();

  /// === 2rd Step: Prepare Geometry ===
  /// Instantiation of a cuboidGeometry with weights
  std::vector<T> extend(3, T());
  extend[0] = lengthX;
  extend[1] = lengthY;
  extend[2] = lengthZ;
  std::vector<T> origin(3, T());
  IndicatorCuboid3D<T> cuboid(extend, origin);

#ifdef PARALLEL_MODE_MPI
  CuboidGeometry3D<T> cuboidGeometry(cuboid, converter.getConversionFactorLength(), singleton::mpi().getSize());
#else
  CuboidGeometry3D<T> cuboidGeometry(cuboid, converter.getConversionFactorLength(), 7);
#endif
  cuboidGeometry.print();

  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);
  SuperGeometry3D<T> superGeometry(cuboidGeometry, loadBalancer, 2);
  prepareGeometry(converter, superGeometry);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice3D<T, DESCRIPTOR> sLattice(superGeometry);

  PorousParticleBGKdynamics<T, DESCRIPTOR, false> designDynamics(converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>());

  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sBoundaryCondition(sLattice);
  createLocalBoundaryCondition3D<T, DESCRIPTOR>(sBoundaryCondition);

  prepareLattice(sLattice, converter, designDynamics, sBoundaryCondition, superGeometry);

  /// === 4th Step: Main Loop with Timer ===
  Timer<double> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel());
  timer.start();

  // Create Particle Dynamics
  ParticleDynamics3D<T, DESCRIPTOR> particleDynamics(sLattice, converter, superGeometry, lengthX, lengthY, lengthZ, externalAcceleration);
  
  // Create Cube Indicator
  T epsilon = 0.5*converter.getConversionFactorLength();
  
  //Cube indicator
  SmoothIndicatorCuboid3D<T, T, true> particleIndicator(cubeCenter, cubeEdgeLength, cubeEdgeLength, cubeEdgeLength, epsilon, cubeOrientation, cubeDensity, cubeVelocity);
  
  //Sphere indicator
  //SmoothIndicatorSphere3D<T, T, true> particleIndicator(cubeCenter, 0.5*cubeEdgeLength, epsilon, cubeDensity, cubeVelocity);
  
  //Cylinder indicator
  //SmoothIndicatorCylinder3D<T, T, true> particleIndicator(cubeCenter, { 1, 0, 0 }, 0.5*cubeEdgeLength, cubeEdgeLength, epsilon, cubeOrientation, cubeDensity, cubeVelocity);

  SuperExternal3D<T,DESCRIPTOR,POROSITY> superExtPorosity(superGeometry, sLattice, sLattice.getOverlap());
  SuperExternal3D<T,DESCRIPTOR,VELOCITY_NUMERATOR> superExtNumerator(superGeometry, sLattice, sLattice.getOverlap());
  SuperExternal3D<T,DESCRIPTOR,VELOCITY_DENOMINATOR> superExtDenominator(superGeometry, sLattice, sLattice.getOverlap());
  particleDynamics.addParticle( particleIndicator );
  particleDynamics.print();

  /// === 5th Step: Definition of Initial and Boundary Conditions ===
  setBoundaryValues(sLattice, converter, 0, superGeometry);

  clout << "MaxIT: " << converter.getLatticeTime(maxPhysT) << std::endl;
  for (int iT = 0; iT < converter.getLatticeTime(maxPhysT)+10; ++iT) {
    particleDynamics.simulateTimestep("verlet");
    getResults(sLattice, converter, iT, superGeometry, timer, particleDynamics);
    sLattice.collideAndStream();
    superExtPorosity.communicate();
    superExtNumerator.communicate();
    superExtDenominator.communicate();
  }

  timer.stop();
  timer.printSummary();
}
