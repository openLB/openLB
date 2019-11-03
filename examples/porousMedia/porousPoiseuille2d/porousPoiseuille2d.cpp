/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2017 Davide Dapelo, Mathias J. Krause
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

/* porousPoiseuille2d.cpp:
 * Poiseuille flow through porous media.
 * This implementation is the reproduction of the Guo and Zhao (2002)'s
 * benchmark example A. The theoretical maximum velocity is calculated
 * as in Equation 21, and the velocity profile as in Equation 23 of
 * the original reference.
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
#define DESCRIPTOR GuoZhaoD2Q9Descriptor
#define DYNAMICS GuoZhaoBGKdynamics
//#define DYNAMICS SmagorinskyGuoZhaoBGKdynamics

/// Functional to calculate velocity profile on pipe with porous media.
template <typename T>
class PorousPipe2D : public AnalyticalF2D<T,T> {
protected:
  std::vector<T> axisPoint;
  std::vector<T> axisDirection;
  T radius, rFactor, u0;

public:
  PorousPipe2D(std::vector<T> axisPoint_, std::vector<T> axisDirection_, T radius_, T rFactor_, T u0_);
  bool operator()(T output[], const T x[]) override;
};

template <typename T>
PorousPipe2D<T>::PorousPipe2D(std::vector<T> axisPoint_, std::vector<T> axisDirection_, T radius_, T rFactor_, T u0_)
  : AnalyticalF2D<T,T>(2)
{
  this->getName() = "PorousPipe2D";
  axisPoint.resize(2);
  axisDirection.resize(2);
  for (int i = 0; i < 2; ++i) {
    axisDirection[i] = axisDirection_[i];
    axisPoint[i] = axisPoint_[i];
  }
  radius = radius_;
  rFactor = rFactor_;
  u0 = u0_;
}

template <typename T>
bool PorousPipe2D<T>::operator()(T output[], const T x[])
{
  output[0] = axisDirection[0]*u0*(cosh(rFactor*radius) - cosh(rFactor*x[1] - rFactor*radius))/(cosh(rFactor*radius) - 1);
  output[1] = axisDirection[1]*u0*(cosh(rFactor*radius) - cosh(rFactor*x[0] - rFactor*radius))/(cosh(rFactor*radius) - 1);

  return true;
}

/// Stores geometry information in form of material numbers
void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter, T lx, T ly,
                     SuperGeometry2D<T>& superGeometry)
{

  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,2);

  std::vector<T> extend(2,T());
  extend[0] = lx;
  extend[1] = ly - 1.8*converter.getPhysLength(1);
  std::vector<T> origin(2,T());
  origin[1] = 0.9*converter.getPhysLength(1);
  IndicatorCuboid2D<T> cuboid2(extend, origin);

  superGeometry.rename(2,1,cuboid2);

  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

/// Set up the geometry of the simulation
void prepareLattice(UnitConverter<T,DESCRIPTOR> const& converter, T lx, T ly,
                    T epsilonIn, T KIn, T bodyForceIn,
                    SuperLattice2D<T, DESCRIPTOR>& sLattice,
                    Dynamics<T, DESCRIPTOR>& bulkDynamics,
                    sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& sBoundaryCondition,
                    SuperGuoZhaoInstantiator2D<T, DESCRIPTOR, DYNAMICS<T, DESCRIPTOR> >& sGuoZhaoInstantiator,
                    SuperGeometry2D<T>& superGeometry )
{

  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  T   omega = converter.getLatticeRelaxationFrequency();

  /// Material=0 -->do nothing
  sLattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>());

  /// Material=1 -->bulk dynamics
  sLattice.defineDynamics(superGeometry, 1, &bulkDynamics);

  /// Material=2 -->bulk dynamics
  sLattice.defineDynamics(superGeometry, 2, &bulkDynamics);

  /// Setting of the boundary conditions
  sBoundaryCondition.addVelocityBoundary(superGeometry, 2, omega);

  /// Initial conditions
  std::vector<T> epsilonValue(1, epsilonIn);
  AnalyticalConst2D<T,T> epsilon(epsilonValue);

  std::vector<T> KValue(1, KIn);
  AnalyticalConst2D<T,T> K(KValue);

  std::vector<T> bodyForceValue (2, (T)0);
  bodyForceValue[0] = bodyForceIn;
  AnalyticalConst2D<T,T> bodyForce(bodyForceValue);

  // Initialize porosity
  sGuoZhaoInstantiator.defineEpsilon(superGeometry, 1, epsilon);
  sGuoZhaoInstantiator.defineEpsilon(superGeometry, 2, epsilon);

  sGuoZhaoInstantiator.defineK(converter, superGeometry, 1, K);
  sGuoZhaoInstantiator.defineK(converter, superGeometry, 2, K);

  sGuoZhaoInstantiator.defineNu(converter, superGeometry, 1);
  sGuoZhaoInstantiator.defineNu(converter, superGeometry, 2);

  sGuoZhaoInstantiator.defineBodyForce(converter, superGeometry, 1, bodyForce);
  sGuoZhaoInstantiator.defineBodyForce(converter, superGeometry, 2, bodyForce);

  /// Make the lattice ready for simulation
  clout << "Ready to initialize the lattice..." << std::endl;
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

/// Compute error norms
void error( SuperGeometry2D<T>& superGeometry,
            SuperLattice2D<T, DESCRIPTOR>& sLattice,
            UnitConverter<T,DESCRIPTOR> const& converter,
            Dynamics<T, DESCRIPTOR>& bulkDynamics,
            AnalyticalF2D<T,T>& uSol) {

  OstreamManager clout( std::cout,"error" );

  int input[1] = { };
  T result[1]  = { };

  SuperLatticePhysVelocity2D<T,DESCRIPTOR> u( sLattice,converter );
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
  clout << "velocity-Linf-error(abs)=" << result[0] << std::endl;
}

/// Output to console and files
void getResults(SuperLattice2D<T,DESCRIPTOR>& sLattice, Dynamics<T, DESCRIPTOR>& bulkDynamics,
                UnitConverter<T,DESCRIPTOR> const& converter, T lx, T ly, T G, T K, T nu, T epsilon, T maxPhysT, int iT, int numOfIterations,
                SuperGeometry2D<T>& superGeometry, Timer<T>& timer, bool hasConverged)
{

  OstreamManager clout(std::cout,"getResults");

  SuperVTMwriter2D<T> vtkWriter("porousPoiseuille2d");
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(sLattice, converter);
//  SuperLatticeEpsilon2D<T, DESCRIPTOR> epsilonVTM(sLattice, converter);
//  SuperLatticePhysK2D<T, DESCRIPTOR> KVTM(sLattice, converter);
//  SuperLatticePhysBodyForce2D<T, DESCRIPTOR> bodyForce(sLattice, converter);
  vtkWriter.addFunctor( velocity );
  vtkWriter.addFunctor( pressure );
//  vtkWriter.addFunctor( epsilonVTM );
//  vtkWriter.addFunctor( KVTM );
//  vtkWriter.addFunctor( bodyForce );

  const int vtkIter  = converter.getLatticeTime(maxPhysT/numOfIterations);
  const int statIter = converter.getLatticeTime(maxPhysT/numOfIterations);

  static Gnuplot<T> gplot_uCentre( "uCentre" );
  static Gnuplot<T> gplot_profile( "profile" );

  if (iT==0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry(sLattice, superGeometry);
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank2D<T, DESCRIPTOR> rank(sLattice);
    superGeometry.rename(0,2);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  /// Writes the vtk files and profile text file
 T Ly = converter.getLatticeLength(ly);
    AnalyticalFfromSuperF2D<T> intpolateVelocity( velocity, true );
    T centre[2] = {converter.getCharPhysLength()/2, converter.getCharPhysLength()/2};
    T uCentre[2];
    intpolateVelocity(uCentre, centre);
    T r = sqrt(epsilon/K);
    T dx = converter.getPhysDeltaX();
    const T radius = ly/2.;
    std::vector<T> axisPoint(2,T());
    axisPoint[0] = lx/2.;
    axisPoint[1] = ly/2.;
    std::vector<T> axisDirection(2,T());
    axisDirection[0] = 1;
    axisDirection[1] = 0;
    PorousPipe2D<T> uSol(axisPoint, axisDirection, radius, r, uCentre[0]);


  if (iT%vtkIter==0 || hasConverged) {
    vtkWriter.write(iT);

    SuperEuklidNorm2D<T, DESCRIPTOR> normVel(velocity);
    BlockReduction2D2D<T> planeReduction( normVel, 600, BlockDataSyncMode::ReduceOnly );
    // write output of velocity as JPEG
    heatmap::write(planeReduction, iT);

    ofstream *ofile = nullptr;
    if (singleton::mpi().isMainProcessor()) {
      ofile = new ofstream((singleton::directories().getLogOutDir()+"centerVel.dat").c_str());
    }
    for (int iY=0; iY<=Ly; ++iY) {
      T point[2]= {T(),T()};
      point[0] = lx/2.;
      point[1] = (T)iY*converter.getPhysLength(1);
      T analytical[2] = {T(),T()};
      uSol(analytical,point);
      SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(sLattice, converter);
      AnalyticalFfromSuperF2D<T> intpolateVelocity(velocity, true);
      T numerical[2] = {T(),T()};
      intpolateVelocity(numerical,point);
      if (singleton::mpi().isMainProcessor()) {
        *ofile << iY*dx << " " << analytical[0]
               << " " << numerical[0] << "\n";
        if ( iT == .8*converter.getLatticeTime( maxPhysT ) ) {
//          if ( iT == converter.numTimeSteps( maxPhysT )-1 ) {
          gplot_profile.setData( point[1], {numerical[0], analytical[0]}, {"Numerical profile", "Analytical profile"}, "bottom right" );
          gplot_profile.writePNG();
        }
      }
    }
    delete ofile;
  }

  /// Writes output on the console
  if (iT%statIter==0 || hasConverged) {
    /// Timer console output
    timer.update(iT);
    timer.printStep();

    /// Lattice statistics console output
    sLattice.getStatistics().print(iT,converter.getPhysTime(iT));

    /// Error norms
    error(superGeometry, sLattice, converter, bulkDynamics, uSol);

    AnalyticalFfromSuperF2D<T> intpolatePressure( pressure, true );
    AnalyticalFfromSuperF2D<T> intpolateVelocity( velocity, true );
    T centre[2] = {converter.getCharPhysLength()/2, converter.getCharPhysLength()/2};
    T uCentre[2];
    intpolateVelocity(uCentre, centre);

    T uMaxMeas = sLattice.getStatistics().getMaxU() / converter.getCharLatticeVelocity();
    T uMaxTheo = G*K/nu*(1.-1./cosh(converter.getCharPhysLength()/2 * sqrt(epsilon/K)));
    clout << "uMaxTheo=" << uMaxTheo
          << "; uMaxMeas=" << uMaxMeas
          << "; uCentre="  << uCentre[0]
          << std::endl;

    gplot_uCentre.setData( converter.getPhysTime( iT ), uCentre[0], "Centre velocity", "bottom right" );
    gplot_uCentre.writePNG( iT, maxPhysT );
  }
}

int main(int argc, char* argv[])
{

  /// === 1st Step: Initialization ===
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout,"main");
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  // Parameters for the simulation setup
  int nx;      // length of the channel
  int ny;      // height of the channel
  int N;       // resolution of the model
  T tau;       // Relaxation time
  T Re;        // Reynolds number
  T Da;        // Darcy number

  T epsilon;        // Porosity (non-dimensional)
  T K;              // Permeability (SI units)

  T maxPhysT; // max. simulation time in s, SI unit
  int numOfIterations; // number of iterations reported in paraview-Gnuplot.
  T residuum;

  string fName("input.xml");
  XMLreader config(fName);

  config["setup"]["nx"].read(nx);
  config["setup"]["ny"].read(ny);
  config["setup"]["tau"].read(tau);
  config["setup"]["N"].read(N);
  config["setup"]["Da"].read(Da);
  config["setup"]["Re"].read(Re);
  config["porous"]["epsilon"].read(epsilon);
  config["porous"]["K"].read(K);
  config["time"]["maxPhysT"].read(maxPhysT);
  config["time"]["numOfIterations"].read(numOfIterations);
  config["convergence"]["residuum"].read(residuum);

  T charL = sqrt(K/Da);
  T lx = nx*charL;
  T ly = ny*charL;

  //T latticeU = Re*(tau-(T).5)/((T)3*N);
  T charU = (T)1.;
  T bodyForceValue = charU*charU*charL / ( Re*K*((T)1. - (T)1./cosh(charL/(T)2 * sqrt(epsilon/K))) );
  T nu = bodyForceValue*K/charU*((T)1. - (T)1./cosh(charL/(T)2 * sqrt(epsilon/K)));

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},     // resolution: number of voxels per charPhysL
    (T)   tau,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   charL, // charPhysLength: reference length of simulation geometry
    (T)   charU, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   nu,    // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0    // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("porousPoiseuille2d");

  /// === 2nd Step: Prepare Geometry ===
  std::vector<T> extend(2,T());
  extend[0] = lx;
  extend[1] = ly;
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extend, origin);

  /// Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 7;
#endif
  CuboidGeometry2D<T> cuboidGeometry(cuboid, converter.getPhysDeltaX(), noOfCuboids);

  /// Periodic boundaries in x-direction
  cuboidGeometry.setPeriodicity(true, false);

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  /// Instantiation of a superGeometry
  SuperGeometry2D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(converter, lx, ly, superGeometry);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice2D<T, DESCRIPTOR> sLattice(superGeometry);

  DYNAMICS<T, DESCRIPTOR> bulkDynamics (
    converter.getLatticeRelaxationFrequency(),
    instances::getBulkMomenta<T,DESCRIPTOR>()
  );

  SuperGuoZhaoInstantiator2D<T, DESCRIPTOR, DYNAMICS<T, DESCRIPTOR> > sGuoZhaoInstantiator(sLattice);

  // choose between local and non-local boundary condition
  sOnLatticeBoundaryCondition2D<T, DESCRIPTOR> sBoundaryCondition(sLattice);
  createInterpBoundaryCondition2D<T, DESCRIPTOR, DYNAMICS<T, DESCRIPTOR> > (sBoundaryCondition);

  prepareLattice(converter, lx, ly, epsilon, K, bodyForceValue, sLattice, bulkDynamics, sBoundaryCondition, sGuoZhaoInstantiator, superGeometry);

  SuperExternal2D<T, DESCRIPTOR, descriptors::FORCE>      externalForce(superGeometry, sLattice, 2);
  SuperExternal2D<T, DESCRIPTOR, descriptors::EPSILON>    externalEpsilon(superGeometry, sLattice, 2);
  SuperExternal2D<T, DESCRIPTOR, descriptors::K>          externalK(superGeometry, sLattice, 2);
  SuperExternal2D<T, DESCRIPTOR, descriptors::NU>         externalNu(superGeometry, sLattice, 2);
  SuperExternal2D<T, DESCRIPTOR, descriptors::BODY_FORCE> externalBodyForce(superGeometry, sLattice, 2);

  /// === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << endl;
  Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( converter.getLatticeTime(maxPhysT/numOfIterations), residuum );
  timer.start();

  for (int iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {
    if ( converge.hasConverged() ) {
       clout << "Simulation converged." << endl;
       getResults(sLattice, bulkDynamics, converter, lx, ly, bodyForceValue, K, converter.getPhysViscosity(), epsilon, maxPhysT, iT, numOfIterations, superGeometry, timer, converge.hasConverged() );
        break;
     }

    /// === 5th Step: Definition of Initial and Boundary Conditions ===
    // in this application no boundary conditions have to be adjusted

    /// === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    externalForce.communicate();
    externalEpsilon.communicate();
    externalK.communicate();
    externalNu.communicate();
    externalBodyForce.communicate();

    /// === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, bulkDynamics, converter, lx, ly, bodyForceValue, K, converter.getPhysViscosity(), epsilon, maxPhysT, iT, numOfIterations, superGeometry, timer, converge.hasConverged() );
  converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
  }

  timer.stop();
  timer.printSummary();
}

