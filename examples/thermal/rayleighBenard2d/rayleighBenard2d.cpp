/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2008 Orestis Malaspinas
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

/* rayleighBenard2d.cpp:
 * Rayleigh-Benard convection rolls in 2D, simulated with
 * the thermal LB model by Z. Guo e.a., between a hot plate at
 * the bottom and a cold plate at the top.
 */


#include "olb2D.h"
#include "olb2D.hh"   // use only generic version!

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

#define NSDESCRIPTOR D2Q9<FORCE>
#define TDESCRIPTOR D2Q5<VELOCITY>

// Parameters for the simulation setup
const T lx  = 2.;          // length of the channel
const T ly  = 1.;          // height of the channel
const int N = 10;         // resolution of the model
const T Ra = 1e4;          // Rayleigh number
const T Pr = 0.71;         // Prandtl number
const T maxPhysT = 1000.; // max. simulation time in s, SI unit
const T epsilon = 1.e-5;   // precision of the convergence (residuum)

const T Thot = 274.15;     // temperature of the lower wall in Kelvin
const T Tcold = 273.15;    // temperature of the fluid in Kelvin
const T Tperturb = 1./5. * Tcold + 4./5. * Thot; // temperature of the perturbation

/// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry2D<T>& superGeometry,
                     ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter)
{

  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,2);
  superGeometry.rename(2,1,0,1);

  std::vector<T> extend( 2, T(0) );
  extend[0] = lx;
  extend[1] = converter.getPhysLength(1);
  std::vector<T> origin( 2, T(0) );
  IndicatorCuboid2D<T> bottom(extend, origin);

  origin[1] = ly-converter.getPhysLength(1);
  IndicatorCuboid2D<T> top(extend, origin);

  origin[0] = lx/2.;
  origin[1] = converter.getPhysLength(1);
  extend[0] = converter.getPhysLength(1);
  extend[1] = converter.getPhysLength(1);
  IndicatorCuboid2D<T> perturbation(extend, origin);

  /// Set material numbers for bottom, top and pertubation
  superGeometry.rename(2,2,1,bottom);
  superGeometry.rename(2,3,1,top);
  superGeometry.rename(1,4,perturbation);

  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                     SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                     SuperLattice2D<T, TDESCRIPTOR>& ADlattice,
                     Dynamics<T, NSDESCRIPTOR> &bulkDynamics,
                     Dynamics<T, TDESCRIPTOR>& advectionDiffusionBulkDynamics,
                     sOnLatticeBoundaryCondition2D<T,NSDESCRIPTOR>& NSboundaryCondition,
                     sOnLatticeBoundaryCondition2D<T,TDESCRIPTOR>& TboundaryCondition,
                     SuperGeometry2D<T>& superGeometry )
{

  OstreamManager clout(std::cout,"prepareLattice");

  T Tomega  = converter.getLatticeThermalRelaxationFrequency();

  /// define lattice Dynamics
  clout << "defining dynamics" << endl;

  ADlattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, TDESCRIPTOR>());
  NSlattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, NSDESCRIPTOR>());

  ADlattice.defineDynamics(superGeometry, 1, &advectionDiffusionBulkDynamics);
  ADlattice.defineDynamics(superGeometry, 2, &advectionDiffusionBulkDynamics);
  ADlattice.defineDynamics(superGeometry, 3, &advectionDiffusionBulkDynamics);
  ADlattice.defineDynamics(superGeometry, 4, &advectionDiffusionBulkDynamics);
  NSlattice.defineDynamics(superGeometry, 1, &bulkDynamics);
  NSlattice.defineDynamics(superGeometry, 2, &instances::getBounceBack<T, NSDESCRIPTOR>());
  NSlattice.defineDynamics(superGeometry, 3, &instances::getBounceBack<T, NSDESCRIPTOR>());
  NSlattice.defineDynamics(superGeometry, 4, &bulkDynamics);

  /// sets boundary
  TboundaryCondition.addTemperatureBoundary(superGeometry, 2, Tomega);
  TboundaryCondition.addTemperatureBoundary(superGeometry, 3, Tomega);

  /// define initial conditions
  AnalyticalConst2D<T,T> rho(1.);
  AnalyticalConst2D<T,T> u0(0.0, 0.0);
  AnalyticalConst2D<T,T> T_cold(converter.getLatticeTemperature(Tcold));
  AnalyticalConst2D<T,T> T_hot(converter.getLatticeTemperature(Thot));
  AnalyticalConst2D<T,T> T_perturb(converter.getLatticeTemperature(Tperturb));

  /// for each material set Rho, U and the Equilibrium
  NSlattice.defineRhoU(superGeometry, 1, rho, u0);
  NSlattice.iniEquilibrium(superGeometry, 1, rho, u0);
  NSlattice.defineRhoU(superGeometry, 2, rho, u0);
  NSlattice.iniEquilibrium(superGeometry, 2, rho, u0);
  NSlattice.defineRhoU(superGeometry, 3, rho, u0);
  NSlattice.iniEquilibrium(superGeometry, 3, rho, u0);
  NSlattice.defineRhoU(superGeometry, 4, rho, u0);
  NSlattice.iniEquilibrium(superGeometry, 4, rho, u0);

  ADlattice.defineRho(superGeometry, 1, T_cold);
  ADlattice.iniEquilibrium(superGeometry, 1, T_cold, u0);
  ADlattice.defineRho(superGeometry, 2, T_hot);
  ADlattice.iniEquilibrium(superGeometry, 2, T_hot, u0);
  ADlattice.defineRho(superGeometry, 3, T_cold);
  ADlattice.iniEquilibrium(superGeometry, 3, T_cold, u0);
  ADlattice.defineRho(superGeometry, 4, T_perturb);
  ADlattice.iniEquilibrium(superGeometry, 4, T_perturb, u0);

  /// Make the lattice ready for simulation
  NSlattice.initialize();
  ADlattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                       SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                       SuperLattice2D<T, TDESCRIPTOR>& ADlattice,
                       int iT, SuperGeometry2D<T>& superGeometry)
{
  // nothing to do here
}

void getResults(ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                SuperLattice2D<T, TDESCRIPTOR>& ADlattice, int iT,
                SuperGeometry2D<T>& superGeometry,
                Timer<T>& timer,
                bool converged)
{

  OstreamManager clout(std::cout,"getResults");

  SuperVTMwriter2D<T> vtkWriter("rayleighBenard2d");
  SuperLatticePhysVelocity2D<T, NSDESCRIPTOR> velocity(NSlattice, converter);
  SuperLatticePhysPressure2D<T, NSDESCRIPTOR> presure(NSlattice, converter);
  SuperLatticePhysTemperature2D<T, NSDESCRIPTOR, TDESCRIPTOR> temperature(ADlattice, converter);
  vtkWriter.addFunctor( presure );
  vtkWriter.addFunctor( velocity );
  vtkWriter.addFunctor( temperature );

  const int saveIter = converter.getLatticeTime(10.0);

  if (iT == 0) {
    /// Writes the converter log file
    // writeLogFile(converter,"rayleighBenard2d");

    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, NSDESCRIPTOR> geometry(NSlattice, superGeometry);
    SuperLatticeCuboid2D<T, NSDESCRIPTOR> cuboid(NSlattice);
    SuperLatticeRank2D<T, NSDESCRIPTOR> rank(NSlattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  /// Writes the VTK files and prints statistics
  if (iT%saveIter == 0 || converged) {
    /// Timer console output
    timer.update(iT);
    timer.printStep();

    /// Lattice statistics console output
    NSlattice.getStatistics().print(iT,converter.getPhysTime(iT));

    vtkWriter.write(iT);

    BlockReduction2D2D<T> planeReduction(temperature, 600, BlockDataSyncMode::ReduceOnly);
    BlockGifWriter<T> gifWriter;
    gifWriter.write(planeReduction, Tcold-0.1, Thot+0.1, iT, "temperature");
  }

}

int main(int argc, char *argv[])
{

  /// === 1st Step: Initialization ===
  OstreamManager clout(std::cout,"main");
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> converter(
    (T) 0.1/N, // physDeltaX
    (T) 0.1 / (1e-5 / 0.1 * sqrt( Ra / Pr)) * 0.1 / N, // physDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX
    (T) 0.1,  // charPhysLength
    (T) 1e-5 / 0.1 * sqrt( Ra / Pr ), // charPhysVelocity
    (T) 1e-5,  // physViscosity
    (T) 1.0, // physDensity
    (T) 0.03, // physThermalConductivity
    (T) Pr * 0.03 / 1e-5 / 1.0,    // physSpecificHeatCapacity
    (T) Ra * 1e-5 * 1e-5 / Pr / 9.81 / (Thot - Tcold) / pow(0.1, 3), // physThermalExpansionCoefficient
    (T) Tcold, // charPhysLowTemperature
    (T) Thot // charPhysHighTemperature
  );
  converter.print();

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

  cuboidGeometry.setPeriodicity(true, false);

  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  SuperGeometry2D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(superGeometry, converter);

  /// === 3rd Step: Prepare Lattice ===

  SuperLattice2D<T, TDESCRIPTOR> ADlattice(superGeometry);
  SuperLattice2D<T, NSDESCRIPTOR> NSlattice(superGeometry);

  sOnLatticeBoundaryCondition2D<T, NSDESCRIPTOR> NSboundaryCondition(NSlattice);
  createLocalBoundaryCondition2D<T, NSDESCRIPTOR>(NSboundaryCondition);

  sOnLatticeBoundaryCondition2D<T, TDESCRIPTOR> TboundaryCondition(ADlattice);
  createAdvectionDiffusionBoundaryCondition2D<T, TDESCRIPTOR>(TboundaryCondition);

  ForcedBGKdynamics<T, NSDESCRIPTOR> NSbulkDynamics(
    converter.getLatticeRelaxationFrequency(),
    instances::getBulkMomenta<T,NSDESCRIPTOR>());

  AdvectionDiffusionBGKdynamics<T, TDESCRIPTOR> TbulkDynamics (
    converter.getLatticeThermalRelaxationFrequency(),
    instances::getAdvectionDiffusionBulkMomenta<T,TDESCRIPTOR>()
  );

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
  // This coupling must be necessarily be put on the Navier-Stokes lattice!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

  std::vector<T> dir{0.0, 1.0};

  T boussinesqForcePrefactor = 9.81 / converter.getConversionFactorVelocity() * converter.getConversionFactorTime() *
    converter.getCharPhysTemperatureDifference() * converter.getPhysThermalExpansionCoefficient();

  NavierStokesAdvectionDiffusionCouplingGenerator2D<T,NSDESCRIPTOR> coupling(0, converter.getLatticeLength(lx), 0, converter.getLatticeLength(ly), boussinesqForcePrefactor, converter.getLatticeTemperature(Tcold), 1., dir);

  NSlattice.addLatticeCoupling(coupling, ADlattice);
  NSlattice.addLatticeCoupling(superGeometry, 2, coupling, ADlattice);
  NSlattice.addLatticeCoupling(superGeometry, 3, coupling, ADlattice);
  NSlattice.addLatticeCoupling(superGeometry, 4, coupling, ADlattice);

  prepareLattice(converter,
                 NSlattice, ADlattice,
                 NSbulkDynamics, TbulkDynamics,
                 NSboundaryCondition, TboundaryCondition, superGeometry );

  /// === 4th Step: Main Loop with Timer ===
  Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  util::ValueTracer<T> converge(converter.getLatticeTime(50.),epsilon);
  for (int iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {

    if (converge.hasConverged()) {
      clout << "Simulation converged." << endl;
      getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, converge.hasConverged());

      clout << "Time " << iT << "." << std::endl;

      break;
    }

    /// === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(converter, NSlattice, ADlattice, iT, superGeometry);

    /// === 6th Step: Collide and Stream Execution ===
    ADlattice.collideAndStream();
    NSlattice.collideAndStream();

    NSlattice.executeCoupling();

    /// === 7th Step: Computation and Output of the Results ===
    getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, converge.hasConverged());
    converge.takeValue(ADlattice.getStatistics().getAverageEnergy(),true);
  }

  timer.stop();
  timer.printSummary();
}
