/* dkt2d.cpp:
 * The case examines the settling of two circles under gravity
 * in a surrounding fluid. The rectangular domain is limited
 * by no-slip boundary conditions.
 * For the calculation of forces a DNS approach is chosen
 * which also leads to a back-coupling of the particle on the fluid,
 * inducing a flow.
 * The simulation is based on the homogenised lattice Boltzmann approach
 * (HLBM) introduced by Krause et al. in "Particle flow simulations
 * with homogenised lattice Boltzmann methods".
 * The drafting-kissing-tumbling benchmark case is e.g. described
 * in "Drafting, kissing and tumbling process of two particles
 * with different sizes" by Wang et al.
 * or "The immersed boundary-lattice Boltzmann method
 * for solving fluid-particles interaction problems" by Feng and Michaelides.
 * The example demonstrates the usage of HLBM in the OpenLB framework
 * as well as the utilisation of the Gnuplot-writer
 * to print simulation results.
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
#define DESCRIPTOR D2Q9<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR>

#define WriteVTK
#define WriteGnuPlot

std::string gnuplotFilename = "gnuplot.dat";

// Parameters for the simulation setup
int N = 1;
int M = N;

T eps = 0.5;      // eps*latticeL: width of transition area

T maxPhysT = 6.;  // max. simulation time in s, SI unit
T iTwrite = 0.125;  //converter.getLatticeTime(.3);

T lengthX = 0.02;
T lengthY = 0.08;

T centerX1 = 0.01;
T centerY1 = 0.068;
Vector<T,2> center1 = {centerX1,centerY1};
T centerX2 = 0.00999;
T centerY2 = 0.072;
Vector<T,2> center2 = {centerX2,centerY2};

T rhoP = 1010.;
T radiusP = 0.001;
Vector<T,2> accExt = {.0, -9.81 * (1. - 1000. / rhoP)};

void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry2D<T>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 2);
  superGeometry.rename(2, 1, 1, 1);

  superGeometry.clean();
  superGeometry.innerClean();

  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(
  SuperLattice2D<T, DESCRIPTOR>& sLattice, UnitConverter<T,DESCRIPTOR> const& converter,
  Dynamics<T, DESCRIPTOR>& designDynamics,
  sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>& sBoundaryCondition,
  SuperGeometry2D<T>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  /// Material=0 -->do nothing
  sLattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>());
  sLattice.defineDynamics(superGeometry, 1, &designDynamics);
  sLattice.defineDynamics(superGeometry, 2, &instances::getBounceBack<T, DESCRIPTOR>());

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(SuperLattice2D<T, DESCRIPTOR>& sLattice,
                       UnitConverter<T,DESCRIPTOR> const& converter,
                       SuperGeometry2D<T>& superGeometry)
{
  OstreamManager clout(std::cout, "setBoundaryValues");

  AnalyticalConst2D<T, T> one(1.);
  sLattice.defineField<POROSITY>(superGeometry.getMaterialIndicator({1,2}), one);
  
  // Set initial condition
  AnalyticalConst2D<T, T> ux(0.);
  AnalyticalConst2D<T, T> uy(0.);
  AnalyticalConst2D<T, T> rho(1.);
  AnalyticalComposed2D<T, T> u(ux, uy);

  //Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU(superGeometry, 1, rho, u);
  sLattice.iniEquilibrium(superGeometry, 1, rho, u);

  // Make the lattice ready for simulation
  sLattice.initialize();
}

void getResults(SuperLattice2D<T, DESCRIPTOR>& sLattice,
                UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                SuperGeometry2D<T>& superGeometry, Timer<double>& timer, SmoothIndicatorF2D<T,T,true> &particle1, SmoothIndicatorF2D<T,T,true> &particle2)
{
  OstreamManager clout(std::cout, "getResults");

#ifdef WriteVTK
  SuperVTMwriter2D<T> vtkWriter("sedimentation");
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(sLattice, converter);
  SuperLatticePhysExternalPorosity2D<T, DESCRIPTOR> externalPor(sLattice, converter);
  vtkWriter.addFunctor(velocity);
  vtkWriter.addFunctor(pressure);
  vtkWriter.addFunctor(externalPor);

  if (iT == 0) {
    converter.write("dkt");
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry(sLattice, superGeometry);
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank2D<T, DESCRIPTOR> rank(sLattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
  }

  if (iT % converter.getLatticeTime(iTwrite) == 0) {
    vtkWriter.write(iT);
  }
#endif

#ifdef WriteGnuPlot
  if (iT % converter.getLatticeTime(iTwrite) == 0) {
    if (singleton::mpi().getRank() == 0) {

      ofstream myfile;
      myfile.open (gnuplotFilename.c_str(), ios::app);
      myfile
          << converter.getPhysTime(iT) << " "
          << std::setprecision(9)
          << particle2.getPos()[1] << " "
          << particle1.getPos()[1] << " "
          << particle2.getPos()[0] << " "
          << particle1.getPos()[0] << endl;
      myfile.close();
    }
  }
#endif

  /// Writes output on the console
  if (iT % converter.getLatticeTime(iTwrite) == 0) {
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }

  return;
}

int main(int argc, char* argv[])
{
  /// === 1st Step: Initialization ===
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  UnitConverter<T,DESCRIPTOR> converter(
    ( T )   0.0001/ N, //physDeltaX
    ( T )   5.e-4/(N*M), //physDeltaT,
    ( T )   .002, //charPhysLength
    ( T )   0.2, //charPhysVelocity
    ( T )   1E-6, //physViscosity
    ( T )   1000. //physDensity
  );
  converter.print();

  /// === 2nd Step: Prepare Geometry ===
  std::vector<T> extend(2, T());
  extend[0] = lengthX;
  extend[1] = lengthY;
  std::vector<T> origin(2, T());
  IndicatorCuboid2D<T> cuboid(extend, origin);

#ifdef PARALLEL_MODE_MPI
  CuboidGeometry2D<T> cuboidGeometry(cuboid, converter.getConversionFactorLength(), singleton::mpi().getSize());
#else
  CuboidGeometry2D<T> cuboidGeometry(cuboid, converter.getConversionFactorLength(), 1);
#endif

  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);
  SuperGeometry2D<T> superGeometry(cuboidGeometry, loadBalancer, 2);
  prepareGeometry(converter, superGeometry);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice2D<T, DESCRIPTOR> sLattice(superGeometry);
  PorousParticleBGKdynamics<T, DESCRIPTOR> designDynamics(converter.getLatticeRelaxationFrequency(), instances::getBulkMomenta<T, DESCRIPTOR>());

  sOnLatticeBoundaryCondition2D<T, DESCRIPTOR> sBoundaryCondition(sLattice);
  createLocalBoundaryCondition2D<T, DESCRIPTOR>(sBoundaryCondition);

  prepareLattice(sLattice, converter, designDynamics, sBoundaryCondition, superGeometry);

  /// === 4th Step: Main Loop with Timer ===
  Timer<double> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel());
  timer.start();

  ParticleDynamics2D<T, DESCRIPTOR> particle(sLattice, converter, superGeometry, lengthX, lengthY, accExt);
  SmoothIndicatorCircle2D<T,T,true> circle2(center1, radiusP, eps*converter.getConversionFactorLength(), rhoP);
  SmoothIndicatorCircle2D<T,T,true> circle1(center2, radiusP, eps*converter.getConversionFactorLength(), rhoP);
  particle.addParticle(circle2);
  particle.addParticle(circle1);

  SuperExternal2D<T,DESCRIPTOR,POROSITY> superExt1(superGeometry, sLattice, sLattice.getOverlap());
  SuperExternal2D<T,DESCRIPTOR,VELOCITY_NUMERATOR> superExt2(superGeometry, sLattice, sLattice.getOverlap());
  SuperExternal2D<T,DESCRIPTOR,VELOCITY_DENOMINATOR> superExt3(superGeometry, sLattice, sLattice.getOverlap());

  /// === 5th Step: Definition of Initial and Boundary Conditions ===
  setBoundaryValues(sLattice, converter, superGeometry);

  clout << "MaxIT: " << converter.getLatticeTime(maxPhysT) << std::endl;
  for (int iT = 0; iT < converter.getLatticeTime(maxPhysT)+10; ++iT) {
    particle.simulateTimestep("verlet");
    getResults(sLattice, converter, iT, superGeometry, timer, circle1, circle2);
    sLattice.collideAndStream();
    superExt1.communicate();
    superExt2.communicate();
    superExt3.communicate();

  }

  // Run Gnuplot
  if (singleton::mpi().getRank() == 0) {
    if (!system(NULL)) {
      exit (EXIT_FAILURE);
    }
    int ret = system("gnuplot dkt.p");
    if (ret == -1) {
      clout << "Writing Gnuplot failed!" << endl;
    }
  }

  timer.stop();
  timer.printSummary();
}
