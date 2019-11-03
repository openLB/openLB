/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2018 Sascha Janz, Marie-Luise Maier
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

/* magneticParticles3d.cpp:
 * High-gradient magnetic separation is a method
 * to separate ferromagnetic particles from a suspension.
 * The simulation shows the deposition of magnetic particles
 * on a single magnetized wire and models
 * the magnetic separation step of the complete process.
 */


#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used
#include "olb3D.hh"     // Include full template code
#endif

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

#define DESCRIPTOR D3Q19<>

#define PARTICLE MagneticParticle3D

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef double T;

int iT = 0;

// simulation geometry dimensions
T geometryLengthX = 3.e-3 ; // length x-dircetion [m]
T geometryLengthY = 2.e-3 ; // length y-dircetion [m]
T geometryLengthZ = 1.e-3 ; // length z-dircetion [m]

// magnetic particle paramters
T pRadius = 1.e-5 ;  // particle radius [m]
T pDensity = 3250.; // particle density [kg / m^3]
T pVolume = 4. / 3. * M_PI * std::pow(pRadius, 3.); // particle volume [m^3]
T pMass = pVolume * pDensity; // particle mass [kg]
T pMag = 67. * pDensity * 0.06; // particle sat. magnetization [A / m]
T elastModulus = 2.e4; // elastic modulus [Pa]
T poissonRatio = 0.3; // poisson's ratio [-]
T shearModulus = elastModulus / (2. * (1. + poissonRatio)); // shear modulus [Pa]

// magnetic wire paramters
T wRadius = 5.e-4;  // wire radius [m]
T wLength = geometryLengthZ; // wire length [m]
T wDensity = 7874.; // wire (iron) density [kg / m^3]
T wMag = 215. * wDensity * 8.5e-3; // wire sat. magnetization [A / m]
// mass-specific sat. magnetization iron = 215. [A m^2 / kg]

// Rayleigh time step
// can be used as an orientation value for particle time step sizes
T eta = 0.8766 + 0.1631 * poissonRatio; // [-]
T t_Rayleigh = (M_PI * pRadius / eta) * std::sqrt(pDensity / shearModulus); // [s]

void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converter, IndicatorF3D<T>& wire,
                     IndicatorF3D<T>& indicator, IndicatorF3D<T>& extendedDomain,
                     SuperGeometry3D<T>& superGeometry)
{

  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 2, extendedDomain);
  superGeometry.rename(2, 1, indicator);
  superGeometry.rename(1, 5, wire);

  T minR = superGeometry.getStatistics().getMinPhysR(2)[0];
  minR -= 0.5 * converter.getConversionFactorLength();
  std::vector < T > center = superGeometry.getStatistics().getCenterPhysR(2);
  std::vector<T> physExtend = superGeometry.getStatistics().getPhysExtend(1);

  Vector<T, 3> origin = { minR, center[1], center[2] };
  Vector<T, 3> extend = { converter.getConversionFactorLength(), physExtend[1], physExtend[2] };

  IndicatorCuboid3D<T> inlet(extend[0], extend[1], extend[2], origin);
  origin[0] += superGeometry.getStatistics().getPhysExtend(2)[0];
  IndicatorCuboid3D<T> outlet(extend[0], extend[1], extend[2], origin);

  // rename the material at the inlet
  superGeometry.rename(2, 3, inlet);
  // rename the material at the outlet
  superGeometry.rename(2, 4, outlet);

  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

/// Set up the geometry of the simulation
void prepareLattice(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                    UnitConverter<T, DESCRIPTOR> const& converter, IndicatorF3D<T>& wire,
                    Dynamics<T, DESCRIPTOR>& bulkDynamics,
                    sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& sOnBC,
                    sOffLatticeBoundaryCondition3D<T, DESCRIPTOR>& offBc,
                    SuperGeometry3D<T>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = (1. / converter.getLatticeRelaxationTime());

  /// Material=0 -->do nothing
  sLattice.defineDynamics(superGeometry, 0,
                          &instances::getNoDynamics<T, DESCRIPTOR>());

  /// Material=1 -->bulk dynamics
  sLattice.defineDynamics(superGeometry, 1, &bulkDynamics);

  /// Material=2 -->bounce back
  sLattice.defineDynamics(superGeometry, 2,
                          &instances::getBounceBack<T, DESCRIPTOR>());

  sLattice.defineDynamics(superGeometry, 3, &bulkDynamics);
  sLattice.defineDynamics(superGeometry, 4, &bulkDynamics);

  /// Material=5 -->do nothing
  sLattice.defineDynamics(superGeometry, 5,
                          &instances::getNoDynamics<T, DESCRIPTOR>());
  offBc.addZeroVelocityBoundary(superGeometry, 5, wire);

  // boundary conditions for fluid

  // inlet
  sOnBC.addVelocityBoundary(superGeometry, 3, omega);
  // outlet
  sOnBC.addPressureBoundary(superGeometry, 4, omega);

  // initialisation
  AnalyticalConst3D<T, T> roh(1.);
  AnalyticalConst3D<T, T> u0(0., 0., 0.);

  sLattice.defineRhoU(superGeometry, 1, roh, u0);
  sLattice.iniEquilibrium(superGeometry, 1, roh, u0);
  sLattice.defineRhoU(superGeometry, 3, roh, u0);
  sLattice.iniEquilibrium(superGeometry, 3, roh, u0);
  sLattice.defineRhoU(superGeometry, 4, roh, u0);
  sLattice.iniEquilibrium(superGeometry, 4, roh, u0);

  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                       UnitConverter<T, DESCRIPTOR> const& converter, SuperParticleSystem3D<T, PARTICLE>& spSys,
                       int iT, int itStartScaleT, SuperGeometry3D<T>& superGeometry, bool outNS, bool outLP)
{

  OstreamManager clout(std::cout, "setBoundaryValues");
  std::vector < T > maxVelocity(3, T());
  T distanceToBoundary = converter.getConversionFactorLength() / 2.;

  if (outNS && iT <= itStartScaleT) {

    SinusStartScale<T, int> startScale(itStartScaleT, T(1));
    int help[1] = { iT };
    T frac[3] = { T() };
    startScale(frac, help);

    maxVelocity[0] = converter.getCharPhysVelocity() * frac[0];

    // outlet
    RectanglePoiseuille3D<T> u3(superGeometry, 3, maxVelocity,
                                distanceToBoundary, distanceToBoundary, distanceToBoundary);

    sLattice.defineU(superGeometry, 3, u3);
  }

  // to reset and adopt boundaries in case of loaded fluid data
  if (outLP && iT == 0) {

    maxVelocity[0] = converter.getCharPhysVelocity();

    // inlet
    RectanglePoiseuille3D<T> u3(superGeometry, 3, maxVelocity,
                                distanceToBoundary, distanceToBoundary, distanceToBoundary);

    clout << "BC for loaded fluid data is reset on mat 3 " << std::endl;

    sLattice.defineU(superGeometry, 3, u3);
  }
}

void getResults(SuperGeometry3D<T>& superGeometry,
                SuperLattice3D<T, DESCRIPTOR>& sLattice,
                UnitConverter<T, DESCRIPTOR> const& converter,
                AnalyticalF3D<T, S>& magForce, AnalyticalF3D<T, S>& magField,
                SuperParticleSystem3D<T, PARTICLE>& spSys,
                SuperParticleSysVtuWriterMag<T>& particleOut, Timer<double>& timer, int& iT,
                int& itConsoleOutputFluid, int& itVtkOutputFluid,
                int& itConsoleOutputMagParticles, int& itVtkOutputMagParticles,
                bool& outNS, bool& outLP)
{

  OstreamManager clout(std::cout, "getResults");

  // start vtm objects data until fluid reaches equilibrium
  SuperVTMwriter3D<T> vtmWriterFluidStart("fluidReachingEquilibrium");

  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticeDensity3D<T, DESCRIPTOR> density(sLattice);
  SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> superMagPForceOne(magForce,
      sLattice);
  SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> superMagPFieldOne (magField, sLattice);

  vtmWriterFluidStart.addFunctor(velocity);
  vtmWriterFluidStart.addFunctor(density);
  vtmWriterFluidStart.addFunctor(superMagPForceOne);

  if (outNS && iT == 0) {
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry(sLattice, superGeometry);
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank3D<T, DESCRIPTOR> rank(sLattice);

    vtmWriterFluidStart.write(geometry);
    vtmWriterFluidStart.write(cuboid);
    vtmWriterFluidStart.write(rank);

    vtmWriterFluidStart.createMasterFile();

    converter.write("magnetics 1. loop");

    clout << "fluid computation" << std::endl;
  }
  if (outLP && iT == 0) {

    converter.write("magnetics 2. loop");
    clout << "magParticles computation" << std::endl;
  }

  // === fluid
  // console output fluid
  if (outNS && iT % itConsoleOutputFluid == 0) {
    timer.update(iT);
    timer.print( iT );
    // output for latticeStatistics
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }
  // vtk output fluid
  if (outNS && iT % itVtkOutputFluid == 0) {
    vtmWriterFluidStart.write(iT);
  }

  // === particles
  // Writes the console output files of the magnetic particles
  if (outLP && iT % itConsoleOutputMagParticles == 0) {
    /// Lattice statistics console output
    timer.print( iT );
  }
  // Writes the pvd files of the magnetic particles
  if (outLP && iT % itVtkOutputMagParticles == 0) {
    particleOut.write(iT);
  }
}

int main(int argc, char* argv[])
{

  /// === 1st Step: Initialization ===
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  string fName("magneticParticles3d.xml");
  XMLreader config(fName);

  // converter contains parameters of fluid simulation
  UnitConverter<T, DESCRIPTOR>* converter = createUnitConverter<T, DESCRIPTOR>(config);

  clout << "fluid converter: ..." << std::endl;
  converter->print(); // gives overview of converter into console

  // particles time step size
  T physDeltaTParticles = t_Rayleigh * 0.4; // [s]
  // particles relaxation time
  T tau_particles = (physDeltaTParticles * 3 * converter->getPhysViscosity() /
                     std::pow(converter->getConversionFactorLength(), 2.) + 0.5);

  // converter contains parameters of particle simulation
  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converterParticles(
  int {converter->getResolution()},    // resolution: number of voxels per charPhysL
  (T) tau_particles, // latticeRelaxationTime: relaxation time, has to be greater than 0.5!
  (T) converter->getCharPhysLength(), // (hydraulic raius [m]) charPhysLength: reference length of simulation geometry
  (T) converter->getCharPhysVelocity(), // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
  (T) converter->getPhysViscosity(),  // physViscosity: physical kinematic viscosity in __m^2 / s__
  (T) converter->getPhysDensity()   // physDensity: physical density in __kg / m^3__
  );

  clout << "particle converter: ..." << std::endl;
  converterParticles.print();

  // name of particle data output file
  string filename = "particleData.txt";
  std::ofstream outputFile;

  bool multiOutput = false;
  config["Output"]["MultiOutput"].read(multiOutput);
  clout.setMultiOutput(multiOutput);

  std::string olbdir, outputdir;
  config["Application"]["OlbDir"].read(olbdir);
  config["Output"]["OutputDir"].read(outputdir);
  singleton::directories().setOlbDir(olbdir);
  singleton::directories().setOutputDir(outputdir);

  /// Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 7;
#endif

  // read physical simulation times from xml-file
  T physFluidNST;   // time for fluid simulation
  config["Application"]["SimulationTimes"]["PhysFluidNST"].read(physFluidNST);

  T physStartScT;    // factor * physFluidNST;
  config["Application"]["SimulationTimes"]["PhysStartScT"].read(physStartScT);

  T physParticleT;
  config["Application"]["SimulationTimes"]["PhysParticleT"].read(physParticleT);

  // convert physical times to lattice time steps
  int itFluidNST = converter->getLatticeTime(physFluidNST);
  int itStartScaleT = converter->getLatticeTime(physStartScT);
  int itParticleT = converterParticles.getLatticeTime(physParticleT);

  int itConsoleOutputFluid = itFluidNST / 10. , itVtkOutputFluid = itFluidNST / 4.;
  int itConsoleOutputMagParticles = itParticleT / 75., itVtkOutputMagParticles = itParticleT / 249.;

  clout << "Fluid: physFluidNST = " << physFluidNST << "s itFluidNST = "
        << itFluidNST << std::endl;

  clout << "Particles: physParticleT = " << physParticleT << "s itParticleT = "
        << itParticleT << std::endl;

  clout << "noOfCuboids = " << noOfCuboids << std::endl;

  /// === 2nd Step: Prepare Geometry ===

  // A: Cuboid as Channel
  std::vector < T > center(3, T());
  std::vector<T> length(3, T());
  length[0] = geometryLengthX;
  length[1] = geometryLengthY;
  length[2] = geometryLengthZ;

  IndicatorCuboid3D<T> cuboid(length, center);
  IndicatorLayer3D<T> extendedDomainCuboid(cuboid, converter->getConversionFactorLength());

  std::vector<T> centerW(3, T());
  centerW[0] = geometryLengthX - 0.675e-3;
  centerW[1] = geometryLengthY / 2;
  centerW[2] = geometryLengthZ / 2.;

  std::vector<T> normalW(3, T());
  normalW[0] = T(0);
  normalW[1] = T(0);
  normalW[2] = T(1);
  wLength += 2 * converter->getConversionFactorLength();
  IndicatorCylinder3D<T> wire(centerW, normalW, wRadius, wLength);

  CuboidGeometry3D<T> cuboidGeometry(extendedDomainCuboid,
                                     converter->getConversionFactorLength(), noOfCuboids);

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  /// Instantiation of a superGeometry
  SuperGeometry3D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(*converter, wire, cuboid, extendedDomainCuboid,
                  superGeometry);

  /// === 3rd Step: Prepare Lattice ===

  SuperLattice3D<T, DESCRIPTOR> sLattice(superGeometry);

  BGKdynamics<T, DESCRIPTOR> bulkDynamics((1. / converter->getLatticeRelaxationTime()),
                                          instances::getBulkMomenta<T, DESCRIPTOR>());

  // choose between local and non-local boundary condition
  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sOnBoundaryCondition(sLattice);
  createInterpBoundaryCondition3D<T, DESCRIPTOR>(sOnBoundaryCondition);

  // for the velocity field around the wire
  sOffLatticeBoundaryCondition3D<T, DESCRIPTOR> sOffBoundaryCondition(sLattice);
  createBouzidiBoundaryCondition3D<T, DESCRIPTOR>(sOffBoundaryCondition);

  // gives dynamics to cells
  prepareLattice(sLattice, *converter, wire, bulkDynamics, sOnBoundaryCondition,
                 sOffBoundaryCondition, superGeometry);

  /// === 4th Step: Prepare Lagrange Particles

  SuperParticleSystem3D<T, PARTICLE> spSys(cuboidGeometry, loadBalancer,
      superGeometry);

  // contact detection
  NanoflannContact<T, PARTICLE> nanoflannContact(
    *(spSys.getParticleSystems()[0]), 6. * pRadius);
  spSys.setContactDetection(nanoflannContact);

  // magnetic force field from wire
  std::vector < T > origin = superGeometry.getStatistics().getCenterPhysR(1);
  origin[0] = geometryLengthX - 0.675e-3;
  origin[2] = superGeometry.getStatistics().getMinPhysR(1)[2];

  std::vector<T> orientation(3, T());
  orientation[0] = 1.;
  T degree = 0.;

  // car2cyl helps to know which cylinder coordinates fit to the cartesion ones
  CartesianToCylinder3D<T, T> car2cyl(origin, degree, orientation);

  /// Forces of Lagrange particles

  // Object of magnetic force to know magnetic value at any place in geometry
  // but needs to be called by operator function.
  // inhomogeneous magnetic field round cylinder
  MagneticForceFromCylinder3D<T, T> magForceFromCylOnParticle(car2cyl, wLength,
      wRadius, 4., wMag, pMag, pRadius);
  MagneticFieldFromCylinder3D<T, T> magFieldFromCylOnParticle(car2cyl, wLength,
      wRadius, 4., wMag);

  // external magnetic force from the wire
  auto magneticPForceOne = make_shared
                           < MagneticForceForMagP3D<T, PARTICLE, DESCRIPTOR>
                           > (magForceFromCylOnParticle, magFieldFromCylOnParticle, 1.e0);
  spSys.addForce(magneticPForceOne);

  // interparticular magnetic force
  auto interpMagF = make_shared
                    < InterpMagForceForMagP3D<T, PARTICLE, DESCRIPTOR>
                    > (1.e0, 1.e0);
  spSys.addForce(interpMagF);

  // dynamic viscosity
  T dynVisc = converter->getPhysViscosity() * converter->getPhysDensity();

  // damping force for the particle rotation
  auto rotDampingF = make_shared
                     < LinearDampingForceForMagDipoleMoment3D<T, PARTICLE, DESCRIPTOR>
                     > (dynVisc, 1.e0);
  spSys.addForce(rotDampingF);

  // mechanic contact force
  auto hertzMindlinContF = make_shared
                           < HertzMindlinDeresiewicz3D<T, PARTICLE, DESCRIPTOR>
                           > (shearModulus, shearModulus, poissonRatio, poissonRatio, 1.e0, 1.e0, false);
  spSys.addForce(hertzMindlinContF);

  // stokes drag force
  SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR> getVel( sLattice, *converter );
  auto stokesDragF = make_shared
                     < StokesDragForce3D<T, PARTICLE, DESCRIPTOR>
                     > (getVel, converterParticles);
  spSys.addForce(stokesDragF);

  /// Boundarys

  T dT = converter->getConversionFactorTime() ;
  std::set<int> wireBMaterial = { 5, 4 };
  std::set<int> reflBMat = { 2, 3 };

  // particle boundary for deposited particles
  auto wireBoundary = make_shared
                      < WireBoundaryForMagP3D<T, PARTICLE>
                      > (superGeometry, wireBMaterial);
  spSys.addBoundary(wireBoundary);

  // particle boundary
  auto materialreflectBoundary = make_shared
                                 < SimpleReflectBoundary3D<T, PARTICLE>
                                 > (dT, superGeometry, reflBMat);
  spSys.addBoundary(materialreflectBoundary);

  // cuboid overlap
  spSys.setOverlap(2.*converter->getConversionFactorLength());

  // particles vtu output
  std::string particleOutputName = "vtuOutMagP";
  std::bitset < 9 > properties;
  properties.set();

  SuperParticleSysVtuWriterMag<T> particleOut(spSys, particleOutputName, properties);

  clout << "Number of Forces: " << spSys.numOfForces()[0] << std::endl;

  /// Paramters for magnetic particle generation
  std::vector<T> pPos = { 0., 0., 0. }; // position
  std::vector<T> pVel = { 0., 0., 0. }; // velocity
  std::vector<T> pAVel = { 0., 0., 0. }; // angular velocity
  std::vector<T> pTrq = { 0., 0., 0. }; // torque
  std::vector<T> pDMoment = { -1., 0., 0. }; // orientation mag. dipole moment
  if(!util::nearZero(util::norm(pDMoment))){util::normalize(pDMoment);}
  else {clout << "Norm of pDMoment near zero!" << endl; exit(0);}

  // magnetic particle magPartTemplate is used as copy template
  PARTICLE<T> magPartTemplate(pPos, pVel, pMass, pRadius, 0, pDMoment, pAVel, pTrq, pMag, 1);

  // intialization geometry particles
  std::vector<T> particlesOrigin(3, T(0));
  particlesOrigin = superGeometry.getStatistics().getMinPhysR(1);
  std::vector<T> particlesExtend(3, T(0));

  particlesOrigin[0] = 1.e-4;
  particlesOrigin[1] = (superGeometry.getStatistics().getMaxPhysR(1)[1]) * (1./2. - 1./6.);
  particlesOrigin[2] = (superGeometry.getStatistics().getMaxPhysR(1)[2]) * (1./2. - 1./6.);
  particlesExtend[0] = (superGeometry.getStatistics().getMaxPhysR(1)[0]) * 1./15. ;
  particlesExtend[1] = (superGeometry.getStatistics().getMaxPhysR(1)[1]) * 2./6. ;
  particlesExtend[2] = (superGeometry.getStatistics().getMaxPhysR(1)[2]) * 2./6. ;
  IndicatorCuboid3D<T> cuboidPart(particlesExtend, particlesOrigin);

  clout << "starting simulation" << endl;

  /// === 5th Step: Main Loop with Timer ===

  // is set on true for the appropriate loop of fluid / particle simulation
  bool outNS = false;
  bool outLP = false;

  // === Calculating Navier-Stokes Fluid

  // if there is no computed fluid data in a external data files
  // compute fluid data in data file, else load it and jump to next loop
  if (!(sLattice.load("water"))) {

    clout << std::endl;
    clout << "Fluid data has to be computed " << std::endl;
    clout << std::endl;

    outNS = true;

    Timer<double> timerFluid(itFluidNST, superGeometry.getStatistics().getNvoxel());
    timerFluid.start();

    for (int iT = 0; iT < itFluidNST; ++iT) {

      setBoundaryValues(sLattice, *converter, spSys, iT, itStartScaleT,
                        superGeometry, outNS, outLP);

      sLattice.collideAndStream();

      getResults(superGeometry, sLattice, *converter,
                 magForceFromCylOnParticle, magFieldFromCylOnParticle, spSys,
                 particleOut, timerFluid, iT, itConsoleOutputFluid, itVtkOutputFluid,
                 itConsoleOutputMagParticles, itVtkOutputMagParticles, outNS, outLP);
    }

    outNS = false;

    // save computed fluid data in external data files
    clout << "save fluid data " << std::endl;

    sLattice.communicate();
    sLattice.save("water");

    timerFluid.stop();
    timerFluid.printSummary();

  } else {
    sLattice.communicate();
    sLattice.collideAndStream();
    sLattice.getStatistics().print(iT, converter->getPhysTime(iT));
  }

  // === Calculating Particles

  clout << std::endl;
  clout << "particles computation" << std::endl;
  clout << std::endl;

  outLP = true;

  // number of Particles in initialization geometry
  int noOfPartX = 4 ;
  int noOfPartY = 12 ;
  int noOfPartZ = 6 ;
  int noOfPart = noOfPartX * noOfPartY * noOfPartZ;

  Timer<double> timerParticles(itParticleT, noOfPart);
  timerParticles.start();

  setBoundaryValues(sLattice, converterParticles, spSys, iT, itStartScaleT,
                    superGeometry, outNS, outLP);

  for (int iT = 0; iT < itParticleT; ++iT) {

    // particles initialization
    if (iT == 0) {

      spSys.addParticleEquallyDistributed(cuboidPart, noOfPartX, noOfPartY, noOfPartZ, magPartTemplate);
      spSys.setParticlesPosRandom(0.75 * pRadius, 0.75 * pRadius, 0.55 * pRadius);
    }

    getResults(superGeometry, sLattice, converterParticles,
               magForceFromCylOnParticle, magFieldFromCylOnParticle, spSys,
               particleOut, timerParticles, iT, itConsoleOutputFluid, itVtkOutputFluid,
               itConsoleOutputMagParticles, itVtkOutputMagParticles,
               outNS, outLP);

    // particles simulation
    spSys.simulate(converterParticles.getConversionFactorTime());

  }

  timerParticles.stop();
  timerParticles.printSummary();

  clout << "End of simulation" << std::endl;
}
