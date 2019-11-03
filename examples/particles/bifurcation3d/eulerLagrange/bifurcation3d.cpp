/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2016 Thomas Henn, Mathias J. Krause,
 *  Marie-Luise Maier
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

/* bifurcation3d.cpp:
 * This example examines a steady particulate flow past a bifurcation. At the inlet,
 * an inflow condition with grad_n u = 0 and rho = 1 is implemented.
 * At both outlets, a Poiseuille profile is imposed on the velocity.
 * After a start time, particles are put into the bifurcation at the
 * inlet and experience a stokes drag force.
 *
 * A publication using the same geometry can be found here:
 * http://link.springer.com/chapter/10.1007/978-3-642-36961-2_5
 *  *
 */

#include "olb3D.h"
#ifndef OLB_PRECOMPILED // Unless precompiled version is used,
#include "olb3D.hh"   // include full template code;
#endif

using namespace std;
using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;

typedef double T;
#define DESCRIPTOR D3Q19<>
#define PARTICLE Particle3D

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const T Re = 50;                   // Reynolds number
const int N = 19;                   // resolution of the model
const T radius = 1.5e-4;           // particles radius
const T partRho = 998.2;           //particles density

const T fluidMaxPhysT = T( 5 );    // max. fluid simulation time in s, SI unit
const T particleMaxPhysT = T( 10 ); // max. particle simulation time in s, SI unit

const int noOfParticles = 1000;    // total number of inserted particles

// center of inflow and outflow regions [m]
Vector<T, 3> inletCenter( T(), T(), 0.0786395 );
Vector<T, 3> outletCenter0( -0.0235929682287551, -0.000052820468762797,
                            -0.021445708949909 );
Vector<T, 3> outletCenter1( 0.0233643529416147, 0.00000212439067050152,
                            -0.0211994104877918 );

// radii of inflow and outflow regions [m]
T inletRadius = 0.00999839;
T outletRadius0 = 0.007927;
T outletRadius1 = 0.00787134;

// normals of inflow and outflow regions
Vector<T, 3> inletNormal( T(), T(), T( -1 ) );
Vector<T, 3> outletNormal0( 0.505126, -0.04177, 0.862034 );
Vector<T, 3> outletNormal1( -0.483331, -0.0102764, 0.875377 );

void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      IndicatorF3D<T>& indicator, STLreader<T>& stlReader,
                      SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout, "prepareGeometry" );

  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0, 2, indicator );
  superGeometry.rename( 2, 1, stlReader );

  superGeometry.clean();

  // rename the material at the inlet
  IndicatorCircle3D<T> inletCircle( inletCenter, inletNormal,
                                    inletRadius );
  IndicatorCylinder3D<T> inlet( inletCircle,
                                2 * converter.getConversionFactorLength() );
  superGeometry.rename( 2, 3, 1, inlet );

  // rename the material at the outlet0
  IndicatorCircle3D<T> outletCircle0( outletCenter0, outletNormal0,
                                      0.95 * outletRadius0 );
  IndicatorCylinder3D<T> outlet0( outletCircle0,
                                  4 * converter.getConversionFactorLength() );
  superGeometry.rename( 2, 4, outlet0 );

  // rename the material at the outlet1
  IndicatorCircle3D<T> outletCircle1( outletCenter1, outletNormal1,
                                      0.95 * outletRadius1 );
  IndicatorCylinder3D<T> outlet1( outletCircle1,
                                  4 * converter.getConversionFactorLength() );
  superGeometry.rename( 2, 5, outlet1 );

  superGeometry.clean();
  superGeometry.innerClean( true );
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}

void prepareLattice( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter, Dynamics<T, DESCRIPTOR>&
                     bulkDynamics,
                     sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& bc,
                     sOffLatticeBoundaryCondition3D<T, DESCRIPTOR>& offBc,
                     STLreader<T>& stlReader, SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout, "prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLattice.defineDynamics( superGeometry, 0,
                           &instances::getNoDynamics<T, DESCRIPTOR>() );

  // Material=1 -->bulk dynamics
  sLattice.defineDynamics( superGeometry, 1, &bulkDynamics );

  // Material=2 -->bounce back
  sLattice.defineDynamics( superGeometry, 2,
                           &instances::getBounceBack<T, DESCRIPTOR>() );

  // Material=3 -->bulk dynamics (inflow)
  sLattice.defineDynamics( superGeometry, 3, &bulkDynamics );

  // Material=4 -->bulk dynamics (outflow)
  sLattice.defineDynamics( superGeometry, 4, &bulkDynamics );
  sLattice.defineDynamics( superGeometry, 5, &bulkDynamics );

  // Setting of the boundary conditions
  bc.addPressureBoundary( superGeometry, 3, omega );
  bc.addVelocityBoundary( superGeometry, 4, omega );
  bc.addVelocityBoundary( superGeometry, 5, omega );

  clout << "Prepare Lattice ... OK" << std::endl;
  return;
}

// Generates a slowly increasing sinuidal inflow for the first iTMax timesteps
void setBoundaryValues( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T,DESCRIPTOR> const& converter, int iT, T maxPhysT,
                        SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout, "setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime( 0.8*maxPhysT );
  int iTperiod = 100;  // amount of timesteps when new boundary conditions are reset

  if ( iT == 0 ) {

    AnalyticalConst3D<T, T> rhoF( 1 );
    std::vector<T> velocity( 3, T() );
    AnalyticalConst3D<T, T> uF( velocity );

    sLattice.iniEquilibrium( superGeometry, 1, rhoF, uF );
    sLattice.iniEquilibrium( superGeometry, 2, rhoF, uF );
    sLattice.iniEquilibrium( superGeometry, 3, rhoF, uF );
    sLattice.iniEquilibrium( superGeometry, 4, rhoF, uF );
    sLattice.iniEquilibrium( superGeometry, 5, rhoF, uF );

    sLattice.defineRhoU( superGeometry, 1, rhoF, uF );
    sLattice.defineRhoU( superGeometry, 2, rhoF, uF );
    sLattice.defineRhoU( superGeometry, 3, rhoF, uF );
    sLattice.defineRhoU( superGeometry, 4, rhoF, uF );
    sLattice.defineRhoU( superGeometry, 5, rhoF, uF );

    // Make the lattice ready for simulation
    sLattice.initialize();
  }

  else if ( iT <= iTmaxStart && iT % iTperiod == 0 ) {
    SinusStartScale<T, int> startScale( iTmaxStart, T( 1 ) );
    int iTvec[1] = { iT };
    T frac[1] = { T( 0 ) };
    startScale( frac, iTvec );
    T maxVelocity = frac[0] * converter.getCharLatticeVelocity() * 3. / 4.
                    * std::pow( inletRadius, 2 ) / std::pow( outletRadius0, 2 );

    CirclePoiseuille3D<T> poiseuilleU4( outletCenter0[0], outletCenter0[1],
                                        outletCenter0[2], outletNormal0[0],
                                        outletNormal0[1], outletNormal0[2],
                                        outletRadius0 * 0.95, -maxVelocity );

    CirclePoiseuille3D<T> poiseuilleU5( outletCenter1[0], outletCenter1[1],
                                        outletCenter1[2], outletNormal1[0],
                                        outletNormal1[1], outletNormal1[2],
                                        outletRadius1 * 0.95, -maxVelocity );

    sLattice.defineU( superGeometry, 4, poiseuilleU4 );
    sLattice.defineU( superGeometry, 5, poiseuilleU5 );
  }
}

// Computes the pressure drop between voxels before and after the cylinder
bool getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT, int iTperiod,
                 SuperGeometry3D<T>& superGeometry,
                 Timer<double>& fluidTimer, STLreader<T>& stlReader,
                 SuperParticleSystem3D<T, PARTICLE>& supParticleSystem,
                 T radii, T partRho, Timer<double>& particleTimer,
                 SuperParticleSysVtuWriter<T, PARTICLE>& supParticleWriter,
                 bool fluidExists)
{

  OstreamManager clout( std::cout, "getResults" );
  SuperVTMwriter3D<T> vtmWriter( "bifurcation3d" );
  SuperVTMwriter3D<T> vtmWriterStartTime( "startingTimeBifurcation3d" );

  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  vtmWriterStartTime.addFunctor( velocity );
  vtmWriterStartTime.addFunctor( pressure );

  int fluidMaxT = converter.getLatticeTime( fluidMaxPhysT );

  if ( iT == 0 ) {
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
    vtmWriterStartTime.createMasterFile();


    // Print some output of the chosen simulation setup
    clout << "N=" << N <<"; maxTimeSteps(fluid)="
          << converter.getLatticeTime( fluidMaxPhysT ) << "; noOfCuboid="
          << superGeometry.getCuboidGeometry().getNc() << "; Re=" << Re
          <<  "; noOfParticles=" << noOfParticles << "; maxTimeSteps(particle)="
          << converter.getLatticeTime( particleMaxPhysT )
          << "; St=" << ( 2.*partRho*radius*radius*converter.getCharPhysVelocity() ) / ( 9.*converter.getPhysViscosity()*converter.getPhysDensity()*converter.getCharPhysLength() ) << std::endl;
  }

  // Writes the vtk and gif files
  if ( iT % iTperiod == 0 ) {
    if ( !fluidExists && iT <= fluidMaxT ) {
      vtmWriterStartTime.write(iT);
      SuperEuklidNorm3D<T, DESCRIPTOR> normVel( velocity );
      BlockReduction3D2D<T> planeReduction( normVel, {0, -1, 0}, 600, BlockDataSyncMode::ReduceOnly );
      // write output as JPEG
      heatmap::write(planeReduction, iT);
    }
    if (iT > fluidMaxT) {
      // only write vtk-files after the fluid calculation is finished
      vtmWriter.write(iT - fluidMaxT);
    }
  }

  // Writes output on the console for the fluid phase
  if (iT < converter.getLatticeTime( fluidMaxPhysT ) && iT%iTperiod == 0 ) {

    // Timer statics
    fluidTimer.update( iT );
    fluidTimer.printStep();

    // Lattice statistics
    sLattice.getStatistics().print( iT, converter.getPhysTime( iT ) );

    // Flux at the inlet and outlet regions
    const std::vector<int> materials = { 1, 3, 4, 5 };

    IndicatorCircle3D<T> inlet(
      inletCenter + 2. * converter.getConversionFactorLength() * inletNormal,
      inletNormal, inletRadius + 2. * converter.getConversionFactorLength() );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxInflow( sLattice, converter, superGeometry, inlet, materials );
    vFluxInflow.print( "inflow", "ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxInflow( sLattice, converter, superGeometry, inlet, materials );
    pFluxInflow.print( "inflow", "N", "Pa" );

    IndicatorCircle3D<T> outlet0(
      outletCenter0 + 2. * converter.getConversionFactorLength() * outletNormal0,
      outletNormal0, outletRadius0 + 2. * converter.getConversionFactorLength() );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow0( sLattice, converter, superGeometry, outlet0, materials );
    vFluxOutflow0.print( "outflow0", "ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow0( sLattice, converter, superGeometry, outlet0, materials );
    pFluxOutflow0.print( "outflow0", "N", "Pa" );

    IndicatorCircle3D<T> outlet1(
      outletCenter1 + 2. * converter.getConversionFactorLength() * outletNormal1,
      outletNormal1, outletRadius1 + 2. * converter.getConversionFactorLength() );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow1( sLattice, converter, superGeometry, outlet1, materials );
    vFluxOutflow1.print( "outflow1", "ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow1( sLattice, converter, superGeometry, outlet1, materials );
    pFluxOutflow1.print( "outflow1", "N", "Pa" );
  }

  // Writes output on the console for the fluid phase
  if ( iT >= converter.getLatticeTime( fluidMaxPhysT ) &&
       (iT%iTperiod == 0 || iT == converter.getLatticeTime( fluidMaxPhysT )) ) {

    particleTimer.print( iT - fluidMaxT );

    // console output number of particles at different material numbers mat
    supParticleSystem.print( {1,2,3,4,5} );
    // console output of escape (E), capture (C) rate for material numbers mat
    supParticleSystem.captureEscapeRate( {4,5} );

    // only write vtk-files after the fluid calculation is finished
    supParticleWriter.write( iT - fluidMaxT );

    // true as long as certain amount of active particles
    if ( supParticleSystem.globalNumOfActiveParticles() < 0.001 * noOfParticles
         && iT > 0.9*converter.getLatticeTime( fluidMaxPhysT + particleMaxPhysT ) ) {
      return false;
    }
  }
  return true;
}

int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===

  olbInit( &argc, &argv );

  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout, "main" );

  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> const converter(
    int {N}, // resolution: number of voxels per charPhysL
    (T)   0.557646,                    // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   inletRadius*2.,              // charPhysLength: reference length of simulation geometry
    (T)   Re*1.5e-5/( inletRadius*2 ), // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1.5e-5,                      // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.225                        // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("bifurcation3d");

  // === 2nd Step: Prepare Geometry ===
  STLreader<T> stlReader( "../bifurcation3d.stl", converter.getConversionFactorLength() );
  IndicatorLayer3D<T> extendedDomain( stlReader,
                                      converter.getConversionFactorLength() );

  // Instantiation of an empty cuboidGeometry
  int noOfCuboids = std::max( 16, 4 * singleton::mpi().getSize() );

  CuboidGeometry3D<T> cuboidGeometry( extendedDomain, converter.getConversionFactorLength(),
                                      noOfCuboids );

  // Instantiation of an empty loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );
  // Default instantiation of superGeometry
  SuperGeometry3D<T> superGeometry( cuboidGeometry, loadBalancer, 2 );

  prepareGeometry( converter, extendedDomain, stlReader, superGeometry );

  // === 3rd Step: Prepare Lattice ===

  SuperLattice3D<T, DESCRIPTOR> sLattice( superGeometry );

  BGKdynamics<T, DESCRIPTOR> bulkDynamics( converter.getLatticeRelaxationFrequency(),
      instances::getBulkMomenta<T, DESCRIPTOR>() );

  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sBoundaryCondition( sLattice );
  createInterpBoundaryCondition3D<T, DESCRIPTOR>( sBoundaryCondition );

  sOffLatticeBoundaryCondition3D<T, DESCRIPTOR> sOffBoundaryCondition(
    sLattice );
  createBouzidiBoundaryCondition3D<T, DESCRIPTOR>( sOffBoundaryCondition );

  prepareLattice( sLattice, converter, bulkDynamics, sBoundaryCondition,
                  sOffBoundaryCondition, stlReader, superGeometry );

  // === 3.1 Step: Particles ===
  clout << "Prepare Particles ..." << std::endl;

  // SuperParticleSystems3D
  SuperParticleSystem3D<T, PARTICLE> supParticleSystem( superGeometry );
  // define which properties are to be written in output data
  SuperParticleSysVtuWriter<T, PARTICLE> supParticleWriter( supParticleSystem,
      "particles", SuperParticleSysVtuWriter<T, PARTICLE>::particleProperties::
      velocity
      | SuperParticleSysVtuWriter<T, PARTICLE>::particleProperties::mass
      | SuperParticleSysVtuWriter<T, PARTICLE>::particleProperties::radius
      | SuperParticleSysVtuWriter<T, PARTICLE>::particleProperties::active );

  SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR> getVel( sLattice, converter );

  auto stokesDragForce = make_shared
                         < StokesDragForce3D<T, PARTICLE, DESCRIPTOR>
                         > ( getVel, converter );

  // material numbers where particles should be reflected
  std::set<int> boundMaterial = { 2, 4, 5};
  auto materialBoundary = make_shared
                          < MaterialBoundary3D<T, PARTICLE>
                          > ( superGeometry, boundMaterial );

  supParticleSystem.addForce( stokesDragForce );
  supParticleSystem.addBoundary( materialBoundary );
  supParticleSystem.setOverlap( 2. * converter.getConversionFactorLength() );

  // particles generation at inlet3
  Vector<T, 3> c( inletCenter );
  c[2] = 0.074;
  IndicatorCircle3D<T> inflowCircle( c, inletNormal, inletRadius -
                                     converter.getConversionFactorLength() * 2.5 );
  IndicatorCylinder3D<T> inletCylinder( inflowCircle, 0.01 *
                                        converter.getConversionFactorLength() );
  supParticleSystem.addParticle( inletCylinder, 4. / 3. * M_PI *
                                 std::pow( radius, 3 ) * partRho, radius, noOfParticles );

  clout << "Prepare Particles ... OK" << std::endl;

  // === 4th Step: Main Loop with Timer ===

  Timer<double> fluidTimer( converter.getLatticeTime( fluidMaxPhysT ),
                            superGeometry.getStatistics().getNvoxel() );

  Timer<double> particleTimer( converter.getLatticeTime( particleMaxPhysT ),
                               noOfParticles );
  fluidTimer.start();

  int iT = 0;
  // amount of timesteps when getResults rewrites data
  int iTperiod = converter.getLatticeTime( .2 );

  bool fluidExists = true;

  // checks whether there is already data of the fluid from an earlier calculation
  if ( !( sLattice.load( "fluidSolution" ) ) ) {

    fluidExists = false;

    // if there is no data available, it is generated
    for ( ; iT <= converter.getLatticeTime( fluidMaxPhysT ); ++iT ) {

      // during run up time boundary values are set, collide and stream step,
      // results of fluid, afterwards only particles are simulated
      setBoundaryValues( sLattice, converter, iT, fluidMaxPhysT, superGeometry );
      sLattice.collideAndStream();

      getResults( sLattice, converter, iT, iTperiod, superGeometry, fluidTimer, stlReader,
                  supParticleSystem, radius, partRho, particleTimer,
                  supParticleWriter, fluidExists );
    }

    fluidTimer.stop();
    fluidTimer.printSummary();

    sLattice.communicate();
    // calculated results are written in a file
    sLattice.save( "fluidSolution" );
  }

  // if there exists already data of the fluid from an earlier calculation, this is used
  else {

    iT = converter.getLatticeTime( fluidMaxPhysT );
    getResults( sLattice, converter, iT,
                iTperiod, superGeometry, fluidTimer, stlReader,
                supParticleSystem, radius, partRho, particleTimer,
                supParticleWriter, fluidExists );

  }

  // after the fluid calculation, particle simulation starts
  supParticleSystem.setVelToFluidVel( getVel );
  particleTimer.start();

  for ( ; iT <= converter.getLatticeTime( fluidMaxPhysT + particleMaxPhysT ); ++iT ) {
    // particles simulation starts after run up time is over
    supParticleSystem.simulate( converter.getConversionFactorTime() );

    if ( !getResults( sLattice, converter, iT,
                      iTperiod, superGeometry, fluidTimer,
                      stlReader, supParticleSystem, radius, partRho,
                      particleTimer, supParticleWriter, fluidExists) ) {
      break;
    }
  }
  particleTimer.stop();
  particleTimer.printSummary();
}
