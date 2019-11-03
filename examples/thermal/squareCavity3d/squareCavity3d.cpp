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

// natural convection of air in a square cavity in 3D

#include "olb3D.h"
#include "olb3D.hh"   // use only generic version!

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

#define NSDESCRIPTOR D3Q19<FORCE>
#define TDESCRIPTOR D3Q7<VELOCITY>

// Parameters for the simulation setup
T Ra = 1e3;  // Rayleigh-Zahl
const T Pr = 0.71; // Prandtl-Zahl

T lx;

int N = 64; // resolution of the model

const T maxPhysT = 1e4;   // max. simulation time in s, SI unit
const T epsilon = 1.e-3;  // precision of the convergence (residuum)

const T Tcold = 275.15;
const T Thot = 285.15;
const T Tmean = (Tcold + Thot) / 2.0;

/// Values from the literature studies from Davis
T LitVelocity3[] = { 3.649, 3.696, 1.013 };
T LitPosition3[] = { 0.813, 0.178 };
T LitVelocity4[] = { 16.178, 19.617, 1.212 };
T LitPosition4[] = { 0.823, 0.119 };
T LitVelocity5[] = { 34.730, 68.590, 1.975 };
T LitPosition5[] = { 0.855, 0.066 };
T LitVelocity6[] = { 64.530, 219.36, 3.400 };
T LitPosition6[] = { 0.850, 0.036 };
T LitNusselt3 = 1.117;
T LitNusselt4 = 2.238;
T LitNusselt5 = 4.509;
T LitNusselt6 = 8.817;

/// Compute the nusselt number at the left wall
T computeNusselt(SuperGeometry3D<T>& superGeometry,
                 SuperLattice3D<T, NSDESCRIPTOR>& NSlattice,
                 SuperLattice3D<T, TDESCRIPTOR>& ADlattice)
{
  int voxel = 0, material = 0;
  T T_x = 0, T_xplus1 = 0, T_xplus2 = 0;
  T q = 0;

  for (int iC = 0; iC < NSlattice.getLoadBalancer().size(); iC++) {
    int ny = NSlattice.getBlockLattice(iC).getNy();

    int iX = 0;
    int iZ = 1;

    for (int iY = 0; iY < ny; ++iY) {
      material = superGeometry.getBlockGeometry(iC).getMaterial(iX,iY,iZ);

      T_x = ADlattice.getBlockLattice(iC).get(iX,iY,iZ).computeRho();
      T_xplus1 = ADlattice.getBlockLattice(iC).get(iX+1,iY,iZ).computeRho();
      T_xplus2 = ADlattice.getBlockLattice(iC).get(iX+2,iY,iZ).computeRho();

      if ( material == 2 ) {
        q += (3.0*T_x - 4.0*T_xplus1 + 1.0*T_xplus2)/2.0*N;
        voxel++;
      }
    }
  }

#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(q, MPI_SUM);
  singleton::mpi().reduceAndBcast(voxel, MPI_SUM);
#endif

  return q / (T)voxel;
}

/// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry3D<T>& superGeometry,
                     ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter)
{

  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,4);

  std::vector<T> extend(3,T());
  extend[0] = lx;
  extend[1] = lx;
  extend[2] = 3.0 * converter.getPhysLength(1);
  std::vector<T> origin(3,T());
  origin[0] = converter.getPhysLength(1);
  origin[1] = 0.5*converter.getPhysLength(1);
  origin[2] = 0.0;
  IndicatorCuboid3D<T> cuboid2(extend, origin);

  superGeometry.rename(4,1,cuboid2);

  std::vector<T> extendwallleft(3,T(0));
  extendwallleft[0] = converter.getPhysLength(1);
  extendwallleft[1] = lx;
  extendwallleft[2] = 0.1;
  std::vector<T> originwallleft(3,T(0));
  originwallleft[0] = 0.0;
  originwallleft[1] = 0.0;
  originwallleft[2] = 0.0;
  IndicatorCuboid3D<T> wallleft(extendwallleft, originwallleft);

  std::vector<T> extendwallright(3,T(0));
  extendwallright[0] = converter.getPhysLength(1);
  extendwallright[1] = lx;
  extendwallright[2] = 0.1;
  std::vector<T> originwallright(3,T(0));
  originwallright[0] = lx+converter.getPhysLength(1);
  originwallright[1] = 0.0;
  originwallright[2] = 0.0;
  IndicatorCuboid3D<T> wallright(extendwallright, originwallright);

  superGeometry.rename(4,2,1,wallleft);
  superGeometry.rename(4,3,1,wallright);


  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                     SuperLattice3D<T, NSDESCRIPTOR>& NSlattice,
                     SuperLattice3D<T, TDESCRIPTOR>& ADlattice,
                     Dynamics<T, NSDESCRIPTOR> &bulkDynamics,
                     Dynamics<T, TDESCRIPTOR>& advectionDiffusionBulkDynamics,
                     sOnLatticeBoundaryCondition3D<T,NSDESCRIPTOR>& NSboundaryCondition,
                     sOnLatticeBoundaryCondition3D<T,TDESCRIPTOR>& TboundaryCondition,
                     SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  T omega  = converter.getLatticeRelaxationFrequency();
  T Tomega  = converter.getLatticeThermalRelaxationFrequency();

  /// define lattice Dynamics
  ADlattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, TDESCRIPTOR>());
  NSlattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, NSDESCRIPTOR>());

  ADlattice.defineDynamics(superGeometry.getMaterialIndicator({1, 2, 3}), &advectionDiffusionBulkDynamics);
  ADlattice.defineDynamics(superGeometry, 4, &instances::getBounceBack<T, TDESCRIPTOR>());

  NSlattice.defineDynamics(superGeometry.getMaterialIndicator({1, 2, 3}), &bulkDynamics);
  NSlattice.defineDynamics(superGeometry, 4, &instances::getBounceBack<T, NSDESCRIPTOR>());

  /// sets boundary
  TboundaryCondition.addTemperatureBoundary(superGeometry.getMaterialIndicator({2, 3}), Tomega);
  NSboundaryCondition.addVelocityBoundary(superGeometry.getMaterialIndicator({2, 3}), omega);

  /// define initial conditions
  AnalyticalConst3D<T,T> rho(1.);
  AnalyticalConst3D<T,T> u0(0.0, 0.0, 0.0);
  AnalyticalConst3D<T,T> T_cold(converter.getLatticeTemperature(Tcold));
  AnalyticalConst3D<T,T> T_hot(converter.getLatticeTemperature(Thot));
  AnalyticalConst3D<T,T> T_mean(converter.getLatticeTemperature(Tmean));

  /// for each material set Rho, U and the Equilibrium
  NSlattice.defineRhoU(superGeometry, 1, rho, u0);
  NSlattice.iniEquilibrium(superGeometry, 1, rho, u0);
  NSlattice.defineRhoU(superGeometry, 2, rho, u0);
  NSlattice.iniEquilibrium(superGeometry, 2, rho, u0);
  NSlattice.defineRhoU(superGeometry, 3, rho, u0);
  NSlattice.iniEquilibrium(superGeometry, 3, rho, u0);

  ADlattice.defineRho(superGeometry, 1, T_mean);
  ADlattice.iniEquilibrium(superGeometry, 1, T_mean, u0);
  ADlattice.defineRho(superGeometry, 2, T_hot);
  ADlattice.iniEquilibrium(superGeometry, 2, T_hot, u0);
  ADlattice.defineRho(superGeometry, 3, T_cold);
  ADlattice.iniEquilibrium(superGeometry, 3, T_cold, u0);

  /// Make the lattice ready for simulation
  NSlattice.initialize();
  ADlattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                       SuperLattice3D<T, NSDESCRIPTOR>& NSlattice,
                       SuperLattice3D<T, TDESCRIPTOR>& ADlattice,
                       int iT, SuperGeometry3D<T>& superGeometry)
{

  // nothing to do here

}

void getResults(ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> &converter,
                SuperLattice3D<T, NSDESCRIPTOR>& NSlattice,
                SuperLattice3D<T, TDESCRIPTOR>& ADlattice, int iT,
                SuperGeometry3D<T>& superGeometry,
                Timer<T>& timer,
                bool converged)
{

  OstreamManager clout(std::cout,"getResults");

  SuperVTMwriter3D<T> vtkWriter("thermalNaturalConvection3D");
  SuperLatticeGeometry3D<T, NSDESCRIPTOR> geometry(NSlattice, superGeometry);
  SuperLatticePhysVelocity3D<T, NSDESCRIPTOR> velocity(NSlattice, converter);
  SuperLatticePhysPressure3D<T, NSDESCRIPTOR> pressure(NSlattice, converter);
  SuperLatticePhysTemperature3D<T, NSDESCRIPTOR,TDESCRIPTOR> temperature(ADlattice, converter);
  vtkWriter.addFunctor( geometry );
  vtkWriter.addFunctor( pressure );
  vtkWriter.addFunctor( velocity );
  vtkWriter.addFunctor( temperature );

  AnalyticalFfromSuperF3D<T> interpolation(velocity, true);

  const int vtkIter = 2000.;

  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, NSDESCRIPTOR> cuboid(NSlattice);
    SuperLatticeRank3D<T, NSDESCRIPTOR> rank(NSlattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  /// Writes the VTK files
  if (iT%vtkIter == 0 || converged) {

    timer.update(iT);
    timer.printStep();

    /// NSLattice statistics console output
    NSlattice.getStatistics().print(iT,converter.getPhysTime(iT));
    /// ADLattice statistics console output
    ADlattice.getStatistics().print(iT,converter.getPhysTime(iT));

    vtkWriter.write(iT);
    const double a[3] = {0, 0, 1.};
    BlockReduction3D2D<T> planeReduction(temperature, a);
    BlockGifWriter<T> gifWriter;
    gifWriter.write(planeReduction, Tcold*0.98, Thot*1.02, iT, "temperature");

    SuperEuklidNorm3D<T, NSDESCRIPTOR> normVel( velocity );
    BlockReduction3D2D<T> planeReduction2(normVel, {0, 0, 1});
    BlockGifWriter<T> gifWriter2;
    gifWriter2.write( planeReduction2, iT, "velocity" );

  }

  if ( converged ) {

    T nusselt = computeNusselt(superGeometry, NSlattice, ADlattice);

    /// Initialize vectors for data output
    T xVelocity[3] = { T() };
    T outputVelX[3] = { T() };
    T yVelocity[3] = { T() };
    T outputVelY[3] = { T() };
    const int outputSize = 512;
    Vector<T, outputSize> velX;
    Vector<T, outputSize> posX;
    Vector<T, outputSize> velY;
    Vector<T, outputSize> posY;

    /// loop for the resolution of the cavity at x = lx/2 in yDirection and vice versa
    for (int n = 0; n < outputSize; ++n) {
      T yPosition[3] = { lx / 2, lx * n / (T) outputSize, lx / N * 2 / 2 };
      T xPosition[3] = { lx * n / (T) outputSize, lx / 2, lx / N * 2 / 2 };

      /// Interpolate xVelocity at x = lx/2 for each yPosition
      interpolation(xVelocity, yPosition);
      interpolation(yVelocity, xPosition);
      /// Store the interpolated values to compare them among each other in order to detect the maximum
      velX[n] = xVelocity[0];
      posY[n] = yPosition[1];
      velY[n] = yVelocity[1];
      posX[n] = xPosition[0];

      /// Initialize output with the corresponding velocities and positions at the origin
      if (n == 0) {
        outputVelX[0] = velX[0];
        outputVelX[1] = posY[0];
        outputVelY[0] = velY[0];
        outputVelY[1] = posX[0];
      }
      /// look for the maximum velocity in xDirection and the corresponding position in yDirection
      if (n > 0 && velX[n] > outputVelX[0]) {
        outputVelX[0] = velX[n];
        outputVelX[1] = posY[n];
      }
      /// look for the maximum velocity in yDirection and the corresponding position in xDirection
      if (n > 0 && velY[n] > outputVelY[0]) {
        outputVelY[0] = velY[n];
        outputVelY[1] = posX[n];
      }
    }

    // compare to De Vahl Davis' benchmark solutions
    clout << "Comparison against De Vahl Davis (1983):" << endl;
    if (Ra == 1e3) {
      clout << "xVelocity in yDir=" <<  outputVelX[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength() << "; error(rel)=" << (T) fabs((LitVelocity3[0] - outputVelX[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength()) / LitVelocity3[0]) << endl;
      clout << "yVelocity in xDir=" <<  outputVelY[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength() << "; error(rel)=" << (T) fabs((LitVelocity3[1] - outputVelY[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength()) / LitVelocity3[1]) << endl;
      clout << "yMaxVel / xMaxVel="  <<  outputVelY[0] / outputVelX[0] << "; error(rel)=" << (T) fabs((LitVelocity3[2] - outputVelY[0] / outputVelX[0])  / LitVelocity3[2]) << endl;
      clout << "yCoord of xMaxVel=" <<  outputVelX[1]/lx << "; error(rel)=" << (T) fabs((LitPosition3[0] - outputVelX[1] / lx) / LitPosition3[0]) << endl;
      clout << "xCoord of yMaxVel=" <<   outputVelY[1]/lx << "; error(rel)=" << (T) fabs((LitPosition3[1] - outputVelY[1] / lx) / LitPosition3[1]) << endl;
      clout << "Nusselt=" <<  nusselt << "; error(rel)=" << (T) fabs((LitNusselt3 - nusselt) / nusselt) << endl;
    }
    else if (Ra == 1e4) {
      clout << "xVelocity in yDir=" <<  outputVelX[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength() << "; error(rel)=" << (T) fabs((LitVelocity4[0] - outputVelX[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength()) / LitVelocity4[0]) << endl;
      clout << "yVelocity in xDir=" <<  outputVelY[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength() << "; error(rel)=" << (T) fabs((LitVelocity4[1] - outputVelY[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength()) / LitVelocity4[1]) << endl;
      clout << "yMaxVel / xMaxVel="  <<  outputVelY[0] / outputVelX[0] << "; error(rel)=" << (T) fabs((LitVelocity4[2] - outputVelY[0] / outputVelX[0])  / LitVelocity4[2]) << endl;
      clout << "yCoord of xMaxVel=" <<  outputVelX[1]/lx << "; error(rel)=" << (T) fabs((LitPosition4[0] - outputVelX[1] / lx) / LitPosition4[0]) << endl;
      clout << "xCoord of yMaxVel=" <<   outputVelY[1]/lx << "; error(rel)=" << (T) fabs((LitPosition4[1] - outputVelY[1] / lx) / LitPosition4[1]) << endl;
      clout << "Nusselt=" <<  nusselt << "; error(rel)=" << (T) fabs((LitNusselt4 - nusselt) / nusselt) << endl;
    }
    else if (Ra == 1e5) {
      clout << "xVelocity in yDir=" <<  outputVelX[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength() << "; error(rel)=" << (T) fabs((LitVelocity5[0] - outputVelX[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength()) / LitVelocity5[0]) << endl;
      clout << "yVelocity in xDir=" <<  outputVelY[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength() << "; error(rel)=" << (T) fabs((LitVelocity5[1] - outputVelY[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength()) / LitVelocity5[1]) << endl;
      clout << "yMaxVel / xMaxVel="  <<  outputVelY[0] / outputVelX[0] << "; error(rel)=" << (T) fabs((LitVelocity5[2] - outputVelY[0] / outputVelX[0])  / LitVelocity5[2]) << endl;
      clout << "yCoord of xMaxVel=" <<  outputVelX[1]/lx << "; error(rel)=" << (T) fabs((LitPosition5[0] - outputVelX[1] / lx) / LitPosition5[0]) << endl;
      clout << "xCoord of yMaxVel=" <<   outputVelY[1]/lx << "; error(rel)=" << (T) fabs((LitPosition5[1] - outputVelY[1] / lx) / LitPosition5[1]) << endl;
      clout << "Nusselt=" <<  nusselt << "; error(rel)=" << (T) fabs((LitNusselt5 - nusselt) / nusselt) << endl;
    }
    else if (Ra == 1e6) {
      clout << "xVelocity in yDir=" <<  outputVelX[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength() << "; error(rel)=" << (T) fabs((LitVelocity6[0] - outputVelX[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength()) / LitVelocity6[0]) << endl;
      clout << "yVelocity in xDir=" <<  outputVelY[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength() << "; error(rel)=" << (T) fabs((LitVelocity6[1] - outputVelY[0] / converter.getPhysThermalDiffusivity() * converter.getCharPhysLength()) / LitVelocity6[1]) << endl;
      clout << "yMaxVel / xMaxVel="  <<  outputVelY[0] / outputVelX[0] << "; error(rel)=" << (T) fabs((LitVelocity6[2] - outputVelY[0] / outputVelX[0])  / LitVelocity6[2]) << endl;
      clout << "yCoord of xMaxVel=" <<  outputVelX[1]/lx << "; error(rel)=" << (T) fabs((LitPosition6[0] - outputVelX[1] / lx) / LitPosition6[0]) << endl;
      clout << "xCoord of yMaxVel=" <<   outputVelY[1]/lx << "; error(rel)=" << (T) fabs((LitPosition6[1] - outputVelY[1] / lx) / LitPosition6[1]) << endl;
      clout << "Nusselt=" <<  nusselt << "; error(rel)=" << (T) fabs((LitNusselt6 - nusselt) / nusselt) << endl;
    }
  }
}

int main(int argc, char *argv[])
{

  /// === 1st Step: Initialization ===
  OstreamManager clout(std::cout,"main");
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  T tau = 0.9;

  if (argc>=2) {
    Ra = atof(argv[1]);
  }

  lx  = pow(Ra * 15.126e-6 * 15.126e-6 / Pr / 9.81 / (Thot - Tcold) / 0.00341, (T) 1/3);  // length of the square
  T charU = 1.0 / lx /( Pr * 25.684e-3 / 15.126e-6 / 1.0 * 1.0 / 25.684e-3);

  if (Ra==1e3) {
    charU *= LitVelocity3[1];
    N = 64;
  }
  if (Ra==1e4) {
    charU *= LitVelocity4[1];
    N = 128;
  }
  if (Ra==1e5) {
    charU *= LitVelocity5[1];
    N = 256;
  }
  if (Ra==1e6) {
    charU *= LitVelocity6[1];
    N = 512;
  }

  ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> converter(
    (T) lx / N,
    (T) (tau - 0.5) / descriptors::invCs2<T,NSDESCRIPTOR>() * pow((lx/N),2) / 15.126e-6,
    (T) lx,
    (T) charU,
    (T) 15.126e-6,
    (T) 1.0,
    (T) 25.684e-3,
    (T) Pr * 25.684e-3 / 15.126e-6 / 1.0,
    (T) 0.00341,
    (T) Tcold,
    (T) Thot
  );
  converter.print();

  /// === 2nd Step: Prepare Geometry ===
  std::vector<T> extend(3,T());
  extend[0] = lx + 2*converter.getPhysLength(1);
  extend[1] = lx + converter.getPhysLength(1);
  extend[2] = 2.*converter.getPhysLength(1);
  std::vector<T> origin(3,T());
  IndicatorCuboid3D<T> cuboid(extend, origin);

  /// Instantiation of an empty cuboidGeometry
  CuboidGeometry3D<T> cuboidGeometry(cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize());
  cuboidGeometry.setPeriodicity(false, false, true);  // x, y, z

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  /// Instantiation of a superGeometry
  SuperGeometry3D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(superGeometry, converter);

  /// === 3rd Step: Prepare Lattice ===

  SuperLattice3D<T, TDESCRIPTOR> ADlattice(superGeometry);
  SuperLattice3D<T, NSDESCRIPTOR> NSlattice(superGeometry);

  sOnLatticeBoundaryCondition3D<T, NSDESCRIPTOR> NSboundaryCondition(NSlattice);
  createLocalBoundaryCondition3D<T, NSDESCRIPTOR>(NSboundaryCondition);

  sOnLatticeBoundaryCondition3D<T, TDESCRIPTOR> TboundaryCondition(ADlattice);
  createAdvectionDiffusionBoundaryCondition3D<T, TDESCRIPTOR>(TboundaryCondition);

  ForcedBGKdynamics<T, NSDESCRIPTOR> NSbulkDynamics(
    converter.getLatticeRelaxationFrequency(),
    instances::getBulkMomenta<T,NSDESCRIPTOR>());

  AdvectionDiffusionBGKdynamics<T, TDESCRIPTOR> TbulkDynamics (
    converter.getLatticeThermalRelaxationFrequency(),
    instances::getAdvectionDiffusionBulkMomenta<T,TDESCRIPTOR>());

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
  // This coupling must be necessarily be put on the Navier-Stokes lattice!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

  std::vector<T> dir{0.0, 1.0, 0.0};

  T boussinesqForcePrefactor = 9.81 / converter.getConversionFactorVelocity() * converter.getConversionFactorTime() *
                               converter.getCharPhysTemperatureDifference() * converter.getPhysThermalExpansionCoefficient();

  NavierStokesAdvectionDiffusionCouplingGenerator3D<T,NSDESCRIPTOR>
  coupling(0, converter.getLatticeLength(lx), 0, converter.getLatticeLength(lx), 0, converter.getLatticeLength(lx),
           boussinesqForcePrefactor, converter.getLatticeTemperature(Tcold), 1., dir);

  NSlattice.addLatticeCoupling(coupling, ADlattice);

  prepareLattice(converter,
                 NSlattice, ADlattice,
                 NSbulkDynamics, TbulkDynamics,
                 NSboundaryCondition, TboundaryCondition, superGeometry );


  /// === 4th Step: Main Loop with Timer ===
  Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  util::ValueTracer<T> converge(6,epsilon);
  for (int iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {

    if (converge.hasConverged()) {
      clout << "Simulation converged." << endl;
      clout << "Time " << iT << "." << std::endl;

      getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, converge.hasConverged());

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
    if (iT % 1000 == 0) {
      converge.takeValue(computeNusselt(superGeometry, NSlattice, ADlattice),true);
    }
  }

  timer.stop();
  timer.printSummary();
}
