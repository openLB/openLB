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

#include "olb3D.h"
#include "olb3D.hh"   // use only generic version!
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


#define NSDESCRIPTOR D3Q19<FORCE>
#define TDESCRIPTOR D3Q7<VELOCITY>




// const int maxIter = 1000000;
// const int saveIter = 5000;

// Parameters for the simulation setup
const T lx = 1.0;      // length of the channel
const T ly = 1.0;      // height of the channel
T lz = 1.0;
int N = 20;        // resolution of the model
T tau = 1.;      // relaxation time
const T Re = 5.;       // Reynolds number
const T Ra = 100.;     // Rayleigh number
const T Pr = 0.71;     // Prandtl number
const T maxPhysT = 1e4; // max. simulation time in s, SI unit
const T epsilon = 1.e-7; // precision of the convergence (residuum)

const T Tcold = 273.15;
const T Thot = 274.15;


template <typename T, typename S>
class AnalyticalVelocityPorousPlate3D : public AnalyticalF3D<T, S> {
private:
  T _Re;
  T _u0;
  T _v0;
  T _ly;
public:
  AnalyticalVelocityPorousPlate3D(T Re, T u0, T v0, T ly) : AnalyticalF3D<T, S>(3),
    _Re(Re), _u0(u0), _v0(v0), _ly(ly)
  {
    this->getName() = "AnalyticalVelocityPorousPlate3D";
  };

  bool operator()(T output[3], const S x[3])
  {
    output[0] = _u0*((exp(_Re* x[1] / _ly) - 1) / (exp(_Re) - 1));
    output[1] = _v0;
    output[2] = 0.0;
    return true;
  };
};

template <typename T, typename S>
class AnalyticalTemperaturePorousPlate3D : public AnalyticalF3D<T, S> {
private:
  T _Re;
  T _Pr;
  T _ly;
  T _T0;
  T _deltaT;
public:
  AnalyticalTemperaturePorousPlate3D(T Re, T u0, T ly, T T0, T deltaT) : AnalyticalF3D<T, S>(1),
    _Re(Re), _Pr(Pr), _ly(ly), _T0(T0), _deltaT(deltaT)
  {
    this->getName() = "AnalyticalTemperaturePorousPlate3D";
  };

  bool operator()(T output[1], const S x[3])
  {
    output[0] = _T0 + _deltaT*((exp(_Pr*_Re*x[1] / _ly) - 1) / (exp(_Pr*_Re) - 1));
    return true;
  };
};

template <typename T, typename S>
class AnalyticalHeatFluxPorousPlate3D : public AnalyticalF3D<T, S> {
private:
  T _Re;
  T _Pr;
  T _deltaT;
  T _ly;
  T _lambda;
public:
  AnalyticalHeatFluxPorousPlate3D(T Re, T Pr, T deltaT, T ly,T lambda) : AnalyticalF3D<T, S>(3),
   _Re(Re), _Pr(Pr), _deltaT(deltaT), _ly(ly), _lambda(lambda)
  {
    this->getName() = "AnalyticalHeatFluxPorousPlate3D";
  };

  bool operator()(T output[3], const S x[3])
  {
    output[0] = 0;
    output[1] = - _lambda * _Re * _Pr * _deltaT / _ly * (exp(_Pr * _Re * x[1] / _ly))/(exp(_Pr * _Re) - 1);
    output[2] = 0;
    return true;
  };
};

void error( SuperGeometry3D<T>& superGeometry,
            SuperLattice3D<T, NSDESCRIPTOR>& NSlattice,
            SuperLattice3D<T, TDESCRIPTOR>& ADlattice,
            ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
            T Re )
{
  OstreamManager clout( std::cout, "error" );

  int input[1] = { };
  T result[1]  = { };

  auto indicatorF = superGeometry.getMaterialIndicator(1);

  T u_Re = Re * converter.getPhysViscosity() / converter.getCharPhysLength();
  AnalyticalVelocityPorousPlate3D<T,T> uSol(Re,converter.getCharPhysVelocity(), u_Re, converter.getCharPhysLength());
  SuperLatticePhysVelocity3D<T,NSDESCRIPTOR> u(NSlattice,converter);

  SuperAbsoluteErrorL2Norm3D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, input);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, input);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  AnalyticalTemperaturePorousPlate3D<T,T> tSol(Re, Pr,  converter.getCharPhysLength(),  converter.getCharPhysLowTemperature(), converter.getCharPhysTemperatureDifference());
  SuperLatticePhysTemperature3D<T, NSDESCRIPTOR, TDESCRIPTOR> t(ADlattice, converter);

  SuperAbsoluteErrorL2Norm3D<T> absTemperatureErrorNormL2(t, tSol, indicatorF);
  absTemperatureErrorNormL2(result, input);
  clout << "temperature-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relTemperatureErrorNormL2(t, tSol, indicatorF);
  relTemperatureErrorNormL2(result, input);
  clout << "; temperature-L2-error(rel)=" << result[0] << std::endl;

  AnalyticalHeatFluxPorousPlate3D<T,T> HeatFluxSol(Re, Pr, converter.getCharPhysTemperatureDifference(), converter.getCharPhysLength(), converter.getThermalConductivity());
  SuperLatticePhysHeatFlux3D<T,NSDESCRIPTOR,TDESCRIPTOR> HeatFlux(ADlattice,converter);

  SuperAbsoluteErrorL2Norm3D<T> absHeatFluxErrorNormL2(HeatFlux, HeatFluxSol, indicatorF);
  absHeatFluxErrorNormL2(result, input);
  clout << "heatFlux-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relHeatFluxErrorNormL2(HeatFlux, HeatFluxSol, indicatorF);
  relHeatFluxErrorNormL2(result, input);
  clout << "; heatFlux-L2-error(rel)=" << result[0] << std::endl;
}

/// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry3D<T>& superGeometry,
                     ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const&converter)
{

  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,2);
  superGeometry.rename(2,1,0,1,0);
  //superGeometry.clean();

  std::vector<T> extend( 3, T(0) );
  extend[0] = lx;
  extend[1] = converter.getPhysLength(1);
  extend[2] = lz;
  std::vector<T> origin( 3, T(0) );
  IndicatorCuboid3D<T> bottom(extend, origin);
  /// Set material number for bottom
  superGeometry.rename(2,3,1,bottom);

  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                     SuperLattice3D<T, NSDESCRIPTOR>& NSlattice,
                     SuperLattice3D<T, TDESCRIPTOR>& ADlattice,
                     Dynamics<T, NSDESCRIPTOR> &bulkDynamics,
                     Dynamics<T, TDESCRIPTOR>& advectionDiffusionBulkDynamics,
                     sOnLatticeBoundaryCondition3D<T,NSDESCRIPTOR>& NSboundaryCondition,
                     sOnLatticeBoundaryCondition3D<T,TDESCRIPTOR>& TboundaryCondition,
                     SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout(std::cout,"prepareLattice");

  T Tomega  = converter.getLatticeThermalRelaxationFrequency();
  T NSomega = converter.getLatticeRelaxationFrequency();

  /// define lattice Dynamics
  clout << "defining dynamics" << endl;

  ADlattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, TDESCRIPTOR>());
  NSlattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, NSDESCRIPTOR>());

  ADlattice.defineDynamics(superGeometry.getMaterialIndicator({1, 2, 3}), &advectionDiffusionBulkDynamics);
  NSlattice.defineDynamics(superGeometry.getMaterialIndicator({1, 2, 3}), &bulkDynamics);

  /// sets boundary
  NSboundaryCondition.addVelocityBoundary(superGeometry.getMaterialIndicator({2, 3}), NSomega);
  TboundaryCondition.addTemperatureBoundary(superGeometry.getMaterialIndicator({2, 3}), Tomega);
}

void setBoundaryValues(ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                       SuperLattice3D<T, NSDESCRIPTOR>& NSlattice,
                       SuperLattice3D<T, TDESCRIPTOR>& ADlattice,
                       int iT, SuperGeometry3D<T>& superGeometry)
{

  if (iT == 0) {

//    typedef advectionDiffusionLbHelpers<T,TDESCRIPTOR> TlbH;

    /// for each material set the defineRhoU and the Equilibrium
    std::vector<T> zero(3,T());
    AnalyticalConst3D<T,T> u(zero);
    AnalyticalConst3D<T,T> rho(1.);
    AnalyticalConst3D<T,T> force(zero);

    T u_Re = converter.getLatticeVelocity( Re * converter.getPhysViscosity() / converter.getCharPhysLength() );
    AnalyticalConst3D<T,T> u_top(converter.getCharLatticeVelocity(), u_Re, 0.0);
    AnalyticalConst3D<T,T> u_bot(0.0, u_Re, 0.0);

    NSlattice.defineRhoU(superGeometry, 1, rho, u);
    NSlattice.iniEquilibrium(superGeometry, 1, rho, u);
    NSlattice.defineField<descriptors::FORCE>(superGeometry, 1, force);
    NSlattice.defineRhoU(superGeometry, 2, rho, u_top);
    NSlattice.iniEquilibrium(superGeometry, 2, rho, u_top);
    NSlattice.defineField<descriptors::FORCE>(superGeometry, 2, force);
    NSlattice.defineRhoU(superGeometry, 3, rho, u_bot);
    NSlattice.iniEquilibrium(superGeometry, 3, rho, u_bot);
    NSlattice.defineField<descriptors::FORCE>(superGeometry, 3, force);

    AnalyticalConst3D<T,T> Cold(converter.getLatticeTemperature(Tcold));
    AnalyticalConst3D<T,T> Hot(converter.getLatticeTemperature(Thot));

    ADlattice.defineRho(superGeometry, 1, Cold);
    ADlattice.iniEquilibrium(superGeometry, 1, Cold, u);
    ADlattice.defineField<descriptors::VELOCITY>(superGeometry, 1, u);
    ADlattice.defineRho(superGeometry, 2, Hot);
    ADlattice.iniEquilibrium(superGeometry, 2, Hot, u);
    ADlattice.defineField<descriptors::VELOCITY>(superGeometry, 2, u);
    ADlattice.defineRho(superGeometry, 3, Cold);
    ADlattice.iniEquilibrium(superGeometry, 3, Cold, u);
    ADlattice.defineField<descriptors::VELOCITY>(superGeometry, 3, u);

    /// Make the lattice ready for simulation
    NSlattice.initialize();
    ADlattice.initialize();
  }
}

void getResults(ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                SuperLattice3D<T, NSDESCRIPTOR>& NSlattice,
                SuperLattice3D<T, TDESCRIPTOR>& ADlattice, int iT,
                SuperGeometry3D<T>& superGeometry,
                Timer<T>& timer,
                bool converged)
{

  OstreamManager clout(std::cout,"getResults");

  SuperVTMwriter3D<T> vtkWriter("thermalPorousPlate3d");

  SuperLatticePhysVelocity3D<T, NSDESCRIPTOR> velocity(NSlattice, converter);
  SuperLatticePhysPressure3D<T, NSDESCRIPTOR> pressure(NSlattice, converter);
  SuperLatticePhysTemperature3D<T, NSDESCRIPTOR, TDESCRIPTOR> temperature(ADlattice, converter);

  AnalyticalHeatFluxPorousPlate3D<T,T> HeatFluxSol(Re, Pr, converter.getCharPhysTemperatureDifference(), converter.getCharPhysLength(), converter.getThermalConductivity());
  SuperLatticePhysHeatFlux3D<T,NSDESCRIPTOR,TDESCRIPTOR> HeatFlux1(ADlattice,converter);
  SuperLatticeFfromAnalyticalF3D<T,TDESCRIPTOR> HeatFluxSolLattice(HeatFluxSol,ADlattice);
  vtkWriter.addFunctor( pressure );
  vtkWriter.addFunctor( velocity );
  vtkWriter.addFunctor( temperature );
  vtkWriter.addFunctor( HeatFlux1 );
  vtkWriter.addFunctor( HeatFluxSolLattice );

  const int vtkIter = converter.getLatticeTime(100.);

  if (iT == 0) {
    /// Writes the converter log file
    //writeLogFile(converter,"thermalPorousPlate3d");

    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, NSDESCRIPTOR> geometry(NSlattice, superGeometry);
    SuperLatticeCuboid3D<T, NSDESCRIPTOR> cuboid(NSlattice);
    SuperLatticeRank3D<T, NSDESCRIPTOR> rank(NSlattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  /// Writes the VTK files
  if (iT%vtkIter == 0 || converged) {
    NSlattice.getStatistics().print(iT,converter.getPhysTime(iT));
    timer.print(iT);
    error(superGeometry, NSlattice, ADlattice, converter, Re);

    vtkWriter.write(iT);

      BlockReduction3D2D<T> planeReduction( temperature, {0,0,1}, 600, BlockDataSyncMode::ReduceOnly );
    // write output as JPEG
    heatmap::plotParam<T> jpeg_Param;
    jpeg_Param.maxValue = Thot;
    jpeg_Param.minValue = Tcold;
    heatmap::write(planeReduction, iT, jpeg_Param);

    /* BlockLatticeReduction3D<T, TDESCRIPTOR> planeReduction(temperature);
     BlockGifWriter<T> gifWriter;
     gifWriter.write(planeReduction, iT, "temperature");*/
  }

}


int main(int argc, char *argv[])
{

  /// === 1st Step: Initialization ===
  OstreamManager clout(std::cout,"main");
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  fstream f;
  if (argc >= 2) {
    N = atoi(argv[1]);
  }
  if (argc == 3) {
    tau = atof(argv[2]);
  }

  ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const converter(
    (T) 1.0 / N, // physDeltaX
    (T) 1.0 / N * 1.0 / 1e-3 * (tau - 0.5) / 3 / N, // physDeltaT
    (T) 1.0, // charPhysLength
    (T) sqrt( 9.81 * Ra * 1e-3 * 1e-3 / Pr / 9.81 / (Thot - Tcold) / pow(1.0, 3) * (Thot - Tcold) * 1.0 ), // charPhysVelocity
    (T) 1e-3, // physViscosity
    (T) 1.0, // physDensity
    (T) 0.03, // physThermalConductivity
    (T) Pr * 0.03 / 1e-3 / 1.0, // physSpecificHeatCapacity
    (T) Ra * 1e-3 * 1e-3 / Pr / 9.81 / (Thot - Tcold) / pow(1.0, 3), // physThermalExpansionCoefficient
    (T) Tcold, // charPhysLowTemperature
    (T) Thot // charPhysHighTemperature
  );
  converter.print();

  /// === 2nd Step: Prepare Geometry ===
  lz = converter.getPhysDeltaX() * 3.;      // depth of the channel
  std::vector<T> extend(3,T());
  extend[0] = lx;
  extend[1] = ly;
  extend[2] = lz;
  std::vector<T> origin(3,T());
  IndicatorCuboid3D<T> cuboid(extend, origin);

  /// Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 7;
#endif
  CuboidGeometry3D<T> cuboidGeometry(cuboid, converter.getPhysDeltaX(), noOfCuboids);
  cuboidGeometry.setPeriodicity(true,false, true);

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
    instances::getAdvectionDiffusionBulkMomenta<T,TDESCRIPTOR>()
  );

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
  // This coupling must be necessarily be put on the Navier-Stokes lattice!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

  /*  int nx = converter.numCells(lx) + 1;
    int ny = converter.numCells(ly) + 1;
    int nz = converter.numCells(lz) + 1;
  */
  std::vector<T> dir{0.0, 1.0, 0.0};

  T boussinesqForcePrefactor = 9.81 / converter.getConversionFactorVelocity() * converter.getConversionFactorTime() *
    converter.getCharPhysTemperatureDifference() * converter.getPhysThermalExpansionCoefficient();

  NavierStokesAdvectionDiffusionCouplingGenerator3D<T,NSDESCRIPTOR>
  coupling(0, converter.getLatticeLength(lx), 0, converter.getLatticeLength(ly), 0, converter.getLatticeLength(lz),
           boussinesqForcePrefactor, converter.getLatticeTemperature(Tcold), 1., dir);

  NSlattice.addLatticeCoupling(coupling, ADlattice);

  prepareLattice(converter,
                 NSlattice, ADlattice,
                 NSbulkDynamics, TbulkDynamics,
                 NSboundaryCondition, TboundaryCondition, superGeometry );


  /// === 4th Step: Main Loop with Timer ===
  Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  util::ValueTracer<T> converge(converter.getLatticeTime(1.0),epsilon);
  for (int iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {

    if (converge.hasConverged()) {
      clout << "Simulation converged." << endl;
      getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, converge.hasConverged());
      break;
    }

    /// === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(converter, NSlattice, ADlattice, iT, superGeometry);

    /// === 6th Step: Collide and Stream Execution ===

    NSlattice.collideAndStream();

    NSlattice.executeCoupling();
    ADlattice.collideAndStream();

    /// === 7th Step: Computation and Output of the Results ===
    getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, converge.hasConverged());
    converge.takeValue(ADlattice.getStatistics().getAverageEnergy());
  }

  timer.stop();
  timer.printSummary();
}
