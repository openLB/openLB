/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn
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

#ifndef SUPERPARTICLESYSTEM_3D_HH
#define SUPERPARTICLESYSTEM_3D_HH

#define shadows

#include <utility>
#include <string>
#include <array>
#include <iomanip>
#include <ios>
#include <iostream>
#include <random>
#include "io/fileName.h"
#include "superParticleSystem3D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
SuperParticleSystem3D<T, PARTICLETYPE>::SuperParticleSystem3D(
  CuboidGeometry3D<T>& cuboidGeometry, LoadBalancer<T>& loadBalancer,
  SuperGeometry3D<T>& superGeometry)
  : SuperStructure3D<T>(cuboidGeometry, loadBalancer,
                        superGeometry.getOverlap()),
    clout(std::cout, "SuperParticleSystem3d"),
    _superGeometry(superGeometry),
    _overlap(0)
{
  init();
}

template<typename T, template<typename U> class PARTICLETYPE>
SuperParticleSystem3D<T, PARTICLETYPE>::SuperParticleSystem3D(
  SuperGeometry3D<T>& superGeometry)
  : SuperStructure3D<T>(superGeometry.getCuboidGeometry(),
                        superGeometry.getLoadBalancer(),
                        superGeometry.getOverlap()),
    clout(std::cout, "SuperParticleSystem3d"),
    _superGeometry(superGeometry),
    _overlap(0)
{
  init();
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::init()
{
  int rank = 0;
  int size = 0;

  _stopSorting = 0;

#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
  size = singleton::mpi().getSize();
#endif
  for (int i = 0; i < this->_cuboidGeometry.getNc(); ++i) {
    //clout << i << std::endl;
    if (this->_loadBalancer.rank(i) == rank) {
      auto dummy = new ParticleSystem3D<T, PARTICLETYPE>(i, _superGeometry);
      this->_cuboidGeometry.get(i).getOrigin().toStdVector()[0];
      std::vector<T> physPos = this->_cuboidGeometry.get(i).getOrigin().toStdVector();
      std::vector<T> physExtend(3, 0);
      T physR = this->_cuboidGeometry.get(i).getDeltaR();
      for (int j = 0; j < 3; j++) {
        physPos[j] -= .5 * physR;
        physExtend[j] = (this->_cuboidGeometry.get(i).getExtend()[j] + 1)
                        * physR;
      }
      dummy->setPosExt(physPos, physExtend);
      _pSystems.push_back(dummy);
    }
  }
  //singleton::exit(0);

#ifdef PARALLEL_MODE_MPI
  for (int i = 0; i < this->_cuboidGeometry.getNc(); ++i) {
    if (this->_loadBalancer.rank(i) == rank) {
      std::vector<int> dummy;
      this->getCuboidGeometry().getNeighbourhood(i, dummy, 3);
      _rankNeighbours.insert(_rankNeighbours.end(), dummy.begin(), dummy.end());
      _cuboidNeighbours.push_back(dummy);
    }
  }
  for (auto& N : _rankNeighbours) {
    N = this->_loadBalancer.rank(N);
  }
#endif

  /* Ein jeder ist sein eigener Nachbar*/
  if (rank == 0) {
    _rankNeighbours.push_back(size - 1);
  }
  if (rank == size - 1) {
    _rankNeighbours.push_back(0);
  }
  _rankNeighbours.push_back(rank);
  _rankNeighbours.sort();
  _rankNeighbours.unique();
}

template<typename T, template<typename U> class PARTICLETYPE>
SuperParticleSystem3D<T, PARTICLETYPE>::SuperParticleSystem3D(
  SuperParticleSystem3D<T, PARTICLETYPE>& spSys)
  : SuperStructure3D<T>(spSys._cuboidGeometry, spSys._loadBalancer,
                        int(spSys._overlap)),
    clout(std::cout, "SuperParticleSystem3d"),
    _pSystems(spSys._pSystems),
    _superGeometry(spSys._superGeometry),
    _rankNeighbours(spSys._rankNeighbours),
    _overlap(spSys._overlap)
{
}

template<typename T, template<typename U> class PARTICLETYPE>
SuperParticleSystem3D<T, PARTICLETYPE>::SuperParticleSystem3D(
  SuperParticleSystem3D<T, PARTICLETYPE> const& spSys)
  : SuperStructure3D<T>(spSys._cuboidGeometry, spSys._loadBalancer,
                        int(spSys._overlap)),
    clout(std::cout, "SuperParticleSystem3d"),
    _pSystems(spSys._pSystems),
    _superGeometry(spSys._superGeometry),
    _rankNeighbours(spSys._rankNeighbours),
    _overlap(spSys._overlap)
{
}

template<typename T, template<typename U> class PARTICLETYPE>
SuperParticleSystem3D<T, PARTICLETYPE>::SuperParticleSystem3D(
  SuperParticleSystem3D<T, PARTICLETYPE> && spSys)
  : SuperStructure3D<T>(spSys._cuboidGeometry, spSys._loadBalancer,
                        int(spSys._overlap)),
    clout(std::cout, "SuperParticleSystem3d"),
    _superGeometry(spSys._superGeometry),
    _rankNeighbours(spSys._rankNeighbours),
    _overlap(spSys._overlap)
{
  _pSystems = std::move(spSys._pSystems);
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::print()
{
  int no = globalNumOfParticles();
  int active = globalNumOfActiveParticles();
  clout << "activeParticles= " << active << " (" << no << ") " << std::endl;
  //  cout << "[SuperParticleSystem3D] " << _pSystems.size()
  //      << " pSystems on rank " << singleton::mpi().getRank() << "\n";
  //  for (auto pS : _pSystems) {
  //    cout << pS->_particles.size() << " ";
  //  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::printDeep(std::string message)
{
  clout << "=========================================================================================================================" << std::endl;
  clout << "printDeep diagnostic tool" << message << std::endl;
  clout << "=========================================================================================================================" << std::endl;
  for (auto pS : _pSystems) {
    pS->printDeep ( std::to_string(_pSystems.size()) + " pSystems on rank " + std::to_string(singleton::mpi().getRank()) + ":  " );
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::print(std::list<int> mat)
{
  std::list<int>::iterator _matIter;
  int no = globalNumOfParticles();
  clout << "globalNumOfParticles=" << no;
  int active = globalNumOfActiveParticles();
  clout << "; activeParticles=" << active;
  for (_matIter = mat.begin(); _matIter != mat.end(); _matIter++) {
    clout << "; material" << *_matIter << "=" << countMaterial(*_matIter);
  }
  clout << std::endl;
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::captureEscapeRate(
  std::list<int> mat)
{
  std::list<int>::iterator _matIter;
  T sum = T();
  // capture rate
  for (_matIter = mat.begin(); _matIter != mat.end(); _matIter++) {
    sum += (T) countMaterial(*_matIter);
  }
  clout << "captureRate=" << 1. - sum / globalNumOfParticles()
        << "; escapeRate=" << sum / globalNumOfParticles() << std::endl;
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::diffEscapeRate(
  std::list<int> mat, int& globalPSum, int& pSumOutlet, T& diffEscapeRate, T& maxDiffEscapeRate,
  int iT, int iTConsole, T genPartPerTimeStep)
{
  std::list<int>::iterator _matIter;
  T pSumOutletNew = T();
  T maxDiffEscapeRateNew = maxDiffEscapeRate;
  T diffERtmp = T(); // temporal differencial escape rate

  // count particles at outlet
  for (_matIter = mat.begin(); _matIter != mat.end(); _matIter++) {
    pSumOutletNew += (T) countMaterial(*_matIter);
  }

  // calculate diff. escape rate
  if (globalNumOfParticles() > globalPSum) {
    T diffERtmp = (T) (pSumOutletNew - pSumOutlet) / (globalNumOfParticles() - globalPSum);
    diffEscapeRate += diffERtmp;
  } else {
    if (genPartPerTimeStep != 0.) {
      diffERtmp = (T) (pSumOutletNew - pSumOutlet) / (genPartPerTimeStep);
      diffEscapeRate += diffERtmp;
    }
  }
  // calculate max. diff. escape rate
  if (diffERtmp > maxDiffEscapeRateNew) {
    maxDiffEscapeRateNew = diffERtmp;
  }
  // console output
  if (iT % iTConsole == 0) {
    diffEscapeRate /= iTConsole;
    if (globalNumOfParticles() > globalPSum) {
      clout << "diffEscapeRate = " << diffEscapeRate << std::endl;
    } else {
      if (genPartPerTimeStep != 0.) {
        clout << "no particles in feedstream, continue calculation of diffEscapeRate with theoretical "
              << genPartPerTimeStep  << " generated particles per phys. time step"
              << " diffEscapeRate = " << diffEscapeRate << std::endl;
      } else {
        clout << "no particles in feedstream, calculation of diffEscapeRate not possible" << std::endl;
      }
    }
    if (maxDiffEscapeRateNew > maxDiffEscapeRate) {
      clout << "maxDiffEscapeRate = " << maxDiffEscapeRateNew << std::endl;
    }
    diffEscapeRate = 0. ;
  }
  maxDiffEscapeRate = maxDiffEscapeRateNew;
  pSumOutlet = pSumOutletNew;
  globalPSum = globalNumOfParticles();
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::diffEscapeRate(
  std::list<int> mat, int& globalPSum, int& pSumOutlet, T& diffEscapeRate, T& maxDiffEscapeRate,
  int iT, int iTConsole, T genPartPerTimeStep,
  T& avDiffEscapeRate, T latticeTimeStart, T latticeTimeEnd)
{
  std::list<int>::iterator _matIter;
  T pSumOutletNew = T();
  T maxDiffEscapeRateNew = maxDiffEscapeRate;
  T avDiffEscapeRateNew = T();

  // count particle at outlet
  for (_matIter = mat.begin(); _matIter != mat.end(); _matIter++) {
    pSumOutletNew += (T) countMaterial(*_matIter);
  }
  // calculate diff. escape rate
  if (globalNumOfParticles() > globalPSum) {
    avDiffEscapeRateNew = (T) (pSumOutletNew - pSumOutlet) / (globalNumOfParticles() - globalPSum);
    diffEscapeRate += avDiffEscapeRateNew;
  } else {
    if (genPartPerTimeStep != 0.) {
      avDiffEscapeRateNew = (T) (pSumOutletNew - pSumOutlet) / (genPartPerTimeStep);
      diffEscapeRate += avDiffEscapeRateNew;
    }
  }
  // calculate max. diff. escape rate
  if (avDiffEscapeRateNew > maxDiffEscapeRateNew) {
    maxDiffEscapeRateNew = avDiffEscapeRateNew;
  }
  // calculate average diff. escape rate between tStart and tEnd
  if ((iT >= latticeTimeStart) && (iT <= latticeTimeEnd)) {
    avDiffEscapeRate += avDiffEscapeRateNew;
  }
  if (iT == latticeTimeEnd) {
    avDiffEscapeRate /= (latticeTimeEnd - latticeTimeStart);
    clout << "average diffEscapeRate between t = " << latticeTimeStart << "s and t = " << latticeTimeEnd << "s : "
          << avDiffEscapeRate << std::endl;
  }
  // console output
  if (iT % iTConsole == 0) {
    diffEscapeRate /= iTConsole;
    if (globalNumOfParticles() > globalPSum) {
      clout << "diffEscapeRate = " << diffEscapeRate << std::endl;
    } else {
      if (genPartPerTimeStep != 0.) {
        clout << "no particles in feedstream, continue calculation of diffEscapeRate with theoretical "
              << genPartPerTimeStep  << " generated particles per phys. time step" << std::endl;
        clout << "diffEscapeRate = " << diffEscapeRate << std::endl;
      } else {
        clout << "no particles in feedstream, calculation of diffEscapeRate not possible" << std::endl;
      }
    }
    if (maxDiffEscapeRateNew > maxDiffEscapeRate) {
      clout << "maxDiffEscapeRate = " << maxDiffEscapeRateNew << std::endl;
    }
    diffEscapeRate = 0. ;
  }
  maxDiffEscapeRate = maxDiffEscapeRateNew;
  pSumOutlet = pSumOutletNew;
  globalPSum = globalNumOfParticles();
}

// Output_txt.file for tracer particles (particles that have an id != 0
/* _filename is filename of output txt-file
 * _file is contact to output file
 *  */
template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::getOutput(std::string _filename,
    int it, T conversionFactorTime, unsigned short _properties )
{

  enum particleProperties
    : unsigned short {position = 1,
                    velocity = 2,
                    radius = 4,
                    mass = 8,
                    force = 16,
                    storeForce = 32,
                   };

  /// initialize Parameters for the desired ParticleOutput:
  /// ID, Radius, Position, Velocity, Force, storeForces
  std::vector < T > id, rad;
  std::vector<std::array<T, 3>> pos, vel, forces, storeForces;
  /// initialize iteration Parameters for the loops
  int k, j, numOfTracerParticles = globalNumOfTracerParticles();
  int m = 0;
  /// setPrecision for the decimal places in the txtOutputFile
  std::setprecision(9);
  /// set size of list in correlation to the number of TracerParticles
  for (k = 0; k < numOfTracerParticles; k++) {
    id.push_back(T());
    rad.push_back(T());
    pos.push_back( { T(), T(), T() });
    vel.push_back( { T(), T(), T() });
    forces.push_back( { T(), T(), T() });
    storeForces.push_back( { T(), T(), T() });
  }

  k = 0;
  /// get Information about each Particle in each particleSystem
  std::vector<ParticleSystem3D<T, PARTICLETYPE>*> _pSystems = getPSystems();
  for (unsigned int pS = 0; pS < _pSystems.size(); ++pS) {
    std::deque<PARTICLETYPE<T>*> particles =
      _pSystems[pS]->getParticlesPointer();
    for (auto p : particles) {
      if (p->getID() != 0) {
        /// check whether it is a tracerParticle(!=0) or not (==0)
        id[k] = p->getID();
        rad[k] = p->getRad();
        for (j = 0; j < 3; j++) { /// loop for each coordinate [x,y,z]
          pos[k][j] = p->getPos()[j];
          vel[k][j] = p->getVel()[j];
          forces[k][j] = p->getForce()[j];
          storeForces[k][j] = p->getStoreForce()[j];
        }
       // p->resetStoreForce();
        k++;
        p++;
      }
    }
  }

#ifdef PARALLEL_MODE_MPI
  /// prepare data vectors for parallel computing (mpi)
  for (m = 0; m < numOfTracerParticles; m++) {
    singleton::mpi().reduceAndBcast(id[m], MPI_SUM);
    singleton::mpi().reduceAndBcast(rad[m], MPI_SUM);
    for (j = 0; j < 3; j++) { /// loop for each coordinate [x,y,z]
      singleton::mpi().reduceAndBcast(pos[m][j], MPI_SUM);
      singleton::mpi().reduceAndBcast(vel[m][j], MPI_SUM);
      singleton::mpi().reduceAndBcast(forces[m][j], MPI_SUM);
      singleton::mpi().reduceAndBcast(storeForces[m][j], MPI_SUM);
    }
  }
#endif

  /// Particle Output to txtFile just on main process ID (==0)
  if (!singleton::mpi().getRank()) {
    std::ofstream _file;

    int i = 0;
    if (it == 0) {
      /// write headers of each column at timeStep zero
      /// after header name, there is the number of the column
      _file.open(_filename, std::ios::out | std::ios::trunc);
      _file << "Timestep" << i + 1 << " " << "physTime" << i + 2;
      i = i + 2;

      for (m = 0; m < numOfTracerParticles; m++) {
        _file << " id" << i + 1;
        i = i + 1;
        _file << " rad" << i + 1;
        i = i + 1;
        if (_properties & particleProperties::position) {
          _file << " pos0_" << i + 1 << " pos1_" << " pos2_" << i + 3;
          i = i + 3;
        }
        if (_properties & particleProperties::force) {
          _file << " forc0_" << i + 1 << " forc1_" << i + 2 << " forc2_"
                << i + 3;
          i = i + 3;
        }
        if (_properties & particleProperties::velocity) {
          _file << " vel0_" << i + 1 << " vel1_" << i + 2 << " vel2_" << i + 3;
          i = i + 3;
        }
        if (_properties & particleProperties::storeForce) {
          _file << " hforc0_" << i + 1 << " hforc1_" << i + 2 << " hforc2_"
                << i + 3;
          i = i + 3;
        }
      }
      _file << "\n";
      _file.close();
    }

    if (it >= 0) {
      /// write the results in regart to the above defined headers
      _file.open(_filename, std::ios::out | std::ios::app);
      _file << it << " " << conversionFactorTime*it << " ";

      for (m = 0; m < numOfTracerParticles; m++) {
        _file << id[m] << " ";
        _file << rad[m] << " ";
        if (_properties & particleProperties::position) {
          _file << pos[m][0] << " " << pos[m][1] << " " << pos[m][2] << " ";
        }
        if (_properties & particleProperties::force) {
          _file << forces[m][0] << " " << forces[m][1] << " " << forces[m][2]
                << " ";
        }
        if (_properties & particleProperties::velocity) {
          _file << vel[m][0] << " " << vel[m][1] << " " << vel[m][2] << " ";
        }
        if (_properties & particleProperties::storeForce) {
          _file << storeForces[m][0] << " " << storeForces[m][1] << " "
                << storeForces[m][2] << " ";
        }
      }
    }
    _file << "\n";
    _file.close();
  }
}


template<typename T, template<typename U> class PARTICLETYPE>
std::vector<ParticleSystem3D<T, PARTICLETYPE>*>& SuperParticleSystem3D<T,
    PARTICLETYPE>::getPSystems()    //+*
{
  return _pSystems;
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::simulate(T dT, bool scale)
{

  //  for (auto pS : _pSystems) {
  //      pS->velocityVerlet1(dT);
  //      if (_overlap > 0) {
  //        pS->makeTree();
  //        cout << "makeTree" << endl;
  //      }
  //  }
  //  updateParticleDistribution();
  //  for (auto pS : _pSystems) {
  //      pS->velocityVerlet2(dT);
  //  }

  //  updateParticleDistribution();
  //  for (auto pS : _pSystems) {
  //    pS->computeForce();
  //    pS->predictorCorrector1(dT);
  //  }
  //
  //  updateParticleDistribution();
  //  for (auto pS : _pSystems) {
  //    pS->computeForce();
  //    pS->predictorCorrector2(dT);
  //  }

  //  for (auto pS : _pSystems) {
  //    pS->rungeKutta4_1(dT);
  //  }
  //  updateParticleDistribution();
  //  for (auto pS : _pSystems) {
  //    pS->rungeKutta4_2(dT);
  //  }
  //  updateParticleDistribution();
  //  for (auto pS : _pSystems) {
  //    pS->rungeKutta4_3(dT);
  //  }
  //  updateParticleDistribution();
  //  for (auto pS : _pSystems) {
  //    pS->rungeKutta4_4(dT);
  //  }
  //  updateParticleDistribution();
  for (auto pS : _pSystems) {
    time_t delta = clock();
    pS->_contactDetection->sort();
    _stopSorting += clock() - delta;
    pS->simulate(dT, scale);
    pS->computeBoundary();

  }
  updateParticleDistribution();
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::simulateWithTwoWayCoupling_Mathias ( T dT,
                                    ForwardCouplingModel<T,PARTICLETYPE>& forwardCoupling,
                                    BackCouplingModel<T,PARTICLETYPE>& backCoupling,
                                    int material, int subSteps, bool resetExternalField, bool scale )
{
  // reset external field
  if (resetExternalField)
    backCoupling.resetExternalField(material);

  for (auto pS : _pSystems) {
    time_t delta = clock();
    pS->_contactDetection->sort();
    _stopSorting += clock() - delta;
    pS->simulateWithTwoWayCoupling_Mathias(dT, forwardCoupling, backCoupling, material, subSteps, scale);
    pS->computeBoundary();

  }
  updateParticleDistribution();
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::simulateWithTwoWayCoupling_Davide ( T dT,
                                    ForwardCouplingModel<T,PARTICLETYPE>& forwardCoupling,
                                    BackCouplingModel<T,PARTICLETYPE>& backCoupling,
                                    int material, int subSteps, bool resetExternalField, bool scale )
{
  // reset external field
  if (resetExternalField)
    backCoupling.resetExternalField(material);

  for (auto pS : _pSystems) {
    time_t delta = clock();
    pS->_contactDetection->sort();
    _stopSorting += clock() - delta;
    pS->simulateWithTwoWayCoupling_Davide(dT, forwardCoupling, backCoupling, material, subSteps, scale);
    pS->computeBoundary();
  }
  updateParticleDistribution();
}

// multiple collision models
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::simulate(double dT, std::set<int> sActivityOfParticle, bool scale)
{
  for (auto pS : _pSystems) {
    time_t delta = clock();
    if (pS->getIGeometry() == singleton::mpi().getRank()) {
      pS->_contactDetection->sort();
    }
    _stopSorting += clock() - delta;
    if (pS->getIGeometry() == singleton::mpi().getRank()) {
      pS->simulate(dT, sActivityOfParticle, scale);
      pS->computeBoundary();
    }
  }
  updateParticleDistribution();
}

template<>
bool SuperParticleSystem3D<double, MagneticParticle3D>::particleSActivityTest(int sActivity)
{
  for (auto pS : _pSystems) {
    for (auto p : pS->_particles) {
      if (p.getSActivity() == sActivity) {
        return false;
      }
    }
  }
  return true;
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::setOverlap(T overlap)
{
  if ( int(overlap) + 1 > _superGeometry.getOverlap() ) {
    clout << "Overlap of SuperParticleSystem3D should be < overlap "
          "of SuperStructure3D" << std::endl;
    exit(-1);
  }
  _overlap = overlap;
}

template<typename T, template<typename U> class PARTICLETYPE>
T SuperParticleSystem3D<T, PARTICLETYPE>::getOverlap()
{
  return _overlap;
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::numOfPSystems()
{
  return _pSystems.size();
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::globalNumOfParticles()
{
#ifdef PARALLEL_MODE_MPI
  // cout << "return1" << endl;
  int buffer = rankNumOfParticles();
  // cout << "return2" << endl;
  singleton::mpi().reduceAndBcast(buffer, MPI_SUM);
  // cout << "return3" << endl;
  return buffer;
#else
  // cout << "return4" << endl;
  return rankNumOfParticles();
#endif
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::globalNumOfShadowParticles()
{
#ifdef PARALLEL_MODE_MPI
  int buffer = rankNumOfShadowParticles();
  singleton::mpi().reduceAndBcast(buffer, MPI_SUM);
  return buffer;
#else
  return rankNumOfShadowParticles();
#endif
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::rankNumOfParticles()
{
  int num = 0;
  for (auto pS : _pSystems) {
    num += pS->size();
  }
  return num;
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::rankNumOfShadowParticles()
{
  int num = 0;
  for (auto pS : _pSystems) {
    num += pS->_shadowParticles.size();
  }
  return num;
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::rankNumOfActiveParticles()
{
  int num = 0;
  for (auto pS : _pSystems) {
    num += pS->numOfActiveParticles();
  }
  return num;
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::countLocMaterial(int mat)
{
  int num = 0;
  for (auto pS : _pSystems) {
    num += pS->countMaterial(mat);
  }
  return num;
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::countMaterial(int mat)
{
#ifdef PARALLEL_MODE_MPI
  int buffer = countLocMaterial(mat);
  singleton::mpi().reduceAndBcast(buffer, MPI_SUM);
  return buffer;
#else
  return countLocMaterial(mat);
#endif
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::globalNumOfActiveParticles()
{
#ifdef PARALLEL_MODE_MPI
  int buffer = rankNumOfActiveParticles();
  singleton::mpi().reduceAndBcast(buffer, MPI_SUM);
  return buffer;
#else
  return rankNumOfActiveParticles();
#endif
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::globalNumOfTracerParticles()
{
#if PARALLEL_MODE_MPI
  int buffer = rankNumOfTracerParticles();
  singleton::mpi().reduceAndBcast(buffer, MPI_SUM);
  return buffer;
#else
  return rankNumOfTracerParticles();
#endif
}

template<typename T, template<typename U> class PARTICLETYPE>
int SuperParticleSystem3D<T, PARTICLETYPE>::rankNumOfTracerParticles()
{
  int num = 0;
  for (auto pS : _pSystems) {
    std::deque<PARTICLETYPE<T>*> particles = pS->getParticlesPointer();
    int pSNum = 0;
    for (auto p : particles) {
      if (p->getID() != 0) {
        pSNum++;
      }
    }
    num += pSNum;
  }
  return num;
}

// TODO class olb::ParticleSystem3D<double, olb::Particle3D>’ has no member named ‘radiusDistribution
/*
 template<typename T, template<typename U> class PARTICLETYPE>
 std::map<T, int> SuperParticleSystem3D<T, PARTICLETYPE>::rankRadiusDistribution()
 {
 std::map<T, int> distrAll;
 typename std::map<T, int>::iterator ita;
 for (auto pS : _pSystems) {
 std::map<T, int> distrPs;
 distrPs = pS->radiusDistribution();
 for (typename std::map<T, int>::iterator itp = distrPs.begin();
 itp != distrPs.end(); ++itp) {
 T radKeyP = itp->first;
 int countP = itp->second;
 if (distrAll.count(radKeyP) > 0) {
 ita = distrAll.find(radKeyP);
 ita->second = ita->second + countP;
 } else {
 distrAll.insert(std::pair<T, int>(radKeyP, countP));
 }
 }
 }
 return distrAll;
 }
 */

template<typename T, template<typename U> class PARTICLETYPE>
std::vector<int> SuperParticleSystem3D<T, PARTICLETYPE>::numOfForces()
{
  std::vector<int> dummy;
  for (auto pS : _pSystems) {
    dummy.push_back(pS->numOfForces());
  }
  return dummy;
}

template<typename T, template<typename U> class PARTICLETYPE>
template<typename DESCRIPTOR>
void SuperParticleSystem3D<T, PARTICLETYPE>::setVelToFluidVel(
  SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& fVel)
{
  for (auto pS : _pSystems) {
    pS->setVelToFluidVel(fVel);
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::setVelToAnalyticalVel(
  AnalyticalConst3D<T, T>& aVel)
{
  for (auto pS : _pSystems) {
    pS->setVelToAnalyticalVel(aVel);
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
bool SuperParticleSystem3D<T, PARTICLETYPE>::findCuboid(PARTICLETYPE<T> &p)
{
  return findCuboid(p, int(_overlap));
}

template<typename T, template<typename U> class PARTICLETYPE>
bool SuperParticleSystem3D<T, PARTICLETYPE>::findCuboid(PARTICLETYPE<T> &p,
    int overlap)
{
  int C = this->_cuboidGeometry.get_iC(p.getPos()[0], p.getPos()[1],
                                       p.getPos()[2], overlap);
  if (C != this->_cuboidGeometry.getNc()) {
    p.setCuboid(C);
    return 1;
  } else {
    clout << "Lost Particle! Pos: " << p.getPos()[0] << " " << p.getPos()[1]
          << " " << p.getPos()[2] << " Vel: " << p.getVel()[0] << " "
          << p.getVel()[1] << " " << p.getVel()[2] << " Cuboid: " << C << std::endl;
    p.setActive(false);
    return 0;
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
bool SuperParticleSystem3D<T, PARTICLETYPE>::checkCuboid(PARTICLETYPE<T>& p,
    T overlap)
{
  return checkCuboid(p, overlap, p.getCuboid());
}

template<typename T, template<typename U> class PARTICLETYPE>
bool SuperParticleSystem3D<T, PARTICLETYPE>::checkCuboid(PARTICLETYPE<T>& p,
    T overlap, int iC)
{
  // Checks whether particle is contained in the cuboid
  // extended with an layer of size overlap*delta
  return this->_cuboidGeometry.get(iC).physCheckPoint(p.getPos()[0],
         p.getPos()[1],
         p.getPos()[2], overlap);
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::updateParticleDistribution()
{
  // Find particles on wrong cuboid, store in relocate and delete
  // maps particles to their new rank
  // std::multimap<int, PARTICLETYPE<T> > relocate;
  _relocate.clear();
#ifdef shadows
  _relocateShadow.clear();
#endif
  for (unsigned int pS = 0; pS < _pSystems.size(); ++pS) {

    _pSystems[pS]->_shadowParticles.clear();
    auto par = _pSystems[pS]->_particles.begin();

    while (par != _pSystems[pS]->_particles.end()) {

      //Check if particle is still in its cuboid
      if (checkCuboid(*par, 0)) {
#ifdef shadows
        // Check if its inside boundary layer of its cuboid
        if (!(checkCuboid(*par, -_overlap))) {

          std::set<int> sendTo;
          // Run through all cuboids to search the one that
          // shares its boundary with the cuboid containing the particle
          for (int iC = 0; iC < this->_cuboidGeometry.getNc(); iC++) {
            // Check if iC is neighbor cuboid in which overlap the particle is
            // located
            if (par->getCuboid() != iC && checkCuboid(*par, _overlap, iC)) {
              // rank is the processing thread of the global cuboid number iC
              int rank = this->_loadBalancer.rank(iC);
              // insert rank into sendTo, if sendTo does not already contain
              // rank for the actual particle
              if (!sendTo.count(rank)) {
                _relocateShadow.insert(std::make_pair(rank, (*par)));
                sendTo.insert(rank);
              }
            }
          }
        }
#endif
        //If not --> find new cuboid
      } else {
        // check to which cuboid particle is located
        // and give new cuboid to par
        findCuboid(*par, 0);
        // gives particle the new cuboid number where it is now located to
        _relocate.insert(
          std::make_pair(this->_loadBalancer.rank(par->getCuboid()), (*par)));
#ifdef shadows
        // If the new particle position is in boundary layer
        // If yes --> check if its inside boundary layer
        if (!(checkCuboid(*par, -_overlap))) {
          // yes, particle is in boundary layer
          std::set<int> sendTo;
          for (int iC = 0; iC < this->_cuboidGeometry.getNc(); iC++) {
            // checks if iC is neighbor cuboid in whose overlap the particle
            // is located
            if (par->getCuboid() != iC && checkCuboid(*par, _overlap, iC)) {
              int rank = this->_loadBalancer.rank(iC);
              if (!sendTo.count(rank)) {
                // shadow particle (copy of original) is inserted for iC
                // but just with information of rank, not of cuboid!
                _relocateShadow.insert(std::make_pair(rank, (*par)));
                sendTo.insert(rank);
              }
            }
          }
        }
#endif
        par = _pSystems[pS]->_particles.erase(par);
        par--;
      }
      par++;
    }
  }
  /* Communicate number of Particles per cuboid*/
#ifdef PARALLEL_MODE_MPI
  singleton::MpiNonBlockingHelper mpiNbHelper;
  //mpiNbHelper.allocate(_rankNeighbours.size());

  int k = 0;

  /* Serialize particles */
  // it is: std::map<int, std::vector<double> > _send_buffer
  _send_buffer.clear();

  // fill std::map<int, std::vector<double> > _send_buffer with
  // particle properties of particles in _relocate

  // buffer size equals number of particle properties
  T buffer[PARTICLETYPE<T>::serialPartSize];
  // insert all regular particles of _relocate into _send_buffer
  for (auto rN : _relocate) {
    // second element in rN is (*par)
    // write particle properties into buffer
    rN.second.serialize(buffer);
    // first element in rN is rank
    // it is: std::map<int, std::vector<double> > _send_buffer;
    // enlarge the _send_buffer element (rank), which is of vector type, with
    // new elements of buffer (begin: buffer, end: buffer + serialPartSize)
    _send_buffer[rN.first].insert(_send_buffer[rN.first].end(),
                                  buffer, buffer + PARTICLETYPE<T>::serialPartSize);
  }

  /*Send Particles */
  k = 0;
  // noSends contains number of neighboring cuboids
  int noSends = 0;
  // it is: std::list<int> _rankNeighbours, rank of neighboring cuboids
  for (auto rN : _rankNeighbours) {
    // if there are particles contained in _send_buffer, increase noSends
    if (_send_buffer[rN].size() > 0) {
      ++noSends;
    }
  }
  //  int noSends = _send_buffer.size();
  if (noSends > 0) {
    mpiNbHelper.allocate(noSends);
    //    cout << mpiNbHelper.get_size() << std::endl;
    for (auto rN : _rankNeighbours) {
      if (_send_buffer[rN].size() > 0) {
        singleton::mpi().iSend<double>(&_send_buffer[rN][0],
                                       _relocate.count(rN)*PARTICLETYPE<T>::serialPartSize,
                                       rN, &mpiNbHelper.get_mpiRequest()[k++], 1);
      }
    }
  }

  /*Receive and add particles*/
  singleton::mpi().barrier();
  k = 0;
  int flag = 0;
  MPI_Iprobe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
  if (flag) {
    for (auto rN : _rankNeighbours) {
      MPI_Status status;
      int tmpFlag = 0;
      MPI_Iprobe(rN, 1, MPI_COMM_WORLD, &tmpFlag, &status);
      if (tmpFlag) {
        int number_amount = 0;
        MPI_Get_count(&status, MPI_DOUBLE, &number_amount);
        T recv_buffer[number_amount];
        singleton::mpi().receive(recv_buffer, number_amount, rN, 1);

        for (int iPar = 0; iPar < number_amount; iPar += PARTICLETYPE<T>::serialPartSize) {
          PARTICLETYPE<T> p;
          p.unserialize(&recv_buffer[iPar]);
          if (singleton::mpi().getRank() == this->_loadBalancer.rank(p.getCuboid())) {
            // particle added to pSystem with local cuboid number on actual thread
            _pSystems[this->_loadBalancer.loc(p.getCuboid())]->addParticle(p);
          }
        }
      }
    }
  }
  if (noSends > 0) {
    singleton::mpi().waitAll(mpiNbHelper);
  }
#ifdef shadows
  /**************************************************************************************************************/
  //Same Again for shadowParticles
  mpiNbHelper.allocate(_rankNeighbours.size());
  k = 0;
  /* Serialize particles */
  _send_buffer.clear();
  //  T buffer[PARTICLETYPE<T>::serialPartSize];
  for (auto rN : _relocateShadow) {
    //    std::vector<T> buffer = rN.second.serialize();
    rN.second.serialize(buffer);
    _send_buffer[rN.first].insert(_send_buffer[rN.first].end(),
                                  buffer, buffer + PARTICLETYPE<T>::serialPartSize);
  }

  /*Send Particles */
  k = 0;
  noSends = _send_buffer.size();
  if (noSends > 0) {
    mpiNbHelper.allocate(noSends);
    for (auto rN : _rankNeighbours) {
      if (_send_buffer[rN].size() > 0) {
        singleton::mpi().iSend<double>(&_send_buffer[rN][0],
                                       _relocateShadow.count(rN)*PARTICLETYPE<T>::serialPartSize,
                                       rN, &mpiNbHelper.get_mpiRequest()[k++], 4);
      }
    }
  }
  /*Receive and add particles*/
  singleton::mpi().barrier();
  k = 0;
  flag = 0;
  MPI_Iprobe(MPI_ANY_SOURCE, 4, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
  if (flag) {
    for (auto rN : _rankNeighbours) {
      MPI_Status status;
      int tmpFlag = 0;
      MPI_Iprobe(rN, 4, MPI_COMM_WORLD, &tmpFlag, &status);
      //      cout << "Message from " << rN << " found on " << singleton::mpi().getRank() << std::endl;
      if (tmpFlag) {
        int number_amount = 0;
        MPI_Get_count(&status, MPI_DOUBLE, &number_amount);
        //        cout << "Contains " << number_amount << " infos" << std::endl;
        T recv_buffer[number_amount];
        singleton::mpi().receive(recv_buffer, number_amount, rN, 4);
        for (int iPar = 0; iPar < number_amount; iPar += PARTICLETYPE<T>::serialPartSize) {
          //          std::cout << "Particle unserialized" << std::endl;
          // par contains information of shadow particle
          // par is already shadow particle !!
          PARTICLETYPE<T> p;
          p.unserialize(&recv_buffer[iPar]);
          addShadowParticle(p);
        }
      }
    }
  }
  /* WARNING:
   * Performance warning: mpi barrier added - Asher 04/01/2017
   * MPI hanging on application /apps/asher/particle/simpleParticleExample
   * apparent communication error, tested with mpirun -np 6 when particle is
   * in overlap regions. This could be due to tag value of 4 coinciding with
   * tags used in fluid parallelization and conditional blocking receives.
   * Unaffected threads could continue on to collision step where these buffers
   * are altered/overwritten. When tested with tag of 4000007 no such error arises,
   * implying such a conflict.
   */
  singleton::mpi().barrier();
  if (noSends > 0) {
    singleton::mpi().waitAll(mpiNbHelper);
  }
#endif
  mpiNbHelper.free();
#else
  for (auto& par : _relocate) {
    // _relocate contains make_pair(rank, (*par))
    // therefore par.second is pointer to particle par
    addParticle(par.second);
  }
  for (auto& par : _relocateShadow) {
    // _relocateShadow contains make_pair(rank, (*par))
    // therefore par.second is pointer to particle par
    addShadowParticle(par.second);
  }
#endif
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addParticle(PARTICLETYPE<T>& p)
{
  if (findCuboid(p)) {
#ifdef PARALLEL_MODE_MPI
    if (singleton::mpi().getRank() == this->_loadBalancer.rank(p.getCuboid())) {
      _pSystems[this->_loadBalancer.loc(p.getCuboid())]->addParticle(p);
    } else {
//      clout << "Particle not found on Cuboid: " << p.getCuboid() << std::endl;
//      clout << "Ppos: " << p.getPos()[0] << " " << p.getPos()[1] << " " << p.getPos()[2]<< std::endl;
    }
#else
    _pSystems[this->_loadBalancer.loc(p.getCuboid())]->addParticle(p);
#endif
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addParticle(
  IndicatorF3D<T>& ind, T mas, T rad, int no, std::vector<T> vel)
{
  srand(clock());
  //  srand(rand());
  std::vector<T> pos(3, 0.);
  bool indic[1] = { false };

  no += globalNumOfParticles();
  while (globalNumOfParticles() < no) {
    pos[0] = ind.getMin()[0]
             + (T) (rand() % 100000) / 100000. * (ind.getMax()[0] - ind.getMin()[0]);
    pos[1] = ind.getMin()[1]
             + (T) (rand() % 100000) / 100000. * (ind.getMax()[1] - ind.getMin()[1]);
    pos[2] = ind.getMin()[2]
             + (T) (rand() % 100000) / 100000. * (ind.getMax()[2] - ind.getMin()[2]);

#ifdef PARALLEL_MODE_MPI
    singleton::mpi().bCast(&*pos.begin(), 3);
#endif

    int x0, y0, z0, C;
    std::vector<int> locLat(4, 0);
    if (this->_cuboidGeometry.getFloorLatticeR(pos, locLat)) {
      C = locLat[0];
      if (this->_loadBalancer.rank(C) == singleton::mpi().getRank()) {
        x0 = locLat[1];
        y0 = locLat[2];
        z0 = locLat[3];
        if (_superGeometry.get(C, x0, y0, z0) == 1
            && _superGeometry.get(C, x0, y0 + 1, z0) == 1
            && _superGeometry.get(C, x0, y0, z0 + 1) == 1
            && _superGeometry.get(C, x0, y0 + 1, z0 + 1) == 1
            && _superGeometry.get(C, x0 + 1, y0, z0) == 1
            && _superGeometry.get(C, x0 + 1, y0 + 1, z0) == 1
            && _superGeometry.get(C, x0 + 1, y0, z0 + 1) == 1
            && _superGeometry.get(C, x0 + 1, y0 + 1, z0 + 1) == 1
            && ind(indic, &pos[0])) {
          PARTICLETYPE<T> p(pos, vel, mas, rad);
          addParticle(p);
        }
      }
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addParticle(std::set<int>
    material, int no, T mas,
    T rad, std::vector<T> vel)
{
  srand(time(nullptr));
  std::vector<T> pos(3, 0.), min(3, std::numeric_limits<T>::max()),
      max(3, std::numeric_limits<T>::min());
  std::vector<T> tmpMin(3, 0.), tmpMax(3, 0.);
  std::set<int>::iterator it = material.begin();
  for (; it != material.end(); ++it) {
    tmpMin = _superGeometry.getStatistics().getMinPhysR(*it);
    tmpMax = _superGeometry.getStatistics().getMaxPhysR(*it);
    max[0] = std::max(max[0], tmpMax[0]);
    max[1] = std::max(max[1], tmpMax[1]);
    max[2] = std::max(max[2], tmpMax[2]);
    min[0] = std::min(min[0], tmpMin[0]);
    min[1] = std::min(min[1], tmpMin[1]);
    min[2] = std::min(min[2], tmpMin[2]);
  }
//  cout << "Min: " << min[0] << " " << min[1] << " " << min[2] << std::endl;
//  cout << "Max: " << max[0] << " " << max[1] << " " << max[2] << std::endl;
  for (int i = 0; i < 3; i++) {
    min[i] -= this->_cuboidGeometry.get(0).getDeltaR();
    max[i] += 2 * this->_cuboidGeometry.get(0).getDeltaR();
  }
  no += globalNumOfParticles();
  while (globalNumOfParticles() < no) {
    pos[0] = min[0] + (T) (rand() % 100000) / 100000. * (max[0] - min[0]);
    pos[1] = min[1] + (T) (rand() % 100000) / 100000. * (max[1] - min[1]);
    pos[2] = min[2] + (T) (rand() % 100000) / 100000. * (max[2] - min[2]);

#ifdef PARALLEL_MODE_MPI
    singleton::mpi().bCast(&*pos.begin(), 3);
#endif

    std::vector<int> locLat(4, 0);
    if (this->_cuboidGeometry.getFloorLatticeR(pos, locLat)) {
      if (this->_loadBalancer.rank(locLat[0]) == singleton::mpi().getRank()) {
        if (material.find(_superGeometry.get(locLat[0], locLat[1], locLat[2], locLat[3]))
            != material.end()
            && material.find(_superGeometry.get(locLat[0], locLat[1], locLat[2] + 1,
                                                locLat[3])) != material.end()
            && material.find(_superGeometry.get(locLat[0], locLat[1], locLat[2],
                                                locLat[3] + 1)) != material.end()
            && material.find(_superGeometry.get(locLat[0], locLat[1], locLat[2] + 1,
                                                locLat[3] + 1)) != material.end()
            && material.find(_superGeometry.get(locLat[0], locLat[1] + 1, locLat[2],
                                                locLat[3])) != material.end()
            && material.find(_superGeometry.get(locLat[0], locLat[1] + 1, locLat[2] + 1,
                                                locLat[3])) != material.end()
            && material.find(_superGeometry.get(locLat[0], locLat[1] + 1, locLat[2],
                                                locLat[3] + 1)) != material.end()
            && material.find(_superGeometry.get(locLat[0], locLat[1] + 1, locLat[2] + 1,
                                                locLat[3] + 1)) != material.end()) {
          PARTICLETYPE<T> p(pos, vel, mas, rad);
          addParticle(p);
        }
      }
    }
  }
  singleton::mpi().barrier();
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addParticle(
  IndicatorF3D<T>& ind, std::set<int> material, T mas, T rad, int no, std::vector<T> vel)
{
  srand(clock());
  std::vector<T> pos(3, 0.);
  bool indic[1] = { false };

  no += globalNumOfParticles();
  while (globalNumOfParticles() < no) {
    pos[0] = ind.getMin()[0]
             + (T) (rand() % 100000) / 100000. * (ind.getMax()[0] - ind.getMin()[0]);
    pos[1] = ind.getMin()[1]
             + (T) (rand() % 100000) / 100000. * (ind.getMax()[1] - ind.getMin()[1]);
    pos[2] = ind.getMin()[2]
             + (T) (rand() % 100000) / 100000. * (ind.getMax()[2] - ind.getMin()[2]);

#ifdef PARALLEL_MODE_MPI
    singleton::mpi().bCast(&*pos.begin(), 3);
#endif

    int x0, y0, z0;
    std::vector<int> locLat(4, 0);
    if (this->_cuboidGeometry.getFloorLatticeR(pos, locLat)) {
      if (this->_loadBalancer.rank(locLat[0]) == singleton::mpi().getRank()) {
        x0 = locLat[1];
        y0 = locLat[2];
        z0 = locLat[3];
        if (_superGeometry.get(locLat[0], x0, y0, z0) == 1
            && _superGeometry.get(locLat[0], x0, y0 + 1, z0) == 1
            && _superGeometry.get(locLat[0], x0, y0, z0 + 1) == 1
            && _superGeometry.get(locLat[0], x0, y0 + 1, z0 + 1) == 1
            && _superGeometry.get(locLat[0], x0 + 1, y0, z0) == 1
            && _superGeometry.get(locLat[0], x0 + 1, y0 + 1, z0) == 1
            && _superGeometry.get(locLat[0], x0 + 1, y0, z0 + 1) == 1
            && _superGeometry.get(locLat[0], x0 + 1, y0 + 1, z0 + 1) == 1
            && ind(indic, &pos[0])) {
          if (material.find(
                _superGeometry.get(locLat[0], locLat[1], locLat[2], locLat[3]))
              != material.end()
              && material.find(
                _superGeometry.get(locLat[0], locLat[1], locLat[2] + 1,
                                   locLat[3])) != material.end()
              && material.find(
                _superGeometry.get(locLat[0], locLat[1], locLat[2],
                                   locLat[3] + 1)) != material.end()
              && material.find(
                _superGeometry.get(locLat[0], locLat[1], locLat[2] + 1,
                                   locLat[3] + 1)) != material.end()
              && material.find(
                _superGeometry.get(locLat[0], locLat[1] + 1, locLat[2],
                                   locLat[3])) != material.end()
              && material.find(
                _superGeometry.get(locLat[0], locLat[1] + 1, locLat[2] + 1,
                                   locLat[3])) != material.end()
              && material.find(
                _superGeometry.get(locLat[0], locLat[1] + 1, locLat[2],
                                   locLat[3] + 1)) != material.end()
              && material.find(
                _superGeometry.get(locLat[0], locLat[1] + 1, locLat[2] + 1,
                                   locLat[3] + 1)) != material.end()) {
            PARTICLETYPE<T> p(pos, vel, mas, rad);
            addParticle(p);
          }
        }
      }
    }
  }
}

template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::addParticle(
  IndicatorF3D<double>& ind, double mas, double rad, int no, int id,
  std::vector<double> vel, std::vector<double> dMoment, std::vector<double> aVel, std::vector<double> torque, double magnetisation,
  int sActivity)
{
  std::vector<double> pos(3, 0.);
  bool indic[1] = { false };

  no += globalNumOfParticles();
  while (globalNumOfParticles() < no) {
    pos[0] = ind.getMin()[0]
             + (double) (rand() % 100000) / 100000. * (ind.getMax()[0] - ind.getMin()[0]);
    pos[1] = ind.getMin()[1]
             + (double) (rand() % 100000) / 100000. * (ind.getMax()[1] - ind.getMin()[1]);
    pos[2] = ind.getMin()[2]
             + (double) (rand() % 100000) / 100000. * (ind.getMax()[2] - ind.getMin()[2]);

#ifdef PARALLEL_MODE_MPI
    singleton::mpi().bCast(&*pos.begin(), 3);
#endif

    int x0, y0, z0, C;
    std::vector<int> locLat(4, 0);
    if (this->_cuboidGeometry.getFloorLatticeR(pos, locLat)) {
      C = locLat[0];
      if (this->_loadBalancer.rank(C) == singleton::mpi().getRank()) {
        x0 = locLat[1];
        y0 = locLat[2];
        z0 = locLat[3];
        if (_superGeometry.get(C, x0, y0, z0) == 1
            && _superGeometry.get(C, x0, y0 + 1, z0) == 1
            && _superGeometry.get(C, x0, y0, z0 + 1) == 1
            && _superGeometry.get(C, x0, y0 + 1, z0 + 1) == 1
            && _superGeometry.get(C, x0 + 1, y0, z0) == 1
            && _superGeometry.get(C, x0 + 1, y0 + 1, z0) == 1
            && _superGeometry.get(C, x0 + 1, y0, z0 + 1) == 1
            && _superGeometry.get(C, x0 + 1, y0 + 1, z0 + 1) == 1
            && ind(indic, &pos[0])) {
          MagneticParticle3D<double> p(pos, vel, mas, rad, id, dMoment, aVel, torque, magnetisation, sActivity);
          id++;
          addParticle(p);
        }
      }
    }
  }
}

template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::addParticle(IndicatorF3D<double>& ind,  std::set<int>  material, double mas, double rad, int no, int id,
    std::vector<double> vel, std::vector<double> dMoment, std::vector<double> aVel, std::vector<double> torque, double magnetisation,
    int sActivity)

{
  std::vector<double> pos(3, 0.);
  bool indic[1] = { false };

  no += globalNumOfParticles();
  while (globalNumOfParticles() < no) {
    pos[0] = ind.getMin()[0]
             + (double) (rand() % 100000) / 100000. * (ind.getMax()[0] - ind.getMin()[0]);
    pos[1] = ind.getMin()[1]
             + (double) (rand() % 100000) / 100000. * (ind.getMax()[1] - ind.getMin()[1]);
    pos[2] = ind.getMin()[2]
             + (double) (rand() % 100000) / 100000. * (ind.getMax()[2] - ind.getMin()[2]);

#ifdef PARALLEL_MODE_MPI
    singleton::mpi().bCast(&*pos.begin(), 3);
#endif

    int x0, y0, z0;
    std::vector<int> locLat(4, 0);
    if (this->_cuboidGeometry.getFloorLatticeR(pos, locLat)) {
      if (this->_loadBalancer.rank(locLat[0]) == singleton::mpi().getRank()) {
        x0 = locLat[1];
        y0 = locLat[2];
        z0 = locLat[3];
        if (_superGeometry.get(locLat[0], x0, y0, z0) == 1
            && _superGeometry.get(locLat[0], x0, y0 + 1, z0) == 1
            && _superGeometry.get(locLat[0], x0, y0, z0 + 1) == 1
            && _superGeometry.get(locLat[0], x0, y0 + 1, z0 + 1) == 1
            && _superGeometry.get(locLat[0], x0 + 1, y0, z0) == 1
            && _superGeometry.get(locLat[0], x0 + 1, y0 + 1, z0) == 1
            && _superGeometry.get(locLat[0], x0 + 1, y0, z0 + 1) == 1
            && _superGeometry.get(locLat[0], x0 + 1, y0 + 1, z0 + 1) == 1
            && ind(indic, &pos[0])) {
          if (material.find(
                _superGeometry.get(locLat[0], locLat[1], locLat[2], locLat[3]))
              != material.end()
              && material.find(
                _superGeometry.get(locLat[0], locLat[1], locLat[2] + 1,
                                   locLat[3])) != material.end()
              && material.find(
                _superGeometry.get(locLat[0], locLat[1], locLat[2],
                                   locLat[3] + 1)) != material.end()
              && material.find(
                _superGeometry.get(locLat[0], locLat[1], locLat[2] + 1,
                                   locLat[3] + 1)) != material.end()
              && material.find(
                _superGeometry.get(locLat[0], locLat[1] + 1, locLat[2],
                                   locLat[3])) != material.end()
              && material.find(
                _superGeometry.get(locLat[0], locLat[1] + 1, locLat[2] + 1,
                                   locLat[3])) != material.end()
              && material.find(
                _superGeometry.get(locLat[0], locLat[1] + 1, locLat[2],
                                   locLat[3] + 1)) != material.end()
              && material.find(
                _superGeometry.get(locLat[0], locLat[1] + 1, locLat[2] + 1,
                                   locLat[3] + 1)) != material.end()) {

            MagneticParticle3D<double> p(pos, vel, mas, rad, id, dMoment, aVel, torque, magnetisation, sActivity);
            id++;
            addParticle(p);
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
template<typename DESCRIPTOR>
void SuperParticleSystem3D<T, PARTICLETYPE>::generateParticlesCircleInletMassConcentration(
  IndicatorCircle3D<T>& indicatorCircle, T particleMassConcentration, T charPhysVelocity,
  T conversionFactorTime, SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& getVel,
  PARTICLETYPE<T>& p, std::set<int> material, int iT, T& particlesPerPhyTimeStep,
  std::vector<T>& inletVec, std::deque<std::vector<T>>& posDeq, int deqSize)
{

  std::vector<T> pos(3, 0.);
  std::vector<T> vel(3, 0.);

  PARTICLETYPE<T> pCopy(p);

  // r can be initialized with arbitrary values
  // r is calculated to lie in the same plane as indicatorCircle
  Vector<T, 3> r(inletVec);
  // s is calculated to lie in the same plane as indicatorCircle
  // s is orthogonal to r
  Vector<T, 3> s(0., 0., 0.);

  T fluidVel[3] = {0., 0., 0.};

  // calculation of r in the first step of iteration
  if (iT == 0) {
    T fluxDensity = charPhysVelocity * M_PI * std::pow(indicatorCircle.getRadius(), 2.) ;
    T particleVolume = 4. / 3.* M_PI * std::pow(p.getRad(), 3);
    T particleDensity = p.getMass() / particleVolume;
    T particleVolumeConcentration = particleMassConcentration / particleDensity;
    T particleVolumeFlux = fluxDensity * particleVolumeConcentration;
    T particleFlux = particleVolumeFlux / particleVolume;
    particlesPerPhyTimeStep =  particleFlux * conversionFactorTime;

    bool b = false;
    for (int i = 0; i < 3; i++) {
      if (indicatorCircle.getNormal()[i] == 0.) {
        b = true;
      }
    }
    if (b == true) {
      for (int i = 0; i < 3; i++) {
        if (indicatorCircle.getNormal()[i] == 0.) {
          r[i] = 1.;
        } else {
          r[i] = 0.;
        }
      }
    } else {
      r[0] = -(indicatorCircle.getNormal()[1] + indicatorCircle.getNormal()[2]) / indicatorCircle.getNormal()[0] ;
      r[1] = 1.;
      r[2] = 1.;
    }
    normalize(r) ;
    for (int i = 0; i <= 2; i++) {
      inletVec[i] = r[i];
    }
  }

  // norm of r
  T r_max = indicatorCircle.getRadius();
  T r_min = -1. * r_max;

  s = crossProduct3D(r, indicatorCircle.getNormal());

  // Non-deterministic random number generator
  std::random_device rd;
  // Pseudo-random number engine: Mersenne Twister 19937 generator
  std::mt19937 engine(rd());
  int id = this->globalNumOfParticles();

  while (this->globalNumOfParticles() < (iT + 1) * particlesPerPhyTimeStep) {

gt_mark:
    normalize(r);
    normalize(s);

    // Random number distribution that produces floating-point values according to a uniform distribution
    std::uniform_real_distribution<T> distR(r_min, r_max);

    // r_norm is between r_min and r_max
    T r_norm = distR(engine);
    r *= r_norm ;

    T s_max = std::sqrt(std::pow(indicatorCircle.getRadius(), 2.) - std::pow(r_norm, 2.)) ;
    T s_min = -1. * s_max ;
    std::uniform_real_distribution<T> distS(s_min, s_max);

    // s_norm is between s_min and s_max
    T s_norm = distS(engine);
    s *= s_norm ;

    std::vector<T> posVecTmp(3, 0.);
    posVecTmp[0] = indicatorCircle.getCenter()[0] + r[0] + s[0];
    posVecTmp[1] = indicatorCircle.getCenter()[1] + r[1] + s[1];
    posVecTmp[2] = indicatorCircle.getCenter()[2] + r[2] + s[2];

    for (auto a : posDeq) {
      T dist = std::sqrt(std::pow(a[0] - posVecTmp[0], 2.)
                         + std::pow(a[1] - posVecTmp[1], 2.)
                         + std::pow(a[2] - posVecTmp[2], 2.));
      if (dist <= 3 * p.getRad()) {
        goto gt_mark;
      }
    }

    pos[0] = posVecTmp[0];
    pos[1] = posVecTmp[1];
    pos[2] = posVecTmp[2];

    std::vector<int> latticeRoundedPos(4, 0);

    if (this->_cuboidGeometry.getFloorLatticeR(pos, latticeRoundedPos)) {
      int globCuboid = latticeRoundedPos[0]; // is global cuboid number
      if (this->_loadBalancer.rank(globCuboid) == singleton::mpi().getRank()) {

        if (material.find(_superGeometry.get(latticeRoundedPos[0], latticeRoundedPos[1], latticeRoundedPos[2], latticeRoundedPos[3]))
            != material.end()
            && material.find(_superGeometry.get(latticeRoundedPos[0], latticeRoundedPos[1], latticeRoundedPos[2] + 1,
                                                latticeRoundedPos[3])) != material.end()
            && material.find(_superGeometry.get(latticeRoundedPos[0], latticeRoundedPos[1], latticeRoundedPos[2],
                                                latticeRoundedPos[3] + 1)) != material.end()
            && material.find(_superGeometry.get(latticeRoundedPos[0], latticeRoundedPos[1], latticeRoundedPos[2] + 1,
                                                latticeRoundedPos[3] + 1)) != material.end()
            && material.find(_superGeometry.get(latticeRoundedPos[0], latticeRoundedPos[1] + 1, latticeRoundedPos[2],
                                                latticeRoundedPos[3])) != material.end()
            && material.find(_superGeometry.get(latticeRoundedPos[0], latticeRoundedPos[1] + 1, latticeRoundedPos[2] + 1,
                                                latticeRoundedPos[3])) != material.end()
            && material.find(_superGeometry.get(latticeRoundedPos[0], latticeRoundedPos[1] + 1, latticeRoundedPos[2],
                                                latticeRoundedPos[3] + 1)) != material.end()
            && material.find(_superGeometry.get(latticeRoundedPos[0], latticeRoundedPos[1] + 1, latticeRoundedPos[2] + 1,
                                                latticeRoundedPos[3] + 1)) != material.end()) {
          getVel(fluidVel, &pos[0], globCuboid);
          vel[0] = fluidVel[0];
          vel[1] = fluidVel[1];
          vel[2] = fluidVel[2];

          pCopy.setPos(pos);
          pCopy.setVel(vel);
          pCopy.setID(id);
          this->addParticle(pCopy);
          id++;

          if (posDeq.size() <= (unsigned)deqSize) {
            posDeq.push_front(posVecTmp);
          }
          else {
            posDeq.push_front(posVecTmp);
            posDeq.pop_back();
          }

        }
      }
    }
  }
  normalize(r);
}


template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::setMagneticParticlesdMomRandom()
{

  for (auto pS : _pSystems) {
    std::deque<MagneticParticle3D<double>*> particles = pS->getParticlesPointer();

    for (auto p : particles) {
      std::vector<double> dMoment = { 0., 0., 0. };
      for (int i = 0; i < 3; i++) {
        dMoment[i] = rand() % (9 - (-9) + 1) + (-9);
      }

      double dMoment_norm = sqrt(pow(dMoment[0], 2.) + pow(dMoment[1], 2.) + pow(dMoment[2], 2.)) ;

      for (int i = 0; i < 3; i++) {
        dMoment[i] /= dMoment_norm ;
      }

      p->setMoment(dMoment);
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::setParticlesVelRandom(T velFactor)
{

  for (auto pS : _pSystems) {
    std::deque<PARTICLETYPE<T>*> particles = pS->getParticlesPointer();

    for (auto p : particles) {
      std::vector<T> vel = { 0., 0., 0. };
      for (int i = 0; i < 3; i++) {
        vel[i] = rand() % (9 - (-9) + 1) + (-9);
      }

      T vel_norm = sqrt(pow(vel[0], 2.) + pow(vel[1], 2.) + pow(vel[2], 2.)) ;

      for (int i = 0; i < 3; i++) {
        vel[i] /= vel_norm ;
        vel[i] *= velFactor ;
      }

      p->setVel(vel);
    }
  }
}


template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::setParticlesPosRandom(T posFactor)
{

  for (auto pS : _pSystems) {
    std::deque<PARTICLETYPE<T>*> particles = pS->getParticlesPointer();

    for (auto p : particles) {
      std::vector<T> pos = { 0., 0., 0. };
      for (int i = 0; i < 3; i++) {
        pos[i] = rand() % (9 - (-9) + 1) + (-9);
      }

      T pos_norm = sqrt(pow(pos[0], 2.) + pow(pos[1], 2.) + pow(pos[2], 2.)) ;

      for (int i = 0; i < 3; i++) {
        pos[i] /= pos_norm ;
        pos[i] *= posFactor ;
      }

      for (int i = 0; i < 3; i++) {
        p->getPos()[0] += pos[0] ;
        p->getPos()[1] += pos[1] ;
        p->getPos()[2] += pos[2] ;
      }
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::setParticlesPosRandom(T posFactorX, T posFactorY, T posFactorZ)
{

  for (auto pS : _pSystems) {
    std::deque<PARTICLETYPE<T>*> particles = pS->getParticlesPointer();

    for (auto p : particles) {
      std::vector<T> pos = { 0., 0., 0. };
      for (int i = 0; i < 3; i++) {
        pos[i] = rand() % (9 - (-9) + 1) + (-9);
      }

      T pos_norm = sqrt(pow(pos[0], 2.) + pow(pos[1], 2.) + pow(pos[2], 2.)) ;

      for (int i = 0; i < 3; i++) {
        pos[i] /= pos_norm ;
      }

      pos[0] *= posFactorX ;
      pos[1] *= posFactorY ;
      pos[2] *= posFactorZ ;

      for (int i = 0; i < 3; i++) {
        p->getPos()[0] += pos[0] ;
        p->getPos()[1] += pos[1] ;
        p->getPos()[2] += pos[2] ;
      }
    }
  }
}

template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::setMagneticParticles(std::vector<double> dMoment,
    std::vector<double> vel, std::vector<double> aVel, std::vector<double> torque, double magnetisation)
{
  int i = 0;
  for (auto pS : _pSystems) {
    std::deque<MagneticParticle3D<double>*> particles = pS->getParticlesPointer();
    
    for (auto p : particles) {

      p->setMoment(dMoment);
      p->setVel(vel);
      p->setAVel(aVel);
      p->setTorque(torque);
      p->setMagnetisation(magnetisation);
      p->setID(i);
      i++;

    }
  }
}

template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::setMagneticParticles(std::vector<double> dMoment,
    std::vector<double> vel, std::vector<double> aVel, std::vector<double> torque, double magnetisation, int sActivity)
{
  int i = 0;
  for (auto pS : _pSystems) {
    std::deque<MagneticParticle3D<double>*> particles = pS->getParticlesPointer();

    for (auto p : particles) {

      p->setMoment(dMoment);
      p->setVel(vel);
      p->setAVel(aVel);
      p->setTorque(torque);
      p->setMagnetisation(magnetisation);
      p->setID(i);
      i++;
      p->setSActivity(sActivity);
    }
  }
}

template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::prepareAgglomerates()
{
  for (auto pS : _pSystems) {
    std::list<MagneticParticle3D<double>*> particlesList;
    pS->_Agglomerates.push_back(particlesList) ;
  }
}

template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::initAggloParticles()
{
  for (auto pS : _pSystems) {
    pS->initAggloParticles() ;
  }
}

template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::findAgglomerates(int iT, int itVtkOutputMagParticles)
{
  int pSi = 0;
  for (auto pS : _pSystems) {
    pS->findAgglomerates() ;

    if (iT % itVtkOutputMagParticles == 0) {

      clout << "Particlesystem number: " << pSi << std::endl;
      clout << "Number of non agglomerated particles" << ": " << pS->_Agglomerates[0].size() << std::endl;
      clout << "Number of agglomerated particles" << ": " << pS->size() - pS->_Agglomerates[0].size() << std::endl;
      clout << "Proportion of agglomeratet particles" << ": "
            << double(pS->size() - pS->_Agglomerates[0].size()) / double(pS->size()) * 100. << "%" << std::endl;
      clout << "Number of agglomerates" << ": " << pS->_Agglomerates.size() - 1 << std::endl;
    }
    pSi++;
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addParticleEquallyDistributed(
  IndicatorCuboid3D<T>& cuboid, T pMass, T pRad,/* number of particles on x, y, z axis*/
  int nox, int noy, int noz, std::vector<T> vel)
{
  std::vector < T > pos(3, 0.);
  Vector<T, 3> minPos(cuboid.getMin());
  bool indic[1] = { false };

  clout << "Number of particles to create: nox*noy*noz = "
            << nox * noy * noz << std::endl;

  T xlength = cuboid.getMax()[0] - cuboid.getMin()[0];
  T ylength = cuboid.getMax()[1] - cuboid.getMin()[1];
  T zlength = cuboid.getMax()[2] - cuboid.getMin()[2];
  int modNox = nox - 1, modNoy = noy - 1, modNoz = noz - 1;
  if (nox == 1) {
    modNox = 1;
  }
  if (noy == 1) {
    modNoy = 1;
  }
  if (noz == 1) {
    modNoz = 1;
  }

  for (int i = 0; i < nox; ++i) {
    pos[0] = minPos[0] + (T) (i) * xlength / modNox;
    for (int j = 0; j < noy; ++j) {
      pos[1] = minPos[1] + (T) (j) * ylength / modNoy;
      for (int k = 0; k < noz; ++k) {
        pos[2] = minPos[2] + (T) (k) * zlength / modNoz;

        if (cuboid(indic, &pos[0])) {
          PARTICLETYPE<T> p(pos, vel, pMass, pRad);
          addParticle(p);
        }
      }
    }
  }
  clout << "Number of created particles = "
            << globalNumOfParticles() << std::endl;
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addParticleEquallyDistributed(
  IndicatorCuboid3D<T>& cuboid, int nox, int noy, int noz, PARTICLETYPE<T>& p)
{
  std::vector < T > pos(3, 0.);
  int id = 0;
  Vector<T, 3> minPos(cuboid.getMin());
  PARTICLETYPE<T> pCopy(p);
  bool indic[1] = { false };

  clout << "Number of particles to create: nox*noy*noz = "
        << nox * noy * noz << std::endl;

  T xlength = cuboid.getMax()[0] - cuboid.getMin()[0];
  T ylength = cuboid.getMax()[1] - cuboid.getMin()[1];
  T zlength = cuboid.getMax()[2] - cuboid.getMin()[2];
  int modNox = nox - 1, modNoy = noy - 1, modNoz = noz - 1;
  if (nox == 1) {
    modNox = 1;
  }
  if (noy == 1) {
    modNoy = 1;
  }
  if (noz == 1) {
    modNoz = 1;
  }

  for (int i = 0; i < nox; ++i) {
    pos[0] = minPos[0] + (T) (i) * xlength / modNox;
    for (int j = 0; j < noy; ++j) {
      pos[1] = minPos[1] + (T) (j) * ylength / modNoy;
      for (int k = 0; k < noz; ++k) {
        pos[2] = minPos[2] + (T) (k) * zlength / modNoz;

        if (cuboid(indic, &pos[0])) {
          pCopy.setPos(pos);
          pCopy.setID(id);
          addParticle(pCopy);
          id++;
        }
      }
    }
  }
  clout << "Number of created particles = "
        << globalNumOfParticles() << std::endl;
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addTracerParticle(
  IndicatorF3D<T>& ind, T idTP, T mas, T rad, int noTP, std::vector<T> vel)
{
  srand(clock());
  std::vector < T > pos(3, 0.);
  bool indic[1] = { false };

  noTP += globalNumOfTracerParticles();
  while (globalNumOfTracerParticles() < noTP) {
    pos[0] = ind.getMin()[0]
             + (T) (rand() % 100000) / 100000. * (ind.getMax()[0] - ind.getMin()[0]);
    pos[1] = ind.getMin()[1]
             + (T) (rand() % 100000) / 100000. * (ind.getMax()[1] - ind.getMin()[1]);
    pos[2] = ind.getMin()[2]
             + (T) (rand() % 100000) / 100000. * (ind.getMax()[2] - ind.getMin()[2]);

#ifdef PARALLEL_MODE_MPI
    singleton::mpi().bCast(&*pos.begin(), 3);
#endif

    int x0, y0, z0, C;
    std::vector<int> locLat(4, 0);
    if (this->_cuboidGeometry.getFloorLatticeR(pos, locLat)) {
      C = locLat[0];
      if (this->_loadBalancer.rank(C) == singleton::mpi().getRank()) {
        x0 = locLat[1];
        y0 = locLat[2];
        z0 = locLat[3];
        if (_superGeometry.get(C, x0, y0, z0) == 1
            && _superGeometry.get(C, x0, y0 + 1, z0) == 1
            && _superGeometry.get(C, x0, y0, z0 + 1) == 1
            && _superGeometry.get(C, x0, y0 + 1, z0 + 1) == 1
            && _superGeometry.get(C, x0 + 1, y0, z0) == 1
            && _superGeometry.get(C, x0 + 1, y0 + 1, z0) == 1
            && _superGeometry.get(C, x0 + 1, y0, z0 + 1) == 1
            && _superGeometry.get(C, x0 + 1, y0 + 1, z0 + 1) == 1
            && ind(indic, &pos[0])
            && !util::nearZero(idTP)) {
          PARTICLETYPE<T> p(pos, vel, mas, rad, idTP);
          addParticle(p);
        }
      }
    }
  }
}


template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addParticleBoxMuller(
  IndicatorF3D<T>& ind, T partRho, T mu, T sigma, int no, T appProb,  std::vector<T> vel)
{

  srand(clock());

  //Randomization for Appearance-Likelyhood
  T rndmApp[1] = {(T) (rand() % 100000) / 100000.};

#ifdef PARALLEL_MODE_MPI
  singleton::mpi().bCast(rndmApp, 1);
#endif

  if (rndmApp[0] <= appProb ) {

    std::vector<T> pos(3, 0.);
    T rad;
    T mas;

    bool indic[1] = { false };

    no += globalNumOfParticles();
    while (globalNumOfParticles() < no) {
      pos[0] = ind.getMin()[0]
               + (T) (rand() % 100000) / 100000. * (ind.getMax()[0] - ind.getMin()[0]);
      pos[1] = ind.getMin()[1]
               + (T) (rand() % 100000) / 100000. * (ind.getMax()[1] - ind.getMin()[1]);
      pos[2] = ind.getMin()[2]
               + (T) (rand() % 100000) / 100000. * (ind.getMax()[2] - ind.getMin()[2]);

      //Normally Distributed Particle Radius (Box-Muller Method)
      T u1 = (T) (rand() % 100000) / 100000.;
      T u2 = (T) (rand() % 100000) / 100000.;

      T x = cos(2 * M_PI * u1) * sqrt(-2 * log(u2));
      rad = mu + x * sigma;
      mas = 4. / 3. * M_PI * std::pow( rad, 3 ) * partRho;

#ifdef PARALLEL_MODE_MPI
      singleton::mpi().bCast(&*pos.begin(), 3);
#endif

      int x0, y0, z0, C;
      std::vector<int> locLat(4, 0);
      if (this->_cuboidGeometry.getFloorLatticeR(pos, locLat)) {
        C = locLat[0];
        if (this->_loadBalancer.rank(C) == singleton::mpi().getRank()) {
          x0 = locLat[1];
          y0 = locLat[2];
          z0 = locLat[3];
          if (_superGeometry.get(C, x0, y0, z0) == 1
              && _superGeometry.get(C, x0, y0 + 1, z0) == 1
              && _superGeometry.get(C, x0, y0, z0 + 1) == 1
              && _superGeometry.get(C, x0, y0 + 1, z0 + 1) == 1
              && _superGeometry.get(C, x0 + 1, y0, z0) == 1
              && _superGeometry.get(C, x0 + 1, y0 + 1, z0) == 1
              && _superGeometry.get(C, x0 + 1, y0, z0 + 1) == 1
              && _superGeometry.get(C, x0 + 1, y0 + 1, z0 + 1) == 1
              && ind(indic, &pos[0])) {
            PARTICLETYPE<T> p(pos, vel, mas, rad);
            addParticle(p);
          }
        }
      }
    }
  }
}


template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addParticleWithDistance(
  IndicatorCuboid3D<T>& ind, T pMass, T pRad, std::vector<T> vel,
  T conc, T minDist, bool checkDist)
{

  srand(clock());
  bool indic[1] = { false };
  std::vector < T > pos(3, 0.);
  int size = T();
  T dist = T();
  Vector<T, 3> diff;
  int C;

  T indicatorVol = (ind.getMax()[0] - ind.getMin()[0])
                   * (ind.getMax()[1] - ind.getMin()[1])
                   * (ind.getMax()[2] - ind.getMin()[2]);

  int noParticles = (int) (conc * indicatorVol / (4. / 3 * M_PI * pow(pRad, 3)));

  if (checkDist == true && (noParticles * pow(minDist, 3) * 4 / 3. * M_PI > indicatorVol) ) {
    std::cout << "Error: minDist too large" << std::endl;
    exit(-1);
  }

  std::cout << " noparticles " << noParticles << std::endl;

  noParticles += globalNumOfParticles();

  while (globalNumOfParticles() < noParticles) {

    pos[0] = ind.getMin()[0]
             + (T) (rand() % 100000) / 100000. * (ind.getMax()[0] - ind.getMin()[0]);
    pos[1] = ind.getMin()[1]
             + (T) (rand() % 100000) / 100000. * (ind.getMax()[1] - ind.getMin()[1]);
    pos[2] = ind.getMin()[2]
             + (T) (rand() % 100000) / 100000. * (ind.getMax()[2] - ind.getMin()[2]);

#ifdef PARALLEL_MODE_MPI
    singleton::mpi().bCast(&*pos.begin(), 3);
#endif

    std::vector<int> locLat(4, 0);
    if (this->_cuboidGeometry.getFloorLatticeR(pos, locLat)) {
      C = locLat[0];
      //related cuboidID of floor lattice position
      if (this->_loadBalancer.rank(C) == singleton::mpi().getRank()) {

        if (ind(indic, &pos[0])) {

          int psno = 0;
          int globIC = 0;
          for (unsigned int pS = 0; pS < _pSystems.size(); ++pS) {
            globIC = this->_loadBalancer.glob(pS);
            if (globIC == C) {
              size = _pSystems[pS]->sizeInclShadow();
              psno = pS;
            }
          }
          if (size == 0) {
            PARTICLETYPE<T> p(pos, vel, pMass, pRad);
            addParticle(p);
          } else {
            for (int j = 0; j < size;) {
              diff[0] = _pSystems[psno]->operator[](j).getPos()[0]
                        - pos[0];
              diff[1] = _pSystems[psno]->operator[](j).getPos()[1]
                        - pos[1];
              diff[2] = _pSystems[psno]->operator[](j).getPos()[2]
                        - pos[2];
              dist = norm(diff);

              if (dist < minDist) {
                goto marke;
              }
              if (j == (size - 1)) {
                PARTICLETYPE<T> p(pos, vel, pMass, pRad);
                addParticle(p);
                j += 1;
              } else {
                j += 1;
              }
            }
          }
marke:
          ;
        }
      }
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addShadowParticle(
  PARTICLETYPE<T>& p)
{
  for (unsigned int pS = 0; pS < _pSystems.size(); ++pS) {
    int globIC = this->_loadBalancer.glob(pS);
    if (globIC != p.getCuboid() && checkCuboid(p, _overlap, globIC)
        && !checkCuboid(p, 0, globIC)) {
      _pSystems[pS]->addShadowParticle(p);
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
std::vector<ParticleSystem3D<T, PARTICLETYPE>*> SuperParticleSystem3D<T,
    PARTICLETYPE>::getParticleSystems()
{
  return _pSystems;
}

template<typename T, template<typename U> class PARTICLETYPE>
ParticleSystem3D<T, PARTICLETYPE>& SuperParticleSystem3D<T, PARTICLETYPE>::operator[](
  int i)
{
  return *(_pSystems[i]);
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addForce(
  std::shared_ptr<Force3D<T, PARTICLETYPE> > f)
{
  for (auto pS : _pSystems) {
    pS->addForce(f);
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addBoundary(
  std::shared_ptr<Boundary3D<T, PARTICLETYPE> > b)
{
  for (auto pS : _pSystems) {
    pS->addBoundary(b);
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::saveToFile(std::string name)
{
  std::string fullName = createFileName(name) + ".particles";

  int rank = 0;

#ifdef PARALLEL_MODE_MPI
  int size = 1;
  size = singleton::mpi().getSize();
  rank = singleton::mpi().getRank();
#endif

  if (rank == 0) {
    std::ofstream fout(fullName.c_str(), std::ios::trunc);
    if (!fout) {
      clout << "Error: could not open " << fullName << std::endl;
    }
    fout.close();
  }
#ifdef PARALLEL_MODE_MPI
  if (rank > 0) {
    int prev = rank - 1;
    int buffer = 0;
    MPI_Status status;
    MPI_Recv(&buffer, 1, MPI_INT, prev, 0, MPI_COMM_WORLD, &status);
  }
#endif
  for (auto pS : _pSystems) {
    pS->saveToFile(fullName);
  }
#ifdef PARALLEL_MODE_MPI
  if (rank < size - 1) {
    int next = rank + 1;
    int buffer = 0;
    MPI_Send(&buffer, 1, MPI_INT, next, 0, MPI_COMM_WORLD);
  }
#endif
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::addParticlesFromFile(
  std::string name, T mass, T radius)
{
  std::string fullName = createFileName(name) + ".particles";
  std::ifstream fin(fullName.c_str());

  std::string line;
  while (std::getline(fin, line)) {
    std::istringstream iss(line);
    T buffer[PARTICLETYPE<T>::serialPartSize];
    for (int i = 0; i < PARTICLETYPE<T>::serialPartSize; i++) {
      iss >> buffer[i];
    }
    PARTICLETYPE<T> p;
    p.unserialize(buffer);
    if ( !util::nearZero(radius) ) {
      p.setRad(radius);
    }
    if ( !util::nearZero(mass) ) {
      p.setMass(mass);
    }
    addParticle(p);
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::clearParticles()
{
  for (auto pS : _pSystems) {
    pS->clearParticles();
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::setContactDetection(
  ContactDetection<T, PARTICLETYPE>& contactDetection)
{
  clout << "Setting ContactDetectionAlgorithm " << contactDetection.getName() << std::endl;

  for (auto pS : _pSystems) {
    pS->setContactDetection(contactDetection);
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
void SuperParticleSystem3D<T, PARTICLETYPE>::setContactDetectionForPSys(
  ContactDetection<T, PARTICLETYPE>& contactDetection, int pSysNr)
{
  clout << "Setting ContactDetectionAlgorithm for pSys: " << pSysNr
        << " = " << contactDetection.getName() << std::endl;

  _pSystems[pSysNr]->setContactDetection(contactDetection);
  // this->getParticleSystems()[pSysNr]->setContactDetection(contactDetection);
}


//template<typename T, template<typename U> class PARTICLETYPE>
//template<typename DESCRIPTOR>
//void SuperParticleSystem3D<T, PARTICLETYPE>::particleOnFluid(
//    SuperLattice3D<T, DESCRIPTOR>& sLattice, T eps,
//    SuperGeometry3D<T>& sGeometry) {
//  for (unsigned int i = 0; i < _pSystems.size(); ++i) {
//    _pSystems[i]->particleOnFluid(
//        //        sLattice.getExtendedBlockLattice(i),
//        sLattice.getBlockLattice(i),
//        sLattice.getCuboidGeometry().get(this->_loadBalancer.glob(i)),
//        sLattice.getOverlap(), eps, sGeometry.getExtendedBlockGeometry(i));
//  }
//}
//
//template<typename T, template<typename U> class PARTICLETYPE>
//template<typename DESCRIPTOR>
//void SuperParticleSystem3D<T, PARTICLETYPE>::resetFluid(
//    SuperLattice3D<T, DESCRIPTOR>& sLattice) {
//  for (unsigned int i = 0; i < _pSystems.size(); ++i) {
//    _pSystems[i]->resetFluid(
//        //        sLattice.getExtendedBlockLattice(i),
//        sLattice.getBlockLattice(i),
//        sLattice.getCuboidGeometry().get(this->_loadBalancer.glob(i)),
//        sLattice.getOverlap());
//  }
//}

//template<typename T>
//class SuperRotParticleSystem3D : public SuperParticleSystem3D<T, RotatingParticle3D> {
// public:
//  SuperRotParticleSystem3D(SuperGeometry3D<T>& sg, UnitConverter<T,DESCRIPTOR>& conv): SuperParticleSystem3D<T, RotatingParticle3D>(sg, conv)
//  {};
//  void simulate(T dT) {
//   for (auto pS : this->_pSystems) {
//     pS->_contactDetection->sort();
//     pS->simulate(dT);
//     pS->computeTorque();
//     pS->computeBoundary();
//   }
//   this->updateParticleDistribution();
// }
//};

}//namespace olb

#endif /* SUPERPARTICLESYSTEM_3D_HH */
