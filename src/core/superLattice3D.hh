/*  This file is part of the OpenLB library
 *  Copyright (C) 2007 Mathias J. Krause
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

/** \file
 * The description of a 3D super lattice -- generic implementation.
 */

#ifndef SUPER_LATTICE_3D_HH
#define SUPER_LATTICE_3D_HH

#include <limits>
#include <numeric>

#include "communication/mpiManager.h"
#include "cell.h"
#include "superLattice3D.h"
#include "io/base64.h"
#include "functors/analytical/analyticalF.h"
#include "functors/lattice/superBaseF3D.h"
#include "functors/lattice/indicator/superIndicatorF3D.h"
#include "io/serializerIO.h"
#include "geometry/superGeometry3D.h"
#include "communication/loadBalancer.h"
#include "geometry/cuboidGeometry3D.h"


namespace olb {

////////////////////// Class SuperLattice3D /////////////////////////

template<typename T, typename DESCRIPTOR>
SuperLattice3D<T, DESCRIPTOR>::SuperLattice3D(SuperGeometry3D<T>& superGeometry)
  : SuperStructure3D<T>(superGeometry.getCuboidGeometry(), superGeometry.getLoadBalancer()),
    _commStream(*this), _commBC(*this)
{
  int overlapBC = superGeometry.getOverlap();
  if (overlapBC >= 1) {
    _commBC_on = true;
    this->_overlap = overlapBC;
  }
  else {
    _commBC_on = false;
    this->_overlap = 1;
  }

  _commStream.init_nh();
  _commStream.add_cells(1);
  _commStream.init();

  this->_communicator.init_nh();
  this->_communicator.add_cells(this->_overlap);
  this->_communicator.init();

  _extendedBlockLattices.reserve(this->_loadBalancer.size());
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    int nX = this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).getNx() + 2 * this->_overlap;
    int nY = this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).getNy() + 2 * this->_overlap;
    int nZ = this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).getNz() + 2 * this->_overlap;

    _extendedBlockLattices.emplace_back(nX, nY, nZ, superGeometry.getExtendedBlockGeometry(iC));
  }
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _blockLattices.emplace_back(
      _extendedBlockLattices[iC],
      this->_overlap, _extendedBlockLattices[iC].getNx() - this->_overlap - 1,
      this->_overlap, _extendedBlockLattices[iC].getNy() - this->_overlap - 1,
      this->_overlap, _extendedBlockLattices[iC].getNz() - this->_overlap - 1
    );
  }
  _statistics = new LatticeStatistics<T> ;
  _statistics_on = true;

  if (_commBC_on) {
    _commBC.init_nh();
  }

  this->_communicationNeeded=true;
}

template<typename T, typename DESCRIPTOR>
SuperLattice3D<T,DESCRIPTOR>::~SuperLattice3D ()
{
  delete _statistics;
}

template<typename T, typename DESCRIPTOR>
bool SuperLattice3D<T, DESCRIPTOR>::set(T iX, T iY, T iZ, Cell<T, DESCRIPTOR> const& cell)
{
  bool found = false;
  int locX, locY, locZ;
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    if (this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).checkPoint(iX, iY, iZ, locX,
        locY, locZ, this->_overlap)) {
      _extendedBlockLattices[iC].get(locX, locY, locZ) = cell;
      found = true;
    }
  }
  this->_communicationNeeded = true;
  return found;
}

// Note: The dynamics of the cell are not communicated here
template<typename T, typename DESCRIPTOR>
bool SuperLattice3D<T, DESCRIPTOR>::get(T iX, T iY, T iZ, Cell<T, DESCRIPTOR>& cell) const
{
  int locX = 0;
  int locY = 0;
  int locZ = 0;
  bool found = false;
  int foundIC = 0;

  T physR[] = {iX, iY, iZ};
  int latticeR[4];
  if (this->_cuboidGeometry.getLatticeR(latticeR,physR) ) {
    found = true;
    foundIC = latticeR[0];
    locX = latticeR[1];
    locY = latticeR[2];
    locZ = latticeR[3];
  }

  /*for (int iC = 0; iC < this->_cuboidGeometry.getNc(); ++iC) {
    if (this->_cuboidGeometry.get(iC).checkPoint(iX, iY, iZ, locX, locY, locZ)) {
      found = true;
      foundIC = iC;
    }
  }*/

#ifdef PARALLEL_MODE_MPI
  const int sizeOfCell = DESCRIPTOR::size();
  T* cellData = new T[sizeOfCell];

  if (found) {
    if (this->_loadBalancer.rank(foundIC)==singleton::mpi().getRank()) {
      _blockLattices[this->_loadBalancer.loc(foundIC)].get(locX,locY,locZ).serialize(cellData);
    }
    singleton::mpi().bCast(cellData, sizeOfCell, this->_loadBalancer.rank(foundIC));
    cell.unSerialize(cellData);
    delete [] cellData;
  }
#else
  if (found) {
    cell = _blockLattices[this->_loadBalancer.loc(foundIC)].get(locX, locY, locZ);
  }
#endif

  return found;
}

template<typename T, typename DESCRIPTOR>
Cell<T,DESCRIPTOR> SuperLattice3D<T,DESCRIPTOR>::get(int iC, int locX, int locY, int locZ) const
{
  Cell<T,DESCRIPTOR> cell;
#ifdef PARALLEL_MODE_MPI
  const int sizeOfCell = DESCRIPTOR::size();
  T* cellData = new T[sizeOfCell];

  if (this->_loadBalancer.rank(iC)==singleton::mpi().getRank()) {
    _blockLattices[this->_loadBalancer.loc(iC)].get(locX,locY,locZ).serialize(cellData);
  }
  singleton::mpi().bCast(cellData, sizeOfCell, this->_loadBalancer.rank(iC));
  cell.unSerialize(cellData);

  delete [] cellData;
#else
  cell = _blockLattices[this->_loadBalancer.loc(iC)].get(locX,locY,locZ);
#endif
  return cell;
}

template<typename T, typename DESCRIPTOR>
Cell<T,DESCRIPTOR> SuperLattice3D<T,DESCRIPTOR>::get(std::vector<int> latticeR) const
{
  return _blockLattices[this->_loadBalancer.loc(latticeR[0])].get(latticeR[1],latticeR[2],latticeR[3]);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::initialize()
{
  if (_commBC_on) {
    _commBC.init();
  }

  //reset_statistics();
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    //AENDERN VON INI in BLOCKLATTICEVIEW!!!!
    //_blockLattices[iC].initialize();
    _blockLattices[iC].postProcess();
  }

  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::defineDynamics(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, Dynamics<T, DESCRIPTOR>* dynamics)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].defineDynamics(
      indicator->getExtendedBlockIndicatorF(iC), dynamics);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::defineDynamics(
  SuperGeometry3D<T>& sGeometry, int material, Dynamics<T, DESCRIPTOR>* dynamics)
{
  defineDynamics(sGeometry.getMaterialIndicator(material), dynamics);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::defineRho(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, AnalyticalF3D<T,T>& rho)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].defineRho(indicator->getExtendedBlockIndicatorF(iC),
                                         rho);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::defineRho(
  SuperGeometry3D<T>& sGeometry, int material, AnalyticalF3D<T,T>& rho)
{
  defineRho(sGeometry.getMaterialIndicator(material), rho);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::defineU(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, AnalyticalF3D<T,T>& u)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].defineU(indicator->getExtendedBlockIndicatorF(iC),
                                       u);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::defineU(
  SuperGeometry3D<T>& sGeometry, int material, AnalyticalF3D<T,T>& u)
{
  defineU(sGeometry.getMaterialIndicator(material), u);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::defineRhoU(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
    AnalyticalF3D<T,T>& rho, AnalyticalF3D<T,T>& u)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].defineRhoU(indicator->getExtendedBlockIndicatorF(iC),
                                          rho, u);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::defineRhoU(SuperGeometry3D<T>& sGeometry, int material,
    AnalyticalF3D<T,T>& rho, AnalyticalF3D<T,T>& u)
{
  defineRhoU(sGeometry.getMaterialIndicator(material), rho, u);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::definePopulations(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, AnalyticalF3D<T,T>& Pop)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].definePopulations(
      indicator->getExtendedBlockIndicatorF(iC), Pop);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::definePopulations(
  SuperGeometry3D<T>& sGeometry, int material, AnalyticalF3D<T,T>& Pop)
{
  definePopulations(sGeometry.getMaterialIndicator(material), Pop);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::definePopulations(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, SuperF3D<T,T>& Pop)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].definePopulations(
      indicator->getExtendedBlockIndicatorF(iC), Pop.getBlockF(iC));
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::definePopulations(
  SuperGeometry3D<T>& sGeometry, int material, SuperF3D<T,T>& Pop)
{
  definePopulations(sGeometry.getMaterialIndicator(material), Pop);
}


template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice3D<T,DESCRIPTOR>::defineField(
  SuperGeometry3D<T>& sGeometry, int material, SuperF3D<T,T>& field)
{
  if (sGeometry.getStatistics().getNvoxel(material)!=0) {
    int overlap = sGeometry.getOverlap();
    for (int iC=0; iC < this->_loadBalancer.size(); ++iC) {
      if (sGeometry.getExtendedBlockGeometry(iC).getStatistics().getNvoxel(material)!=0) {
        const int x0 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMinLatticeR(material)[0];
        const int y0 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMinLatticeR(material)[1];
        const int z0 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMinLatticeR(material)[2];
        const int x1 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMaxLatticeR(material)[0];
        const int y1 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMaxLatticeR(material)[1];
        const int z1 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMaxLatticeR(material)[2];
        for (int iX=x0; iX<=x1; ++iX) {
          for (int iY=y0; iY<=y1; ++iY) {
            for (int iZ=z0; iZ<=z1; ++iZ) {
              if (sGeometry.getExtendedBlockGeometry(iC).getMaterial(iX,iY,iZ) == material) {
                T fieldTmp[DESCRIPTOR::template size<FIELD>()];
                int inputTmp[4]= {this->_loadBalancer.glob(iC),iX-overlap,iY-overlap,iZ-overlap};
                field(fieldTmp,inputTmp);
                _extendedBlockLattices[iC].get(iX,iY,iZ).template defineField<FIELD>(fieldTmp);
              }
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice3D<T,DESCRIPTOR>::defineField(SuperGeometry3D<T>& sGeometry, IndicatorF3D<T>& indicator,
    AnalyticalF3D<T,T>& field)
{
  SuperIndicatorFfromIndicatorF3D<T> indicatorF(indicator, sGeometry);
  defineField<FIELD>(indicatorF, field);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::iniEquilibrium(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
  AnalyticalF3D<T,T>& rho, AnalyticalF3D<T,T>& u)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].iniEquilibrium(
      indicator->getExtendedBlockIndicatorF(iC), rho, u);
  }
  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::iniEquilibrium(
  SuperGeometry3D<T>& sGeometry, int material,
  AnalyticalF3D<T,T>& rho, AnalyticalF3D<T,T>& u)
{
  iniEquilibrium(sGeometry.getMaterialIndicator(material), rho, u);
}

#ifndef OLB_PRECOMPILED
template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::setExternalParticleField(SuperGeometry3D<T>& sGeometry,
    AnalyticalF3D<T,T>& velocity, SmoothIndicatorF3D<T,T,true>& sIndicator)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].setExternalParticleField(
      sGeometry.getExtendedBlockGeometry(iC), velocity, sIndicator);
  }
}
#endif

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::collide()
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _blockLattices[iC].collide();
  }
  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::collide(T x0, T x1, T y0, T y1, T z0, T z1)
{
  int locX0, locX1, locY0, locY1, locZ0, locZ1;
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    if (this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).checkInters(x0, x1, y0, y1,
        z0, z1, locX0, locX1, locY0, locY1, locZ0, locZ1)) {
      _blockLattices[iC].collide(locX0, locX1, locY0, locY1, locZ0, locZ1);
    }
  }
  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::stream()
{
  _commStream.send();
  _commStream.receive();
  _commStream.write();

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].stream(this->_overlap - 1, _extendedBlockLattices[iC].getNx()
                                      - this->_overlap, this->_overlap - 1,
                                      _extendedBlockLattices[iC].getNy() - this->_overlap, this->_overlap - 1,
                                      _extendedBlockLattices[iC].getNz() - this->_overlap);
  }
  if (_commBC_on) {
    _commBC.send();
    _commBC.receive();
    _commBC.write();
  }

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _blockLattices[iC].postProcess();
  }
  if (_statistics_on) {
    reset_statistics();
  }
  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::stream(T x0, T x1, T y0, T y1, T z0, T z1)
{
  _commStream.send();
  _commStream.receive();
  _commStream.write();

  int locX0, locX1, locY0, locY1, locZ0, locZ1;
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    if (this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).checkInters(x0, x1, y0, y1,
        z0, z1, locX0, locX1, locY0, locY1, locZ0, locZ1, this->_overlap)) {
      _extendedBlockLattices[iC].stream(locX0, locX1, locY0, locY1, locZ0, locZ1);
    }
  }
  if (_commBC_on) {
    _commBC.send();
    _commBC.receive();
    _commBC.write();
  }

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _blockLattices[iC].postProcess();
  }
  if (_statistics_on) {
    reset_statistics();
  }
  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::collideAndStream()
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    int x1 = _blockLattices[iC].getNx() - 1;
    int y1 = _blockLattices[iC].getNy() - 1;
    int z1 = _blockLattices[iC].getNz() - 1;

    _blockLattices[iC].collide(0, x1, 0, 0, 0, z1);
    _blockLattices[iC].collide(0, x1, y1, y1, 0, z1);
    if (y1>1) {
      _blockLattices[iC].collide(0, 0, 1, y1 - 1, 0, z1);
      _blockLattices[iC].collide(x1, x1, 1, y1 - 1, 0, z1);
      if (x1>1) {
        _blockLattices[iC].collide(1, x1 - 1, 1, y1 - 1, 0, 0);
        _blockLattices[iC].collide(1, x1 - 1, 1, y1 - 1, z1, z1);
      }
    }
  }

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    if (2*this->_overlap + 2 < _extendedBlockLattices[iC].getNx() &&
        2*this->_overlap + 2 < _extendedBlockLattices[iC].getNy() &&
        2*this->_overlap + 2 < _extendedBlockLattices[iC].getNz() )
      _extendedBlockLattices[iC].bulkCollideAndStream(this->_overlap + 1,
          _extendedBlockLattices[iC].getNx() - this->_overlap - 2, this->_overlap + 1,
          _extendedBlockLattices[iC].getNy() - this->_overlap - 2, this->_overlap + 1,
          _extendedBlockLattices[iC].getNz() - this->_overlap - 2);
  }

  _commStream.send();
  _commStream.receive();
  _commStream.write();

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    int x1 = _extendedBlockLattices[iC].getNx() - 1;
    int y1 = _extendedBlockLattices[iC].getNy() - 1;
    int z1 = _extendedBlockLattices[iC].getNz() - 1;

    if (2*this->_overlap-3 < x1 &&
        2*this->_overlap-3 < z1) {
      _extendedBlockLattices[iC].boundaryStream(0, x1, 0, y1, 0, z1,
          this->_overlap - 1, x1 - this->_overlap + 1,
          this->_overlap - 1, this->_overlap,
          this->_overlap - 1, z1 - this->_overlap + 1);

      _extendedBlockLattices[iC].boundaryStream(0, x1, 0, y1, 0, z1,
          this->_overlap - 1, x1 - this->_overlap + 1,
          y1 - this->_overlap, y1 - this->_overlap + 1,
          this->_overlap - 1, z1 - this->_overlap + 1);
    }

    if (2*this->_overlap+1 < y1 &&
        2*this->_overlap-3 < z1) {
      _extendedBlockLattices[iC].boundaryStream(0, x1, 0, y1, 0, z1,
          this->_overlap - 1, this->_overlap,
          this->_overlap + 1, y1 - this->_overlap - 1,
          this->_overlap - 1, z1 - this->_overlap + 1);

      _extendedBlockLattices[iC].boundaryStream(0, x1, 0, y1, 0, z1,
          x1 - this->_overlap, x1 - this->_overlap + 1,
          this->_overlap + 1, y1 - this->_overlap - 1,
          this->_overlap - 1, z1 - this->_overlap + 1);
    }

    if (2*this->_overlap+1 < x1 &&
        2*this->_overlap+1 < y1) {
      _extendedBlockLattices[iC].boundaryStream(0, x1, 0, y1, 0, z1,
          this->_overlap + 1, x1 - this->_overlap - 1,
          this->_overlap + 1, y1 - this->_overlap - 1,
          this->_overlap - 1, this->_overlap);
      _extendedBlockLattices[iC].boundaryStream(0, x1, 0, y1, 0, z1,
          this->_overlap + 1, x1 - this->_overlap - 1,
          this->_overlap + 1, y1 - this->_overlap - 1,
          z1 - this->_overlap, z1 - this->_overlap + 1);
    }
  }

  this->_communicationNeeded = true;

  if (_commBC_on) {
    _commBC.send();
    _commBC.receive();
    _commBC.write();
    //SuperLattice3D<T, DESCRIPTOR>::communicate();
  }

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _blockLattices[iC].postProcess();
  }
  if (_statistics_on) {
    reset_statistics();
  }
  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::stripeOffDensityOffset ( int x0, int x1, int y0,
    int y1, int z0, int z1, T offset )
{

  int locX0, locX1, locY0, locY1, locZ0, locZ1;
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    if (this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).checkInters(x0, x1, y0, y1,
        z0, z1, locX0, locX1, locY0, locY1, locZ0, locZ1, this->_overlap)) {
      _extendedBlockLattices[iC].stripeOffDensityOffset(locX0, locX1, locY0, locY1,
          locZ0, locZ1, offset);
    }
  }
  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::stripeOffDensityOffset(T offset)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].stripeOffDensityOffset(offset);
  }
  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T>& SuperLattice3D<T, DESCRIPTOR>::getStatistics()
{
  return *_statistics;
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T> const& SuperLattice3D<T, DESCRIPTOR>::getStatistics() const
{
  return *_statistics;
}


template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::reset_statistics()
{
  T weight;
  T sum_weight = 0;
  T average_rho = 0;
  T average_energy = 0;
  T maxU = 0;
  T delta = 0;

  getStatistics().reset();

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    delta = this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).getDeltaR();
    weight = _extendedBlockLattices[iC].getStatistics().getNumCells() * delta
             * delta * delta;
    sum_weight += weight;
    average_rho += _extendedBlockLattices[iC].getStatistics().getAverageRho()
                   * weight;
    average_energy += _extendedBlockLattices[iC].getStatistics().getAverageEnergy()
                      * weight;
    if (maxU < _extendedBlockLattices[iC].getStatistics().getMaxU()) {
      maxU = _extendedBlockLattices[iC].getStatistics().getMaxU();
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(sum_weight, MPI_SUM);
  singleton::mpi().reduceAndBcast(average_rho, MPI_SUM);
  singleton::mpi().reduceAndBcast(average_energy, MPI_SUM);
  singleton::mpi().reduceAndBcast(maxU, MPI_MAX);
#endif

  average_rho = average_rho / sum_weight;
  average_energy = average_energy / sum_weight;

  getStatistics().reset(average_rho, average_energy, maxU, (int) sum_weight);
  getStatistics().incrementTime();
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    delta = this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).getDeltaR();
    _extendedBlockLattices[iC].getStatistics().reset(average_rho, average_energy,
        maxU, (int) sum_weight);
    _extendedBlockLattices[iC].getStatistics().incrementTime();
  }
}


template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice3D<T, DESCRIPTOR>::addLatticeCoupling( LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen, SuperLattice3D<T,Slattice>& partnerLattice)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    std::vector< SpatiallyExtendedObject3D* > partnerOne;
    partnerOne.push_back(&partnerLattice.getExtendedBlockLattice(iC));

    int nx = getExtendedBlockLattice(iC).getNx();
    int ny = getExtendedBlockLattice(iC).getNy();
    int nz = getExtendedBlockLattice(iC).getNz();

    LatticeCouplingGenerator3D<T, DESCRIPTOR> *extractedLcGen = lcGen.clone();
    extractedLcGen->reset(1, nx-2, 1, ny-2, 1, nz-2);
    getExtendedBlockLattice(iC).addLatticeCoupling(*extractedLcGen, partnerOne);

    delete extractedLcGen;
  }
  return;
}


template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice3D<T, DESCRIPTOR>::addLatticeCoupling(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
  LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen,
  SuperLattice3D<T,Slattice>& partnerLattice)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    std::vector<SpatiallyExtendedObject3D*> partnerOne;
    partnerOne.push_back(&partnerLattice.getExtendedBlockLattice(iC));

    for (int iX = 1; iX < _extendedBlockLattices[iC].getNx()-1; ++iX) {
      for (int iY = 1; iY < _extendedBlockLattices[iC].getNy()-1; ++iY) {
        for (int iZ = 1; iZ < _extendedBlockLattices[iC].getNz()-1; ++iZ) {
          std::unique_ptr<LatticeCouplingGenerator3D<T, DESCRIPTOR>> extractedLcGen{
            lcGen.clone() };
          //TODO done quick and dirty
          if (extractedLcGen->extract(0, 0, 0, 0, 0, 0) ) {
            if (indicator->getExtendedBlockIndicatorF(iC)(iX, iY, iZ)) {
              extractedLcGen->shift(iX, iY, iZ, this->_loadBalancer.glob(iC));
              _extendedBlockLattices[iC].addLatticeCoupling(*extractedLcGen, partnerOne);
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice3D<T, DESCRIPTOR>::addLatticeCoupling(
  SuperGeometry3D<T>& sGeometry, int material,
  LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen,
  SuperLattice3D<T,Slattice>& partnerLattice)
{
  addLatticeCoupling(sGeometry.getMaterialIndicator(material),
                     lcGen, partnerLattice);
}

template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice3D<T, DESCRIPTOR>::addLatticeCoupling(
  LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen,
  std::vector<SuperLattice3D<T,Slattice>*> partnerLattices)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    std::vector< SpatiallyExtendedObject3D* > partners;
    for (auto partnerLattice: partnerLattices) {
      partners.push_back(&partnerLattice->getExtendedBlockLattice(iC));
    }

    int nx = getExtendedBlockLattice(iC).getNx();
    int ny = getExtendedBlockLattice(iC).getNy();
    int nz = getExtendedBlockLattice(iC).getNz();

    LatticeCouplingGenerator3D<T, DESCRIPTOR> *extractedLcGen = lcGen.clone();
    extractedLcGen->reset(1, nx-2, 1, ny-2, 1, nz-2);
    getExtendedBlockLattice(iC).addLatticeCoupling(*extractedLcGen, partners);

    delete extractedLcGen;
  }
  return;
}


template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice3D<T, DESCRIPTOR>::addLatticeCoupling(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
  LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen,
  std::vector<SuperLattice3D<T,Slattice>*> partnerLattices)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    std::vector< SpatiallyExtendedObject3D* > partners;
    for (auto partnerLattice: partnerLattices) {
      partners.push_back(&partnerLattice->getExtendedBlockLattice(iC));
    }

    for (int iX = 1; iX < _extendedBlockLattices[iC].getNx()-1; ++iX) {
      for (int iY = 1; iY < _extendedBlockLattices[iC].getNy()-1; ++iY) {
        for (int iZ = 1; iZ < _extendedBlockLattices[iC].getNz()-1; ++iZ) {
          std::unique_ptr<LatticeCouplingGenerator3D<T, DESCRIPTOR>> extractedLcGen{
            lcGen.clone() };
          //TODO done quick and dirty
          if (extractedLcGen->extract(0, 0, 0, 0, 0, 0) ) {
            if (indicator->getExtendedBlockIndicatorF(iC)(iX, iY, iZ)) {
              extractedLcGen->shift(iX, iY, iZ, this->_loadBalancer.glob(iC));
              _extendedBlockLattices[iC].addLatticeCoupling(*extractedLcGen, partners);
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice3D<T, DESCRIPTOR>::addLatticeCoupling(
  SuperGeometry3D<T>& sGeometry, int material,
  LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen,
  std::vector<SuperLattice3D<T,Slattice>*> partnerLattices)
{
  addLatticeCoupling(sGeometry.getMaterialIndicator(material),
                     lcGen, partnerLattices);
}


template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::executeCoupling()
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].executeCoupling();
  }
  this->_communicationNeeded = true;
  return;
}


template<typename T, typename DESCRIPTOR>
std::size_t SuperLattice3D<T,DESCRIPTOR>::getNblock() const
{
  return std::accumulate(_extendedBlockLattices.begin(), _extendedBlockLattices.end(), size_t(0),
                         Serializable::sumNblock());
}


template<typename T, typename DESCRIPTOR>
std::size_t SuperLattice3D<T,DESCRIPTOR>::getSerializableSize() const
{
  return std::accumulate(_extendedBlockLattices.begin(), _extendedBlockLattices.end(), size_t(0),
                         Serializable::sumSerializableSize());
}

template<typename T, typename DESCRIPTOR>
bool* SuperLattice3D<T,DESCRIPTOR>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  // NOTE: _extendedBlockLattices is correctly sized after constructing SuperLattice, so no resize should be done!
  for ( auto& bLattice : _extendedBlockLattices ) {
    registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, bLattice, loadingMode);
  }

  return dataPtr;
}

} // namespace olb

#endif
