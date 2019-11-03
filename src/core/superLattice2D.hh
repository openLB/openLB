/*  This file is part of the OpenLB library
 *
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
 * The description of a 2D super lattice -- generic implementation.
 */


#ifndef SUPER_LATTICE_2D_HH
#define SUPER_LATTICE_2D_HH

#include <limits>
#include <numeric>

#include "communication/mpiManager.h"
#include "blockLattice2D.h"
#include "cell.h"
#include "geometry/cuboidGeometry2D.h"
#include "geometry/superGeometry2D.h"
#include "communication/loadBalancer.h"
#include "superLattice2D.h"
#include "io/base64.h"
#include "functors/lattice/superBaseF2D.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"
#include "io/serializerIO.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
SuperLattice2D<T, DESCRIPTOR>::SuperLattice2D(SuperGeometry2D<T>& superGeometry,
    int overlapRefinement)
  : SuperStructure2D<T>(superGeometry.getCuboidGeometry(), superGeometry.getLoadBalancer()),
    _overlapRefinement(overlapRefinement), _commStream(*this), _commBC(*this)
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

    _extendedBlockLattices.emplace_back(nX, nY);
  }
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    BlockLatticeView2D<T, DESCRIPTOR> lattice(_extendedBlockLattices[iC], this->_overlap,
                                           _extendedBlockLattices[iC].getNx() - this->_overlap - 1, this->_overlap,
                                           _extendedBlockLattices[iC].getNy() - this->_overlap - 1);
    _blockLattices.push_back(lattice);
  }
  _statistics = new LatticeStatistics<T> ;
  _statistics_on = true;

  if (_commBC_on) {
    _commBC.init_nh();
  }

  this->_communicationNeeded=true;
}


template<typename T, typename DESCRIPTOR>
SuperLattice2D<T,DESCRIPTOR>::~SuperLattice2D ()
{
  delete _statistics;
}


template<typename T, typename DESCRIPTOR>
bool SuperLattice2D<T,DESCRIPTOR>::set(T iX, T iY, Cell<T,DESCRIPTOR> const& cell)
{

  bool found = false;
  int locX, locY;
  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    if (this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).checkPoint(iX, iY, locX, locY, this->_overlap)) {
      _extendedBlockLattices[iC].get(locX,locY) = cell;
      found = true;
    }
  }
  return found;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::set(Vector<int,3> pos, const Cell<T,DESCRIPTOR>& cell)
{
#ifdef PARALLEL_MODE_MPI
  if (this->_loadBalancer.isLocal(pos[0])) {
    _blockLattices[this->_loadBalancer.loc(pos[0])].get(pos[1],pos[2]) = cell;
  }
#else
  _blockLattices[this->_loadBalancer.loc(pos[0])].get(pos[1],pos[2]) = cell;
#endif
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::get(Vector<int,3> pos, Cell<T,DESCRIPTOR>& cell) const
{
#ifdef PARALLEL_MODE_MPI
  const int sizeOfCell = DESCRIPTOR::size();
  T* cellData = new T[sizeOfCell];

  if (this->_loadBalancer.isLocal(pos[0])) {
    _blockLattices[this->_loadBalancer.loc(pos[0])].get(pos[1],pos[2]).serialize(cellData);
  }

  singleton::mpi().bCast(cellData, sizeOfCell, this->_loadBalancer.rank(pos[0]));

  if (!this->_loadBalancer.isLocal(pos[0])) {
    cell.unSerialize(cellData);
  }
  else {
    cell = _blockLattices[this->_loadBalancer.loc(pos[0])].get(pos[1],pos[2]);
  }

  delete [] cellData;
#else
  cell = _blockLattices[this->_loadBalancer.loc(pos[0])].get(pos[1],pos[2]);
#endif
}

template<typename T, typename DESCRIPTOR>
bool SuperLattice2D<T,DESCRIPTOR>::get(T iX, T iY, Cell<T,DESCRIPTOR>& cell) const
{
  int locX, locY;
  bool found = false;
  int foundIC = 0;

  for (int iC=0; iC<this->_cuboidGeometry.getNc(); ++iC) {
    if (this->_cuboidGeometry.get(iC).checkPoint(iX, iY, locX, locY)) {
      found = true;
      foundIC = iC;
      break;
    }
  }

#ifdef PARALLEL_MODE_MPI
  const int sizeOfCell = DESCRIPTOR::size();
  T* cellData = new T[sizeOfCell];

  if (found) {
    if (this->_loadBalancer.rank(foundIC)==singleton::mpi().getRank()) {
      _blockLattices[this->_loadBalancer.loc(foundIC)].get(locX,locY).serialize(cellData);
    }
    singleton::mpi().bCast(cellData, sizeOfCell, this->_loadBalancer.rank(foundIC));
    cell.unSerialize(cellData);
    delete [] cellData;
  }
#else
  if (found) {
    cell = _blockLattices[this->_loadBalancer.loc(foundIC)].get(locX,locY);
  }
#endif

  return found;
}

template<typename T, typename DESCRIPTOR>
Cell<T,DESCRIPTOR> SuperLattice2D<T,DESCRIPTOR>::get(int iC, int locX, int locY) const
{
  Cell<T,DESCRIPTOR> cell;
#ifdef PARALLEL_MODE_MPI
  const int sizeOfCell = DESCRIPTOR::size();
  T* cellData = new T[sizeOfCell];

  if (this->_loadBalancer.rank(iC)==singleton::mpi().getRank()) {
    _blockLattices[this->_loadBalancer.loc(iC)].get(locX,locY).serialize(cellData);
  }
  singleton::mpi().bCast(cellData, sizeOfCell, this->_loadBalancer.rank(iC));
  cell.unSerialize(cellData);

  delete [] cellData;
#else
  cell = _blockLattices[this->_loadBalancer.loc(iC)].get(locX,locY);
#endif
  return cell;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::initialize()
{
  if (_commBC_on) {
    _commBC.init();
  }
  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    //AENDERN VON INI in BLOCKLATTICEVIEW!!!!
    //_blockLattices[iC].initialize();
    _blockLattices[iC].postProcess();
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T, DESCRIPTOR>::defineDynamics(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator, Dynamics<T, DESCRIPTOR>* dynamics)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].defineDynamics(
      indicator->getExtendedBlockIndicatorF(iC), dynamics);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T, DESCRIPTOR>::defineDynamics(
  SuperGeometry2D<T>& superGeometry, int material, Dynamics<T, DESCRIPTOR>* dynamics)
{
  defineDynamics(superGeometry.getMaterialIndicator(material), dynamics);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::defineRhoU(T x0, T x1, T y0, T y1, T rho, const T u[DESCRIPTOR::d] )
{
  int locX0, locX1, locY0, locY1;
  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    if (this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).checkInters(x0, x1, y0, y1, locX0, locX1, locY0, locY1, this->_overlap)) {
      for (int iX=locX0; iX<=locX1; ++iX) {
        for (int iY=locY0; iY<=locY1; ++iY) {
          _extendedBlockLattices[iC].get(iX,iY).defineRhoU(rho, u);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::defineRhoU(FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
    AnalyticalF2D<T,T>& rho, AnalyticalF2D<T,T>& u)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].defineRhoU(indicator->getExtendedBlockIndicatorF(iC),
                                          rho, u);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::defineRhoU(SuperGeometry2D<T>& sGeometry, int material,
    AnalyticalF2D<T,T>& rho, AnalyticalF2D<T,T>& u)
{
  defineRhoU(sGeometry.getMaterialIndicator(material), rho, u);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::defineRho ( T x0, T x1, T y0, T y1, T rho )
{

  int locX0, locX1, locY0, locY1;
  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    if (this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).checkInters( x0, x1, y0, y1,
        locX0, locX1, locY0, locY1, this->_overlap)) {
      for (int iX=locX0; iX<=locX1; ++iX) {
        for (int iY=locY0; iY<=locY1; ++iY) {
          _extendedBlockLattices[iC].get(iX,iY).defineRho(rho);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::defineRho(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator, AnalyticalF2D<T,T>& rho)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].defineRho(indicator->getExtendedBlockIndicatorF(iC),
                                         rho);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::defineRho(
  SuperGeometry2D<T>& sGeometry, int material, AnalyticalF2D<T,T>& rho)
{
  defineRho(sGeometry.getMaterialIndicator(material), rho);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::defineU( T x0, T x1, T y0, T y1, const T u[DESCRIPTOR::d] )
{
  int locX0, locX1, locY0, locY1;
  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    if (this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).checkInters( x0, x1, y0, y1,
        locX0, locX1, locY0, locY1, this->_overlap)) {
      for (int iX=locX0; iX<=locX1; ++iX) {
        for (int iY=locY0; iY<=locY1; ++iY) {
          _extendedBlockLattices[iC].get(iX,iY).defineU(u);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::defineU(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator, AnalyticalF2D<T,T>& u)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].defineU(indicator->getExtendedBlockIndicatorF(iC),
                                       u);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::defineU(SuperGeometry2D<T>& sGeometry, int material,
                                        AnalyticalF2D<T,T>& u)
{
  defineU(sGeometry.getMaterialIndicator(material), u);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::definePopulations(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator, AnalyticalF2D<T,T>& Pop)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].definePopulations(
      indicator->getExtendedBlockIndicatorF(iC), Pop);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::definePopulations(SuperGeometry2D<T>& sGeometry,
    int material, AnalyticalF2D<T,T>& Pop)
{
  definePopulations(sGeometry.getMaterialIndicator(material), Pop);
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice2D<T,DESCRIPTOR>::defineField( T x0, T x1, T y0, T y1, T* field )
{
  int locX0, locX1, locY0, locY1;
  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    if (this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).checkInters(x0, x1, y0, y1, locX0, locX1, locY0, locY1, this->_overlap)) {
      for (int iX=locX0; iX<=locX1; ++iX) {
        for (int iY=locY0; iY<=locY1; ++iY) {
          _extendedBlockLattices[iC].get(iX,iY).template defineField<FIELD>(field);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice2D<T,DESCRIPTOR>::defineField(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator, SuperLatticeF2D<T,DESCRIPTOR>& field)
{
  const int overlap = indicator->getSuperGeometry().getOverlap();
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    BlockIndicatorF2D<T>& blockIndicator = indicator->getExtendedBlockIndicatorF(iC);
    if ( !blockIndicator.isEmpty() ) {
      const Vector<int,2> min = blockIndicator.getMin();
      const Vector<int,2> max = blockIndicator.getMax();
      for (int iX = min[0]; iX <= max[0]; ++iX) {
        for (int iY = min[1]; iY <= max[1]; ++iY) {
          if (blockIndicator(iX,iY)) {
            T fieldTmp[DESCRIPTOR::template size<FIELD>()];
            int inputTmp[3]= { this->_loadBalancer.glob(iC), iX-overlap, iY-overlap };
            field(fieldTmp,inputTmp);
            _extendedBlockLattices[iC].get(iX,iY).template defineField<FIELD>(fieldTmp);
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice2D<T,DESCRIPTOR>::defineField(
  SuperGeometry2D<T>& sGeometry, int material, SuperLatticeF2D<T,DESCRIPTOR>& field)
{
  defineField<FIELD>(sGeometry.getMaterialIndicator(material), field);
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice2D<T,DESCRIPTOR>::defineField(SuperGeometry2D<T>& sGeometry, IndicatorF2D<T>& indicator,
    AnalyticalF2D<T,T>& field)
{
  SuperIndicatorFfromIndicatorF2D<T> indicatorF(indicator, sGeometry);
  defineField<FIELD>(indicatorF, field);
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice2D<T,DESCRIPTOR>::addField(SuperGeometry2D<T>& sGeometry,
    IndicatorF2D<T>& indicator, AnalyticalF2D<T,T>& field)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].template addField<FIELD>(sGeometry.getExtendedBlockGeometry(iC), indicator, field);
  }
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice2D<T,DESCRIPTOR>::addField( SuperGeometry2D<T>& sGeometry,
    IndicatorF2D<T>& indicator,
    AnalyticalF2D<T,T>& field, AnalyticalF2D<T,T>& porous)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].template addField<FIELD>(sGeometry.getExtendedBlockGeometry(iC), indicator, field, porous);
  }
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice2D<T,DESCRIPTOR>::multiplyField(SuperGeometry2D<T>& sGeometry,
    IndicatorF2D<T>& indicator, AnalyticalF2D<T,T>& field)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].template multiplyField<FIELD>(sGeometry.getExtendedBlockGeometry(iC), indicator, field);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::iniEquilibrium ( T x0, T x1, T y0, T y1, T rho,
    const T u[DESCRIPTOR::d] )
{
  int locX0, locX1, locY0, locY1;
  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    if (this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).checkInters(x0, x1, y0, y1, locX0, locX1, locY0, locY1, this->_overlap)) {
      for (int iX=locX0; iX<=locX1; ++iX) {
        for (int iY=locY0; iY<=locY1; ++iY) {
          _extendedBlockLattices[iC].get(iX,iY).iniEquilibrium(rho, u);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::iniEquilibrium(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
  AnalyticalF2D<T,T>& rho, AnalyticalF2D<T,T>& u)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].iniEquilibrium(
      indicator->getExtendedBlockIndicatorF(iC), rho, u);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::iniEquilibrium(
  SuperGeometry2D<T>& sGeometry, int material,
  AnalyticalF2D<T,T>& rho, AnalyticalF2D<T,T>& u)
{
  iniEquilibrium(sGeometry.getMaterialIndicator(material), rho, u);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::collide ()
{
  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    _blockLattices[iC].collide();
  }

  this->_communicationNeeded=true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::collide (T x0, T x1, T y0, T y1)
{
  int locX0, locX1, locY0, locY1;
  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    if (this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).checkInters(x0, x1, y0, y1, locX0, locX1, locY0, locY1, this->_overlap)) {
      _blockLattices[iC].collide(locX0, locX1, locY0, locY1);
    }
  }

  this->_communicationNeeded=true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::stream ()
{
  _commStream.send();
  _commStream.receive();
  _commStream.write();

  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].stream(this->_overlap-1, _extendedBlockLattices[iC].getNx() - this->_overlap,
                                      this->_overlap-1, _extendedBlockLattices[iC].getNy() - this->_overlap);
  }
  if (_commBC_on) {
    //_commBC.send();
    //_commBC.receive();
    //_commBC.write();
    _commStream.send();
    _commStream.receive();
    _commStream.write();
  }

  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    _blockLattices[iC].postProcess();
  }
  if (_statistics_on) {
    reset_statistics();
  }

  this->_communicationNeeded=true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::stream (T x0, T x1, T y0, T y1)
{
  _commStream.send();
  _commStream.receive();
  _commStream.write();

  int locX0, locX1, locY0, locY1;
  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    if (this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).checkInters(x0, x1, y0, y1, locX0, locX1, locY0, locY1, this->_overlap)) {
      _extendedBlockLattices[iC].stream(locX0, locX1, locY0, locY1);
    }
  }
  if (_commBC_on) {
    //_commBC.send();
    //_commBC.receive();
    //_commBC.write();
    _commStream.send();
    _commStream.receive();
    _commStream.write();
  }

  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    _blockLattices[iC].postProcess();
  }
  if (_statistics_on) {
    reset_statistics();
  }

  this->_communicationNeeded=true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::collideAndStream ()
{
  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    int x1 = _blockLattices[iC].getNx() - 1;
    int y1 = _blockLattices[iC].getNy() - 1;
    _blockLattices[iC].collide(0,x1,0,0);
    _blockLattices[iC].collide(0,x1,y1,y1);
    if (y1>1) {
      _blockLattices[iC].collide(0,0,1,y1-1);
      _blockLattices[iC].collide(x1,x1,1,y1-1);
    }
  }

  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    if (2*this->_overlap + 2 < _extendedBlockLattices[iC].getNx() &&
        2*this->_overlap + 2 < _extendedBlockLattices[iC].getNy() )
      _extendedBlockLattices[iC].bulkCollideAndStream(this->_overlap+1, _extendedBlockLattices[iC].getNx() - this->_overlap-2,
          this->_overlap+1, _extendedBlockLattices[iC].getNy() - this->_overlap-2);
  }

  _commStream.send();
  _commStream.receive();
  _commStream.write();

  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    int x1 = _extendedBlockLattices[iC].getNx() - 1;
    int y1 = _extendedBlockLattices[iC].getNy() - 1;
    if (2*this->_overlap-3<x1) {
      _extendedBlockLattices[iC].boundaryStream(0,x1,0,y1, this->_overlap-1,x1-this->_overlap+1,  this->_overlap-1,this->_overlap);
      _extendedBlockLattices[iC].boundaryStream(0,x1,0,y1, this->_overlap-1,x1-this->_overlap+1,  y1-this->_overlap,y1-this->_overlap+1);
    }

    if (2*this->_overlap+1<y1) {
      _extendedBlockLattices[iC].boundaryStream(0,x1,0,y1, this->_overlap-1, this->_overlap,      this->_overlap+1, y1-this->_overlap-1);
      _extendedBlockLattices[iC].boundaryStream(0,x1,0,y1, x1-this->_overlap,x1-this->_overlap+1, this->_overlap+1, y1-this->_overlap-1);
    }
  }


  if (_commBC_on) {
    _commStream.send();
    _commStream.receive();
    _commStream.write();
    //_commBC.send();
    //_commBC.receive();
    //_commBC.write();
  }

  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    _blockLattices[iC].postProcess();
  }
  if (_statistics_on) {
    reset_statistics();
  }

  this->_communicationNeeded=true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::stripeOffDensityOffset ( int x0, int x1, int y0,
    int y1, T offset )
{
  int locX0, locX1, locY0, locY1;
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    if (this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).checkInters(x0, x1, y0, y1,
        locX0, locX1, locY0, locY1, this->_overlap)) {
      _extendedBlockLattices[iC].stripeOffDensityOffset(locX0, locX1, locY0, locY1,
          offset);
    }
  }

  this->_communicationNeeded=true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::stripeOffDensityOffset(T offset)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].stripeOffDensityOffset(offset);
  }
}

template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice2D<T, DESCRIPTOR>::addLatticeCoupling(
  LatticeCouplingGenerator2D<T, DESCRIPTOR> const& lcGen,
  SuperLattice2D<T,Slattice>& partnerLattice)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    std::vector< SpatiallyExtendedObject2D* > partnerOne;
    partnerOne.push_back(&partnerLattice.getExtendedBlockLattice(iC));

    int nx = getExtendedBlockLattice(iC).getNx();
    int ny = getExtendedBlockLattice(iC).getNy();

    LatticeCouplingGenerator2D<T, DESCRIPTOR> *extractedLcGen = lcGen.clone();
    extractedLcGen->reset(1, nx-2, 1, ny-2);
    getExtendedBlockLattice(iC).addLatticeCoupling(*extractedLcGen, partnerOne);

    delete extractedLcGen;
  }
  return;
}


template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice2D<T, DESCRIPTOR>::addLatticeCoupling(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
  LatticeCouplingGenerator2D<T, DESCRIPTOR> const& lcGen,
  SuperLattice2D<T,Slattice>& partnerLattice)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    std::vector<SpatiallyExtendedObject2D*> partnerOne;
    partnerOne.push_back(&partnerLattice.getExtendedBlockLattice(iC));

    for (int iX = 1; iX < _extendedBlockLattices[iC].getNx()-1; ++iX) {
      for (int iY = 1; iY < _extendedBlockLattices[iC].getNy()-1; ++iY) {
        std::unique_ptr<LatticeCouplingGenerator2D<T, DESCRIPTOR>> extractedLcGen{
          lcGen.clone() };
        //TODO done quick and dirty
        if (extractedLcGen->extract(0, 0, 0, 0) ) {
          if (indicator->getExtendedBlockIndicatorF(iC)(iX, iY)) {
            extractedLcGen->shift(iX, iY);
            _extendedBlockLattices[iC].addLatticeCoupling(*extractedLcGen, partnerOne);
          }
        }
      }
    }
  }

  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice2D<T, DESCRIPTOR>::addLatticeCoupling(
  SuperGeometry2D<T>& sGeometry, int material,
  LatticeCouplingGenerator2D<T, DESCRIPTOR> const& lcGen,
  SuperLattice2D<T,Slattice>& partnerLattice)
{
  addLatticeCoupling(sGeometry.getMaterialIndicator(material),
                     lcGen, partnerLattice);
}


template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice2D<T, DESCRIPTOR>::addLatticeCoupling(
  LatticeCouplingGenerator2D<T, DESCRIPTOR> const& lcGen,
  std::vector<SuperLattice2D<T,Slattice>*> partnerLattices)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    std::vector< SpatiallyExtendedObject2D* > partners;
    for (auto partnerLattice: partnerLattices) {
      partners.push_back(&partnerLattice->getExtendedBlockLattice(iC));
    }

    int nx = getExtendedBlockLattice(iC).getNx();
    int ny = getExtendedBlockLattice(iC).getNy();

    LatticeCouplingGenerator2D<T, DESCRIPTOR> *extractedLcGen = lcGen.clone();
    extractedLcGen->reset(1, nx-2, 1, ny-2);
    getExtendedBlockLattice(iC).addLatticeCoupling(*extractedLcGen, partners);

    delete extractedLcGen;
  }
  return;
}


template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice2D<T, DESCRIPTOR>::addLatticeCoupling(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
  LatticeCouplingGenerator2D<T, DESCRIPTOR> const& lcGen,
  std::vector<SuperLattice2D<T,Slattice>*> partnerLattices)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    std::vector<SpatiallyExtendedObject2D*> partners;
    for (auto partnerLattice: partnerLattices) {
      partners.push_back(&partnerLattice->getExtendedBlockLattice(iC));
    }

    for (int iX = 1; iX < _extendedBlockLattices[iC].getNx()-1; ++iX) {
      for (int iY = 1; iY < _extendedBlockLattices[iC].getNy()-1; ++iY) {
        std::unique_ptr<LatticeCouplingGenerator2D<T, DESCRIPTOR>> extractedLcGen{
          lcGen.clone() };
        //TODO done quick and dirty
        if (extractedLcGen->extract(0, 0, 0, 0) ) {
          if (indicator->getExtendedBlockIndicatorF(iC)(iX, iY)) {
            extractedLcGen->shift(iX, iY);
            _extendedBlockLattices[iC].addLatticeCoupling(*extractedLcGen, partners);
          }
        }
      }
    }
  }

  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice2D<T, DESCRIPTOR>::addLatticeCoupling(
  SuperGeometry2D<T>& sGeometry, int material,
  LatticeCouplingGenerator2D<T, DESCRIPTOR> const& lcGen,
  std::vector<SuperLattice2D<T,Slattice>*> partnerLattices)
{
  addLatticeCoupling(sGeometry.getMaterialIndicator(material),
                     lcGen, partnerLattices);
}


template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T, DESCRIPTOR>::executeCoupling()
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].executeCoupling();
  }
  this->_communicationNeeded = true;
  return;
}


template<typename T, typename DESCRIPTOR>
std::size_t SuperLattice2D<T,DESCRIPTOR>::getNblock() const
{
  return std::accumulate(_extendedBlockLattices.begin(), _extendedBlockLattices.end(), size_t(0),
                         Serializable::sumNblock());
}


template<typename T, typename DESCRIPTOR>
std::size_t SuperLattice2D<T,DESCRIPTOR>::getSerializableSize() const
{
  return std::accumulate(_extendedBlockLattices.begin(), _extendedBlockLattices.end(), size_t(0),
                         Serializable::sumSerializableSize());
}

template<typename T, typename DESCRIPTOR>
bool* SuperLattice2D<T,DESCRIPTOR>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  // NOTE: _extendedBlockLattices is correctly sized after constructing SuperLattice, so no resize should be done!
  for ( auto& bLattice : _extendedBlockLattices ) {
    registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, bLattice, loadingMode);
  }

  return dataPtr;
}


template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::reset_statistics ()
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
LatticeStatistics<T>& SuperLattice2D<T,DESCRIPTOR>::getStatistics()
{
  return *_statistics;
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T> const&
SuperLattice2D<T,DESCRIPTOR>::getStatistics() const
{
  return *_statistics;
}
/*
template<typename T, typename DESCRIPTOR>
void SuperLattice2D<T,DESCRIPTOR>::communicate(bool verbose)
{
  _commStream.send();
  _commStream.receive();
  _commStream.write();
}*/



////////// FREE FUNCTIONS //////////


template<typename T, typename DESCRIPTOR>
void setSuperExternalParticleField( SuperGeometry2D<T>& sGeometry, AnalyticalF2D<T,T>& velocity,
                                    SmoothIndicatorF2D<T,T,true>& sIndicator,
                                    SuperLattice2D<T, DESCRIPTOR>& sLattice )
{
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    setBlockExternalParticleField( sGeometry.getExtendedBlockGeometry(iC), velocity, sIndicator,
    sLattice.getExtendedBlockLattice(iC));
  }
}


} // namespace olb

#endif
