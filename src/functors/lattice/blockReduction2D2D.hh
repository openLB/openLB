/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#ifndef BLOCK_REDUCTION_2D2D_HH
#define BLOCK_REDUCTION_2D2D_HH

#include <limits>
#include <cmath>

#include "blockReduction2D2D.h"
#include "utilities/vectorHelpers.h"
#include "functors/analytical/interpolationF2D.h"
#include "communication/mpiManager.h"
#include "utilities/functorPtr.hh"
#include "utilities/hyperplane3D.h"
#include "utilities/hyperplaneLattice3D.h"
#include "utilities/hyperplane3D.hh"
#include "utilities/hyperplaneLattice3D.hh"

namespace olb {


template <typename T>
void BlockReduction2D2D<T>::updateToWantedResolution(int resolution)
{
  if (resolution>0) {
    if (_nx > _ny) {
      T newH = _nx*_h/(T)resolution;
      _nx = resolution;
      _ny = (int)(_ny*_h/newH) + 1;
      _h = newH;
    }
    else {
      T newH = _ny*_h/(T)resolution;
      _ny = resolution;
      _nx = (int)(_nx*_h/newH) + 1;
      _h = newH;
    }
    _blockDataMemory.reset(new BlockData2D<T,T>(_nx, _ny, _f->getTargetDim()));
  }
}

template <typename T>
BlockReduction2D2D<T>::BlockReduction2D2D(
  FunctorPtr<SuperF2D<T>>&& f, int resolution, BlockDataSyncMode mode)
  : BlockDataF2D<T,T>(1, 1, f->getTargetDim()),
    _f(std::move(f)),
    _origin(_f->getSuperStructure().getCuboidGeometry().getMotherCuboid().getOrigin()),
    _h(_f->getSuperStructure().getCuboidGeometry().getMinDeltaR()),
    _nx(_f->getSuperStructure().getCuboidGeometry().getMotherCuboid().getNx()),
    _ny(_f->getSuperStructure().getCuboidGeometry().getMotherCuboid().getNy()),
    _blockDataMemory(new BlockData2D<T,T>(_nx, _ny, _f->getTargetDim())),
    _syncMode(mode)
{
  this->getName() = "planeReduction(" + _f->getName() + ")";

  _origin[0] -= 10*std::numeric_limits<T>::epsilon();
  _origin[1] -= 10*std::numeric_limits<T>::epsilon();

  // changes _h, _nx, _ny to match resolution
  updateToWantedResolution(resolution);
  // expose block data fields
  this->_blockData = *_blockDataMemory;
  // intialize list of relevant rank local points
  initialize();
  // first update of data
  update();
}

template <typename T>
Vector<T,2> BlockReduction2D2D<T>::getPhysR(const int& iX, const int& iY) const
{
  return Vector<T,2>{
    _origin[0] + double(iX) * _h,
    _origin[1] + double(iY) * _h
  };
}

template <typename T>
HyperplaneLattice3D<T> BlockReduction2D2D<T>::getPlaneDiscretizationIn3D() const
{
  return HyperplaneLattice3D<T>(
      Hyperplane3D<T>()
      .originAt({_origin[0], _origin[1], 0})
      .spannedBy({1,0,0}, {0,1,0}),
      _h, _nx, _ny);
}

template <typename T>
void BlockReduction2D2D<T>::initialize()
{
  const CuboidGeometry2D<T>& geometry = _f->getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>&           load     = _f->getSuperStructure().getLoadBalancer();

  _rankLocalSubplane.clear();

  for ( int iX = 0; iX < _nx; ++iX ) {
    for ( int iY = 0; iY < _ny; ++iY ) {
      const Vector<T,2> physR = getPhysR(iX, iY);

      // Schedule plane point for storage if its physical position intersects the
      // mother cuboid and the cuboid of the nearest lattice position is local to
      // the current rank:
      int iC;
      if ( geometry.getC(physR, iC) ) {
        if ( load.isLocal(iC) ) {
          _rankLocalSubplane.emplace_back(iX, iY, iC);
        }
      }
    }
  }
}

template <typename T>
void BlockReduction2D2D<T>::update()
{
  _f->getSuperStructure().communicate();

  AnalyticalFfromSuperF2D<T> analyticalF(*_f);

#ifdef PARALLEL_MODE_MPI
  BlockData2D<T,T> localBlockData(_nx, _ny, this->getTargetDim());
#endif

  for ( std::tuple<int,int,int>& pos : _rankLocalSubplane ) {
    const int& iX = std::get<0>(pos);
    const int& iY = std::get<1>(pos);
    const Vector<T,2> physR = getPhysR(iX, iY);

    for ( int iSize = 0; iSize < _f->getTargetDim(); ++iSize ) {
      this->_blockData.get(iX, iY, iSize) = T();
    }

    T output[_f->getTargetDim()];
    const T input[2] { physR[0], physR[1] };

    if (analyticalF(output, input)) {
      for ( int iSize = 0; iSize < _f->getTargetDim(); ++iSize ) {
#ifdef PARALLEL_MODE_MPI
        localBlockData.get(iX, iY, iSize) += output[iSize];
#else
        this->_blockData.get(iX, iY, iSize) += output[iSize];
#endif
      }
    }
  }

#ifdef PARALLEL_MODE_MPI
  switch ( _syncMode ) {
  case BlockDataSyncMode::ReduceAndBcast:
    singleton::mpi().reduce(localBlockData, this->getBlockData(), MPI_SUM);
    singleton::mpi().bCast(this->getBlockData());
    break;
  case BlockDataSyncMode::ReduceOnly:
    singleton::mpi().reduce(localBlockData, this->getBlockData(), MPI_SUM);
    break;
  case BlockDataSyncMode::None:
    this->_blockData.swap(localBlockData);
    break;
  }
#endif
}

template <typename T>
BlockStructure2D& BlockReduction2D2D<T>::getBlockStructure()
{
  return this->_blockData;
}


} // end namespace olb

#endif
