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

#ifndef BLOCK_REDUCTION_2D1D_HH
#define BLOCK_REDUCTION_2D1D_HH

#include "blockReduction2D1D.h"

#include <limits>
#include <cmath>

#include "utilities/vectorHelpers.h"
#include "functors/analytical/interpolationF2D.h"
#include "communication/mpiManager.h"
#include "utilities/functorPtr.hh"


namespace olb {


template <typename T>
void BlockReduction2D1D<T>::updateBlockAnalytical(BlockData2D<T,T>& block)
{
  AnalyticalFfromSuperF2D<T> analyticalF(*_f);

  for ( std::tuple<int,int>& pos : _rankLocalSubplane ) {
    const int& i = std::get<0>(pos);
    const Vector<T,2> physR = this->getPhysR(i);

    for ( int iSize = 0; iSize < _f->getTargetDim(); ++iSize ) {
      block.get(i, 0, iSize) = T();
    }

    T output[_f->getTargetDim()];
    const T input[2] { physR[0], physR[1] };

    if (analyticalF(output, input)) {
      for ( int iSize = 0; iSize < _f->getTargetDim(); ++iSize ) {
        block.get(i, 0, iSize) += output[iSize];
      }
    }
  }
}

template <typename T>
void BlockReduction2D1D<T>::updateBlockDiscrete(BlockData2D<T,T>& block)
{
  CuboidGeometry2D<T>& geometry = _f->getSuperStructure().getCuboidGeometry();

  for ( std::tuple<int,int>& pos : _rankLocalSubplane ) {
    const int& i  = std::get<0>(pos);
    const int& iC = std::get<1>(pos);
    const Vector<T,2> physR = this->getPhysR(i);

    for ( int iSize = 0; iSize < _f->getTargetDim(); ++iSize ) {
      block.get(i, 0, iSize) = T();
    }

    T output[_f->getTargetDim()];
    int input[3] { iC, 0, 0 };
    geometry.get(iC).getLatticeR(&input[1], physR);

    if (_f(output, input)) {
      for ( int iSize = 0; iSize < _f->getTargetDim(); ++iSize ) {
        block.get(i, 0, iSize) += output[iSize];
      }
    }
  }
}

template <typename T>
BlockReduction2D1D<T>::BlockReduction2D1D(
  FunctorPtr<SuperF2D<T>>&& f,
  const HyperplaneLattice2D<T>& lattice,
  BlockDataSyncMode      syncMode,
  BlockDataReductionMode reductionMode)
  : HyperplaneLattice2D<T>(lattice),
    BlockDataF2D<T,T>(1, 1, f->getTargetDim()),
    _blockDataMemory(new BlockData2D<T,T>(lattice.getN(), 1, f->getTargetDim())),
    _f(std::move(f)),
    _syncMode(syncMode),
    _reductionMode(reductionMode)
{
  this->getName() = "lineReduction(" + _f->getName() + ")";

  if ( _reductionMode == BlockDataReductionMode::Discrete ) {
    const CuboidGeometry2D<T>& geometry = _f->getSuperStructure().getCuboidGeometry();
    const Hyperplane2D<T>& hyperplane   = this->getHyperplane();
    const bool spansAxisPlane = hyperplane.isParallelToX() ||
                                hyperplane.isParallelToY();
    // verify axes alignment and spacing of hyperplane parametrization
    if ( !spansAxisPlane ||
         lattice.getPhysSpacing() != geometry.getMinDeltaR() ) {
      // hyperplane lattice doesn't describe a trivially discretizable plane
      OstreamManager clerr(std::cerr, "BlockReduction2D1D");
      clerr << "Given hyperplane is not trivially discretizable. "
            << "Use BlockDataReductionMode::Analytical instead."
            << std::endl;
      exit(-1);
    }
  }

  // expose block data fields
  this->_blockData = *_blockDataMemory;
  // intialize list of relevant rank local points making up the reduced line
  initialize();
  // first update of data
  update();
}

template <typename T>
BlockReduction2D1D<T>::BlockReduction2D1D(
  FunctorPtr<SuperF2D<T>>&& f,
  const Hyperplane2D<T>& hyperplane,
  BlockDataSyncMode      syncMode,
  BlockDataReductionMode reductionMode)
  : BlockReduction2D1D(
      std::forward<decltype(f)>(f),
      HyperplaneLattice2D<T>(f->getSuperStructure().getCuboidGeometry(),
                             hyperplane),
      syncMode,
      reductionMode)
{ }

template <typename T>
BlockReduction2D1D<T>::BlockReduction2D1D(
  FunctorPtr<SuperF2D<T>>&& f,
  const Hyperplane2D<T>& hyperplane,
  int resolution, BlockDataSyncMode mode)
  : BlockReduction2D1D(
      std::forward<decltype(f)>(f),
      HyperplaneLattice2D<T>(f->getSuperStructure().getCuboidGeometry(),
                             hyperplane, resolution),
      mode)
{ }

template <typename T>
BlockReduction2D1D<T>::BlockReduction2D1D(
  FunctorPtr<SuperF2D<T>>&& f,
  const Vector<T,2>& origin, const Vector<T,2>& direction,
  int resolution, BlockDataSyncMode mode)
  : BlockReduction2D1D(
      std::forward<decltype(f)>(f),
      Hyperplane2D<T>().originAt(origin).parallelTo(direction),
      resolution, mode) { }

template <typename T>
bool BlockReduction2D1D<T>::operator()(T output[], int i)
{
  const int input[2] = { i, 0 };
  return static_cast<BlockDataF2D<T,T>*>(this)->operator()(output, input);
}

template <typename T>
void BlockReduction2D1D<T>::initialize()
{
  const CuboidGeometry2D<T>& geometry = _f->getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>&           load     = _f->getSuperStructure().getLoadBalancer();

  _rankLocalSubplane.clear();

  for ( int i = 0; i < this->getN(); ++i ) {
    const Vector<T,2> physR = this->getPhysR(i);

    // Schedule line point for storage if its physical position intersects the
    // mother cuboid and the cuboid of the nearest lattice position is local to
    // the current rank:
    int iC;
    if ( geometry.getC(physR, iC) ) {
      if ( load.isLocal(iC) ) {
        _rankLocalSubplane.emplace_back(i, iC);
      }
    }
  }
}

template <typename T>
void BlockReduction2D1D<T>::update()
{
  _f->getSuperStructure().communicate();

#ifdef PARALLEL_MODE_MPI
  BlockData2D<T,T> localBlockData(this->getN(), 1, _f->getTargetDim());

  switch ( _reductionMode ) {
  case BlockDataReductionMode::Analytical:
    updateBlockAnalytical(localBlockData);
    break;
  case BlockDataReductionMode::Discrete:
    updateBlockDiscrete(localBlockData);
    break;
  }

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
#else
  switch ( _reductionMode ) {
  case BlockDataReductionMode::Analytical:
    updateBlockAnalytical(this->_blockData);
    break;
  case BlockDataReductionMode::Discrete:
    updateBlockDiscrete(this->_blockData);
    break;
  }
#endif
}

template <typename T>
BlockStructure2D& BlockReduction2D1D<T>::getBlockStructure()
{
  return this->_blockData;
}

template <typename T>
const std::vector<std::tuple<int,int>>& BlockReduction2D1D<T>::getRankLocalSubplane() const
{
  return this->_rankLocalSubplane;
}


} // end namespace olb

#endif
