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

#ifndef SUPER_LOCAL_AVERAGE_2D_HH
#define SUPER_LOCAL_AVERAGE_2D_HH

#include "superLocalAverage2D.h"
#include "blockLocalAverage2D.h"
#include "indicator/superIndicatorF2D.h"

namespace olb {


template<typename T, typename W>
SuperLocalAverage2D<T,W>::SuperLocalAverage2D(
  FunctorPtr<SuperF2D<T>>&&          f,
  FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF,
  T radius)
  : SuperF2D<T,W>(f->getSuperStructure(), f->getTargetDim()),
    _f(std::move(f)),
    _indicatorF(std::move(indicatorF)),
    _radius(radius)
{
  this->getName() = "LocalAverage(" + _f->getName() + ")";

  LoadBalancer<T>& load = _f->getSuperStructure().getLoadBalancer();

  if ( _f->getBlockFSize()          == load.size() &&
       _indicatorF->getBlockFSize() == load.size() ) {
    for (int iC = 0; iC < load.size(); ++iC) {
      this->_blockF.emplace_back(
        new BlockLocalAverage2D<T,W>(_f->getBlockF(iC),
                                     _indicatorF->getBlockIndicatorF(iC),
                                     _radius)
      );
    }
  }
}

template<typename T, typename W>
bool SuperLocalAverage2D<T,W>::operator() (W output[], const int input[])
{
  const auto& geometry = this->getSuperStructure().getCuboidGeometry();
  const auto& load     = this->getSuperStructure().getLoadBalancer();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = 0.;
  }

  if (!_indicatorF(input)) {
    return true;
  }

  T centerOfCircle[2];
  geometry.getPhysR(centerOfCircle, input);
  IndicatorCircle2D<T> analyticalCircle(centerOfCircle, _radius);
  SuperIndicatorFfromIndicatorF2D<T> latticeCircle(
    analyticalCircle,
    _indicatorF->getSuperGeometry());

  std::size_t voxels(0);
  int inputTmp[3];

  for (int iC = 0; iC < load.size(); ++iC) {
    inputTmp[0] = load.glob(iC);
    const auto& cuboid = geometry.get(inputTmp[0]);

    for (inputTmp[1] = 0; inputTmp[1] < cuboid.getNx(); ++inputTmp[1]) {
      for (inputTmp[2] = 0; inputTmp[2] < cuboid.getNy(); ++inputTmp[2]) {
        if (latticeCircle(inputTmp) && _indicatorF(inputTmp)) {
          T outputTmp[_f->getTargetDim()];
          _f(outputTmp, inputTmp);
          for (int i = 0; i < this->getTargetDim(); ++i) {
            output[i] += outputTmp[i];
          }
          voxels += 1;
        }
      }
    }
  }

#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(voxels, MPI_SUM);
#endif

  if (voxels > 0) {
    for (int i = 0; i < this->getTargetDim(); ++i) {
#ifdef PARALLEL_MODE_MPI
      singleton::mpi().reduceAndBcast(output[i], MPI_SUM);
#endif
      output[i] /= voxels;
    }
  }

  return true;
}


}

#endif
