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

#ifndef BLOCK_LOCAL_AVERAGE_3D_HH
#define BLOCK_LOCAL_AVERAGE_3D_HH

#include "blockLocalAverage3D.h"
#include "indicator/blockIndicatorF3D.h"

namespace olb {


template<typename T, typename W>
BlockLocalAverage3D<T,W>::BlockLocalAverage3D(
  BlockF3D<W>&          f,
  BlockIndicatorF3D<T>& indicatorF,
  T radius)
  : BlockF3D<W>(f.getBlockStructure(), f.getTargetDim()),
    _f(f),
    _indicatorF(indicatorF),
    _radius(radius)
{
  this->getName() = "BlockLocalAverage(" + _f.getName() + ")";
}

template<typename T, typename W>
bool BlockLocalAverage3D<T,W>::operator() (W output[], const int input[])
{
  const auto& geometry = _indicatorF.getBlockGeometryStructure();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = 0.;
  }

  if (!_indicatorF(input)) {
    return true;
  }

  T centerOfSphere[3];
  geometry.getPhysR(centerOfSphere, input);
  IndicatorSphere3D<T> analyticalSphere(centerOfSphere, _radius);
  BlockIndicatorFfromIndicatorF3D<T> latticeSphere(
    analyticalSphere,
    _indicatorF.getBlockGeometryStructure());

  std::size_t voxels(0);
  int inputTmp[3];

  for (inputTmp[0] = 0; inputTmp[0] < geometry.getNx(); ++inputTmp[0]) {
    for (inputTmp[1] = 0; inputTmp[1] < geometry.getNy(); ++inputTmp[1]) {
      for (inputTmp[2] = 0; inputTmp[2] < geometry.getNz(); ++inputTmp[2]) {
        if (latticeSphere(inputTmp) && _indicatorF(inputTmp)) {
          T outputTmp[_f.getTargetDim()];
          _f(outputTmp, inputTmp);
          for (int i = 0; i < this->getTargetDim(); ++i) {
            output[i] += outputTmp[i];
          }
          voxels += 1;
        }
      }
    }
  }

  if (voxels > 0) {
    for (int i = 0; i < this->getTargetDim(); ++i) {
      output[i] /= voxels;
    }
  }

  return true;
}


}

#endif
