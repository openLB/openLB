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

#ifndef BLOCK_LOCAL_AVERAGE_2D_HH
#define BLOCK_LOCAL_AVERAGE_2D_HH

#include "blockLocalAverage2D.h"
#include "indicator/blockIndicatorF2D.h"

namespace olb {


template<typename T, typename W>
BlockLocalAverage2D<T,W>::BlockLocalAverage2D(
  BlockF2D<W>&          f,
  BlockIndicatorF2D<T>& indicatorF,
  T radius)
  : BlockF2D<W>(f.getBlockStructure(), f.getTargetDim()),
    _f(f),
    _indicatorF(indicatorF),
    _radius(radius)
{
  this->getName() = "BlockLocalAverage(" + _f.getName() + ")";
}

template<typename T, typename W>
bool BlockLocalAverage2D<T,W>::operator() (W output[], const int input[])
{
  const auto& geometry = _indicatorF.getBlockGeometryStructure();

  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = 0.;
  }

  if (!_indicatorF(input)) {
    return true;
  }

  T centerOfCircle[2];
  geometry.getPhysR(centerOfCircle, input);
  IndicatorCircle2D<T> analyticalCircle(centerOfCircle, _radius);
  BlockIndicatorFfromIndicatorF2D<T> latticeCircle(
    analyticalCircle,
    _indicatorF.getBlockGeometryStructure());

  std::size_t voxels(0);
  int inputTmp[2];

  for (inputTmp[0] = 0; inputTmp[0] < geometry.getNx(); ++inputTmp[0]) {
    for (inputTmp[1] = 0; inputTmp[1] < geometry.getNy(); ++inputTmp[1]) {
      if (latticeCircle(inputTmp) && _indicatorF(inputTmp)) {
        T outputTmp[_f.getTargetDim()];
        _f(outputTmp, inputTmp);
        for (int i = 0; i < this->getTargetDim(); ++i) {
          output[i] += outputTmp[i];
        }
        voxels += 1;
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
