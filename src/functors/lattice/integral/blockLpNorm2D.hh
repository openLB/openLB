/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerlaender
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

#ifndef BLOCK_LP_NORM_2D_HH
#define BLOCK_LP_NORM_2D_HH

#include "blockLpNorm2D.h"
#include "functors/lattice/indicator/blockIndicatorBaseF2D.h"
#include "geometry/cuboid2D.h"
#include "latticeIntegralCommon.h"

namespace olb {

template <typename T, typename W, int P>
BlockLpNorm2D<T,W,P>::BlockLpNorm2D(BlockF2D<W>&          f,
                                    BlockIndicatorF2D<T>& indicatorF)
  : BlockF2D<W>(f.getBlockStructure(), f.getTargetDim()),
    _f(f),
    _indicatorF(indicatorF)
{
  OLB_ASSERT(_f.getSourceDim() == _indicatorF.getSourceDim(),
             "functor source dimension equals indicator source dimension");

  this->getName() = "BlockL" + std::to_string(P) + "Norm(" + _f.getName() + ")";
}

template <typename T, typename W, int P>
bool BlockLpNorm2D<T,W,P>::operator()(W output[], const int input[])
{
  const auto& blockGeometry = _indicatorF.getBlockGeometryStructure();
  const int nX = blockGeometry.getNx();
  const int nY = blockGeometry.getNy();
  const T weight = pow(blockGeometry.getDeltaR(), 2);

  output[0] = W(0);
  W outputTmp[_f.getTargetDim()];
  int inputTmp[_f.getSourceDim()];

  for (inputTmp[0] = 0; inputTmp[0] < nX; ++inputTmp[0]) {
    for (inputTmp[1] = 0; inputTmp[1] < nY; ++inputTmp[1]) {
      if (_indicatorF(inputTmp)) {
        _f(outputTmp, inputTmp);
        for (int iDim = 0; iDim < _f.getTargetDim(); ++iDim) {
          output[0] = LpNormImpl<T,W,P>()(output[0], outputTmp[iDim], weight);
        }
      }
    }
  }

  output[0] = LpNormImpl<T,W,P>().enclose(output[0]);

  return true;
}

}

#endif
