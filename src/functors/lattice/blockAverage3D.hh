/*  This file is part of the OpenLB library
 *
 *  Copyright (C) Adrian Kummerlaender
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

#ifndef BLOCK_AVERAGE_3D_HH
#define BLOCK_AVERAGE_3D_HH

#include "blockAverage3D.h"

namespace olb {


template <typename T, typename W>
BlockAverage3D<T,W>::BlockAverage3D(BlockF3D<W>&          f,
                                    BlockIndicatorF3D<T>& indicatorF)
  : BlockSum3D<T,W>(f, indicatorF)
{
  this->getName() = "BlockAverage("+f.getName()+")";
}

template <typename T, typename W>
bool BlockAverage3D<T,W>::operator() (W output[], const int input[])
{
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = W(0);
  }

  const bool result = BlockSum3D<T,W>::operator()(output, input);

  for (int i = 0; i < this->getTargetDim()-1; ++i) {
    output[i] /= output[this->getTargetDim()-1];
  }

  return result;
}


}

#endif
