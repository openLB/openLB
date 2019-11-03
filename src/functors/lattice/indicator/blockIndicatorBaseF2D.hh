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

#ifndef BLOCK_INDICATOR_BASE_F_2D_HH
#define BLOCK_INDICATOR_BASE_F_2D_HH

#include "blockIndicatorBaseF2D.h"
#include "geometry/blockGeometry2D.h"

namespace olb {

template <typename T>
BlockIndicatorF2D<T>::BlockIndicatorF2D(BlockGeometryStructure2D<T>& geometry)
  : BlockF2D<bool>(geometry.getBlockStructure(), 1),
    _blockGeometryStructure(geometry),
    _cachedData{nullptr}
{ }

template <typename T>
BlockGeometryStructure2D<T>& BlockIndicatorF2D<T>::getBlockGeometryStructure()
{
  return _blockGeometryStructure;
}

template <typename T>
bool BlockIndicatorF2D<T>::operator() (const int input[])
{
  bool output{};
  if (_cachedData == nullptr) {
    this->operator()(&output, input);
  }
  else {
    _cachedData->get(input[0], input[1]);
  }
  return output;
}

template <typename T>
bool BlockIndicatorF2D<T>::operator() (int iX, int iY)
{
  bool output{};
  if (_cachedData == nullptr) {
    this->operator()(&output, iX, iY);
  }
  else {
    _cachedData->get(iX, iY);
  }
  return output;
}

template <typename T>
void BlockIndicatorF2D<T>::setCache(const BlockData2D<T,bool>& cache)
{
  _cachedData = &cache;
}

template <typename T>
bool BlockIndicatorF2D<T>::isEmpty()
{
  // There is no way to determine domain emptyness in a fashion that is both
  // generic and efficient.
  return false;
}

} // namespace olb

#endif
