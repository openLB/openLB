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

#ifndef SUPER_INDICATOR_BASE_F_2D_HH
#define SUPER_INDICATOR_BASE_F_2D_HH

#include "superIndicatorBaseF2D.h"
#include "geometry/superGeometry2D.h"

namespace olb {


template <typename T>
SuperIndicatorF2D<T>::SuperIndicatorF2D(SuperGeometry2D<T>& geometry)
  : SuperF2D<T, bool>(geometry, 1),
    _superGeometry(geometry)
{ }

template <typename T>
BlockIndicatorF2D<T>& SuperIndicatorF2D<T>::getBlockIndicatorF(int iCloc)
{
  OLB_ASSERT(size_t(iCloc) < this->_blockF.size() && iCloc >= 0,
             "block functor index within bounds");
  // Note: The type system doesn't guarantee this operation to be valid
  //       as blockF may contain any implementation of the BlockF2D interface.
  return *static_cast<BlockIndicatorF2D<T>*>(this->_blockF[iCloc].get());
}

template <typename T>
BlockIndicatorF2D<T>& SuperIndicatorF2D<T>::getExtendedBlockIndicatorF(int iCloc)
{
  OLB_ASSERT(iCloc < int(this->_extendedBlockF.size()) && iCloc >= 0,
             "block functor index within bounds");
  return *(this->_extendedBlockF[iCloc]);
}

template <typename T>
SuperGeometry2D<T>& SuperIndicatorF2D<T>::getSuperGeometry()
{
  return this->_superGeometry;
}

template <typename T>
bool SuperIndicatorF2D<T>::operator() (const int input[])
{
  bool output;
  if (_cachedData) {
    output = _cachedData->get(input[0], input[1], input[2]);
  }
  else {
    this->operator()(&output, input);
  }
  return output;
}

template <typename T>
bool SuperIndicatorF2D<T>::operator() (int iC, int iX, int iY)
{
  bool output;
  if (_cachedData) {
    output = _cachedData->get(iC, iX, iY);
  }
  else {
    this->operator()(&output, iC, iX, iY);
  }
  return output;
}

template <typename T>
void SuperIndicatorF2D<T>::cache()
{
  _cachedData = std::unique_ptr<SuperData2D<T,bool>>(
                  new SuperData2D<T,bool>(*this));
  for (unsigned iC = 0; iC < this->_extendedBlockF.size(); ++iC) {
    getExtendedBlockIndicatorF(iC).setCache(_cachedData->get(iC));
  }
}


} // namespace olb

#endif
