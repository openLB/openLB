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

#ifndef SUPER_CONST_2D_HH
#define SUPER_CONST_2D_HH

#include "superConst2D.h"
#include "core/olbDebug.h"

namespace olb {


template <typename T, typename W>
SuperConst2D<T,W>::SuperConst2D(SuperStructure2D<T>& superStructure, std::vector<W> v)
  : SuperF2D<T,W>(superStructure, v.size()),
    _c{std::move(v)}
{
  this->getName() = "const(" + std::to_string(_c.size()) + ")";
}

template <typename T, typename W>
SuperConst2D<T,W>::SuperConst2D(SuperStructure2D<T>& superStructure, W scalar)
  : SuperConst2D(superStructure, std::vector<W>(1, scalar))
{ }

template <typename T, typename W>
bool SuperConst2D<T,W>::operator() (W output[], const int input[])
{
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = _c[i];
  }
  return true;
}


}

#endif
