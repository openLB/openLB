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

#ifndef HYPERPLANE_2D_HH
#define HYPERPLANE_2D_HH

#include "hyperplane2D.h"
#include "core/olbDebug.h"

namespace olb {


template <typename T>
Hyperplane2D<T>& Hyperplane2D<T>::originAt(const Vector<T,2>& o)
{
  origin[0] = o[0] - 2*std::numeric_limits<T>::epsilon()*fabs(o[0]);
  origin[1] = o[1] - 2*std::numeric_limits<T>::epsilon()*fabs(o[1]);

  return *this;
}

template <typename T>
Hyperplane2D<T>& Hyperplane2D<T>::centeredIn(const Cuboid2D<T>& cuboid)
{
  const Vector<T,2>& cuboidOrigin = cuboid.getOrigin();
  const Vector<int,2>& extend     = cuboid.getExtend();
  const T deltaR = cuboid.getDeltaR();

  origin[0] = (cuboidOrigin[0] + 0.5 * deltaR * extend[0]);
  origin[1] = (cuboidOrigin[1] + 0.5 * deltaR * extend[1]);
  origin[0] -= 2*std::numeric_limits<T>::epsilon()*fabs(origin[0]);
  origin[1] -= 2*std::numeric_limits<T>::epsilon()*fabs(origin[1]);

  return *this;
}

template <typename T>
Hyperplane2D<T>& Hyperplane2D<T>::parallelTo(const Vector<T,2>& direction)
{
  u = direction;
  normal = { u[1], -u[0] };
  normal.normalize();

  OLB_POSTCONDITION(util::nearZero(util::dotProduct2D(u,normal)));

  return *this;
}

template <typename T>
Hyperplane2D<T>& Hyperplane2D<T>::normalTo(const Vector<T,2>& n)
{
  normal = n;

  if ( util::nearZero(normal[0]*normal[1]) ) {
    if ( util::nearZero(normal[0]) ) {
      u = {T(1), T()};
    }
    else if ( util::nearZero(normal[1]) ) {
      u = {T(), T(1)};
    }
  }
  else {
    u = {normal[1], -normal[0]};
  }

  u.normalize();
  normal.normalize();

  OLB_POSTCONDITION(util::nearZero(util::dotProduct2D(u,normal)));

  return *this;
}

template <typename T>
bool Hyperplane2D<T>::isParallelToX() const
{
  return util::nearZero(util::dotProduct2D(normal, {1,0}));
}

template <typename T>
bool Hyperplane2D<T>::isParallelToY() const
{
  return util::nearZero(util::dotProduct2D(normal, {0,1}));
}


}

#endif
