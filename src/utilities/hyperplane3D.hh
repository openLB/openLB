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

#ifndef HYPERPLANE_3D_HH
#define HYPERPLANE_3D_HH

#include "hyperplane3D.h"
#include "core/olbDebug.h"

namespace olb {


template <typename T>
Hyperplane3D<T>& Hyperplane3D<T>::spannedBy(const Vector<T,3>& a, const Vector<T,3>& b)
{
  u = a;
  v = b;
  normal = crossProduct3D(a,b);
  normal.normalize();

  OLB_POSTCONDITION(util::nearZero(util::dotProduct3D(u,v)));
  OLB_POSTCONDITION(util::nearZero(util::dotProduct3D(u,normal)));
  OLB_POSTCONDITION(util::nearZero(util::dotProduct3D(v,normal)));

  return *this;
}

template <typename T>
Hyperplane3D<T>& Hyperplane3D<T>::normalTo(const Vector<T,3>& n)
{
  normal = n;

  if ( util::nearZero(normal[0]*normal[1]*normal[2]) ) {
    if ( util::nearZero(normal[0]) ) {
      u = {T(1), T(), T()};
    }
    else if ( util::nearZero(normal[1]) ) {
      u = {T(), T(1), T()};
    }
    else if ( util::nearZero(normal[2]) ) {
      u = {T(), T(), T(1)};
    }
  }
  else {
    u = {normal[2], T(), -normal[0]};
  }

  v = crossProduct3D(normal,u);
  u.normalize();
  v.normalize();
  normal.normalize();

  OLB_POSTCONDITION(util::nearZero(util::dotProduct3D(u,v)));
  OLB_POSTCONDITION(util::nearZero(util::dotProduct3D(u,normal)));
  OLB_POSTCONDITION(util::nearZero(util::dotProduct3D(v,normal)));

  return *this;
}

template <typename T>
Hyperplane3D<T>& Hyperplane3D<T>::originAt(const Vector<T,3>& o)
{
  origin[0] = o[0] - 2*std::numeric_limits<T>::epsilon()*fabs(o[0]);
  origin[1] = o[1] - 2*std::numeric_limits<T>::epsilon()*fabs(o[1]);
  origin[2] = o[2] - 2*std::numeric_limits<T>::epsilon()*fabs(o[2]);

  return *this;
}

template <typename T>
Hyperplane3D<T>& Hyperplane3D<T>::centeredIn(const Cuboid3D<T>& cuboid)
{
  const Vector<T,3>& cuboidOrigin = cuboid.getOrigin();
  const Vector<int,3>& extend     = cuboid.getExtend();
  const T deltaR = cuboid.getDeltaR();

  origin[0] = (cuboidOrigin[0] + 0.5 * deltaR * extend[0]);
  origin[1] = (cuboidOrigin[1] + 0.5 * deltaR * extend[1]);
  origin[2] = (cuboidOrigin[2] + 0.5 * deltaR * extend[2]);
  origin[0] -= 2*std::numeric_limits<T>::epsilon()*fabs(origin[0]);
  origin[1] -= 2*std::numeric_limits<T>::epsilon()*fabs(origin[1]);
  origin[2] -= 2*std::numeric_limits<T>::epsilon()*fabs(origin[2]);

  return *this;
}

template <typename T>
Hyperplane3D<T>& Hyperplane3D<T>::applyMatrixToSpan(
  const Vector<T,3>& row0,
  const Vector<T,3>& row1,
  const Vector<T,3>& row2)
{
  const auto u_prime = u;
  const auto v_prime = v;

  u[0] = row0 * u_prime;
  u[1] = row1 * u_prime;
  u[2] = row2 * u_prime;

  v[0] = row0 * v_prime;
  v[1] = row1 * v_prime;
  v[2] = row2 * v_prime;

  return *this;
}

template <typename T>
Hyperplane3D<T>& Hyperplane3D<T>::rotateSpanAroundX(T r)
{
  return applyMatrixToSpan(
  {1, 0,       0     },
  {0, cos(r), -sin(r)},
  {0, sin(r),  cos(r)}
  );

}

template <typename T>
Hyperplane3D<T>& Hyperplane3D<T>::rotateSpanAroundY(T r)
{
  return applyMatrixToSpan(
  { cos(r), 0, sin(r)},
  { 0,      1, 0     },
  {-sin(r), 0, cos(r)}
  );
}

template <typename T>
Hyperplane3D<T>& Hyperplane3D<T>::rotateSpanAroundZ(T r)
{
  return applyMatrixToSpan(
  {cos(r), -sin(r), 0},
  {sin(r),  cos(r), 0},
  {0,       0,      1}
  );
}

template <typename T>
bool Hyperplane3D<T>::isXYPlane() const
{
  return util::nearZero(util::dotProduct3D(normal, {1,0,0})) &&
         util::nearZero(util::dotProduct3D(normal, {0,1,0}));
}

template <typename T>
bool Hyperplane3D<T>::isXZPlane() const
{
  return util::nearZero(util::dotProduct3D(normal, {1,0,0})) &&
         util::nearZero(util::dotProduct3D(normal, {0,0,1}));
}

template <typename T>
bool Hyperplane3D<T>::isYZPlane() const
{
  return util::nearZero(util::dotProduct3D(normal, {0,1,0})) &&
         util::nearZero(util::dotProduct3D(normal, {0,0,1}));
}


}

#endif
