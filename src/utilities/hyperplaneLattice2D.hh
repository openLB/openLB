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

#ifndef HYPERPLANE_LATTICE_2D_HH
#define HYPERPLANE_LATTICE_2D_HH

#include "hyperplaneLattice2D.h"

namespace olb {

template <typename T>
int HyperplaneLattice2D<T>::computeMaxLatticeDistance() const
{
  const Cuboid2D<T>&  cuboid = _geometry.getMotherCuboid();
  const Vector<T,2>   origin = cuboid.getOrigin();
  const Vector<int,2> extend = cuboid.getExtend();
  const T             deltaR = cuboid.getDeltaR();

  T maxPhysDistance = T();
  T tmp;
  Vector<T,2> tmpO;
  Vector<T,2> tmpE;

  for(int iDim=0; iDim<2; ++iDim){
  tmpO[iDim] = origin[iDim] - _origin[iDim];
  tmpE[iDim] = origin[iDim] + extend[iDim]*deltaR - _origin[iDim];
  }
  tmp = sqrt(tmpO[0]*tmpO[0] + tmpO[1]*tmpO[1]);
  if (maxPhysDistance < tmp) {
      maxPhysDistance = tmp;
  }
  tmp = sqrt((tmpE[0]*tmpE[0] + tmpO[1]*tmpO[1]));
  if (maxPhysDistance < tmp) {
      maxPhysDistance = tmp;
  }
  tmp = sqrt(tmpO[0]*tmpO[0] + tmpE[1]*tmpE[1]);
  if (maxPhysDistance < tmp) {
      maxPhysDistance = tmp;
  }
  tmp = sqrt(tmpE[0]*tmpE[0] + tmpE[1]*tmpE[1]);
  if (maxPhysDistance < tmp) {
      maxPhysDistance = tmp;
  }

  return int(maxPhysDistance/_h) + 1;
}

template <typename T>
void HyperplaneLattice2D<T>::constructCuboid(int maxLatticeDistance)
{
  int iC;
  int min = -maxLatticeDistance;
  int max = maxLatticeDistance;

  for ( int i = -maxLatticeDistance; i < maxLatticeDistance; ++i ) {
    if ( _geometry.getC(getPhysR(i), iC) ) {
      min = i;
      break;
    }
  }
  for ( int i = maxLatticeDistance; i > -maxLatticeDistance; --i ) {
    if ( _geometry.getC(getPhysR(i), iC) ) {
      max = i;
      break;
    }
  }

  _n = max - min + 1;
  _origin = _origin + double(min)*_u;
}

template <typename T>
void HyperplaneLattice2D<T>::setToResolution(int resolution)
{
  T newH = _n*_h/(T) resolution;
  _n = resolution;
  _h = newH;
  _u.normalize(_h);
}

template<typename T>
HyperplaneLattice2D<T>::HyperplaneLattice2D(
  CuboidGeometry2D<T>& geometry, Hyperplane2D<T> hyperplane)
  : _geometry(geometry),
    _hyperplane(hyperplane),
    _origin(hyperplane.origin),
    _u(hyperplane.u),
    _h(geometry.getMinDeltaR())
{
  _u.normalize(_h);

  const int maxLatticeDistance = computeMaxLatticeDistance();
  // compute _hyperplane.origin, _nx, _ny so that the cuboid is right inside the geometry
  constructCuboid(maxLatticeDistance);
}

template<typename T>
HyperplaneLattice2D<T>::HyperplaneLattice2D(
  CuboidGeometry2D<T>& geometry, Hyperplane2D<T> hyperplane, int resolution)
  : _geometry(geometry),
    _hyperplane(hyperplane),
    _origin(hyperplane.origin),
    _u(hyperplane.u),
    _h(geometry.getMinDeltaR())
{
  _u.normalize(_h);

  const int maxLatticeDistance = computeMaxLatticeDistance();
  // compute _hyperplane.origin, _nx, _ny so that the cuboid is right inside the geometry
  constructCuboid(maxLatticeDistance);

  if ( resolution > 0 ) {
    setToResolution(resolution);
  }
}

template<typename T>
HyperplaneLattice2D<T>::HyperplaneLattice2D(
  CuboidGeometry2D<T>& geometry, Hyperplane2D<T> hyperplane, T h)
  : _geometry(geometry),
    _hyperplane(hyperplane),
    _origin(hyperplane.origin),
    _u(hyperplane.u),
    _h(h)
{
  if ( util::nearZero(_h) ) {
    _h = _geometry.getMinDeltaR();
  }

  _u.normalize(_h);

  const int maxLatticeDistance = computeMaxLatticeDistance();
  // compute _hyperplane.origin, _nx, _ny so that the cuboid is right inside the geometry
  constructCuboid(maxLatticeDistance);
}

template <typename T>
const Hyperplane2D<T>& HyperplaneLattice2D<T>::getHyperplane() const
{
  return _hyperplane;
}

template <typename T>
Vector<T,2> HyperplaneLattice2D<T>::getPhysR(const int& n) const
{
  return Vector<T,2> {
    _origin[0] + double(n)*_u[0],
    _origin[1] + double(n)*_u[1]
  };
}

template <typename T>
int HyperplaneLattice2D<T>::getN() const
{
  return _n;
}

template <typename T>
T HyperplaneLattice2D<T>::getPhysSpacing() const
{
  return _h;
}

template <typename T>
Vector<T,2> HyperplaneLattice2D<T>::getPhysOrigin() const
{
  return _origin;
}

template <typename T>
Vector<T,2> HyperplaneLattice2D<T>::getVectorU() const
{
  return _u;
}


}

#endif
