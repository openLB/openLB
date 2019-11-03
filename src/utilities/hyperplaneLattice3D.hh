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

#ifndef HYPERPLANE_LATTICE_3D_HH
#define HYPERPLANE_LATTICE_3D_HH

#include "hyperplaneLattice3D.h"

namespace olb {

template <typename T>
int HyperplaneLattice3D<T>::computeMaxLatticeDistance(Cuboid3D<T>&& cuboid) const
{
  const Vector<T,3>   origin = cuboid.getOrigin();
  const Vector<int,3> extend = cuboid.getExtend();
  const T             deltaR = cuboid.getDeltaR();

  T maxPhysDistance = T();
  T tmp;
  Vector<T,3> tmpO;
  Vector<T,3> tmpE;

  for(int iDim=0; iDim<3; ++iDim){
  tmpO[iDim] = origin[iDim] - _origin[iDim];
  tmpE[iDim] = origin[iDim] + extend[iDim]*deltaR - _origin[iDim];
  }
  tmp = sqrt(tmpO[0]*tmpO[0] + tmpO[1]*tmpO[1] + tmpO[2]*tmpO[2]);
  if (maxPhysDistance < tmp) {
      maxPhysDistance = tmp;
  }
  tmp = sqrt((tmpE[0]*tmpE[0] + tmpO[1]*tmpO[1] + tmpO[2]*tmpO[2]));
  if (maxPhysDistance < tmp) {
      maxPhysDistance = tmp;
  }
  tmp = sqrt(tmpO[0]*tmpO[0] + tmpE[1]*tmpE[1] + tmpO[2]*tmpO[2]);
  if (maxPhysDistance < tmp) {
      maxPhysDistance = tmp;
  }
  tmp = sqrt(tmpO[0]*tmpO[0] + tmpO[1]*tmpO[1] + tmpE[2]*tmpE[2]);
  if (maxPhysDistance < tmp) {
      maxPhysDistance = tmp;
  }
  tmp = sqrt(tmpO[0]*tmpO[0] + tmpE[1]*tmpE[1] + tmpE[2]*tmpE[2]);
  if (maxPhysDistance < tmp) {
      maxPhysDistance = tmp;
  }
  tmp = sqrt(tmpE[0]*tmpE[0] + tmpO[1]*tmpO[1] + tmpE[2]*tmpE[2]);
  if (maxPhysDistance < tmp) {
      maxPhysDistance = tmp;
  }
  tmp = sqrt(tmpE[0]*tmpE[0] + tmpE[1]*tmpE[1] + tmpO[2]*tmpO[2]);
  if (maxPhysDistance < tmp) {
      maxPhysDistance = tmp;
  }
  tmp = sqrt(tmpE[0]*tmpE[0] + tmpE[1]*tmpE[1] + tmpE[2]*tmpE[2]);
  if (maxPhysDistance < tmp) {
      maxPhysDistance = tmp;
  }

  return int(maxPhysDistance/_h) + 1;
}

template <typename T>
void HyperplaneLattice3D<T>::constructCuboid(
  CuboidGeometry3D<T>& geometry, int maxLatticeDistance)
{
  int iC;
  int minX = -maxLatticeDistance;
  int maxX = maxLatticeDistance;
  int minY = -maxLatticeDistance;
  int maxY = maxLatticeDistance;
  bool found = false;

  for ( int iX = -maxLatticeDistance; iX < maxLatticeDistance; ++iX ) {
    for ( int iY = -maxLatticeDistance; iY < maxLatticeDistance; ++iY ) {
      if ( geometry.getC(getPhysR(iX, iY), iC) ) {
        minX = iX;
        found = true;
        break;
      }
    }
    if (found) {
      break;
    }
  }
  found = false;
  for ( int iX = maxLatticeDistance; iX > -maxLatticeDistance; --iX ) {
    for ( int iY = -maxLatticeDistance; iY < maxLatticeDistance; ++iY ) {
      if ( geometry.getC(getPhysR(iX, iY), iC) ) {
        maxX = iX;
        found = true;
        break;
      }
    }
    if (found) {
      break;
    }
  }
  found = false;
  for ( int iY = -maxLatticeDistance; iY < maxLatticeDistance; ++iY ) {
    for ( int iX = -maxLatticeDistance; iX < maxLatticeDistance; ++iX ) {
      if ( geometry.getC(getPhysR(iX, iY), iC) ) {
        minY = iY;
        found = true;
        break;
      }
    }
    if (found) {
      break;
    }
  }
  found = false;
  for ( int iY = maxLatticeDistance; iY > -maxLatticeDistance; --iY ) {
    for ( int iX = -maxLatticeDistance; iX < maxLatticeDistance; ++iX ) {
      if ( geometry.getC(getPhysR(iX, iY), iC) ) {
        maxY = iY;
        found = true;
        break;
      }
    }
    if (found) {
      break;
    }
  }

  _nx = maxX - minX + 1;
  _ny = maxY - minY + 1;

  _origin = _origin + double(minX)*_u + double(minY)*_v;
}

template <typename T>
void HyperplaneLattice3D<T>::setToResolution(int resolution)
{
  if (_nx > _ny) {
    T newH = _nx*_h/(T) resolution;
    _nx = resolution;
    _ny = (int) (_ny*_h/newH) + 1;
    _h = newH;
  }
  else {
    T newH = _ny*_h/(T) resolution;
    _ny = resolution;
    _nx = (int) (_nx*_h/newH) + 1;
    _h = newH;
  }
  _u.normalize(_h);
  _v.normalize(_h);
}

template<typename T>
HyperplaneLattice3D<T>::HyperplaneLattice3D(
  CuboidGeometry3D<T>& geometry, Hyperplane3D<T> hyperplane)
  : _hyperplane(hyperplane),
    _origin(hyperplane.origin),
    _u(hyperplane.u),
    _v(hyperplane.v),
    _h(geometry.getMinDeltaR())
{
  _u.normalize(_h);
  _v.normalize(_h);

  const int maxLatticeDistance = computeMaxLatticeDistance(geometry.getMotherCuboid());
  // compute _hyperplane.origin, _nx, _ny so that the cuboid is right inside the geometry
  constructCuboid(geometry, maxLatticeDistance);
}

template<typename T>
HyperplaneLattice3D<T>::HyperplaneLattice3D(
  CuboidGeometry3D<T>& geometry, Hyperplane3D<T> hyperplane, int resolution)
  : _hyperplane(hyperplane),
    _origin(hyperplane.origin),
    _u(hyperplane.u),
    _v(hyperplane.v),
    _h(geometry.getMinDeltaR())
{
  _u.normalize(_h);
  _v.normalize(_h);

  const int maxLatticeDistance = computeMaxLatticeDistance(geometry.getMotherCuboid());
  // compute _hyperplane.origin, _nx, _ny so that the cuboid is right inside the geometry
  constructCuboid(geometry, maxLatticeDistance);

  if ( resolution > 0 ) {
    setToResolution(resolution);
  }
}

template<typename T>
HyperplaneLattice3D<T>::HyperplaneLattice3D(
  CuboidGeometry3D<T>& geometry, Hyperplane3D<T> hyperplane, T h)
  : _hyperplane(hyperplane),
    _origin(hyperplane.origin),
    _u(hyperplane.u),
    _v(hyperplane.v),
    _h(h)
{
  if ( util::nearZero(_h) ) {
    _h = geometry.getMinDeltaR();
  }

  _u.normalize(_h);
  _v.normalize(_h);

  const int maxLatticeDistance = computeMaxLatticeDistance(geometry.getMotherCuboid());
  // compute _hyperplane.origin, _nx, _ny so that the cuboid is right inside the geometry
  constructCuboid(geometry, maxLatticeDistance);
}

template<typename T>
HyperplaneLattice3D<T>::HyperplaneLattice3D(
    Hyperplane3D<T> hyperplane,
    T h, int nx, int ny)
  : _hyperplane(hyperplane),
    _origin(hyperplane.origin),
    _u(hyperplane.u),
    _v(hyperplane.v),
    _h(h),
    _nx(nx),
    _ny(ny)
{
  _u.normalize(_h);
  _v.normalize(_h);
}

template <typename T>
const Hyperplane3D<T>& HyperplaneLattice3D<T>::getHyperplane() const
{
  return _hyperplane;
}

template <typename T>
Vector<T,3> HyperplaneLattice3D<T>::getPhysR(const int& planeX, const int& planeY) const
{
  return Vector<T,3> {
    _origin[0] + double(planeX)*_u[0] + double(planeY)*_v[0],
    _origin[1] + double(planeX)*_u[1] + double(planeY)*_v[1],
    _origin[2] + double(planeX)*_u[2] + double(planeY)*_v[2]
  };
}

template <typename T>
int HyperplaneLattice3D<T>::getNx() const
{
  return _nx;
}

template <typename T>
int HyperplaneLattice3D<T>::getNy() const
{
  return _ny;
}

template <typename T>
T HyperplaneLattice3D<T>::getPhysSpacing() const
{
  return _h;
}

template <typename T>
Vector<T,3> HyperplaneLattice3D<T>::getPhysOrigin() const
{
  return _origin;
}

template <typename T>
Vector<T,3> HyperplaneLattice3D<T>::getVectorU() const
{
  return _u;
}

template <typename T>
Vector<T,3> HyperplaneLattice3D<T>::getVectorV() const
{
  return _v;
}

}

#endif
