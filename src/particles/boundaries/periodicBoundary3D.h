/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Thomas Henn, Mathias J. Krause
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

#ifndef PERIODICBOUNDARY3D_H_
#define PERIODICBOUNDARY3D_H_

#include <math.h>
#include <vector>

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D;

/*
 * Particle boundary based on a cube around the area with material number 1.
 * Only applicable to rectangles since if a particle leaves the area with
 * material number 1 it is moved to the opposing side of the area by
 * newPosition = oldPosition +/- extend(MaterialNumber=1).
 **/

template<typename T, template<typename U> class PARTICLETYPE>
class PeriodicBoundary3D : public Boundary3D<T, PARTICLETYPE> {
public:
  PeriodicBoundary3D(SuperGeometry3D<T> sg, bool x,
                     bool y, bool z);
  PeriodicBoundary3D(PeriodicBoundary3D<T, PARTICLETYPE>& f);
  virtual ~PeriodicBoundary3D() { };
  virtual void applyBoundary(typename std::deque<PARTICLETYPE<T> >::iterator& p, ParticleSystem3D<T, PARTICLETYPE>& psSys);
  /// Returns number of particles that moved through the periodic boundary
  /// Order: x+, x-, y+, y-, z+, z-
  unsigned int* getJumper();
private:
  //cube extents with origin (0,0,0)
  std::vector<T> _minPhys, _maxPhys, _extend;
  bool _x, _y, _z;
  unsigned int _jumper[6];
  CuboidGeometry3D<T>& _cuboidGeometry;
  T _overlap;
};

template<typename T, template<typename U> class PARTICLETYPE>
PeriodicBoundary3D<T, PARTICLETYPE>::PeriodicBoundary3D(
  SuperGeometry3D<T> sg, bool x, bool y, bool z) : Boundary3D<T, PARTICLETYPE>(),
  _minPhys(3, T()), _maxPhys(3, T()), _extend(3, T()),
  _x(x),
  _y(y),
  _z(z), _cuboidGeometry(sg.getCuboidGeometry())
{
  _minPhys = sg.getStatistics().getMinPhysR(1);
  _maxPhys = sg.getStatistics().getMaxPhysR(1);
  _extend[0] = _maxPhys[0] - _minPhys[0];
  _extend[1] = _maxPhys[1] - _minPhys[1];
  _extend[2] = _maxPhys[2] - _minPhys[2];
  for (int i=0; i<6; ++i) {
    _jumper[i] = 0;
  }
  _overlap = sg.getOverlap();
}

template<typename T, template<typename U> class PARTICLETYPE>
void PeriodicBoundary3D<T, PARTICLETYPE>::applyBoundary(
  typename std::deque<PARTICLETYPE<T> >::iterator& p,
  ParticleSystem3D<T, PARTICLETYPE>& psSys)
{
  if (_x) {
    if (p->getPos()[0] > _maxPhys[0]) {
      p->getPos()[0] -= _extend[0];
      ++_jumper[0];
      int C = this->_cuboidGeometry.get_iC(p->getPos()[0], p->getPos()[1], p->getPos()[2], _overlap);
      p->setCuboid(C);
    } else if (p->getPos()[0] < _minPhys[0]) {
      p->getPos()[0] += _extend[0];
      ++_jumper[1];
      int C = this->_cuboidGeometry.get_iC(p->getPos()[0], p->getPos()[1], p->getPos()[2], _overlap);
      p->setCuboid(C);
    }
  }
  if (_y) {
    if (p->getPos()[1] > _maxPhys[1]) {
      p->getPos()[1] -= _extend[1];
      ++_jumper[2];
      int C = this->_cuboidGeometry.get_iC(p->getPos()[0], p->getPos()[1], p->getPos()[2], _overlap);
      p->setCuboid(C);
    } else if (p->getPos()[1] < _minPhys[1]) {
      p->getPos()[1] += _extend[1];
      ++_jumper[3];
      int C = this->_cuboidGeometry.get_iC(p->getPos()[0], p->getPos()[1], p->getPos()[2], _overlap);
      p->setCuboid(C);
    }
  }
  if (_z) {
    if (p->getPos()[2] > _maxPhys[2]) {
      p->getPos()[2] -= _extend[2];
      ++_jumper[4];
      int C = this->_cuboidGeometry.get_iC(p->getPos()[0], p->getPos()[1], p->getPos()[2], _overlap);
      p->setCuboid(C);
    } else if (p->getPos()[2] < _minPhys[2]) {
      p->getPos()[2] += _extend[2];
      ++_jumper[5];
      int C = this->_cuboidGeometry.get_iC(p->getPos()[0], p->getPos()[1], p->getPos()[2], _overlap);
      p->setCuboid(C);
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
unsigned int* PeriodicBoundary3D<T, PARTICLETYPE>::getJumper()
{
  return _jumper;
}


}

#endif
