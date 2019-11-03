/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Robin Trunk, Mathias J. Krause
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

#ifndef BOUNDARYSIMPLEREFLECTION3D_HH
#define BOUNDARYSIMPLEREFLECTION3D_HH

#include <set>
#include "boundarySimpleReflection3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
SimpleReflectBoundary3D<T, PARTICLETYPE>::SimpleReflectBoundary3D(T dT, SuperGeometry3D<T>& sg, std::set<int> materials) :
  Boundary3D<T, PARTICLETYPE>(), _dT(dT), _sg(sg), _materials(materials) {
  _matIter = _materials.begin();
}

template<typename T, template<typename U> class PARTICLETYPE>
void SimpleReflectBoundary3D<T, PARTICLETYPE>::applyBoundary(typename std::deque<PARTICLETYPE<T> >::iterator& p, ParticleSystem3D<T, PARTICLETYPE>& psSys) {
  std::vector<T> oldPos(3,T());
  oldPos[0] = p->getPos()[0] - _dT * p->getVel()[0];
  oldPos[1] = p->getPos()[1] - _dT * p->getVel()[1];
  oldPos[2] = p->getPos()[2] - _dT * p->getVel()[2];

  std::vector<T> line(3,T());
  line[0] = p->getPos()[0] - oldPos[0];
  line[1] = p->getPos()[1] - oldPos[1];
  line[2] = p->getPos()[2] - oldPos[2];

  std::vector<T> reflection(3, T());
  int latticeR[4] = { 0,0,0,0 };
  bool outer = false;

  // get material number (Trigger)
  // TODO: take particle radius into account for collision detection
  latticeR[0] = p->getCuboid();
  _sg.getCuboidGeometry().get(latticeR[0]).getLatticeR(&(latticeR[1]),&p->getPos()[0]);
  int mat = _sg.get(latticeR);

  for (_matIter = _materials.begin(); _matIter != _materials.end();_matIter++) {
    if (mat == *_matIter) {
      // compute discrete normal
      // TODO: rearrange the if to make computation more efficient
      std::vector<int> normal(3, T());
      if (_sg.get(latticeR[0], latticeR[1] - 1, latticeR[2], latticeR[3]) == 1)   normal[0] = -1;
      if (_sg.get(latticeR[0], latticeR[1] + 1, latticeR[2], latticeR[3]) == 1)   normal[0] = 1;
      if (_sg.get(latticeR[0], latticeR[1], latticeR[2] - 1, latticeR[3]) == 1)   normal[1] = -1;
      if (_sg.get(latticeR[0], latticeR[1], latticeR[2] + 1, latticeR[3]) == 1)   normal[1] = 1;
      if (_sg.get(latticeR[0], latticeR[1], latticeR[2], latticeR[3] - 1) == 1)   normal[2] = -1;
      if (_sg.get(latticeR[0], latticeR[1], latticeR[2], latticeR[3] + 1) == 1)   normal[2] = 1;

      if (normal[0]==0 && normal[1]==0 && normal[2]==0) {
        outer = true;
        // check for outer edge
        if (_sg.get(latticeR[0], latticeR[1] - 1, latticeR[2] - 1, latticeR[3]) == 1)   {
          normal[0] = -1;
          normal[1] = -1;
        }
        else if (_sg.get(latticeR[0], latticeR[1] - 1, latticeR[2] + 1, latticeR[3]) == 1)   {
          normal[0] = -1;
          normal[1] = 1;
        }
        else if (_sg.get(latticeR[0], latticeR[1] + 1, latticeR[2] - 1, latticeR[3]) == 1)   {
          normal[0] = 1;
          normal[1] = -1;
        }
        else if (_sg.get(latticeR[0], latticeR[1] + 1, latticeR[2] + 1, latticeR[3]) == 1)   {
          normal[0] = 1;
          normal[1] = 1;
        }
        else if (_sg.get(latticeR[0], latticeR[1] - 1, latticeR[2], latticeR[3] - 1) == 1)   {
          normal[0] = -1;
          normal[2] = -1;
        }
        else if (_sg.get(latticeR[0], latticeR[1] - 1, latticeR[2], latticeR[3] + 1) == 1)   {
          normal[0] = -1;
          normal[2] = 1;
        }
        else if (_sg.get(latticeR[0], latticeR[1] + 1, latticeR[2], latticeR[3] - 1) == 1)   {
          normal[0] = 1;
          normal[2] = -1;
        }
        else if (_sg.get(latticeR[0], latticeR[1] + 1, latticeR[2], latticeR[3] + 1) == 1)   {
          normal[0] = 1;
          normal[2] = 1;
        }
        else if (_sg.get(latticeR[0], latticeR[1], latticeR[2] - 1, latticeR[3] - 1) == 1)   {
          normal[1] = -1;
          normal[2] = -1;
        }
        else if (_sg.get(latticeR[0], latticeR[1], latticeR[2] - 1, latticeR[3] + 1) == 1)   {
          normal[1] = -1;
          normal[2] = 1;
        }
        else if (_sg.get(latticeR[0], latticeR[1], latticeR[2] + 1, latticeR[3] - 1) == 1)   {
          normal[1] = +1;
          normal[2] = -1;
        }
        else if (_sg.get(latticeR[0], latticeR[1], latticeR[2] + 1, latticeR[3] + 1) == 1)   {
          normal[1] = +1;
          normal[2] = 1;
        }
        // check for outer corner
        else if (_sg.get(latticeR[0], latticeR[1] - 1, latticeR[2] - 1, latticeR[3] - 1) == 1)   {
          normal[0] = -1;
          normal[1] = -1;
          normal[2] = -1;
        }
        else if (_sg.get(latticeR[0], latticeR[1] - 1, latticeR[2] - 1, latticeR[3] + 1) == 1)   {
          normal[0] = -1;
          normal[1] = -1;
          normal[2] = 1;
        }
        else if (_sg.get(latticeR[0], latticeR[1] - 1, latticeR[2] + 1, latticeR[3] - 1) == 1)   {
          normal[0] = -1;
          normal[1] = +1;
          normal[2] = -1;
        }
        else if (_sg.get(latticeR[0], latticeR[1] - 1, latticeR[2] + 1, latticeR[3] + 1) == 1)   {
          normal[0] = -1;
          normal[1] = +1;
          normal[2] = 1;
        }
        else if (_sg.get(latticeR[0], latticeR[1] + 1, latticeR[2] - 1, latticeR[3] - 1) == 1)   {
          normal[0] = +1;
          normal[1] = -1;
          normal[2] = -1;
        }
        else if (_sg.get(latticeR[0], latticeR[1] + 1, latticeR[2] - 1, latticeR[3] + 1) == 1)   {
          normal[0] = +1;
          normal[1] = -1;
          normal[2] = 1;
        }
        else if (_sg.get(latticeR[0], latticeR[1] + 1, latticeR[2] + 1, latticeR[3] - 1) == 1)   {
          normal[0] = +1;
          normal[1] = +1;
          normal[2] = -1;
        }
        else if (_sg.get(latticeR[0], latticeR[1] + 1, latticeR[2] + 1, latticeR[3] + 1) == 1)   {
          normal[0] = +1;
          normal[1] = +1;
          normal[2] = 1;
        }
        // Error since normal is still zero
        else {
          std::cout << "----->>>>> ERROR Normal is ZERO" << std::endl;
          outer=false;
        }
      }

      T prod = line[0]*normal[0] + line[1]*normal[1] + line[2]*normal[2];
      // if particle has obtuse angle to normal do nothing
      if (!outer) {
        if (prod >= 0) {
          return;
        }
      }

      //  if particle has acute angle to normal reverse normal component of velocity
      T invNormNormal = 1. / util::norm2(normal);
      reflection[0] = line[0] - 2.*prod*normal[0]*invNormNormal;
      reflection[1] = line[1] - 2.*prod*normal[1]*invNormNormal;
      reflection[2] = line[2] - 2.*prod*normal[2]*invNormNormal;

      // compute new velocity vector
      T vel = util::norm(p->getVel());
      reflection = util::normalize(reflection);
      reflection[0] = reflection[0]*vel;
      reflection[1] = reflection[1]*vel;
      reflection[2] = reflection[2]*vel;

      // TODO: take distance to impact point into account, not just 0.5
      oldPos[0] += 0.5 * (line[0] + reflection[0]) * _dT;
      oldPos[1] += 0.5 * (line[1] + reflection[1]) * _dT;
      oldPos[2] += 0.5 * (line[2] + reflection[2]) * _dT;

      p->getVel() = reflection;
      p->setPos(oldPos);

      return;
    }
  }
}

}
#endif
