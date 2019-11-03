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

#ifndef WireBoundaryForMagP3D_HH
#define WireBoundaryForMagP3D_HH

#include <set>
#include "wireBoundaryForMagP3D.h"

namespace olb {


template<typename T, template<typename U> class PARTICLETYPE>
WireBoundaryForMagP3D<T, PARTICLETYPE>::WireBoundaryForMagP3D(
  SuperGeometry3D<T>& sg, std::set<int> materials)
  : Boundary3D<T, PARTICLETYPE>(),
    _sg(sg),
    _materials(materials.begin(), materials.end())
{
}

template<typename T, template<typename U> class PARTICLETYPE>
void WireBoundaryForMagP3D<T, PARTICLETYPE>::applyBoundary(
  typename std::deque<PARTICLETYPE<T> >::iterator& p,
  ParticleSystem3D<T, PARTICLETYPE>& psSys)
{

  int latticeR[3] = { 0 };
  _sg.getCuboidGeometry().get(p->getCuboid()).getFloorLatticeR(latticeR, &p->getPos()[0]);

  // Read only access to the material numbers of nodes around particle position
  const BlockGeometryStructure3D<T>& bg = _sg.getExtendedBlockGeometry(
      _sg.getLoadBalancer().loc(p->getCuboid()));

  // + overlap is because of lower boundaries, latticeR has to be shifted up
  int iX = latticeR[0] + _sg.getOverlap();
  int iY = latticeR[1] + _sg.getOverlap();
  int iZ = latticeR[2] + _sg.getOverlap();
  for (_matIter = _materials.begin(); _matIter != _materials.end(); _matIter++) {

    if (bg.get(iX, iY, iZ) == *_matIter ||
        bg.get(iX, iY + 1, iZ) == *_matIter ||
        bg.get(iX, iY, iZ + 1) == *_matIter ||
        bg.get(iX, iY + 1, iZ + 1) == *_matIter ||
        bg.get(iX + 1, iY, iZ) == *_matIter ||
        bg.get(iX + 1, iY + 1, iZ) == *_matIter ||
        bg.get(iX + 1, iY, iZ + 1) == *_matIter ||
        bg.get(iX + 1, iY + 1, iZ + 1) == *_matIter) {

      if ((*_matIter == 5) && (p->getSActivity() != 3)) {

        std::vector<T> vel(3, T()) ;
        p->setVel(vel) ;
        p->setSActivity(3) ;
        break;
      }
      if ((*_matIter == 4) && (p->getActive() != false)) {

        p->setActive(false) ;
        break;
      }
    }
  }
  return;
}
}

#endif /* WireBoundaryForMagP3D_HH */
