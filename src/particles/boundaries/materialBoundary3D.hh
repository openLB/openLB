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

#ifndef MATERIALBOUNDARY3D_HH
#define MATERIALBOUNDARY3D_HH

#include <set>
#include "materialBoundary3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
MaterialBoundary3D<T, PARTICLETYPE>::MaterialBoundary3D(
  SuperGeometry3D<T>& sg)
  : Boundary3D<T, PARTICLETYPE>(),
    _sg(sg)
{
  _matIter = _materials.begin();
}

template<typename T, template<typename U> class PARTICLETYPE>
MaterialBoundary3D<T, PARTICLETYPE>::MaterialBoundary3D(
  SuperGeometry3D<T>& sg,
  std::set<int> materials)
  : Boundary3D<T, PARTICLETYPE>(),
    _sg(sg),
    _materials(materials.begin(),materials.end())
{
}

template<typename T, template<typename U> class PARTICLETYPE>
void MaterialBoundary3D<T, PARTICLETYPE>::applyBoundary(
  typename std::deque<PARTICLETYPE<T> >::iterator& p,
  ParticleSystem3D<T, PARTICLETYPE>& psSys)
{
  int latticeR[3] = { 0 };
  _sg.getCuboidGeometry().get(p->getCuboid()).getFloorLatticeR(latticeR, &p->getPos()[0]);
  // Read only access to the material numbers of nodes around particle position
  const BlockGeometryStructure3D<T>& bg = _sg.getExtendedBlockGeometry(
      _sg.getLoadBalancer().loc(p->getCuboid()));
  // + overlap is because of lower boundaries, latticeR has to be shifted up
  int iX = latticeR[0]+_sg.getOverlap();
  int iY = latticeR[1]+_sg.getOverlap();
  int iZ = latticeR[2]+_sg.getOverlap();
  for (_matIter = _materials.begin(); _matIter != _materials.end(); _matIter++) {
    if (bg.get(iX, iY, iZ) == *_matIter ||
        bg.get(iX, iY+1, iZ) == *_matIter ||
        bg.get(iX, iY, iZ+1) == *_matIter ||
        bg.get(iX, iY+1, iZ+1) == *_matIter ||
        bg.get(iX+1, iY, iZ) == *_matIter ||
        bg.get(iX+1, iY+1, iZ) == *_matIter ||
        bg.get(iX+1, iY, iZ+1) == *_matIter ||
        bg.get(iX+1, iY+1, iZ+1) == *_matIter
       ) {
      p->setActive(false);
      return;
    }
  }
}

}

#endif /* MATERIALBOUNDARY3D_HH */
