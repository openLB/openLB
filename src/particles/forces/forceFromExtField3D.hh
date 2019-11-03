/*
 *  Copyright (C) 2018 Marie-Luise Maier
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

#ifndef ForceFromExtField3D_HH
#define ForceFromExtField3D_HH

#include "forceFromExtField3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
ForceFromExtField3D<T, PARTICLETYPE, DESCRIPTOR>::ForceFromExtField3D(
//    SuperLattice3D<T, DESCRIPTOR>& sLattice,
    //    SuperLatticeField3D<T, DESCRIPTOR>& sLatticeForceField,
    AnalyticalFfromSuperF3D<T>& analyticalExternalField,
    T scale) :
    Force3D<T, PARTICLETYPE>(),
//    _sLattice(sLattice),
    //_sLatticeForceField(sLatticeForceField),
    _analyticalExternalField(analyticalExternalField),
    _scale(scale)
{ }

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
void ForceFromExtField3D<T, PARTICLETYPE, DESCRIPTOR>::applyForce(
    typename std::deque<PARTICLETYPE<T> >::iterator p, int pInt,
    ParticleSystem3D<T, PARTICLETYPE>& pSys) {

//  int glob = p->getCuboid();
//  int latticePos[3] = { 0, 0, 0 };
//
//  _sLattice.getCuboidGeometry().get(p->getCuboid()).getLatticeR(latticePos,
//      &p->getPos()[0]);
//
//  int X[4] = { glob, latticePos[0], latticePos[1], latticePos[2] };
//
  T F[3] = { T(), T(), T() };
//  _sLatticeForceField(F, X);

  _analyticalExternalField(F, &p->getPos()[0]);
  p->getForce()[0] += F[0] * _scale;
  p->getForce()[1] += F[1] * _scale;
  p->getForce()[2] += F[2] * _scale;
}

}
#endif
