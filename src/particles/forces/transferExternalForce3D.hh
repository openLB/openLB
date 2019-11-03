/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Marie-Luise Maier
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

#ifndef TRANSFEREXTERNALFORCE3D_HH
#define TRANSFEREXTERNALFORCE3D_HH

#include<cmath>
#include "transferExternalForce3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D;

template<typename T, template<typename U> class PARTICLETYPE>
TransferExternalForce3D<T, PARTICLETYPE>::TransferExternalForce3D()
  : Force3D<T, PARTICLETYPE>()
{
}

template<typename T, template<typename U> class PARTICLETYPE>
void TransferExternalForce3D<T, PARTICLETYPE>::applyForce(
  typename std::deque<PARTICLETYPE<T> >::iterator p, int pInt,
  ParticleSystem3D<T, PARTICLETYPE>& psSys)
{
  p->getForce()[0] += p->getStoreForce()[0];
  p->getForce()[1] += p->getStoreForce()[1];
  p->getForce()[2] += p->getStoreForce()[2];

  p->resetStoreForce();
}

}
#endif /* TRANSFEREXTERNALFORCE3D_HH */
