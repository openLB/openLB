/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn
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

#ifndef BUOYANCYFORCE_3D_HH
#define BUOYANCYFORCE_3D_HH

#include <cmath>
#include "buoyancyForce3D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
BuoyancyForce3D<T, PARTICLETYPE, DESCRIPTOR>::BuoyancyForce3D(
  UnitConverter<T,DESCRIPTOR> const& converter,
  std::vector<T> direction, T g) :
  Force3D<T, PARTICLETYPE>(),
  _direction(direction), _g(g), _physDensity(converter.getPhysDensity() )
{
  T directionNorm = sqrt(
                      pow(_direction[0], 2.) + pow(_direction[1], 2.)
                      + pow(_direction[2], 2.));
  for (int i = 0; i < 3; ++i) {
    _direction[i] /= directionNorm;
  }
}

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
void BuoyancyForce3D<T, PARTICLETYPE, DESCRIPTOR>::applyForce(
  typename std::deque<PARTICLETYPE<T> >::iterator p, int pInt,
  ParticleSystem3D<T, PARTICLETYPE>& psSys)
{
  // weight force of fluid in similar volume like particle that replaces the fluid
  // acts in direction opposite to weight force of particle

  T factor = 4. / 3. * M_PI * std::pow(p->getRad(), 3) * _g
             * _physDensity;
  for (int j = 0; j < 3; ++j) {
    p->getForce()[j] -= factor * _direction[j];
  }
}

}
#endif /* BUOYANCYFORCE_3D_HH */
