/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Marie-Luise Maier, Mathias J. Krause, Sascha Janz
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

/// Linear damping Force for magnetic dipolemoment in magnetic Field

#ifndef LinearDampingForceForMagDipoleMoment3D_HH
#define LinearDampingForceForMagDipoleMoment3D_HH

#include <cmath>
#include <vector>
#include "linearDampingForceForMagDipoleMoment3D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
LinearDampingForceForMagDipoleMoment3D<T, PARTICLETYPE, DESCRIPTOR>::LinearDampingForceForMagDipoleMoment3D(T dynVisc, T frictionFac):
  Force3D<T, PARTICLETYPE>(), _dynVisc(dynVisc), _frictionFac(frictionFac)
{
//  this->_name = "LinearDampingForceForMagDipoleMoment3D";
}

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
void LinearDampingForceForMagDipoleMoment3D<T, PARTICLETYPE, DESCRIPTOR>::applyForce(
  typename std::deque<PARTICLETYPE<T> >::iterator p, int pInt,
  ParticleSystem3D<T, PARTICLETYPE>& pSys)
{

  Vector<T, 3> dampingForce = {T(0), T(0), T(0)} ;

  for (int i = 0; i < 3; i++) {

    dampingForce[i] = - 8 * M_PI * pow(p->getRad(), 3.) * _dynVisc * p->getAVel()[i] * _frictionFac ;
    p->getTorque()[i] += dampingForce[i] ;
  }
}

}

#endif // LinearDampingForceForMagDipoleMoment3D_HH
