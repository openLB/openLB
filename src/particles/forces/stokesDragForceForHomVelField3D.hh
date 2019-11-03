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
#ifndef StokesDragForceForHomVelField3D_HH
#define StokesDragForceForHomVelField3D_HH

#include <cmath>
#include "stokesDragForceForHomVelField3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
StokesDragForceForHomVelField3D<T, PARTICLETYPE, DESCRIPTOR>::StokesDragForceForHomVelField3D(T dynVisc, Vector<T, 3> fluidVel):
  Force3D<T, PARTICLETYPE>(), _dynVisc(dynVisc), _fluidVel(fluidVel)
{
//  this->_name = "StokesDragForceForHomVelField3D";
}

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
void StokesDragForceForHomVelField3D<T, PARTICLETYPE, DESCRIPTOR>::applyForce(
  typename std::deque<PARTICLETYPE<T> >::iterator p, int pInt,
  ParticleSystem3D<T, PARTICLETYPE>& pSys)
{

  Vector<T, 3> force = {T(0), T(0), T(0)} ;

  for (int i = 0; i < 3; i++) {

    force[i] = -1 * 6 * M_PI * p->getRad() * (p->getVel()[i] - _fluidVel[i]) * _dynVisc ;
    p->getForce()[i] += force[i] ;
  }

}
}

#endif // StokesDragForceForHomVelField3D
