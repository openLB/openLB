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

#ifndef STOKESDRAGFORCE_3D_HH
#define STOKESDRAGFORCE_3D_HH

#include<cmath>
#include "stokesDragForce3D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
StokesDragForce3D<T, PARTICLETYPE, DESCRIPTOR>::StokesDragForce3D(SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& getVel, T dT, T mu)
  : Force3D<T, PARTICLETYPE>(),
    _getVel(getVel),
    _mu(mu)
{
  _C1 = 6. * M_PI * _mu * dT ;
  _dTinv = 1. / dT;
  _scaleFactor = 1.;
}

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
StokesDragForce3D<T, PARTICLETYPE, DESCRIPTOR>::StokesDragForce3D(SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& getVel, UnitConverter<T, DESCRIPTOR> const& converter)
  : Force3D<T, PARTICLETYPE>(),
    _getVel(getVel)
{
  //implicit formulation
  _C1 = 6. * M_PI * converter.getPhysViscosity() * converter.getPhysDensity() * converter.getConversionFactorTime();
  _mu = converter.getPhysViscosity() * converter.getPhysDensity();
  // explicit formulation
  //  _C1 = 6. * M_PI * converter.getDynamicViscosity();
  _dTinv = 1. / converter.getConversionFactorTime();
  _scaleFactor = 1. ;
}

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
StokesDragForce3D<T, PARTICLETYPE, DESCRIPTOR>::StokesDragForce3D(SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& getVel, UnitConverter<T, DESCRIPTOR> const& converter, T scaleFactor)
  : Force3D<T, PARTICLETYPE>(),
    _getVel(getVel),
    _scaleFactor(scaleFactor)
{
  //implicit formulation
  _C1 = 6. * M_PI * converter.getPhysViscosity() * converter.getPhysDensity() * converter.getConversionFactorTime();
  _mu = converter.getPhysViscosity() * converter.getPhysDensity();
  // explicit formulation
  //  _C1 = 6. * M_PI * converter.getDynamicViscosity();
  _dTinv = 1. / converter.getConversionFactorTime();
}

/// 6 Pi r mu (u_f-u_p)
template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
void StokesDragForce3D<T, PARTICLETYPE, DESCRIPTOR>::applyForce(
  typename std::deque<PARTICLETYPE<T> >::iterator p, int pInt,
  ParticleSystem3D<T, PARTICLETYPE>& psSys)
{
  T fluidVel[3] = {0., 0., 0.};

  //implicit formulation
  _getVel(fluidVel, &p->getPos()[0], p->getCuboid());
  fluidVel[0] *= _scaleFactor;
  fluidVel[1] *= _scaleFactor;
  fluidVel[2] *= _scaleFactor;

  T c = _C1 * p->getRad() * p->getInvMass();
  T C2 = 1. / (1. + c);

  // p->getVel() is particle velocity of the last time step
  // formulation of new force with particle velocity of the last time step with
  // implicit Euler
  p->getForce()[0] += p->getMass() * _dTinv
                      * ((c * fluidVel[0] + p->getVel()[0]) * C2 - p->getVel()[0]);
  p->getForce()[1] += p->getMass() * _dTinv
                      * ((c * fluidVel[1] + p->getVel()[1]) * C2 - p->getVel()[1]);
  p->getForce()[2] += p->getMass() * _dTinv
                      * ((c * fluidVel[2] + p->getVel()[2]) * C2 - p->getVel()[2]);

  // explicit formulation
  //  T cex = 6. * M_PI * _mu * p->getRad();
  //  p->getForce()[0] += cex * (fluidVel[0]-p->getVel()[0]);
  //  p->getForce()[1] += cex * (fluidVel[1]-p->getVel()[1]);
  //  p->getForce()[2] += cex * (fluidVel[2]-p->getVel()[2]);
}

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
void StokesDragForce3D<T, PARTICLETYPE, DESCRIPTOR>::computeForce(
  int pInt, ParticleSystem3D<T, PARTICLETYPE>* psSys, T force[3])
{
  T fluidVel[3] = {0., 0., 0.};

  _getVel(fluidVel, &psSys->operator[](pInt).getPos()[0], psSys->operator[](pInt).getCuboid());

  T c = _C1 * psSys->operator[](pInt).getRad() * psSys->operator[](pInt).getInvMass();
  T C2 = 1. / (1. + c);
  T mass = psSys->operator[](pInt).getMass();
  std::vector<T> vel = psSys->operator[](pInt).getVel();
  force[0] = mass * _dTinv * ((c * fluidVel[0] + vel[0]) * C2 - vel[0]);
  force[1] = mass * _dTinv * ((c * fluidVel[1] + vel[1]) * C2 - vel[1]);
  force[2] = mass * _dTinv * ((c * fluidVel[2] + vel[2]) * C2 - vel[2]);
}

}
#endif /* STOKESDRAGFORCE_3D_HH */
