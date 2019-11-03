/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Jonas Latt
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

/** \file
 * BGK Dynamics with adjustable speed of sound -- generic implementation.
 */
#ifndef CHOPARD_DYNAMICS_HH
#define CHOPARD_DYNAMICS_HH

#include "chopardDynamics.h"
#include "core/util.h"

namespace olb {

////////////////////// Class ChopardDynamics //////////////////////////

/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
ChopardDynamics<T,DESCRIPTOR>::ChopardDynamics (
  T vs2_, T omega_, Momenta<T,DESCRIPTOR>& momenta_ )
  : BasicDynamics<T,DESCRIPTOR>(momenta_),
    vs2(vs2_),
    omega(omega_)
{ }

/** With this constructor, the speed of sound is vs2 = cs2
 *  \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
ChopardDynamics<T,DESCRIPTOR>::ChopardDynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_ )
  : BasicDynamics<T,DESCRIPTOR>(momenta_),
    vs2((T)1/descriptors::invCs2<T,DESCRIPTOR>()),
    omega(omega_)
{ }

template<typename T, typename DESCRIPTOR>
T ChopardDynamics<T,DESCRIPTOR>::computeEquilibrium
(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return chopardEquilibrium(iPop, rho, u, uSqr, vs2);
}

template<typename T, typename DESCRIPTOR>
void ChopardDynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T uSqr = chopardBgkCollision(cell, rho, u, vs2, omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T ChopardDynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void ChopardDynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  omega = omega_;
}

template<typename T, typename DESCRIPTOR>
T ChopardDynamics<T,DESCRIPTOR>::getVs2() const
{
  return vs2;
}

template<typename T, typename DESCRIPTOR>
void ChopardDynamics<T,DESCRIPTOR>::setVs2(T vs2_)
{
  vs2 = vs2_;
}

template<typename T, typename DESCRIPTOR>
T ChopardDynamics<T,DESCRIPTOR>::chopardBgkCollision (
  Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T vs2, T omega)
{
  const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] *= (T)1-omega;
    cell[iPop] += omega * chopardEquilibrium(iPop, rho, u, uSqr, vs2);
  }
  return uSqr;
}

template<typename T, typename DESCRIPTOR>
T ChopardDynamics<T,DESCRIPTOR>::chopardEquilibrium (
  int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr, T vs2)
{
  if (iPop==0) {
    return rho*( (T)1 -
                 vs2*descriptors::invCs2<T,DESCRIPTOR>()*((T)1-descriptors::t<T,DESCRIPTOR>(0)) -
                 descriptors::t<T,DESCRIPTOR>(0)/(T)2*descriptors::invCs2<T,DESCRIPTOR>()*uSqr ) - descriptors::t<T,DESCRIPTOR>(0);
  } else {
    T c_u = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
    }
    return rho * descriptors::t<T,DESCRIPTOR>(iPop) * descriptors::invCs2<T,DESCRIPTOR>()* (
             vs2 + c_u +
             descriptors::invCs2<T,DESCRIPTOR>() / (T)2 * c_u*c_u -
             uSqr / (T)2
           ) - descriptors::t<T,DESCRIPTOR>(iPop);
  }
}


}

#endif
