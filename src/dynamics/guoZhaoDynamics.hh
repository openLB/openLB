/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016-2017 Davide Dapelo, Mathias J. Krause
 *  OpenLB e-mail contact: info@openlb.net
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
 * Specific dynamics classes for Guo and Zhao (2002) porous model, with
 * which a Cell object can be instantiated -- generic implementation.
 */
#ifndef LB_GUOZHAO_DYNAMICS_HH
#define LB_GUOZHAO_DYNAMICS_HH

#include <algorithm>
#include <limits>
#include "dynamics/dynamics.h"
#include "core/cell.h"
#include "dynamics/guoZhaoLbHelpers.h"
#include "dynamics/firstOrderLbHelpers.h"
#include "dynamics/guoZhaoDynamics.h"

namespace olb {

////////////////////// Class GuoZhaoBGKdynamics /////////////////////////

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
GuoZhaoBGKdynamics<T,DESCRIPTOR>::GuoZhaoBGKdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta )
  : BasicDynamics<T,DESCRIPTOR>(momenta),
    _omega(omega)
{
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::FORCE>() );

  _epsilon = (T)1.0; // This to avoid a NaN error at the first timestep.
}

template<typename T, typename DESCRIPTOR>
T GuoZhaoBGKdynamics<T,DESCRIPTOR>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return GuoZhaoLbHelpers<T,DESCRIPTOR>::equilibrium(iPop, _epsilon, rho, u, uSqr);
}

template<typename T, typename DESCRIPTOR>
void GuoZhaoBGKdynamics<T,DESCRIPTOR>::computeU (Cell<T,DESCRIPTOR> const& cell, T u[DESCRIPTOR::d] ) const
{
  T rho;
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += cell.template getFieldPointer<descriptors::FORCE>()[iVel] / (T)2.;
  }
}

template<typename T, typename DESCRIPTOR>
void GuoZhaoBGKdynamics<T,DESCRIPTOR>::computeRhoU (Cell<T,DESCRIPTOR> const& cell, T& rho, T u[DESCRIPTOR::d] ) const
{
  this->_momenta.computeRhoU(cell, rho, u);

  T epsilon = cell.template getFieldPointer<descriptors::EPSILON>()[0];
  T nu = cell.template getFieldPointer<descriptors::NU>()[0];
  T k = cell.template getFieldPointer<descriptors::K>()[0];

  T Fe = 0.; //1.75/sqrt(150.*pow(epsilon,3));

  T c_0 = 0.5*(1 + 0.5*epsilon*nu/k);
  T c_1 = 0.5*epsilon*Fe/sqrt(k);

  T uMag = 0.;
  for (int iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
    uMag += u[iDim]*u[iDim];
  }
  uMag = sqrt(uMag);

  for (int iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
    u[iDim] /= (c_0 + sqrt(c_0*c_0 + c_1*uMag));
    u[iDim] += 0.5*epsilon*cell.template getFieldPointer<descriptors::BODY_FORCE>()[iDim];
  }
}

template<typename T, typename DESCRIPTOR>
void GuoZhaoBGKdynamics<T,DESCRIPTOR>::updateEpsilon (Cell<T,DESCRIPTOR>& cell)
{
  _epsilon = *cell.template getFieldPointer<descriptors::EPSILON>(); //Copying epsilon from
  // external to member variable to provide access for computeEquilibrium.
}


template<typename T, typename DESCRIPTOR>
void GuoZhaoBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  // Copying epsilon from
  // external to member variable to provide access for computeEquilibrium.
  updateEpsilon(cell);
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  GuoZhaoLbHelpers<T,DESCRIPTOR>::updateGuoZhaoForce(cell, u);
  T uSqr = GuoZhaoLbHelpers<T,DESCRIPTOR>::bgkCollision(cell, _epsilon, rho, u, _omega);
  GuoZhaoLbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, _omega, rho, _epsilon);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T GuoZhaoBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
T GuoZhaoBGKdynamics<T,DESCRIPTOR>::getEpsilon()
{
  return _epsilon;
}

template<typename T, typename DESCRIPTOR>
void GuoZhaoBGKdynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  _omega = omega;
}

}

#endif
