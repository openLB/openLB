/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Asher Zarth, Thomas Henn, Mathias J. Krause, Jonas Latt
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
 * BGK Dynamics for porous -- generic implementation.
 */
#ifndef POROUS_FORCED_BGK_DYNAMICS_HH
#define POROUS_FORCED_BGK_DYNAMICS_HH

#include "porousForcedBGKDynamics.h"
#include "core/cell.h"
#include "dynamics.h"
#include "core/util.h"
#include "lbHelpers.h"
#include "math.h"

namespace olb {

////////////////////// Class PorousForcedBGKdynamics //////////////////////////

template<typename T, typename DESCRIPTOR>
PorousForcedBGKdynamics<T,DESCRIPTOR>::PorousForcedBGKdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta)
  : BasicDynamics<T,DESCRIPTOR>(momenta),
    _omega(omega)
{
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::FORCE>() );
}

template<typename T, typename DESCRIPTOR>
void PorousForcedBGKdynamics<T,DESCRIPTOR>::computeU (Cell<T,DESCRIPTOR> const& cell, T u[DESCRIPTOR::d] ) const
{
  T rho;
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += cell.template getFieldPointer<descriptors::FORCE>()[iVel] / (T)2.;
  }
}

template<typename T, typename DESCRIPTOR>
void PorousForcedBGKdynamics<T,DESCRIPTOR>::computeRhoU (Cell<T,DESCRIPTOR> const& cell, T& rho, T u[DESCRIPTOR::d] ) const
{
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += cell.template getFieldPointer<descriptors::FORCE>()[iVel] / (T)2.;
  }
}


template<typename T, typename DESCRIPTOR>
void PorousForcedBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  T* porosity = cell.template getFieldPointer<descriptors::POROSITY>();
  for (int i=0; i<DESCRIPTOR::d; i++)  {
    u[i] *= porosity[0];
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, _omega);
  lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, _omega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T PorousForcedBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
void PorousForcedBGKdynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  _omega = omega;
}

} // olb

#endif
