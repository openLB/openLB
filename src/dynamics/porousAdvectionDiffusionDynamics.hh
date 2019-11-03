/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef POROUS_ADVECTION_DIFFUSION_DYNAMICS_HH
#define POROUS_ADVECTION_DIFFUSION_DYNAMICS_HH

#include <algorithm>
#include <limits>

#include "porousAdvectionDiffusionDynamics.h"
#include "core/util.h"
#include "lbHelpers.h"

namespace olb {


//==================================================================//
//========== BGK Model for porous Advection diffusion=======//
//==================================================================//

template<typename T, typename DESCRIPTOR>
PorousAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::PorousAdvectionDiffusionBGKdynamics (
  T omega, Momenta<T, DESCRIPTOR>& momenta, T tSolid )
  : BasicDynamics<T, DESCRIPTOR>( momenta ),
    _omega(omega), _tSolid(tSolid)
{ }

template<typename T, typename DESCRIPTOR>
PorousAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::PorousAdvectionDiffusionBGKdynamics (
  const UnitConverter<T,DESCRIPTOR>& converter, Momenta<T, DESCRIPTOR>& momenta, T tSolid )
  : BasicDynamics<T, DESCRIPTOR>( momenta ),
    _omega(converter.getLatticeRelaxationFrequency()), _tSolid(tSolid)
{ }

template<typename T, typename DESCRIPTOR>
T PorousAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::computeEquilibrium( int iPop, T rho,
    const T u[DESCRIPTOR::d], T uSqr ) const
{
  // does temperature need to be considered here?
  return lbHelpers<T, DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
}


template<typename T, typename DESCRIPTOR>
void PorousAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::collide( Cell<T, DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics ) 
{
  T temperature = this->_momenta.computeRho( cell );

  // apply temperature scaling
  T* porosity = cell.template getFieldPointer<descriptors::POROSITY>();
  temperature = scaleTemp(temperature, porosity[0]);

  const T* u = cell.template getFieldPointer<descriptors::VELOCITY>();

  T uSqr = lbHelpers<T, DESCRIPTOR>::
           bgkCollision( cell, temperature, u, _omega );

  statistics.incrementStats( temperature, uSqr );
}

template<typename T, typename DESCRIPTOR>
T PorousAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
void PorousAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::setOmega( T omega )
{
  _omega = omega;
}

template<typename T, typename DESCRIPTOR>
T PorousAdvectionDiffusionBGKdynamics<T,DESCRIPTOR>::computeRho(const Cell<T,DESCRIPTOR>& cell) const
{
  //T rho = this->_momenta.computeRho(cell);
  //const T* porosity = cell.template getFieldPointer<descriptors::POROSITY>();
  //return scaleTemp(rho, porosity[0]);
  return this->_momenta.computeRho(cell);
}

template<typename T, typename DESCRIPTOR>
T PorousAdvectionDiffusionBGKdynamics<T,DESCRIPTOR>::scaleTemp(const T temp, const T porosity) const
{
  T rho = porosity * temp + ( T(1) - porosity ) * _tSolid;
  return rho;
}


} // namespace olb

#endif
