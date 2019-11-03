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
 * can be instantiated -- header file.
 */
#ifndef POROUS_ADVECTION_DIFFUSION_DYNAMICS_H
#define POROUS_ADVECTION_DIFFUSION_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "core/cell.h"

namespace olb {


// ===== the porous BGK advection diffusion dynamics =====//
/// This approach contains a slight error in the diffusion term.
template<typename T, typename DESCRIPTOR>
class PorousAdvectionDiffusionBGKdynamics : public BasicDynamics<T, DESCRIPTOR> {
public:
  /// Constructor
  PorousAdvectionDiffusionBGKdynamics( T omega, Momenta<T, DESCRIPTOR>& momenta, T tSolid );
  PorousAdvectionDiffusionBGKdynamics( const UnitConverter<T,DESCRIPTOR>& converter, Momenta<T, DESCRIPTOR>& momenta, T tSolid );
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr ) const override;
  /// Collision step
  void collide( Cell<T, DESCRIPTOR>& cell, LatticeStatistics<T>& statistics ) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega( T omega ) override;
  /// Set local relaxation parameter of the dynamics
  T computeRho( const Cell<T, DESCRIPTOR>& cell ) const override;
private:
  T scaleTemp(const T rho, const T porosity) const; // scales temperature relative to porosity
  T _omega;  ///< relaxation parameter
  T _tSolid; // temperature in lattice units of material with porosity 0;
};

} // namespace olb

#endif
