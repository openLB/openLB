/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Albert Mink, Christopher McHardy
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
 * A collection of radiative transport dynamics classes -- header file.
 */

#ifndef RTLBM_DYNAMICS_H
#define RTLBM_DYNAMICS_H

#include "dynamics/dynamics.h"

namespace olb {



/**
 * Solves RTE according Christopher McHardy et al 2016.
 * absorption and scattering coefficient:
 * \f$ \sigma_a \f$ and \f$ \sigma_s \f$
 *
 * \param omega             change into beta the extinction coefficient
 * \param singleScatAlbedo  is the single scattering albedo, given by \f$ \frac{\sigma_s}{\sigma_a + \sigma_s} \f$
 */
template<typename T, typename DESCRIPTOR>
class RTLBMdynamicsMcHardy : public BasicDynamics<T, DESCRIPTOR> {
public:
  /// Constructor
  RTLBMdynamicsMcHardy( Momenta<T,DESCRIPTOR>& momenta, T latticeAbsorption, T latticeScattering, std::array<std::array<T,DESCRIPTOR::q>, DESCRIPTOR::q>& anisoMatrix );
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr ) const override;
  /// Collision step
  void collide( Cell<T, DESCRIPTOR>& cell, LatticeStatistics<T>& statistics ) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega( T omega ) override;
  T getSink() const;

protected:
  T _absorption;
  T _scattering;
  std::array<std::array<T,DESCRIPTOR::q>, DESCRIPTOR::q>& _anisoMatrix;
};

template<typename T, typename DESCRIPTOR>
class RTLBMdynamicsMcHardyRK : public BasicDynamics<T, DESCRIPTOR> {
public:
  /// Constructor
  RTLBMdynamicsMcHardyRK( Momenta<T,DESCRIPTOR>& momenta, T latticeAbsorption, T latticeScattering, std::array<std::array<T,DESCRIPTOR::q>, DESCRIPTOR::q>& anisoMatrix );
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr ) const override;
  /// Collision step
  void collide( Cell<T, DESCRIPTOR>& cell, LatticeStatistics<T>& statistics ) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega( T omega ) override;
private:
  void computeEquilibriumAniso( Cell<T,DESCRIPTOR>& cell, std::array<T,DESCRIPTOR::q>& feq );
  std::array<T,DESCRIPTOR::q> doCollision( Cell<T,DESCRIPTOR>& cell, std::array<T,DESCRIPTOR::q>& feq );
  T _absorption;
  T _scattering;
  std::array<std::array<T,DESCRIPTOR::q>, DESCRIPTOR::q>& _anisoMatrix;
};


}  // namespace olb

#endif
