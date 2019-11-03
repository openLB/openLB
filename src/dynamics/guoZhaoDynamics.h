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
 * which a Cell object can be instantiated -- header file.
 */
#ifndef LB_GUOZHAO_DYNAMICS_H
#define LB_GUOZHAO_DYNAMICS_H

#include "dynamics/guoZhaoLatticeDescriptors.h"
#include "core/util.h"
#include "core/postProcessing.h"
#include "core/latticeStatistics.h"
#include "dynamics/dynamics.h"

namespace olb {

/// Implementation of the BGK collision step with porous force according to
/// Guo and Zhao (2012), described as an external force
template<typename T, typename DESCRIPTOR>
class GuoZhaoBGKdynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  GuoZhaoBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  ///  Compute fluid velocity on the cell.
  void computeU (
    Cell<T,DESCRIPTOR> const& cell,
    T u[DESCRIPTOR::d] ) const override;
  /// Compute fluid velocity and particle density on the cell.
  void computeRhoU (
    Cell<T,DESCRIPTOR> const& cell,
    T& rho, T u[DESCRIPTOR::d]) const override;
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega_) override;
  /// Get local porosity as per from the dynamics class member variable.
  T getEpsilon();
protected:
  /// Copies epsilon from external to member variable to provide access to computeEquilibrium.
  void updateEpsilon(Cell<T,DESCRIPTOR>& cell);
  T _omega;  ///< relaxation parameter
  T _epsilon; ///< porosity. Must be re-declared as a member variable to allow
};

}

#endif
