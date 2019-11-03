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
 * BGK Dynamics with adjustable speed of sound -- header file.
 */
#ifndef CHOPARD_DYNAMICS_H
#define CHOPARD_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "core/cell.h"

namespace olb {

/// Implementation of the BGK collision step
template<typename T, typename DESCRIPTOR>
class ChopardDynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ChopardDynamics(T vs2_, T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  ChopardDynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega_) override;
  /// Set local speed of sound
  void setVs2(T vs2_);
  /// Get local speed of sound
  T    getVs2() const;
public:
  static T chopardBgkCollision (
    Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T vs2, T omega);
  static T chopardEquilibrium (
    int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr, T vs2 );
private:
  T vs2;    ///< speed of sound
  T omega;  ///< relaxation parameter
};

}

#endif
