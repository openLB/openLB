/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Robin Trunk, Sam Avis
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

#ifndef FREE_ENERGY_DYNAMICS_H
#define FREE_ENERGY_DYNAMICS_H

#include "dynamics/dynamics.h"

/** \file
 * In this file the dynamic calls for the free energy model is implemented. It
 * is used for the second (and third) lattices, as for the first one a BGK collision with
 * Guo forcing is applied (see ForcedBGKdynamcs).
 */

namespace olb {

template<typename T, typename DESCRIPTOR>
class FreeEnergyBGKdynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// This dynamics describes the propagation of density(fluid1) - density(fluid2). And is
  /// used for the second (and third) lattices in the free energy model.
  /// \param[in] omega_ - lattice relaxation frequency [lattice units]
  /// \param[in] gamma_ - tunable parameter for the equilibrium distribution [lattice units]
  /// \param[in] momenta_ - momenta object describing the calculation of macroscopic values (e.g. rho and u).
  ///                       Usually "BulkMomenta" are used.
  FreeEnergyBGKdynamics(T omega_, T gamma_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Compute equilibrium distribution function.
  /// This should contain an additional term that depends upon the chemical potential. However the external field
  /// cannot be accessed when iniEquilibrium is called and so this has been neglected.
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega_) override;
  /// Compute fluid velocity and particle density on the cell.
  void computeRhoU(
    Cell<T,DESCRIPTOR> const& cell,
    T& rho, T u[DESCRIPTOR::d]) const override;
  /// Compute fluid velocity on the cell.
  void computeU(Cell<T,DESCRIPTOR> const& cell, T u[DESCRIPTOR::d]) const override;
private:
  T omega; /// relaxation parameter
  T gamma; /// tunable parameter
};

template<typename T, typename DESCRIPTOR>
class FreeEnergyWallDynamics : public BounceBack<T,DESCRIPTOR> {
public:
  /// This dynamics is used for the second (and third) lattices in the free energy model at wall boundaries.
  /// It is neccessary for returning the correct equilibrium distributions when iniEquilibrium is called.
  FreeEnergyWallDynamics();
  /// Compute equilibrium distribution function.
  /// This should contain an additional term that depends upon the chemical potential. However the external field
  /// cannot be accessed when iniEquilibrium is called and so this has been neglected.
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega_) override;
};

template<typename T, typename DESCRIPTOR, int direction, int orientation>
class FreeEnergyInletOutletDynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// This dynamics is used for the second (and third) lattices in the free energy model at inlets.
  /// It first defines the missing distribution functions and then performs the normal collision.
  /// \param[in] omega_ - lattice relaxation frequency [lattice units]
  /// \param[in] momenta_ - momenta object describing the calculation of macroscopic values (not including rho).
  ///                       Usually "BulkMomenta" are used.
  FreeEnergyInletOutletDynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics_) override;
  /// Compute equilibrium distribution function.
  /// This should contain an additional term that depends upon the chemical potential. However the external field
  /// cannot be accessed when iniEquilibrium is called and so this has been neglected.
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega_) override;
  /// Compute particle density on the cell.
  T computeRho(Cell<T,DESCRIPTOR> const& cell) const override;
  /// Compute fluid velocity on the cell.
  void computeU(Cell<T,DESCRIPTOR > const &cell, T u[DESCRIPTOR::d]) const override;
  /// Compute fluid velocity and particle density on the cell. 
  void computeRhoU (Cell<T,DESCRIPTOR> const& cell, T& rho, T u[DESCRIPTOR::d]) const override;
  /// Set particle density on the cell.
  void defineRho(Cell<T,DESCRIPTOR>& cell, T rho) override;
  /// Set fluid velocity on the cell.
  void defineU(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d]) override;
private:
  T omega;
};

}

#endif
