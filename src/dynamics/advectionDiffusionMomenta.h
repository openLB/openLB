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
#ifndef ADVECTION_DIFFUSION_MOMENTA_H
#define ADVECTION_DIFFUSION_MOMENTA_H

#include "latticeDescriptors.h"
#include "dynamics.h"

namespace olb {

/// Standard computation of velocity momenta in the bulk
template<typename T, typename DESCRIPTOR>
struct AdvectionDiffusionBulkMomenta : public Momenta<T,DESCRIPTOR> {
  /// Compute particle density on the cell.
  T computeRho(Cell<T,DESCRIPTOR> const& cell) const override;
  /// Compute fluid velocity on the cell.
  void computeU (
    Cell<T,DESCRIPTOR> const& cell,
    T u[DESCRIPTOR::d] ) const override;
  /// Compute fluid momentum on the cell.
  void computeJ (
    Cell<T,DESCRIPTOR> const& cell,
    T j[DESCRIPTOR::d] ) const override;
  /// Compute components of the stress tensor on the cell.
  void computeStress (
    Cell<T,DESCRIPTOR> const& cell,
    T rho, const T u[DESCRIPTOR::d],
    T pi[util::TensorVal<DESCRIPTOR >::n] ) const override;
  /// Compute fluid velocity and particle density on the cell.
  void computeRhoU (
    Cell<T,DESCRIPTOR> const& cell,
    T& rho, T u[DESCRIPTOR::d]) const override;
  /// Compute all momenta on the cell, up to second order.
  void computeAllMomenta (
    Cell<T,DESCRIPTOR> const& cell,
    T& rho, T u[DESCRIPTOR::d],
    T pi[util::TensorVal<DESCRIPTOR >::n] ) const override;
  /// Set particle density on the cell.
  void defineRho(Cell<T,DESCRIPTOR>& cell, T rho) override;
  /// Set fluid velocity on the cell.
  void defineU(Cell<T,DESCRIPTOR>& cell,
                       const T u[DESCRIPTOR::d]) override;
  /// Define fluid velocity and particle density on the cell.
  void defineRhoU (
    Cell<T,DESCRIPTOR>& cell,
    T rho, const T u[DESCRIPTOR::d]) override;
  /// Define all momenta on the cell, up to second order.
  void defineAllMomenta (
    Cell<T,DESCRIPTOR>& cell,
    T rho, const T u[DESCRIPTOR::d],
    const T pi[util::TensorVal<DESCRIPTOR >::n] ) override;
};

namespace instances {

template<typename T, typename DESCRIPTOR>
AdvectionDiffusionBulkMomenta<T,DESCRIPTOR>& getAdvectionDiffusionBulkMomenta();

}

}

#endif
