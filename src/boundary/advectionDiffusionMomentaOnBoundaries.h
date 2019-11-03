/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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
 * Implementation of boundary cell dynamics -- header file.
 */
#ifndef ADE_MOMENTA_ON_BOUNDARIES_H
#define ADE_MOMENTA_ON_BOUNDARIES_H

#include "dynamics/dynamics.h"
#include "dynamics/advectionDiffusionMomenta.h"
#include "core/util.h"
#include "core/cell.h"

namespace olb {

/// Computation of velocity momenta on a velocity boundary
template<typename T, typename DESCRIPTOR>
class AdvectionDiffusionBoundaryMomenta : virtual public AdvectionDiffusionBulkMomenta<T,DESCRIPTOR> {
public:
  virtual void computeJneq( Cell<T,DESCRIPTOR> const& cell, T jNeq[DESCRIPTOR::d] ) const =0;
};

/// Computation of velocity momenta on a velocity boundary
template<typename T, typename DESCRIPTOR, int direction, int orientation>
class RegularizedTemperatureBM : public AdvectionDiffusionBoundaryMomenta<T,DESCRIPTOR> {
public:
  /// Default Constructor: initialization to u=0, rho=1
  RegularizedTemperatureBM(const T temperature = 0.5);

  T computeRho(Cell<T,DESCRIPTOR> const& cell) const override;
  void computeJ( Cell<T,DESCRIPTOR> const& cell, T j[DESCRIPTOR::d] ) const override;
  void computeJneq( Cell<T,DESCRIPTOR> const& cell, T jNeq[DESCRIPTOR::d] ) const override;
  void defineRho(Cell<T,DESCRIPTOR>& cell, T rho) override;
  void defineRhoU (Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d]) override;
private:
  T _temperature;
};


/// Computation of velocity momenta on a velocity boundary
template<typename T, typename DESCRIPTOR, int direction, int orientation>
class RegularizedHeatFluxBM : public AdvectionDiffusionBoundaryMomenta<T,DESCRIPTOR> {
public:
  /// Default Constructor: initialization to u=0, rho=1
  RegularizedHeatFluxBM(T *heatFlux = nullptr);

  T computeRho(Cell<T,DESCRIPTOR> const& cell) const override;
  void computeJ( Cell<T,DESCRIPTOR> const& cell, T j[DESCRIPTOR::d] ) const override;
  void computeJneq( Cell<T,DESCRIPTOR> const& cell, T jNeq[DESCRIPTOR::d] ) const override;
private:
  T _heatFlux[DESCRIPTOR::d];
};


}  // namespace olb


#endif

