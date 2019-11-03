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
#ifndef ADVECTION_DIFFUSION_MOMENTA_HH
#define ADVECTION_DIFFUSION_MOMENTA_HH

#include "advectionDiffusionMomenta.h"
#include "lbHelpers.h"

namespace olb {

////////////////////// Class AdvectionDiffusionBulkMomenta //////////////////////////

template<typename T, typename DESCRIPTOR>
T AdvectionDiffusionBulkMomenta<T,DESCRIPTOR>::computeRho(Cell<T,DESCRIPTOR> const& cell) const
{
  return lbHelpers<T,DESCRIPTOR>::computeRho(cell);
}

template<typename T, typename DESCRIPTOR>
void AdvectionDiffusionBulkMomenta<T,DESCRIPTOR>::computeU(Cell<T,DESCRIPTOR> const& cell, T u[DESCRIPTOR::d]) const
{
  const T *u_ = cell.template getFieldPointer<descriptors::VELOCITY>();
  for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
    u[iD] = u_[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void AdvectionDiffusionBulkMomenta<T,DESCRIPTOR>::computeJ(Cell<T,DESCRIPTOR> const& cell, T j[DESCRIPTOR::d]) const
{
  lbHelpers<T,DESCRIPTOR>::computeJ(cell, j);
}

template<typename T, typename DESCRIPTOR>
void AdvectionDiffusionBulkMomenta<T,DESCRIPTOR>::computeStress (
  Cell<T,DESCRIPTOR> const& cell,
  T rho, const T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
}

template<typename T, typename DESCRIPTOR>
void AdvectionDiffusionBulkMomenta<T,DESCRIPTOR>::computeRhoU (
  Cell<T,DESCRIPTOR> const& cell,
  T& rho, T u[DESCRIPTOR::d] ) const
{
  rho = cell.computeRho();
  const T *u_ = cell.template getFieldPointer<descriptors::VELOCITY>();
  for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
    u[iD] = u_[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void AdvectionDiffusionBulkMomenta<T,DESCRIPTOR>::computeAllMomenta (
  Cell<T,DESCRIPTOR> const& cell,
  T& rho, T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  rho = cell.computeRho();
  const T *u_ = cell.template getFieldPointer<descriptors::VELOCITY>();
  for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
    u[iD] = u_[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void AdvectionDiffusionBulkMomenta<T,DESCRIPTOR>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{
  T *u = cell.template getFieldPointer<descriptors::VELOCITY>();
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = lbHelpers<T, DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
  }
}

template<typename T, typename DESCRIPTOR>
void AdvectionDiffusionBulkMomenta<T,DESCRIPTOR>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d])
{
  T rho = cell.computeRho();
  T *u_ = cell.template getFieldPointer<descriptors::VELOCITY>();
  for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
    u_[iD] = u[iD];
  }
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = lbHelpers<T, DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
  }
}

template<typename T, typename DESCRIPTOR>
void AdvectionDiffusionBulkMomenta<T,DESCRIPTOR>::defineRhoU (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d])
{
  T *u_ = cell.template getFieldPointer<descriptors::VELOCITY>();
  for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
    u_[iD] = u[iD];
  }
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = lbHelpers<T, DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
  }
}

template<typename T, typename DESCRIPTOR>
void AdvectionDiffusionBulkMomenta<T,DESCRIPTOR>::defineAllMomenta (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  const T pi[util::TensorVal<DESCRIPTOR >::n] )
{
  T *u_ = cell.template getFieldPointer<descriptors::VELOCITY>();
  for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
    u_[iD] = u[iD];
  }
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = lbHelpers<T, DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
  }
}

/////////////// Singletons //////////////////////////////////

namespace instances {

template<typename T, typename DESCRIPTOR>
AdvectionDiffusionBulkMomenta<T,DESCRIPTOR>& getAdvectionDiffusionBulkMomenta()
{
  static AdvectionDiffusionBulkMomenta<T,DESCRIPTOR> bulkMomentaSingleton;
  return bulkMomentaSingleton;
}
}

}

#endif
