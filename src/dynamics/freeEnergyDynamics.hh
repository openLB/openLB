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

#ifndef FREE_ENERGY_DYNAMICS_HH
#define FREE_ENERGY_DYNAMICS_HH

#include "core/cell.h"
#include "dynamics/firstOrderLbHelpers.h"
#include "dynamics/freeEnergyDynamics.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
FreeEnergyBGKdynamics<T,DESCRIPTOR>::FreeEnergyBGKdynamics (
  T omega_, T gamma_, Momenta<T,DESCRIPTOR>& momenta )
  : BasicDynamics<T,DESCRIPTOR>(momenta),
    omega(omega_), gamma(gamma_)
{ }

template<typename T, typename DESCRIPTOR>
void FreeEnergyBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  for(int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
    u[iDim] = cell.template getFieldPointer<descriptors::FORCE>()[iDim];
  }
  rho = this->_momenta.computeRho(cell);
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  statistics.incrementStats(rho, uSqr);
  T tmp = gamma * descriptors::invCs2<T,DESCRIPTOR>() * cell.template getFieldPointer<descriptors::CHEM_POTENTIAL>()[0];
  for(int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] -= omega * descriptors::t<T,DESCRIPTOR>(iPop) * (rho - tmp);
  }
  cell[0] += omega * (1.-descriptors::t<T,DESCRIPTOR>(0)) * (rho - tmp);
}

template<typename T, typename DESCRIPTOR>
T FreeEnergyBGKdynamics<T,DESCRIPTOR>::computeEquilibrium(
  int iPop, T rho, const T u[DESCRIPTOR::d], T ) const
{
  T uSqr = 0.;
  T equilibrium = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
  T tmp = 0; //gamma * descriptors::invCs2<T,DESCRIPTOR>() * cell.template getFieldPointer<descriptors::CHEM_POTENTIAL>()[0];
  if (iPop == 0) {
    equilibrium += (1.-descriptors::t<T,DESCRIPTOR>(0)) * (rho - tmp);
  } else {
    equilibrium -= descriptors::t<T,DESCRIPTOR>(iPop) * (rho - tmp);
  }
  return equilibrium;
}

template<typename T, typename DESCRIPTOR>
T FreeEnergyBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyBGKdynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  omega = omega_;
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyBGKdynamics<T,DESCRIPTOR>::computeRhoU (
  Cell<T,DESCRIPTOR> const& cell, T& rho, T u[DESCRIPTOR::d]) const
{
  rho = this->_momenta.computeRho(cell);
  computeU(cell, u);
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyBGKdynamics<T,DESCRIPTOR>::computeU (
  Cell<T,DESCRIPTOR> const& cell, T u[DESCRIPTOR::d]) const
{
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] = cell.template getFieldPointer<descriptors::FORCE>()[iVel];
  }
}


/////////////// FreeEnergyWallDynamics ////////////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyWallDynamics<T,DESCRIPTOR>::FreeEnergyWallDynamics ()
  : BounceBack<T,DESCRIPTOR> ()
{ }

template<typename T, typename DESCRIPTOR>
T FreeEnergyWallDynamics<T,DESCRIPTOR>::computeEquilibrium(
  int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  T equilibrium = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
  if (iPop == 0) {
    equilibrium += (1.-descriptors::t<T,DESCRIPTOR>(0)) * rho;
  } else {
    equilibrium -= descriptors::t<T,DESCRIPTOR>(iPop) * rho;
  }
  return equilibrium;
}

template<typename T, typename DESCRIPTOR>
T FreeEnergyWallDynamics<T,DESCRIPTOR>::getOmega() const
{
  return 0.;
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyWallDynamics<T,DESCRIPTOR>::setOmega(T omega_)
{ }


/////////////// FreeEnergyInletOutletDynamics ////////////////////////////////////////

template<typename T, typename DESCRIPTOR, int direction, int orientation>
FreeEnergyInletOutletDynamics<T,DESCRIPTOR,direction,orientation>::FreeEnergyInletOutletDynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta )
  : BasicDynamics<T,DESCRIPTOR>(momenta), omega(omega_)
{ }

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void FreeEnergyInletOutletDynamics<T,DESCRIPTOR,direction,orientation>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  typedef DESCRIPTOR L;

  // Do a standard collision neglecting the chemical potential term.
  T rho, u[DESCRIPTOR::d];
  computeRhoU(cell, rho, u);
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  statistics.incrementStats(rho, uSqr);
  for(int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] -= omega * descriptors::t<T,DESCRIPTOR>(iPop) * rho;
  }
  cell[0] += omega * (1.-descriptors::t<T,DESCRIPTOR>(0)) * rho;

  // Distribute the missing density to the unknown distribution functions.
  std::vector<int> missingIndices = util::subIndexOutgoing<L,direction,orientation>();
  T missingRho = rho - 1.;
  T missingWeightSum = 0;
  for (int iPop = 0; iPop < L::q; ++iPop) {
    if ( std::find(missingIndices.begin(), missingIndices.end(), iPop) != missingIndices.end() ) {
      missingWeightSum += descriptors::t<T,L>(iPop);
    } else {
      missingRho -= cell[iPop];
    }
  }

  for (unsigned iPop = 0; iPop < missingIndices.size(); ++iPop) {
    cell[missingIndices[iPop]] = missingRho * descriptors::t<T,L>(missingIndices[iPop]) / missingWeightSum;
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
T FreeEnergyInletOutletDynamics<T,DESCRIPTOR,direction,orientation>::computeEquilibrium(
  int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  T equilibrium = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
  if (iPop == 0) {
    equilibrium += (1.-descriptors::t<T,DESCRIPTOR>(0)) * rho;
  } else {
    equilibrium -= descriptors::t<T,DESCRIPTOR>(iPop) * rho;
  }
  return equilibrium;
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
T FreeEnergyInletOutletDynamics<T,DESCRIPTOR,direction,orientation>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void FreeEnergyInletOutletDynamics<T,DESCRIPTOR,direction,orientation>::setOmega(T omega_)
{
  omega = omega_;
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
T FreeEnergyInletOutletDynamics<T,DESCRIPTOR,direction,orientation>::computeRho(
  Cell<T,DESCRIPTOR> const& cell) const
{
  return cell.template getField<descriptors::FORCE>()[0];
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void FreeEnergyInletOutletDynamics<T,DESCRIPTOR,direction,orientation>::computeU (
  Cell<T,DESCRIPTOR> const& cell, T u[DESCRIPTOR::d]) const
{
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] = 0.;
  }
  u[direction] = cell.template getFieldPointer<descriptors::FORCE>()[1];
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void FreeEnergyInletOutletDynamics<T,DESCRIPTOR,direction,orientation>::computeRhoU (
  Cell<T,DESCRIPTOR> const& cell, T& rho, T u[DESCRIPTOR::d]) const
{
  rho = computeRho(cell);
  computeU(cell, u);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void FreeEnergyInletOutletDynamics<T,DESCRIPTOR,direction,orientation>::defineRho(
  Cell<T,DESCRIPTOR>& cell, T rho)
{
  cell.template getFieldPointer<descriptors::FORCE>()[0] = rho;
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void FreeEnergyInletOutletDynamics<T,DESCRIPTOR,direction,orientation>::defineU(
  Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d])
{
  cell.template getFieldPointer<descriptors::FORCE>()[1] = u[direction];
}

}

#endif
