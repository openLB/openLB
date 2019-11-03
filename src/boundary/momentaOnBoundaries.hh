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
 * Implementation of boundary cell dynamics -- generic implementation.
 */
#ifndef MOMENTA_ON_BOUNDARIES_HH
#define MOMENTA_ON_BOUNDARIES_HH

#include "momentaOnBoundaries.h"
#include "core/util.h"
#include "dynamics/lbHelpers.h"

namespace olb {

////////////////////// Class EquilibriumBM //////////////////////

template<typename T, typename DESCRIPTOR>
EquilibriumBM<T,DESCRIPTOR>::EquilibriumBM()
{
  _rho = (T)1;
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = T();
  }
}

template<typename T, typename DESCRIPTOR>
EquilibriumBM<T,DESCRIPTOR>::EquilibriumBM(T rho, const T u[DESCRIPTOR::d])
{
  _rho = rho;
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
T EquilibriumBM<T,DESCRIPTOR>::computeRho( Cell<T,DESCRIPTOR> const& cell ) const
{
  return _rho;
}

template<typename T, typename DESCRIPTOR>
void EquilibriumBM<T,DESCRIPTOR>::computeU(
  Cell<T,DESCRIPTOR> const& cell, T u[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = _u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void EquilibriumBM<T,DESCRIPTOR>::computeJ (
  Cell<T,DESCRIPTOR> const& cell, T j[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    j[iD] = _u[iD]*_rho;
  }
}

template<typename T, typename DESCRIPTOR>
void EquilibriumBM<T,DESCRIPTOR>::computeStress( Cell<T,DESCRIPTOR> const& cell,
    T rho, const T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  for (int iPi=0; iPi<util::TensorVal<DESCRIPTOR >::n; ++iPi) {
    pi[iPi] = T();
  }
}


template<typename T, typename DESCRIPTOR>
void EquilibriumBM<T,DESCRIPTOR>::defineRho( Cell<T,DESCRIPTOR>& cell, T rho )
{
  _rho = rho;
}

template<typename T, typename DESCRIPTOR>
void EquilibriumBM<T,DESCRIPTOR>::defineU(
  Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d])
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void EquilibriumBM<T,DESCRIPTOR>::defineAllMomenta( Cell<T,DESCRIPTOR>& cell,
    T rho, const T u[DESCRIPTOR::d], const T pi[util::TensorVal<DESCRIPTOR >::n] )
{
  defineRho(cell, rho);
  defineU(cell, u);
}



////////////////////// Class VelocityBM //////////////////////

template<typename T, typename DESCRIPTOR, int direction, int orientation>
VelocityBM<T,DESCRIPTOR,direction,orientation>::VelocityBM()
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = T();
  }
}

/** It takes as argument the value of the velocity to be imposed on the
 * boundary.
 */
template<typename T, typename DESCRIPTOR, int direction, int orientation>
VelocityBM<T,DESCRIPTOR,direction,orientation>::VelocityBM(const T u[DESCRIPTOR::d])
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
T VelocityBM<T,DESCRIPTOR,direction,orientation>::computeRho( Cell<T,DESCRIPTOR> const& cell ) const
{
  return velocityBMRho<T, DESCRIPTOR, direction, orientation>(cell, _u);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void VelocityBM<T,DESCRIPTOR,direction,orientation>::computeU (
  Cell<T,DESCRIPTOR> const& cell, T u[DESCRIPTOR::d]) const
{
  computeU(u);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void VelocityBM<T,DESCRIPTOR,direction,orientation>::computeJ (
  Cell<T,DESCRIPTOR> const& cell, T j[DESCRIPTOR::d]) const
{
  T rho = computeRho(cell);
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    j[iD] = _u[iD]*rho;
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void VelocityBM<T,DESCRIPTOR,direction,orientation>::computeU( T u[DESCRIPTOR::d] ) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = _u[iD];
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void VelocityBM<T,DESCRIPTOR,direction,orientation>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho )
{
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void VelocityBM<T,DESCRIPTOR,direction,orientation>::defineU (
  Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d])
{
  defineU(u);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void VelocityBM<T,DESCRIPTOR,direction,orientation>::defineU(const T u[DESCRIPTOR::d])
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void VelocityBM<T,DESCRIPTOR,direction,orientation>::defineAllMomenta (
  Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], const T pi[util::TensorVal<DESCRIPTOR >::n] )
{
  this->defineRhoU(cell, rho, u);
}


////////////////////// Class PressureBM //////////////////////

template<typename T, typename DESCRIPTOR, int direction, int orientation>
PressureBM<T,DESCRIPTOR,direction,orientation>::PressureBM()
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _values[iD] = T();
  }
  _values[direction] = 1.;
}

/** It takes as argument the value of the tangential velocity
 * components, and the pressure, to be imposed on the boundary.
 */
template<typename T, typename DESCRIPTOR, int direction, int orientation>
PressureBM<T,DESCRIPTOR,direction,orientation>::PressureBM(const T values[DESCRIPTOR::d])
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _values[iD] = values[iD];
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
T PressureBM<T,DESCRIPTOR,direction,orientation>::computeRho( Cell<T,DESCRIPTOR> const& cell ) const
{
  return computeRho();
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
T PressureBM<T,DESCRIPTOR,direction,orientation>::computeRho() const
{
  return _values[direction];
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void PressureBM<T,DESCRIPTOR,direction,orientation>::computeU (
  Cell<T,DESCRIPTOR> const& cell, T u[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = _values[iD];
  }
  T rho = _values[direction];

  std::vector<int> const& onWallIndices
    = util::subIndex<DESCRIPTOR, direction, 0>();

  std::vector<int> const& normalIndices
    = util::subIndex<DESCRIPTOR, direction, orientation>();

  T rhoOnWall = T();
  for (auto & e : onWallIndices) {
    rhoOnWall += cell[e];
  }

  T rhoNormal = T();
  for (auto & e : normalIndices) {
    rhoNormal += cell[e];
  }

  u[direction] = (T)orientation*( ((T)2*rhoNormal+rhoOnWall+(T)1 ) / rho-(T)1 );
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void PressureBM<T,DESCRIPTOR,direction,orientation>::computeJ (
  Cell<T,DESCRIPTOR> const& cell, T j[DESCRIPTOR::d]) const
{
  computeU(cell, j);
  T rho = computeRho(cell);
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    j[iD] *= rho;
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void PressureBM<T,DESCRIPTOR,direction,orientation>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{
  defineRho(rho);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void PressureBM<T,DESCRIPTOR,direction,orientation>::defineRho(T rho )
{
  _values[direction] = rho;
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void PressureBM<T,DESCRIPTOR,direction,orientation>::defineU(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d])
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    if (iD != direction) {
      _values[iD] = u[iD];
    }
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void PressureBM<T,DESCRIPTOR,direction,orientation>::defineAllMomenta(
  Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], const T pi[util::TensorVal<DESCRIPTOR >::n] )
{
  this->defineRhoU(cell, rho, u);
}


////////  FreeStressBM //////////////////////////////////////////////

template<typename T, typename DESCRIPTOR>
void FreeStressBM<T,DESCRIPTOR>::computeStress (
  Cell<T,DESCRIPTOR> const& cell, T rho, const T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  lbHelpers<T,DESCRIPTOR>::computeStress(cell, rho, u, pi);
}



////////////////////// Class RegularizedBM //////////////////////

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void RegularizedBM<T,DESCRIPTOR,direction,orientation>::computeStress (
  Cell<T,DESCRIPTOR> const& cell, T rho, const T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  BoundaryHelpers<T,DESCRIPTOR,direction,orientation>::computeStress(cell, rho, u, pi);
}

////////////////////// Class FixedVelocityBM //////////////////////////

template<typename T, typename DESCRIPTOR>
T FixedVelocityBM<T,DESCRIPTOR>::computeRho(Cell<T,DESCRIPTOR> const& cell) const
{
  return _basicMomenta.computeRho(cell);
}

template<typename T, typename DESCRIPTOR>
void FixedVelocityBM<T,DESCRIPTOR>::computeU(Cell<T,DESCRIPTOR> const& cell, T u[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = _fixU[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void FixedVelocityBM<T,DESCRIPTOR>::computeJ(Cell<T,DESCRIPTOR> const& cell, T j[DESCRIPTOR::d]) const
{
  T rho = computeRho(cell);
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    j[iD] = _fixU[iD]*rho;
  }
}

template<typename T, typename DESCRIPTOR>
void FixedVelocityBM<T,DESCRIPTOR>::computeStress (
  Cell<T,DESCRIPTOR> const& cell, T rho, const T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  _basicMomenta.computeStress(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void FixedVelocityBM<T,DESCRIPTOR>::computeRhoU (
  Cell<T,DESCRIPTOR> const& cell, T& rho, T u[DESCRIPTOR::d] ) const
{
  rho = computeRho(cell);
  computeU(cell,u);
}

template<typename T, typename DESCRIPTOR>
void FixedVelocityBM<T,DESCRIPTOR>::computeAllMomenta (
  Cell<T,DESCRIPTOR> const& cell, T& rho, T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  _basicMomenta.computeAllMomenta(cell, rho, u, pi);
  computeU(cell, u);
}

template<typename T, typename DESCRIPTOR>
void FixedVelocityBM<T,DESCRIPTOR>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{
  _basicMomenta.defineRho(cell, rho);
}

template<typename T, typename DESCRIPTOR>
void FixedVelocityBM<T,DESCRIPTOR>::defineU(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d])
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _fixU[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void FixedVelocityBM<T,DESCRIPTOR>::defineRhoU(Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d])
{
  defineRho(cell,rho);
  defineU(cell,u);
}

template<typename T, typename DESCRIPTOR>
void FixedVelocityBM<T,DESCRIPTOR>::defineAllMomenta( Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d],
    const T pi[util::TensorVal<DESCRIPTOR >::n] )
{
  _basicMomenta.defineAllMomenta(cell, rho, u, pi);
  defineU(cell,u);
}

////////////////////// Extracted helper functions //////////////////////////

template<typename T, typename DESCRIPTOR, int direction, int orientation>
T velocityBMRho( Cell<T,DESCRIPTOR> const& cell, const T* u )
{
  std::vector<int> const& onWallIndices
    = util::subIndex<DESCRIPTOR, direction, 0>();

  std::vector<int> const& normalIndices
    = util::subIndex<DESCRIPTOR, direction, orientation>();

  T rhoOnWall = T();
  for (auto & e : onWallIndices) {
    rhoOnWall += cell[e];
  }

  T rhoNormal = T();
  for (auto & e : normalIndices) {
    rhoNormal += cell[e];
  }

  T rho =((T)2*rhoNormal+rhoOnWall+(T)1) /
         ((T)1+(T)orientation*u[direction]);

  return rho;
}



}  // namespace olb

#endif
