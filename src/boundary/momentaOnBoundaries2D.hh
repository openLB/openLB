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
#ifndef MOMENTA_ON_BOUNDARIES_2D_HH
#define MOMENTA_ON_BOUNDARIES_2D_HH

#include "momentaOnBoundaries2D.h"
#include "dynamics/lbHelpers.h"
#include "dynamics/firstOrderLbHelpers.h"

namespace olb {

////////////////////// Class InnerCornerVelBM2D ///////////////

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY>
InnerCornerVelBM2D<T,DESCRIPTOR,normalX,normalY>::InnerCornerVelBM2D()
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = T();
  }
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY>
InnerCornerVelBM2D<T,DESCRIPTOR,normalX,normalY>::InnerCornerVelBM2D (
  const T u_[DESCRIPTOR::d])
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u_[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY>
T InnerCornerVelBM2D<T,DESCRIPTOR,normalX,normalY>::computeRho (
  Cell<T,DESCRIPTOR> const& cell ) const
{
  T rhoX = velocityBMRho<T, DESCRIPTOR, 0, normalX>(cell, _u);
  T rhoY = velocityBMRho<T, DESCRIPTOR, 1, normalY>(cell, _u);
  return (rhoX + rhoY) / (T)2;
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY>
void InnerCornerVelBM2D<T,DESCRIPTOR,normalX,normalY>::computeU (
  Cell<T,DESCRIPTOR> const& cell,
  T u[DESCRIPTOR::d] ) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = _u[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY>
void InnerCornerVelBM2D<T,DESCRIPTOR,normalX,normalY>::computeJ (
  Cell<T,DESCRIPTOR> const& cell,
  T j[DESCRIPTOR::d] ) const
{
  computeU(cell, j);
  T rho = computeRho(cell);
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    j[iD] *= rho;
  }
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY>
void InnerCornerVelBM2D<T,DESCRIPTOR,normalX,normalY>::computeU (
  T u[DESCRIPTOR::d] ) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = _u[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY>
void InnerCornerVelBM2D<T,DESCRIPTOR,normalX,normalY>::defineRho (
  Cell<T,DESCRIPTOR>& cell, T rho )
{ }

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY>
void InnerCornerVelBM2D<T,DESCRIPTOR,normalX,normalY>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d] )
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY>
void InnerCornerVelBM2D<T,DESCRIPTOR,normalX,normalY>::defineU (
  const T u[DESCRIPTOR::d] )
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY>
void InnerCornerVelBM2D<T,DESCRIPTOR,normalX,normalY>::defineAllMomenta (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  const T pi[util::TensorVal<DESCRIPTOR >::n] )
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY>
void InnerCornerVelBM2D<T,DESCRIPTOR,normalX,normalY>::computeStress (
  Cell<T,DESCRIPTOR> const& cell,
  T rho, const T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  typedef lbHelpers<T,DESCRIPTOR> lbH;
  Cell<T,DESCRIPTOR> newCell(cell);
  int v[DESCRIPTOR::d] = { -normalX, -normalY };
  int unknownF  = util::findVelocity<DESCRIPTOR >(v);

  if (unknownF != DESCRIPTOR::q) {
    int oppositeF = util::opposite<DESCRIPTOR >(unknownF);

    T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);

    newCell[unknownF] = newCell[oppositeF]
                        - lbH::equilibrium(oppositeF, rho, u, uSqr)
                        + lbH::equilibrium(unknownF, rho, u, uSqr);
  }

  lbH::computeStress(newCell, rho, u, pi);
}


}  // namespace olb

#endif
