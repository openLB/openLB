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
 * Local boundary cell 3D dynamics -- generic implementation.
 */
#ifndef MOMENTA_ON_BOUNDARIES_3D_HH
#define MOMENTA_ON_BOUNDARIES_3D_HH

#include "momentaOnBoundaries3D.h"
#include "dynamics/lbHelpers.h"
#include "dynamics/firstOrderLbHelpers.h"

namespace olb {

////////////////////// Class InnerEdgeVelBM3D ///////////////

template<typename T, typename DESCRIPTOR,
         int plane, int normal1, int normal2>
InnerEdgeVelBM3D<T,DESCRIPTOR,plane,normal1,normal2>::
InnerEdgeVelBM3D()
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = T();
  }
}

template<typename T, typename DESCRIPTOR,
         int plane, int normal1, int normal2>
InnerEdgeVelBM3D<T,DESCRIPTOR,plane,normal1,normal2>::
InnerEdgeVelBM3D(const T u_[DESCRIPTOR::d])
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u_[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int plane, int normal1, int normal2>
T InnerEdgeVelBM3D<T,DESCRIPTOR,plane,normal1,normal2>::computeRho (
  Cell<T,DESCRIPTOR> const& cell ) const
{
  T rho1 = velocityBMRho<T, DESCRIPTOR, direction1, normal1>(cell, _u);
  T rho2 = velocityBMRho<T, DESCRIPTOR, direction2, normal2>(cell, _u);
  return (rho1 + rho2) / (T)2;
}

template<typename T, typename DESCRIPTOR,
         int plane, int normal1, int normal2>
void InnerEdgeVelBM3D<T,DESCRIPTOR,plane,normal1,normal2>::computeU (
  Cell<T,DESCRIPTOR> const& cell,
  T u[DESCRIPTOR::d] ) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = _u[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int plane, int normal1, int normal2>
void InnerEdgeVelBM3D<T,DESCRIPTOR,plane,normal1,normal2>::computeJ (
  Cell<T,DESCRIPTOR> const& cell,
  T j[DESCRIPTOR::d] ) const
{
  T rho = computeRho(cell);
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    j[iD] = _u[iD]*rho;
  }
}

template<typename T, typename DESCRIPTOR,
         int plane, int normal1, int normal2>
void InnerEdgeVelBM3D<T,DESCRIPTOR,plane,normal1,normal2>::computeU (
  T u[DESCRIPTOR::d] ) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = _u[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int plane, int normal1, int normal2>
void InnerEdgeVelBM3D<T,DESCRIPTOR,plane,normal1,normal2>::defineRho (
  Cell<T,DESCRIPTOR>& cell, T rho )
{ }

template<typename T, typename DESCRIPTOR,
         int plane, int normal1, int normal2>
void InnerEdgeVelBM3D<T,DESCRIPTOR,plane,normal1,normal2>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d] )
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int plane, int normal1, int normal2>
void InnerEdgeVelBM3D<T,DESCRIPTOR,plane,normal1,normal2>::defineU (
  const T u[DESCRIPTOR::d] )
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int plane, int normal1, int normal2>
void InnerEdgeVelBM3D<T,DESCRIPTOR,plane,normal1,normal2>::
defineAllMomenta (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  const T pi[util::TensorVal<DESCRIPTOR >::n] )
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int plane, int normal1, int normal2>
void InnerEdgeVelBM3D<T,DESCRIPTOR,plane,normal1,normal2>::
computeStress (
  Cell<T,DESCRIPTOR> const& cell,
  T rho, const T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  typedef lbHelpers<T,DESCRIPTOR> lbH;

  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);

  Cell<T,DESCRIPTOR> newCell(cell);
  for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
    if ( (descriptors::c<DESCRIPTOR>(iPop,direction1) == -normal1) &&
         (descriptors::c<DESCRIPTOR>(iPop,direction2) == -normal2) ) {
      int opp = util::opposite<DESCRIPTOR >(iPop);
      newCell[iPop] = newCell[opp]
                      - lbH::equilibrium(opp, rho, u, uSqr)
                      + lbH::equilibrium(iPop, rho, u, uSqr);
    }
  }
  lbH::computeStress(newCell, rho, u, pi);
}


////////////////////// Class InnerCornerVelBM3D ///////////////

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY, int normalZ>
InnerCornerVelBM3D<T,DESCRIPTOR,normalX,normalY,normalZ>::InnerCornerVelBM3D()
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = T();
  }
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY, int normalZ>
InnerCornerVelBM3D<T,DESCRIPTOR,normalX,normalY,normalZ>::InnerCornerVelBM3D (
  const T u_[DESCRIPTOR::d])
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u_[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY, int normalZ>
T InnerCornerVelBM3D<T,DESCRIPTOR,normalX,normalY,normalZ>::computeRho (
  Cell<T,DESCRIPTOR> const& cell ) const
{
  T rhoX = velocityBMRho<T, DESCRIPTOR, 0, normalX>(cell, _u);
  T rhoY = velocityBMRho<T, DESCRIPTOR, 1, normalY>(cell, _u);
  T rhoZ = velocityBMRho<T, DESCRIPTOR, 2, normalZ>(cell, _u);
  return (rhoX + rhoY + rhoZ) / (T)3;
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY, int normalZ>
void InnerCornerVelBM3D<T,DESCRIPTOR,normalX,normalY,normalZ>::computeU (
  Cell<T,DESCRIPTOR> const& cell,
  T u[DESCRIPTOR::d] ) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = _u[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY, int normalZ>
void InnerCornerVelBM3D<T,DESCRIPTOR,normalX,normalY,normalZ>::computeJ (
  Cell<T,DESCRIPTOR> const& cell,
  T j[DESCRIPTOR::d] ) const
{
  T rho = computeRho(cell);
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    j[iD] = _u[iD]*rho;
  }
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY, int normalZ>
void InnerCornerVelBM3D<T,DESCRIPTOR,normalX,normalY,normalZ>::computeU (
  T u[DESCRIPTOR::d] ) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = _u[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY, int normalZ>
void InnerCornerVelBM3D<T,DESCRIPTOR,normalX,normalY,normalZ>::defineRho (
  Cell<T,DESCRIPTOR>& cell, T rho )
{ }

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY, int normalZ>
void InnerCornerVelBM3D<T,DESCRIPTOR,normalX,normalY,normalZ>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d] )
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY, int normalZ>
void InnerCornerVelBM3D<T,DESCRIPTOR,normalX,normalY,normalZ>::defineU (
  const T u[DESCRIPTOR::d] )
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY, int normalZ>
void InnerCornerVelBM3D<T,DESCRIPTOR,normalX,normalY,normalZ>::defineAllMomenta (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  const T pi[util::TensorVal<DESCRIPTOR >::n] )
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR,
         int normalX, int normalY, int normalZ>
void InnerCornerVelBM3D<T,DESCRIPTOR,normalX,normalY,normalZ>::computeStress (
  Cell<T,DESCRIPTOR> const& cell,
  T rho, const T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  typedef lbHelpers<T,DESCRIPTOR> lbH;
  Cell<T,DESCRIPTOR> newCell(cell);
  int v[DESCRIPTOR::d] = { -normalX, -normalY, -normalZ };
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
