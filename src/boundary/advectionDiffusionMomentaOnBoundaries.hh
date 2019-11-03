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
#ifndef ADE_MOMENTA_ON_BOUNDARIES_HH
#define ADE_MOMENTA_ON_BOUNDARIES_HH

#include <limits>
#include "advectionDiffusionMomentaOnBoundaries.h"
#include "dynamics/lbHelpers.h"
#include "dynamics/firstOrderLbHelpers.h"

namespace olb {

////////////////////// Class AdvectionDiffusionBM //////////////////////

template<typename T, typename DESCRIPTOR, int direction, int orientation>
RegularizedTemperatureBM<T,DESCRIPTOR,direction,orientation>::RegularizedTemperatureBM(const T temperature) :
  AdvectionDiffusionBulkMomenta<T,DESCRIPTOR>(),
  _temperature(temperature)
{
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
T RegularizedTemperatureBM<T,DESCRIPTOR,direction,orientation>::computeRho( Cell<T,DESCRIPTOR> const& cell ) const
{
  return _temperature;
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void RegularizedTemperatureBM<T,DESCRIPTOR,direction,orientation>::computeJ( Cell<T,DESCRIPTOR> const& cell, T j[DESCRIPTOR::d] ) const
{
  const T* u = cell.template getFieldPointer<descriptors::VELOCITY>();
  computeJneq( cell, j );

  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    j[iD] += _temperature * u[iD];
  }
}


template<typename T, typename DESCRIPTOR, int direction, int orientation>
void RegularizedTemperatureBM<T,DESCRIPTOR,direction,orientation>::computeJneq( Cell<T,DESCRIPTOR> const& cell, T jNeq[DESCRIPTOR::d] ) const
{
  std::vector<int> const& onWallIndices = util::subIndex<DESCRIPTOR, direction, 0>();
  std::vector<int> const& normalIndices = util::subIndex<DESCRIPTOR, direction, orientation>();

  const T* u = cell.template getFieldPointer<descriptors::VELOCITY>();

  T jNeqOnWall[DESCRIPTOR::d], jNeqNormal[DESCRIPTOR::d];
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    jNeqOnWall[iD] = T();
    jNeqNormal[iD] = T();
  }

  for (unsigned fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      jNeqOnWall[iD] += (cell[onWallIndices[fIndex]] - lbHelpers<T,DESCRIPTOR>::equilibriumFirstOrder(onWallIndices[fIndex],_temperature,u)) * descriptors::c<DESCRIPTOR>(onWallIndices[fIndex],iD);
    }
  }

  for (unsigned fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      jNeqNormal[iD] += (cell[normalIndices[fIndex]]-lbHelpers<T,DESCRIPTOR>::equilibriumFirstOrder(normalIndices[fIndex],_temperature,u)) * descriptors::c<DESCRIPTOR>(normalIndices[fIndex],iD);
    }
  }

  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    jNeq[iD] = jNeqOnWall[iD] + (T)2 * jNeqNormal[iD];
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void RegularizedTemperatureBM<T,DESCRIPTOR,direction,orientation>::defineRho( Cell<T,DESCRIPTOR>& cell, T rho )
{
  _temperature = rho;
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void RegularizedTemperatureBM<T,DESCRIPTOR,direction,orientation>::defineRhoU( Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d])
{
  _temperature = rho;
  T *u_ = cell.template getFieldPointer<descriptors::VELOCITY>();
  for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
    u_[iD] = u[iD];
  }
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = lbHelpers<T, DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
  }
}




template<typename T, typename DESCRIPTOR, int direction, int orientation>
RegularizedHeatFluxBM<T,DESCRIPTOR,direction,orientation>::RegularizedHeatFluxBM(T *heatFlux) :
  AdvectionDiffusionBulkMomenta<T,DESCRIPTOR>()
{
  for (int iDim = 0; iDim < DESCRIPTOR::d; iDim++) {
    if (heatFlux == nullptr)
      _heatFlux[iDim] = T();
    else
      _heatFlux[iDim] = heatFlux[iDim];
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
T RegularizedHeatFluxBM<T,DESCRIPTOR,direction,orientation>::computeRho( Cell<T,DESCRIPTOR> const& cell ) const
{
  std::vector<int> const& onWallIndices = util::subIndex<DESCRIPTOR, direction, 0>();
  std::vector<int> const& normalIndices = util::subIndex<DESCRIPTOR, direction, orientation>();
  const T* u = cell.template getFieldPointer<descriptors::VELOCITY>();

  T rhoOnWall = T();
  for (unsigned fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
    rhoOnWall += cell[onWallIndices[fIndex]];
  }

  T rhoNormal = T();
  for (unsigned fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
    rhoNormal += cell[normalIndices[fIndex]];
  }

  T rho =((T)2*rhoNormal+rhoOnWall+(T)1 - (T)orientation*_heatFlux[direction]) /
         ((T)1+(T)orientation*u[direction]);

  return rho;
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void RegularizedHeatFluxBM<T,DESCRIPTOR,direction,orientation>::computeJ( Cell<T,DESCRIPTOR> const& cell, T j[DESCRIPTOR::d] ) const
{
  T temperature = computeRho(cell);
  const T* u = cell.template getFieldPointer<descriptors::VELOCITY>();

  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    j[iD] = _heatFlux[iD] + temperature * u[iD];
  }
}


template<typename T, typename DESCRIPTOR, int direction, int orientation>
void RegularizedHeatFluxBM<T,DESCRIPTOR,direction,orientation>::computeJneq( Cell<T,DESCRIPTOR> const& cell, T jNeq[DESCRIPTOR::d] ) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    jNeq[iD] = _heatFlux[iD];
  }
}


}  // namespace olb

#endif
