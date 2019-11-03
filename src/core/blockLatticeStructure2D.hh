/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2018 Jonas Latt, Adrian Kummerlaender
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
 * Dynamics for a generic 2D block structure -- header file.
 */
#ifndef BLOCK_LATTICE_STRUCTURE_2D_HH
#define BLOCK_LATTICE_STRUCTURE_2D_HH

#include "blockLatticeStructure2D.h"
#include "functors/lattice/indicator/blockIndicatorBaseF2D.hh"
#include "functors/lattice/indicator/blockIndicatorF2D.hh"

namespace olb {


template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineRho(
  BlockIndicatorF2D<T>& indicator, AnalyticalF2D<T,T>& rho)
{
  T physR[2] = { };
  T rhoTmp = T();
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      if (indicator(iX, iY)) {
        indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY);
        rho(&rhoTmp, physR);
        get(iX, iY).defineRho(rhoTmp);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineRho(
  BlockGeometryStructure2D<T>& blockGeometry, int material, AnalyticalF2D<T,T>& rho)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometry, material);
  defineRho(indicator, rho);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineU(
  BlockIndicatorF2D<T>& indicator, AnalyticalF2D<T,T>& u)
{
  T physR[2] = { };
  T uTmp[2] = { };
  const auto& geometry = indicator.getBlockGeometryStructure();
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      if (indicator(iX, iY)) {
        geometry.getPhysR(physR, iX, iY);
        u(uTmp, physR);
        get(iX, iY).defineU(uTmp);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineU(
  BlockGeometryStructure2D<T>& blockGeometry, int material, AnalyticalF2D<T,T>& u)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometry, material);
  defineU(indicator, u);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineRhoU(
  BlockIndicatorF2D<T>& indicator,
  AnalyticalF2D<T,T>& rho, AnalyticalF2D<T,T>& u)
{
  T physR[2] = { };
  T uTmp[2] = { };
  T rhoTmp = T();
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      if (indicator(iX, iY)) {
        indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY);
        rho(&rhoTmp, physR);
        u(uTmp, physR);
        get(iX, iY).defineRhoU(rhoTmp, uTmp);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineRhoU(
  BlockGeometryStructure2D<T>& blockGeometry, int material,
  AnalyticalF2D<T,T>& rho, AnalyticalF2D<T,T>& u)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometry, material);
  defineRhoU(indicator, rho, u);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::definePopulations(
  BlockIndicatorF2D<T>& indicator, AnalyticalF2D<T,T>& Pop)
{
  T physR[2] = { };
  T PopTmp[DESCRIPTOR::q];
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      if (indicator(iX, iY)) {
        indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY);
        Pop(PopTmp, physR);
        get(iX, iY).definePopulations(PopTmp);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::definePopulations(
  BlockGeometryStructure2D<T>& blockGeometry, int material, AnalyticalF2D<T,T>& Pop)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometry, material);
  definePopulations(indicator, Pop);
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineField(
  BlockIndicatorF2D<T>& indicator, AnalyticalF2D<T,T>& field)
{
  T* fieldTmp = new T[DESCRIPTOR::template size<FIELD>()];
  T physR[2] = { };
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      if (indicator(iX, iY)) {
        indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY);
        field(fieldTmp, physR);
        get(iX, iY).template defineField<FIELD>(fieldTmp);
      }
    }
  }
  delete[] fieldTmp;
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineField(
  BlockGeometryStructure2D<T>& blockGeometry, int material, AnalyticalF2D<T,T>& field)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometry, material);
  defineField<FIELD>(indicator, field);
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLatticeStructure2D<T,DESCRIPTOR>::defineField(
  BlockGeometryStructure2D<T>& blockGeometry, IndicatorF2D<T>& indicatorF,
  AnalyticalF2D<T,T>& field)
{
  BlockIndicatorFfromIndicatorF2D<T> indicator(indicatorF, blockGeometry);
  defineField<FIELD>(indicator, field);
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLatticeStructure2D<T,DESCRIPTOR>::addField(
  BlockGeometryStructure2D<T>& blockGeometry, IndicatorF2D<T>& indicator,
  AnalyticalF2D<T,T>& field)
{
  T* fieldTmp = new T[DESCRIPTOR::template size<FIELD>()];
  T physR[2] = { };
  bool inside;
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      blockGeometry.getPhysR(physR, iX, iY);
      indicator(&inside, &(physR[0]));
      if (inside) {
        field(fieldTmp,physR);
        get(iX, iY).template addField<FIELD>(fieldTmp);
      }
    }
  }
  delete[] fieldTmp;
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLatticeStructure2D<T,DESCRIPTOR>::addField(
  BlockGeometryStructure2D<T>& blockGeometry, IndicatorF2D<T>& indicator,
  AnalyticalF2D<T,T>& field, AnalyticalF2D<T,T>& porous)
{
  T* fieldTmp = new T[DESCRIPTOR::template size<FIELD>()];
  bool inside;
  T physR[2] = { };
  T porousA[1] = { };
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      blockGeometry.getPhysR(physR, iX, iY);
      indicator(&inside, physR);
      if (inside) {
        porous(porousA, physR);
        field(fieldTmp,physR);
        for (int i = 0; i < DESCRIPTOR::template size<FIELD>(); ++i) {
          fieldTmp[i] *= porousA[0];
        }
        get(iX, iY).template addField<FIELD>(fieldTmp);
      }
    }
  }
  delete[] fieldTmp;
}


template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLatticeStructure2D<T,DESCRIPTOR>::multiplyField(
  BlockGeometryStructure2D<T>& blockGeometry, IndicatorF2D<T>& indicator,
  AnalyticalF2D<T,T>& field)
{
  T* fieldTmp = new T [DESCRIPTOR::template size<FIELD>()];
  bool inside;
  T physR[3] = { };
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      blockGeometry.getPhysR(physR, iX, iY);
      indicator(&inside, &(physR[0]));
      if (inside) {
        field(fieldTmp, physR);
        get(iX, iY).template multiplyField<FIELD>(fieldTmp);
      }
    }
  }
  delete[] fieldTmp;
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::iniEquilibrium(
  BlockIndicatorF2D<T>& indicator,
  AnalyticalF2D<T,T>& rho, AnalyticalF2D<T,T>& u)
{
  T physR[2] = { };
  T uTmp[2] = { };
  T rhoTmp = T();
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      if (indicator(iX, iY)) {
        indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY);
        u(uTmp, physR);
        rho(&rhoTmp, physR);
        get(iX, iY).iniEquilibrium(rhoTmp, uTmp);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure2D<T,DESCRIPTOR>::iniEquilibrium(
  BlockGeometryStructure2D<T>& blockGeometry, int material,
  AnalyticalF2D<T,T>& rho, AnalyticalF2D<T,T>& u)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometry, material);
  iniEquilibrium(indicator, rho, u);
}

////////// FREE FUNCTIONS //////////

template<typename T, typename DESCRIPTOR>
void setBlockExternalParticleField(BlockGeometryStructure2D<T>& blockGeometry, AnalyticalF2D<T,T>& velocity, SmoothIndicatorF2D<T,T,true>& sIndicator, BlockLattice2D<T,DESCRIPTOR>& extendedBlockLattice)
{

  int start[2] = {0};
  int end[2] = {0};
  // check for intersection of cuboid and indicator
  Cuboid2D<T> tmpCuboid(blockGeometry.getOrigin()[0], blockGeometry.getOrigin()[1], blockGeometry.getDeltaR(), blockGeometry.getNx(), blockGeometry.getNy());
  T posXmin = sIndicator.getPos()[0] - sIndicator.getCircumRadius();
  T posXmax = sIndicator.getPos()[0] + sIndicator.getCircumRadius();
  T posYmin = sIndicator.getPos()[1] - sIndicator.getCircumRadius();
  T posYmax = sIndicator.getPos()[1] + sIndicator.getCircumRadius();
  if(tmpCuboid.checkInters(posXmin, posXmax, posYmin, posYmax, start[0], end[0], start[1], end[1])) {

    for (int k=0; k<2; k++) {
      start[k] -= 1;
      if(start[k] < 0) start[k] = 0;
      end[k] += 2;
      if(end[k] > blockGeometry.getExtend()[k]) end[k] = blockGeometry.getExtend()[k];
    }

    T foo[3] = { }; /// Contains foo[0]=vel0; foo[1]=vel1; foo[2]=porosity
    T physR[2]= { };
    T porosity[1] = { };
    for (int iX = start[0]; iX < end[0]; ++iX) {
      for (int iY = start[1]; iY < end[1]; ++iY) {
        blockGeometry.getPhysR(physR, iX, iY);
        sIndicator(porosity, physR);
        if (!util::nearZero(porosity[0])) {
          // TODO: Check / adapt to use descriptor fields
          velocity(foo,physR);
          foo[0] *= porosity[0];
          foo[1] *= porosity[0];
          foo[2] = porosity[0];
          extendedBlockLattice.get(iX, iY).template addField<descriptors::VELOCITY_NUMERATOR>(foo);
          extendedBlockLattice.get(iX, iY).template addField<descriptors::VELOCITY_DENOMINATOR>(&foo[2]);
          porosity[0] = 1.-porosity[0];
          *(extendedBlockLattice.get(iX, iY).template getFieldPointer<descriptors::POROSITY>()) *= porosity[0];
        }
      }
    }
  }
}


}  // namespace olb

#endif
