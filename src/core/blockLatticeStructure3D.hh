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
 * Dynamics for a generic 3D block lattice view -- generic implementation.
 */
#ifndef BLOCK_LATTICE_STRUCTURE_3D_HH
#define BLOCK_LATTICE_STRUCTURE_3D_HH

#include "blockLatticeStructure3D.h"
#include "functors/lattice/indicator/blockIndicatorBaseF3D.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure3D<T,DESCRIPTOR>::defineRho(
  BlockIndicatorF3D<T>& indicator, AnalyticalF3D<T,T>& rho)
{
  T physR[3] = { };
  T rhoTmp = T();
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      for (int iZ = 0; iZ < getNz(); ++iZ) {
        if (indicator(iX, iY, iZ)) {
          indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY, iZ);
          rho(&rhoTmp,physR);
          get(iX, iY, iZ).defineRho(rhoTmp);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure3D<T,DESCRIPTOR>::defineRho(
  BlockGeometryStructure3D<T>& blockGeometry, int material, AnalyticalF3D<T,T>& rho)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometry, material);
  defineRho(indicator, rho);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure3D<T,DESCRIPTOR>::defineU(
  BlockIndicatorF3D<T>& indicator, AnalyticalF3D<T,T>& u)
{
  T physR[3] = { };
  T uTmp[3] = { };
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      for (int iZ = 0; iZ < getNz(); ++iZ) {
        if (indicator(iX, iY, iZ)) {
          indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY, iZ);
          u(uTmp,physR);
          get(iX, iY, iZ).defineU(uTmp);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure3D<T,DESCRIPTOR>::defineU(
  BlockGeometryStructure3D<T>& blockGeometry, int material, AnalyticalF3D<T,T>& u)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometry, material);
  defineU(indicator, u);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure3D<T,DESCRIPTOR>::defineRhoU(
  BlockIndicatorF3D<T>& indicator,
  AnalyticalF3D<T,T>& rho, AnalyticalF3D<T,T>& u)
{
  T physR[3] = { };
  T uTmp[3] = { };
  T rhoTmp = T();
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      for (int iZ = 0; iZ < getNz(); ++iZ) {
        if (indicator(iX, iY, iZ)) {
          indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY, iZ);
          rho(&rhoTmp,physR);
          u(uTmp,physR);
          get(iX, iY, iZ).defineRhoU(rhoTmp,uTmp);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure3D<T,DESCRIPTOR>::defineRhoU(
  BlockGeometryStructure3D<T>& blockGeometry, int material,
  AnalyticalF3D<T,T>& rho, AnalyticalF3D<T,T>& u)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometry, material);
  defineRhoU(indicator, rho, u);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure3D<T,DESCRIPTOR>::definePopulations(
  BlockIndicatorF3D<T>& indicator, AnalyticalF3D<T,T>& Pop)
{
  T physR[3] = { };
  T PopTmp[DESCRIPTOR::q];
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      for (int iZ = 0; iZ < getNz(); ++iZ) {
        if (indicator(iX, iY, iZ)) {
          indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY, iZ);
          Pop(PopTmp,physR);
          get(iX, iY, iZ).definePopulations(PopTmp);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure3D<T,DESCRIPTOR>::definePopulations(
  BlockGeometryStructure3D<T>& blockGeometry, int material, AnalyticalF3D<T,T>& Pop)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometry, material);
  definePopulations(indicator, Pop);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure3D<T,DESCRIPTOR>::definePopulations(
  BlockIndicatorF3D<T>& indicator, BlockF3D<T>& Pop)
{
  int latticeR[3];
  T PopTmp[DESCRIPTOR::q];
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      for (int iZ = 0; iZ < getNz(); ++iZ) {
        if (indicator(iX, iY, iZ)) {
          latticeR[0] = iX;
          latticeR[1] = iY;
          latticeR[2] = iZ;
          Pop(PopTmp,latticeR);
          get(iX, iY, iZ).definePopulations(PopTmp);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure3D<T,DESCRIPTOR>::definePopulations(
  BlockGeometryStructure3D<T>& blockGeometry, int material, BlockF3D<T>& Pop)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometry, material);
  definePopulations(indicator, Pop);
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLatticeStructure3D<T,DESCRIPTOR>::defineField(
  BlockIndicatorF3D<T>& indicator, AnalyticalF3D<T,T>& field)
{
  T* fieldTmp = new T[DESCRIPTOR::template size<FIELD>()];
  T physR[3] = { };
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      for (int iZ = 0; iZ < getNz(); ++iZ) {
        if (indicator(iX, iY, iZ)) {
          indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY, iZ);
          field(fieldTmp, physR);
          get(iX, iY, iZ).template defineField<FIELD>(fieldTmp);
        }
      }
    }
  }
  delete[] fieldTmp;
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLatticeStructure3D<T,DESCRIPTOR>::defineField(
  BlockGeometryStructure3D<T>& blockGeometry, int material,
  AnalyticalF3D<T,T>& field)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometry, material);
  defineField<FIELD>(indicator, field);
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void BlockLatticeStructure3D<T,DESCRIPTOR>::defineField(
  BlockGeometryStructure3D<T>& blockGeometry, IndicatorF3D<T>& indicatorF,
  AnalyticalF3D<T,T>& field)
{
  BlockIndicatorFfromIndicatorF3D<T> indicator(indicatorF, blockGeometry);
  defineField<FIELD>(indicator, field);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure3D<T,DESCRIPTOR>::iniEquilibrium(
  BlockIndicatorF3D<T>& indicator,
  AnalyticalF3D<T,T>& rho, AnalyticalF3D<T,T>& u)
{
  T physR[3] = { };
  T uTmp[3] = { };
  T rhoTmp = T();
  for (int iX = 0; iX < getNx(); ++iX) {
    for (int iY = 0; iY < getNy(); ++iY) {
      for (int iZ = 0; iZ < getNz(); ++iZ) {
        if (indicator(iX, iY, iZ)) {
          indicator.getBlockGeometryStructure().getPhysR(physR, iX, iY, iZ);
          u(uTmp, physR);
          rho(&rhoTmp, physR);
          get(iX, iY, iZ).iniEquilibrium(rhoTmp, uTmp);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeStructure3D<T,DESCRIPTOR>::iniEquilibrium(
  BlockGeometryStructure3D<T>& blockGeometry, int material,
  AnalyticalF3D<T,T>& rho, AnalyticalF3D<T,T>& u)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometry, material);
  iniEquilibrium(indicator, rho, u);
}

////////// FREE FUNCTIONS //////////

template<typename T, typename DESCRIPTOR>
void setBlockExternalParticleField(BlockGeometryStructure3D<T>& blockGeometry, AnalyticalF3D<T,T>& velocity, SmoothIndicatorF3D<T,T,true>& sIndicator, BlockLattice3D<T,DESCRIPTOR>& extendedBlockLattice)
{

  int start[3] = {0};
  int end[3] = {0};
  // check for intersection of cuboid and indicator
  Cuboid3D<T> tmpCuboid(blockGeometry.getOrigin()[0], blockGeometry.getOrigin()[1], blockGeometry.getOrigin()[2], blockGeometry.getDeltaR(), blockGeometry.getNx(), blockGeometry.getNy(), blockGeometry.getNz());
  T posXmin = sIndicator.getPos()[0] - sIndicator.getCircumRadius();
  T posXmax = sIndicator.getPos()[0] + sIndicator.getCircumRadius();
  T posYmin = sIndicator.getPos()[1] - sIndicator.getCircumRadius();
  T posYmax = sIndicator.getPos()[1] + sIndicator.getCircumRadius();
  T posZmin = sIndicator.getPos()[2] - sIndicator.getCircumRadius();
  T posZmax = sIndicator.getPos()[2] + sIndicator.getCircumRadius();
  if(tmpCuboid.checkInters(posXmin, posXmax, posYmin, posYmax, posZmin, posZmax,
                           start[0], end[0], start[1], end[1], start[2], end[2]))
  {
    for (int k=0; k<3; k++) {
      start[k] -= 1;
      if(start[k] < 0) start[k] = 0;
      end[k] += 2;
      if(end[k] > blockGeometry.getExtend()[k]) end[k] = blockGeometry.getExtend()[k];
    }

    T foo[4] = { }; /// Contains foo[0]=vel0; foo[1]=vel1; foo[2]=vel2; foo[3]=porosity
    T physR[3] = { };
    T porosity[1] = { };
   
    for (int iX = start[0]; iX < end[0]; ++iX) {
      for (int iY = start[1]; iY < end[1]; ++iY) {
        for (int iZ = start[2]; iZ < end[2]; ++iZ) {
          blockGeometry.getPhysR(physR, iX, iY, iZ);
          sIndicator(porosity, physR);
          if (!util::nearZero(porosity[0])) {
            velocity(foo, physR);
            foo[0] *= porosity[0];
            foo[1] *= porosity[0];
            foo[2] *= porosity[0];
            foo[3] = porosity[0];
            extendedBlockLattice.get(iX, iY, iZ).template addField<descriptors::VELOCITY_NUMERATOR>(foo);
            extendedBlockLattice.get(iX, iY, iZ).template addField<descriptors::VELOCITY_DENOMINATOR>(&foo[3]);
            porosity[0] = 1. - porosity[0];
            *(extendedBlockLattice.get(iX, iY, iZ).template getFieldPointer<descriptors::POROSITY>()) *= porosity[0];
          }
        }
      }
    }
  }
}


}  // namespace olb

#endif
