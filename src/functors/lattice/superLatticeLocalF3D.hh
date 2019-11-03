/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2014 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#ifndef SUPER_LATTICE_LOCAL_F_3D_HH
#define SUPER_LATTICE_LOCAL_F_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm

#include "superBaseF3D.h"
#include "superLatticeLocalF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry3D.h"

namespace olb {



template<typename T, typename DESCRIPTOR>
SuperLatticeFpop3D<T, DESCRIPTOR>::SuperLatticeFpop3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, DESCRIPTOR::q)
{
  this->getName() = "fPop";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeFpop3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeFpop3D<T, DESCRIPTOR>::operator()( T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, typename DESCRIPTOR>
SuperLatticeDissipation3D<T, DESCRIPTOR>::SuperLatticeDissipation3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1), _converter(converter)
{
  this->getName() = "dissipation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeDissipation3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),this->_converter));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeDissipation3D<T, DESCRIPTOR>::operator()( T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, typename DESCRIPTOR>
SuperLatticePhysDissipation3D<T, DESCRIPTOR>::SuperLatticePhysDissipation3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physDissipation";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysDissipation3D<T, DESCRIPTOR>(
        this->_sLattice.getExtendedBlockLattice(iC),
        this->_sLattice.getOverlap(),
        this->_converter)
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysDissipation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    const int loc = this->_sLattice.getLoadBalancer().loc(input[0]);
    return this->getBlockF(loc)(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, typename DESCRIPTOR>
SuperLatticeEffevtiveDissipation3D<T, DESCRIPTOR>::SuperLatticeEffevtiveDissipation3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst, LESDynamics<T, DESCRIPTOR>& LESdynamics)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1), _converter(converter)
{
  this->getName() = "EffevtiveDissipation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeEffevtiveDissipation3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),
                               this->_converter, smagoConst, LESdynamics));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeEffevtiveDissipation3D<T, DESCRIPTOR>::operator()( T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, typename DESCRIPTOR>
SuperLatticePhysEffevtiveDissipation3D<T, DESCRIPTOR>::SuperLatticePhysEffevtiveDissipation3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst, LESDynamics<T, DESCRIPTOR>& LESdynamics)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physEffevtiveDissipation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysEffevtiveDissipation3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),
                               this->_converter, smagoConst, LESdynamics));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysEffevtiveDissipation3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticeDensity3D<T, DESCRIPTOR>::SuperLatticeDensity3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice) : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "density";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeDensity3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeDensity3D<T, DESCRIPTOR>::operator()( T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, typename DESCRIPTOR>
SuperLatticeVelocity3D<T, DESCRIPTOR>::SuperLatticeVelocity3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
{
  this->getName() = "velocity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeVelocity3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeVelocity3D<T, DESCRIPTOR>::operator()( T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticeExternalVelocity3D<T, DESCRIPTOR>::SuperLatticeExternalVelocity3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
{
  this->getName() = "externalVelocity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeExternalVelocity3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeExternalVelocity3D<T, DESCRIPTOR>::operator()( T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticeFlux3D<T, DESCRIPTOR>::SuperLatticeFlux3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
{
  this->getName() = "flux";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeFlux3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeFlux3D<T, DESCRIPTOR>::operator()( T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticeStrainRate3D<T, DESCRIPTOR>::SuperLatticeStrainRate3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 9), _converter(converter)
{
  this->getName() = "strainRate";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeStrainRate3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),this->_converter));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeStrainRate3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, typename DESCRIPTOR>
SuperLatticePhysStrainRate3D<T, DESCRIPTOR>::SuperLatticePhysStrainRate3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 9)
{
  this->getName() = "physStrainRate";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysStrainRate3D<T, DESCRIPTOR>(
        this->_sLattice.getExtendedBlockLattice(iC),
        this->_sLattice.getOverlap(),
        this->_converter)
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysStrainRate3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    const int loc = this->_sLattice.getLoadBalancer().loc(input[0]);
    return this->getBlockF(loc)(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, typename DESCRIPTOR>
SuperLatticeGeometry3D<T, DESCRIPTOR>::SuperLatticeGeometry3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1), _superGeometry(superGeometry),
    _material(material)
{
  this->getName() = "geometry";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new  BlockLatticeGeometry3D<T, DESCRIPTOR>(
                                 this->_sLattice.getBlockLattice(iC),
                                 this->_superGeometry.getBlockGeometry(iC),
                                 _material) );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeGeometry3D<T, DESCRIPTOR>::operator()( T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template<typename T, typename DESCRIPTOR>
SuperLatticeRank3D<T, DESCRIPTOR>::SuperLatticeRank3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice) : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "rank";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeRank3D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)) );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeRank3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
    return true;
  } else {
    return false;
  }
}


template<typename T, typename DESCRIPTOR>
SuperLatticeCuboid3D<T, DESCRIPTOR>::SuperLatticeCuboid3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice) : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "cuboid";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeCuboid3D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),
                                this->_sLattice.getLoadBalancer().glob(iC)) );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeCuboid3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
    return true;
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysPressure3D<T, DESCRIPTOR>::SuperLatticePhysPressure3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physPressure";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysPressure3D<T, DESCRIPTOR>(
        this->_sLattice.getExtendedBlockLattice(iC),
        this->_sLattice.getOverlap(),
        this->_converter)
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysPressure3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    const int loc = this->_sLattice.getLoadBalancer().loc(input[0]);
    return this->getBlockF(loc)(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysVelocity3D<T, DESCRIPTOR>::SuperLatticePhysVelocity3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter, bool print)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3), _print(print)
{
  this->getName() = "physVelocity";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysVelocity3D<T, DESCRIPTOR>(
        this->_sLattice.getExtendedBlockLattice(iC),
        this->_sLattice.getOverlap(),
        this->_converter,
        _print)
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysVelocity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    const int loc = this->_sLattice.getLoadBalancer().loc(input[0]);
    return this->getBlockF(loc)(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysExternal3D<T, DESCRIPTOR>::SuperLatticePhysExternal3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, T convFactorToPhysUnits,
  int offset, int size)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
{
  this->getName() = "physExtField";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysExternal3D<T, DESCRIPTOR>(
        this->_sLattice.getBlockLattice(iC), convFactorToPhysUnits,
        offset, size)
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysExternal3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    const int loc = this->_sLattice.getLoadBalancer().loc(input[0]);
    return this->getBlockF(loc)(output,&input[1]);
  } else {
    return false;
  }
}

template <typename T, typename DESCRIPTOR>
SuperLatticePhysExternalPorosity3D<T,DESCRIPTOR>::SuperLatticePhysExternalPorosity3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,1)
{
  this->getName() = "ExtPorosityField";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysExternalPorosity3D<T, DESCRIPTOR>(
        this->_sLattice.getExtendedBlockLattice(iC),
        this->_sLattice.getOverlap(),
        this->_converter)
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysExternalPorosity3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    const int loc = this->_sLattice.getLoadBalancer().loc(input[0]);
    return this->getBlockF(loc)(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysExternalVelocity3D<T, DESCRIPTOR>::SuperLatticePhysExternalVelocity3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3)
{
  this->getName() = "physVelExtField";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysExternalVelocity3D<T, DESCRIPTOR>(
        this->_sLattice.getBlockLattice(iC), this->_converter)
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysExternalVelocity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    const int loc = this->_sLattice.getLoadBalancer().loc(input[0]);
    return this->getBlockF(loc)(output, &input[1]);
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysExternalParticleVelocity3D<T, DESCRIPTOR>::SuperLatticePhysExternalParticleVelocity3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "ExtPartVelField";

  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticePhysExternalParticleVelocity3D<T, DESCRIPTOR>(
        sLattice.getExtendedBlockLattice(iC),
        converter)
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysExternalParticleVelocity3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  auto& load = this->_sLattice.getLoadBalancer();
  const int& globIC = input[0];

  if (load.rank(globIC) == singleton::mpi().getRank()) {
    const int overlap = this->_sLattice.getOverlap();

    int inputLocal[3] = { };
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;
    inputLocal[2] = input[3] + overlap;

    return this->getBlockF(load.loc(globIC))(output, inputLocal);
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysBoundaryForce3D<T, DESCRIPTOR>::SuperLatticePhysBoundaryForce3D(
  SuperLattice3D<T, DESCRIPTOR>&     sLattice,
  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "physBoundaryForce";
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticePhysBoundaryForce3D<T, DESCRIPTOR>(
        this->_sLattice.getBlockLattice(iC),
        _indicatorF->getBlockIndicatorF(iC),
        this->_converter));
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysBoundaryForce3D<T, DESCRIPTOR>::SuperLatticePhysBoundaryForce3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice,
  SuperGeometry3D<T>& superGeometry, const int material,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysBoundaryForce3D(sLattice,
                                    superGeometry.getMaterialIndicator(material),
                                    converter)
{ }

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysBoundaryForce3D<T, DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]))(output, &input[1]);
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePSMPhysForce3D<T, DESCRIPTOR>::SuperLatticePSMPhysForce3D(
  SuperLattice3D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  int mode_)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3)
{
  this->getName() = "PSMPhysForce";
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticePSMPhysForce3D<T, DESCRIPTOR>(
        this->_sLattice.getBlockLattice(iC),
        this->_converter, mode_));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePSMPhysForce3D<T, DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]))(output, &input[1]);
  } else {
    return false;
  }
}


template<typename T, typename DESCRIPTOR>
SuperLatticePhysWallShearStress3D<T, DESCRIPTOR>::SuperLatticePhysWallShearStress3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter,
  IndicatorF3D<T>& indicator)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 1),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physWallShearStress";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysWallShearStress3D<T, DESCRIPTOR>(
        this->_sLattice.getExtendedBlockLattice(iC),
        _superGeometry.getExtendedBlockGeometry(iC),
        this->_sLattice.getOverlap(),
        _material,
        this->_converter,
        indicator)
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysWallShearStress3D<T, DESCRIPTOR>::operator() (T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]))(output, &input[1]);
  } else {
    return false;
  }
}


template<typename T, typename DESCRIPTOR>
SuperLatticePhysCorrBoundaryForce3D<T, DESCRIPTOR>::SuperLatticePhysCorrBoundaryForce3D(
  SuperLattice3D<T, DESCRIPTOR>&     sLattice,
  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "physCorrBoundaryForce";
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticePhysCorrBoundaryForce3D<T, DESCRIPTOR>(
        this->_sLattice.getBlockLattice(iC),
        _indicatorF->getBlockIndicatorF(iC),
        this->_converter));
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysCorrBoundaryForce3D<T, DESCRIPTOR>::SuperLatticePhysCorrBoundaryForce3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice,
  SuperGeometry3D<T>& superGeometry, const int material,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysCorrBoundaryForce3D(sLattice,
                                        superGeometry.getMaterialIndicator(material),
                                        converter)
{ }

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysCorrBoundaryForce3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]))(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR, typename FIELD>
SuperLatticeField3D<T,DESCRIPTOR,FIELD>::SuperLatticeField3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, DESCRIPTOR::template size<FIELD>())
{
  this->getName() = "externalField";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; ++iC) {
    this->_blockF.emplace_back(
        new BlockLatticeField3D<T,DESCRIPTOR,FIELD>(this->_sLattice.getBlockLattice(iC)));
  }
}

template<typename T, typename DESCRIPTOR, typename FIELD>
bool SuperLatticeField3D<T,DESCRIPTOR,FIELD>::operator()(
  T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().isLocal(input[0])) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]))(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePorosity3D<T, DESCRIPTOR>::SuperLatticePorosity3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "porosity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePorosity3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePorosity3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticeVolumeFractionApproximation3D<T, DESCRIPTOR>::SuperLatticeVolumeFractionApproximation3D(
  SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  IndicatorF3D<T>& indicator, int refinementLevel,
  const UnitConverter<T,DESCRIPTOR>& converter, bool insideOut)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "volumeFractionApproximation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeVolumeFractionApproximation3D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),
                               superGeometry.getBlockGeometry(iC),
                               indicator, refinementLevel,
                               converter, insideOut));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeVolumeFractionApproximation3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysPermeability3D<T, DESCRIPTOR>::SuperLatticePhysPermeability3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "permeability";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticePhysPermeability3D<T, DESCRIPTOR>(
                                  this->_sLattice.getBlockLattice(iC), this->getConverter() ) );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysPermeability3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output, &input[1]);
  } else {
    return false;
  }
}


template<typename T, typename DESCRIPTOR>
SuperLatticePhysCroppedPermeability3D<T, DESCRIPTOR>::SuperLatticePhysCroppedPermeability3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "cropped_permeability";
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); iC++ ) {
    this->_blockF.emplace_back( new BlockLatticePhysCroppedPermeability3D<T, DESCRIPTOR>(
                                  this->_sLattice.getBlockLattice(iC), this->getConverter() ) );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysCroppedPermeability3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output, &input[1]);
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysDarcyForce3D<T, DESCRIPTOR>::SuperLatticePhysDarcyForce3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "alphaU";
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysDarcyForce3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  SuperLatticePhysPermeability3D<T, DESCRIPTOR> permeability(this->_sLattice, this->_converter);
  SuperLatticeVelocity3D<T, DESCRIPTOR> velocity(this->_sLattice);

  T nu = this->_converter.getPhysViscosity();
  T K;
  T u[velocity.getTargetDim()];
  permeability(&K,input);
  velocity(u,input);

  output[0] = -nu / K * u[0];
  output[1] = -nu / K * u[1];
  output[2] = -nu / K * u[2];

  return true;
}


template<typename T, typename DESCRIPTOR>
SuperEuklidNorm3D<T, DESCRIPTOR>::SuperEuklidNorm3D(
  SuperLatticeF3D<T, DESCRIPTOR>& f)
  : SuperLatticeF3D<T, DESCRIPTOR>(f.getSuperLattice(), 1), _f(f)
{
  this->getName() = "EuklidNorm(" + _f.getName() + ")";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockEuklidNorm3D<T, DESCRIPTOR>(f.getBlockF(iC)));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperEuklidNorm3D<T, DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]))(output,&input[1]);
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>::SuperLatticeInterpPhysVelocity3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, UnitConverter<T,DESCRIPTOR> const& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3)
{
  this->getName() = "InterpVelocity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int lociC = 0; lociC < maxC; lociC++) {
    int globiC = this->_sLattice.getLoadBalancer().glob(lociC);

    this->_blockF.emplace_back(
      new BlockLatticeInterpPhysVelocity3D<T, DESCRIPTOR>(
        sLattice.getExtendedBlockLattice(lociC),
        converter,
        &sLattice.getCuboidGeometry().get(globiC),
        sLattice.getOverlap())
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  return false;
}

template<typename T, typename DESCRIPTOR>
void SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>::operator()(T output[],
    const T input[], const int globiC)
{
  if (this->_sLattice.getLoadBalancer().isLocal(globiC)) {
    static_cast<BlockLatticeInterpPhysVelocity3D<T, DESCRIPTOR>*>(
      this->_blockF[this->_sLattice.getLoadBalancer().loc(globiC)].get()
    )->operator()(output, input);
  }
}


template <typename T, typename DESCRIPTOR>
SuperLatticePorousMomentumLossForce3D<T,DESCRIPTOR>::SuperLatticePorousMomentumLossForce3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
 std::vector<SmoothIndicatorF3D<T,T,true>* >& indicator, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T,DESCRIPTOR>(sLattice,converter,7*indicator.size())
{
  this->getName() = "physPorousMomentumLossForce";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticePorousMomentumLossForce3D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC), superGeometry.getBlockGeometry(iC), indicator, converter));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticePorousMomentumLossForce3D<T,DESCRIPTOR>::operator() (T output[],
    const int input[])
{
  for (int i=0; i<this->getTargetDim(); i++) {
    output[i] = 0.;
  }
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    int globiC = this->_sLattice.getLoadBalancer().glob(iC);
    if ( this->_sLattice.getLoadBalancer().rank(globiC) == singleton::mpi().getRank() ) {
      this->getBlockF(iC)(output,&input[1]);
    }
  }

#ifdef PARALLEL_MODE_MPI
  for (int i = 0; i < this->getTargetDim(); ++i) {
    singleton::mpi().reduceAndBcast(output[i], MPI_SUM);
  }
#endif
  return true;

}


template <typename T, typename DESCRIPTOR>
SuperLatticePhysBoundaryDistance3D<T,DESCRIPTOR>::SuperLatticePhysBoundaryDistance3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
 XMLreader const& xmlReader)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1),
    _superGeometry(superGeometry)
{
  this->getName() = "physBoundaryDistance";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticePhysBoundaryDistance3D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC), this->_superGeometry.getBlockGeometry(iC), xmlReader));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysBoundaryDistance3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
SuperLatticePhysTemperature3D<T,DESCRIPTOR,TDESCRIPTOR>::SuperLatticePhysTemperature3D(
  SuperLattice3D<T,TDESCRIPTOR>& sLattice, ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter)
  : SuperLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physTemperature";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysTemperature3D<T,DESCRIPTOR,TDESCRIPTOR>(this->_sLattice.getBlockLattice(iC), this->_converter));
  }
}

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
bool SuperLatticePhysTemperature3D<T,DESCRIPTOR,TDESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output, &input[1]);
  } else {
    return false;
  }
}

template<typename T,typename DESCRIPTOR, typename TDESCRIPTOR>
SuperLatticePhysHeatFlux3D<T,DESCRIPTOR,TDESCRIPTOR>::SuperLatticePhysHeatFlux3D(
  SuperLattice3D<T,TDESCRIPTOR>& sLattice, const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter)
  : SuperLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR>(sLattice, converter, 3)
{
  this->getName() = "physHeatFlux";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysHeatFlux3D<T,DESCRIPTOR,TDESCRIPTOR>(this->_sLattice.getBlockLattice(iC),
                               this->_converter));
  }
}

template<typename T,typename DESCRIPTOR,  typename TDESCRIPTOR>
bool SuperLatticePhysHeatFlux3D<T,DESCRIPTOR,TDESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}



template <typename T, typename DESCRIPTOR>
SuperLatticePhysPoreSizeDistribution3D<T,DESCRIPTOR>::SuperLatticePhysPoreSizeDistribution3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry, int material,
 XMLreader const& xmlReader)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1),
    _superGeometry(superGeometry)
{
  this->getName() = "physPoreSizeDistribution";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticePhysPoreSizeDistribution3D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC), this->_superGeometry.getBlockGeometry(iC), material, xmlReader));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysPoreSizeDistribution3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}



template<typename T,typename DESCRIPTOR, typename TDESCRIPTOR>
SuperLatticePhysTauFromBoundaryDistance3D<T,DESCRIPTOR,TDESCRIPTOR>::SuperLatticePhysTauFromBoundaryDistance3D(
  SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& sGeometry, XMLreader const& xmlReader, const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter, const T p, const T T_avg, const T c_p, const T beta, const T lambda_0, const T sigma, const T p_0, const T n_0)
  : SuperLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physTauFromBoundaryDistance";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysTauFromBoundaryDistance3D<T,DESCRIPTOR,TDESCRIPTOR>(this->_sLattice.getBlockLattice(iC), sGeometry.getBlockGeometry(iC), xmlReader, this->_converter, p, T_avg, c_p, beta, lambda_0, sigma, p_0, n_0));
  }
}

template<typename T,typename DESCRIPTOR,  typename TDESCRIPTOR>
bool SuperLatticePhysTauFromBoundaryDistance3D<T,DESCRIPTOR,TDESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  } else {
    return false;
  }
}


template <typename T, typename DESCRIPTOR, bool HLBM>
SuperLatticeIndicatorSmoothIndicatorIntersection3D<T,DESCRIPTOR,HLBM>::SuperLatticeIndicatorSmoothIndicatorIntersection3D (
  SuperLattice3D<T,DESCRIPTOR>& sLattice,
  SuperGeometry3D<T>& superGeometry,
  IndicatorF3D<T>& normalInd, SmoothIndicatorF3D<T,T,HLBM>& smoothInd )
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "Indicator-SmoothIndicator Intersection";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeIndicatorSmoothIndicatorIntersection3D<T,DESCRIPTOR,HLBM>(this->_sLattice.getExtendedBlockLattice(iC), superGeometry.getBlockGeometry(iC), normalInd, smoothInd));
  }
}

template <typename T, typename DESCRIPTOR, bool HLBM>
bool SuperLatticeIndicatorSmoothIndicatorIntersection3D<T,DESCRIPTOR,HLBM>::operator() (T output[], const int input[])
{
  output[0] = 0.;
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    int globiC = this->_sLattice.getLoadBalancer().glob(iC);
    if ( this->_sLattice.getLoadBalancer().rank(globiC) == singleton::mpi().getRank() ) {
      this->getBlockF(iC)(output,&input[1]);
    }
  }

#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(output[0], MPI_MAX);
#endif
  return true;

}

/////////// SuperLatticeGuoZhaoEpsilon3D /////////////////////////////////////////////

template<typename T, typename DESCRIPTOR>
SuperLatticeGuoZhaoEpsilon3D<T, DESCRIPTOR>::SuperLatticeGuoZhaoEpsilon3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "epsilon";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticeGuoZhaoEpsilon3D<T, DESCRIPTOR>(
        this->_sLattice.getExtendedBlockLattice(iC),
        this->_sLattice.getOverlap() )
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeGuoZhaoEpsilon3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    const int loc = this->_sLattice.getLoadBalancer().loc(input[0]);
    return this->getBlockF(loc)(output,&input[1]);
  } else {
    return false;
  }
}

/////////// SuperLatticeGuoZhaoPhysK3D /////////////////////////////////////////////

template<typename T, typename DESCRIPTOR>
SuperLatticeGuoZhaoPhysK3D<T, DESCRIPTOR>::SuperLatticeGuoZhaoPhysK3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physK";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticeGuoZhaoPhysK3D<T, DESCRIPTOR>(
        this->_sLattice.getExtendedBlockLattice(iC),
        this->_sLattice.getOverlap(),
        this->_converter)
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeGuoZhaoPhysK3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    const int loc = this->_sLattice.getLoadBalancer().loc(input[0]);
    return this->getBlockF(loc)(output,&input[1]);
  } else {
    return false;
  }
}


/////////// SuperLatticeGuoZhaoPhysBodyForce3D /////////////////////////////////////////////

template<typename T, typename DESCRIPTOR>
SuperLatticeGuoZhaoPhysBodyForce3D<T, DESCRIPTOR>::SuperLatticeGuoZhaoPhysBodyForce3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3)
{
  this->getName() = "physBodyForce";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticeGuoZhaoPhysBodyForce3D<T, DESCRIPTOR>(
        this->_sLattice.getExtendedBlockLattice(iC),
        this->_sLattice.getOverlap(),
        this->_converter)
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeGuoZhaoPhysBodyForce3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    const int loc = this->_sLattice.getLoadBalancer().loc(input[0]);
    return this->getBlockF(loc)(output,&input[1]);
  } else {
    return false;
  }
}

/////////// SuperLatticeTimeStepScale3D /////////////////////////////////////////////

template<typename T, typename DESCRIPTOR>
SuperLatticeTimeStepScale3D<T, DESCRIPTOR>::SuperLatticeTimeStepScale3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice,
  T oldTau, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, DESCRIPTOR::q)
{
  this->getName() = "latticeTimeStepScale";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticeTimeStepScale3D<T, DESCRIPTOR>(
        this->_sLattice.getExtendedBlockLattice(iC),
        oldTau,
        converter)
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeTimeStepScale3D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    const int loc = this->_sLattice.getLoadBalancer().loc(input[0]);
    return this->getBlockF(loc)(output,&input[1]);
  } else {
    return false;
  }
}

}  // end namespace olb

#endif
