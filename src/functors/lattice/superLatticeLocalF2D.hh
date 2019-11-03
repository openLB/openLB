/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause
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

#ifndef SUPER_LATTICE_LOCAL_F_2D_HH
#define SUPER_LATTICE_LOCAL_F_2D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<limits>

#include "superLatticeLocalF2D.h"
#include "blockLatticeLocalF2D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry2D.h"
#include "indicator/superIndicatorF2D.h"

namespace olb {

template<typename T,typename DESCRIPTOR>
SuperLatticeDissipation2D<T,DESCRIPTOR>::SuperLatticeDissipation2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1), _converter(converter)
{
  this->getName() = "dissipation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeDissipation2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),this->_converter));
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticeDissipation2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{

  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  ////  int lociz= input[3];
  //  if ( this->sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
  //    // local coords are given, fetch local cell and compute value(s)
  //    T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  //    int overlap = this->sLattice.getOverlap();
  //    this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap/*, lociz+overlap*/).computeAllMomenta(rho, uTemp, pi);

  //    T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
  //    if (util::TensorVal<DESCRIPTOR >::n == 6)
  //      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];

  //    T nuLattice = converter.getLatticeNu();
  //    T omega = converter.getOmega();
  //    T finalResult = PiNeqNormSqr*nuLattice*pow(omega*descriptors::invCs2<T,DESCRIPTOR>(),2)/rho/2.;

  //    return std::vector<T>(1,finalResult);
  //  } else {
  //    return std::vector<T>(); // empty vector
  //  }
  return false;
}

template<typename T,typename DESCRIPTOR>
SuperLatticePhysDissipation2D<T,DESCRIPTOR>::SuperLatticePhysDissipation2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physDissipation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysDissipation2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),this->_converter));
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticePhysDissipation2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  }
  else {
    return false;
  }
}

template<typename T,typename DESCRIPTOR>
SuperLatticeDensity2D<T,DESCRIPTOR>::SuperLatticeDensity2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "density";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeDensity2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)) );
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticeDensity2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  }
  else {
    return false;
  }
}

template<typename T,typename DESCRIPTOR>
SuperLatticeVelocity2D<T,DESCRIPTOR>::SuperLatticeVelocity2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 2)
{
  this->getName() = "velocity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeVelocity2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)));
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticeVelocity2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  }
  else {
    return false;
  }
}

template<typename T,typename DESCRIPTOR>
SuperLatticePhysStrainRate2D<T,DESCRIPTOR>::SuperLatticePhysStrainRate2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 4)
{
  this->getName() = "physStrainRate";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysStrainRate2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),this->_converter));
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticePhysStrainRate2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{

  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  }
  else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysWallShearStress2D<T, DESCRIPTOR>::SuperLatticePhysWallShearStress2D(
  SuperLattice2D<T, DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter,
  IndicatorF2D<T>& indicator)
  : SuperLatticePhysF2D<T, DESCRIPTOR>(sLattice, converter, 1),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physWallShearStress";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysWallShearStress2D<T, DESCRIPTOR>(
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
bool SuperLatticePhysWallShearStress2D<T, DESCRIPTOR>::operator() (T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  }
  else {
    return false;
  }
}

template<typename T,typename DESCRIPTOR>
SuperLatticeGeometry2D<T,DESCRIPTOR>::SuperLatticeGeometry2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1), _superGeometry(superGeometry),
    _material(material)
{
  this->getName() = "geometry";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new  BlockLatticeGeometry2D<T,DESCRIPTOR>(
                                  this->_sLattice.getBlockLattice(iC),
                                  this->_superGeometry.getBlockGeometry(iC),
                                  _material) );
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticeGeometry2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  }
  else {
    return false;
  }
}

template<typename T,typename DESCRIPTOR>
SuperLatticeRank2D<T,DESCRIPTOR>::SuperLatticeRank2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "rank";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeRank2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)) );
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticeRank2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
    return true;
  }
  else {
    return false;
  }
}

template<typename T,typename DESCRIPTOR>
SuperLatticeCuboid2D<T,DESCRIPTOR>::SuperLatticeCuboid2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "cuboid";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeCuboid2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),
                                this->_sLattice.getLoadBalancer().glob(iC)) );
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticeCuboid2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
    return true;
  }
  else {
    return false;
  }
}

template<typename T,typename DESCRIPTOR>
SuperLatticePhysPressure2D<T,DESCRIPTOR>::SuperLatticePhysPressure2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physPressure";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysPressure2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC), this->_converter));
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticePhysPressure2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output, &input[1]);
  }
  else {
    return false;
  }
}

template<typename T,typename DESCRIPTOR>
SuperLatticePhysVelocity2D<T,DESCRIPTOR>::SuperLatticePhysVelocity2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "physVelocity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysVelocity2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),
                               this->_converter));
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticePhysVelocity2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  }
  else {
    return false;
  }


}


template<typename T, typename DESCRIPTOR, typename FIELD>
SuperLatticeField2D<T,DESCRIPTOR,FIELD>::SuperLatticeField2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, DESCRIPTOR::template size<FIELD>())
{
  this->getName() = "ExtField";
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); iC++ ) {
    this->_blockF.emplace_back(
      new BlockLatticeField2D<T,DESCRIPTOR,FIELD>(this->_sLattice.getBlockLattice(iC)));
  }
}

template<typename T, typename DESCRIPTOR, typename FIELD>
bool SuperLatticeField2D<T,DESCRIPTOR,FIELD>::operator()(
  T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().isLocal(input[0])) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]))(output,&input[1]);
  }
  else {
    return false;
  }
}

template <typename T,typename DESCRIPTOR>
SuperLatticePhysExternalPorosity2D<T,DESCRIPTOR>::SuperLatticePhysExternalPorosity2D
(SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "ExtPorosityField";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysExternalPorosity2D<T, DESCRIPTOR>(
        this->_sLattice.getExtendedBlockLattice(iC),
        this->_sLattice.getOverlap(),
        this->_converter)
    );
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticePhysExternalPorosity2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    const int loc = this->_sLattice.getLoadBalancer().loc(input[0]);
    return this->getBlockF(loc)(output, &input[1]);
  }
  else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticeVolumeFractionApproximation2D<T, DESCRIPTOR>::SuperLatticeVolumeFractionApproximation2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  IndicatorF2D<T>& indicator, int refinementLevel,
  const UnitConverter<T,DESCRIPTOR>& converter, bool insideOut)
  : SuperLatticeF2D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "volumeFractionApproximation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeVolumeFractionApproximation2D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),
                               superGeometry.getBlockGeometry(iC),
                               indicator, refinementLevel,
                               converter, insideOut));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeVolumeFractionApproximation2D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  }
  else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticeVolumeFractionPolygonApproximation2D<T, DESCRIPTOR>::SuperLatticeVolumeFractionPolygonApproximation2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  IndicatorF2D<T>& indicator, const UnitConverter<T,DESCRIPTOR>& converter, bool insideOut)
  : SuperLatticeF2D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "volumeFractionPolygonApproximation";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeVolumeFractionPolygonApproximation2D<T, DESCRIPTOR>(this->_sLattice.getBlockLattice(iC),
                               superGeometry.getBlockGeometry(iC),
                               indicator, converter, insideOut));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeVolumeFractionPolygonApproximation2D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  }
  else {
    return false;
  }
}


template<typename T,typename DESCRIPTOR>
SuperLatticePhysExternalVelocity2D<T,DESCRIPTOR>::SuperLatticePhysExternalVelocity2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "ExtVelocityField";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysExternalVelocity2D<T, DESCRIPTOR>(
        this->_sLattice.getExtendedBlockLattice(iC),
        this->_sLattice.getOverlap(),
        this->_converter)
    );
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticePhysExternalVelocity2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    const int loc = this->_sLattice.getLoadBalancer().loc(input[0]);
    return this->getBlockF(loc)(output, &input[1]);
  }
  else {
    return false;
  }
}

template<typename T,typename DESCRIPTOR>
SuperLatticePhysExternalParticleVelocity2D<T,DESCRIPTOR>::SuperLatticePhysExternalParticleVelocity2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "ExtPartVelField";

  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticePhysExternalParticleVelocity2D<T, DESCRIPTOR>(
        sLattice.getExtendedBlockLattice(iC),
        converter)
    );
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticePhysExternalParticleVelocity2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  auto& load = this->_sLattice.getLoadBalancer();
  const int& globIC = input[0];

  if (load.rank(globIC) == singleton::mpi().getRank()) {
    const int overlap = this->_sLattice.getOverlap();

    int inputLocal[2] = { };
    inputLocal[0] = input[1] + overlap;
    inputLocal[1] = input[2] + overlap;

    return this->getBlockF(load.loc(globIC))(output, inputLocal);
  }
  else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysBoundaryForce2D<T, DESCRIPTOR>::SuperLatticePhysBoundaryForce2D(
  SuperLattice2D<T, DESCRIPTOR>&     sLattice,
  FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T, DESCRIPTOR>(sLattice, converter, 2),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "physBoundaryForce";
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticePhysBoundaryForce2D<T, DESCRIPTOR>(
        this->_sLattice.getBlockLattice(iC),
        _indicatorF->getBlockIndicatorF(iC),
        this->_converter));
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysBoundaryForce2D<T, DESCRIPTOR>::SuperLatticePhysBoundaryForce2D(
  SuperLattice2D<T, DESCRIPTOR>& sLattice,
  SuperGeometry2D<T>& superGeometry, const int material,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysBoundaryForce2D(sLattice,
                                    superGeometry.getMaterialIndicator(material),
                                    converter)
{ }

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysBoundaryForce2D<T, DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]))(output, &input[1]);
  }
  else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePSMPhysForce2D<T, DESCRIPTOR>::SuperLatticePSMPhysForce2D(
  SuperLattice2D<T, DESCRIPTOR>&     sLattice,
  const UnitConverter<T,DESCRIPTOR>& converter,
  int mode_)
  : SuperLatticePhysF2D<T, DESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "PSMPhysForce";
  for (int iC = 0; iC < this->_sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticePSMPhysForce2D<T, DESCRIPTOR>(
        this->_sLattice.getBlockLattice(iC),
        this->_converter,
        mode_));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePSMPhysForce2D<T, DESCRIPTOR>::operator() (T output[], const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]))(output, &input[1]);
  }
  else {
    return false;
  }
}

/*****************************NEW*****************************/

template<typename T,typename DESCRIPTOR>
SuperLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR>::SuperLatticePhysCorrBoundaryForce2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 2),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "physCorrBoundaryForce";
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  //  int lociz= input[3];

  //  std::vector<T> force(3, T());
  //  if ( this->sLattice.getLoadBalancer().rank(globIC) == singleton::mpi().getRank() ) {
  //    int globX = (int)this->sLattice.getCuboidGeometry().get(globIC).get_globPosX() + locix;
  //    int globY = (int)this->sLattice.getCuboidGeometry().get(globIC).get_globPosY() + lociy;
  //    int globZ = (int)this->sLattice.getCuboidGeometry().get(globIC).get_globPosZ() + lociz;
  //    if(superGeometry.getMaterial(globX, globY, globZ) == material) {
  //      for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
  //        // Get direction
  //        const int* c = DESCRIPTOR::c[iPop];
  //        // Get next cell located in the current direction
  //        // Check if the next cell is a fluid node
  //        if (superGeometry.getMaterial(globX + c[0], globY + c[1], globZ + c[2]) == 1) {
  //          int overlap = this->sLattice.getOverlap();
  //          // Get f_q of next fluid cell where l = opposite(q)
  //          T f = this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap + c[0], lociy+overlap + c[1], lociz+overlap + c[2])[iPop];
  //          // Get f_l of the boundary cell
  //          // Add f_q and f_opp
  //          f += this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap)[util::opposite<DESCRIPTOR >(iPop)];
  //          // Update force
  //          force[0] -= c[0]*(f-2.*descriptors::t<T,DESCRIPTOR>(iPop));
  //          force[1] -= c[1]*(f-2.*descriptors::t<T,DESCRIPTOR>(iPop));
  //          force[2] -= c[2]*(f-2.*descriptors::t<T,DESCRIPTOR>(iPop));
  //        }
  //      }
  //      force[0]=this->converter.physForce(force[0]);
  //      force[1]=this->converter.physForce(force[1]);
  //      force[2]=this->converter.physForce(force[2]);
  //      return force;
  //    }
  //    else {
  //      return force;
  //    }
  //  }
  //  else {
  //    return force;
  //  }
  return false;
}

template<typename T,typename DESCRIPTOR>
SuperLatticePorosity2D<T,DESCRIPTOR>::SuperLatticePorosity2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 1),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "porosity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePorosity2D<T,DESCRIPTOR>(
                                 this->_sLattice.getBlockLattice(iC),
                                 this->_superGeometry.getBlockGeometry(iC),
                                 _material,
                                 converter));
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticePorosity2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  }
  else {
    return false;
  }
}

template<typename T,typename DESCRIPTOR>
SuperLatticePhysPermeability2D<T,DESCRIPTOR>::SuperLatticePhysPermeability2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 1),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "permeability";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticePhysPermeability2D<T,DESCRIPTOR>(
                                  this->_sLattice.getBlockLattice(iC),
                                  this->_superGeometry.getBlockGeometry(iC),
                                  _material,
                                  this->getConverter() ) );
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticePhysPermeability2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  ////  int lociz= input[3];

  //  T* value = new T[1];
  //  int overlap = this->sLattice.getOverlap();
  //  this->sLattice.getExtendedBlockLattice(this->sLattice.getLoadBalancer().loc(globIC)).get(locix+overlap, lociy+overlap/*, lociz+overlap*/).computeField(0,1,value);
  //  std::vector<T> result(1,this->converter.physPermeability(value[0]));
  //  delete value;
  //  if (!(result[0]<42)&&!(result[0]>42)&&!(result[0]==42)) result[0] = 999999;
  //  if (isinf(result[0])) result[0] = 1e6;
  //  return result;
  return false;
}

template<typename T,typename DESCRIPTOR>
SuperLatticePhysDarcyForce2D<T,DESCRIPTOR>::SuperLatticePhysDarcyForce2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 2),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "alphaU";
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticePhysDarcyForce2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  //  SuperLatticePhysPermeability2D<T,DESCRIPTOR> permeability(this->sLattice,this->superGeometry,this->material,this->converter);
  //  SuperLatticeVelocity2D<T,DESCRIPTOR> velocity(this->sLattice);

  //  T nu = this->converter.getCharNu();
  //  T K = permeability(input)[0];
  //  std::vector<T> u = velocity(input);

  //  std::vector<T> result(2,T());
  //  result[0] = -nu/K*u[0];
  //  result[1] = -nu/K*u[1];
  ////  result[2] = -nu/K*u[2];

  //  return result;
  return false;
}

template<typename T,typename DESCRIPTOR>
SuperEuklidNorm2D<T,DESCRIPTOR>::SuperEuklidNorm2D(
  SuperLatticeF2D<T,DESCRIPTOR>& f)
  : SuperLatticeF2D<T,DESCRIPTOR>(f.getSuperLattice(), 1), _f(f)
{
  this->getName() = "l2(" + _f.getName() + ")";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockEuklidNorm2D<T,DESCRIPTOR>(f.getBlockF(iC)));
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperEuklidNorm2D<T,DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  }
  else {
    return false;
  }
}

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
SuperLatticePhysTemperature2D<T,DESCRIPTOR,TDESCRIPTOR>::SuperLatticePhysTemperature2D(
  SuperLattice2D<T,TDESCRIPTOR>& sLattice, ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter)
  : SuperLatticeThermalPhysF2D<T,DESCRIPTOR,TDESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physTemperature";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysTemperature2D<T,DESCRIPTOR,TDESCRIPTOR>(this->_sLattice.getBlockLattice(iC), this->_converter));
  }
}

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
bool SuperLatticePhysTemperature2D<T,DESCRIPTOR,TDESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output, &input[1]);
  }
  else {
    return false;
  }
}

template<typename T,typename DESCRIPTOR, typename TDESCRIPTOR>
SuperLatticePhysHeatFlux2D<T,DESCRIPTOR,TDESCRIPTOR>::SuperLatticePhysHeatFlux2D(
  SuperLattice2D<T,TDESCRIPTOR>& sLattice, const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter)
  : SuperLatticeThermalPhysF2D<T,DESCRIPTOR,TDESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "physHeatFlux";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysHeatFlux2D<T,DESCRIPTOR,TDESCRIPTOR>(this->_sLattice.getBlockLattice(iC),
                               this->_converter));
  }
}

template<typename T,typename DESCRIPTOR,  typename TDESCRIPTOR>
bool SuperLatticePhysHeatFlux2D<T,DESCRIPTOR,TDESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(this->_sLattice.getLoadBalancer().loc(input[0]) )(output,&input[1]);
  }
  else {
    return false;
  }
}

template <typename T, typename DESCRIPTOR>
SuperLatticePorousMomentumLossForce2D<T,DESCRIPTOR>::SuperLatticePorousMomentumLossForce2D
(SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
 std::vector<SmoothIndicatorF2D<T,T,true>* >& indicator, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice,converter,4*indicator.size())
{
  this->getName() = "physPorousMomentumLossForce";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticePorousMomentumLossForce2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC), superGeometry.getBlockGeometry(iC), indicator, converter));
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticePorousMomentumLossForce2D<T,DESCRIPTOR>::operator() (T output[],
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

template <typename T, typename DESCRIPTOR, bool HLBM>
SuperLatticeIndicatorSmoothIndicatorIntersection2D<T,DESCRIPTOR,HLBM>::SuperLatticeIndicatorSmoothIndicatorIntersection2D (
  SuperLattice2D<T,DESCRIPTOR>& sLattice,
  SuperGeometry2D<T>& superGeometry,
  IndicatorF2D<T>& normalInd, SmoothIndicatorF2D<T,T,HLBM>& smoothInd )
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "Indicator-SmoothIndicator Intersection";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeIndicatorSmoothIndicatorIntersection2D<T,DESCRIPTOR,HLBM>(this->_sLattice.getExtendedBlockLattice(iC), superGeometry.getBlockGeometry(iC), normalInd, smoothInd));
  }
}

template <typename T, typename DESCRIPTOR, bool HLBM>
bool SuperLatticeIndicatorSmoothIndicatorIntersection2D<T,DESCRIPTOR,HLBM>::operator() (T output[], const int input[])
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

/////////// SuperLatticeGuoZhaoEpsilon2D /////////////////////////////////////////////

template<typename T, typename DESCRIPTOR>
SuperLatticeGuoZhaoEpsilon2D<T, DESCRIPTOR>::SuperLatticeGuoZhaoEpsilon2D(
  SuperLattice2D<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T, DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "epsilon";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeGuoZhaoEpsilon2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC)));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeGuoZhaoEpsilon2D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output, &input[1]);
  }
  else {
    return false;
  }
}

/////////// SuperLatticeGuoZhaoPhysK2D /////////////////////////////////////////////

template<typename T, typename DESCRIPTOR>
SuperLatticeGuoZhaoPhysK2D<T, DESCRIPTOR>::SuperLatticeGuoZhaoPhysK2D(
  SuperLattice2D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T, DESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physK";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeGuoZhaoPhysK2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC), this->_converter));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeGuoZhaoPhysK2D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output, &input[1]);
  }
  else {
    return false;
  }
}

/////////// SuperLatticeGuoZhaoPhysBodyForce2D /////////////////////////////////////////////

template<typename T, typename DESCRIPTOR>
SuperLatticeGuoZhaoPhysBodyForce2D<T, DESCRIPTOR>::SuperLatticeGuoZhaoPhysBodyForce2D(
  SuperLattice2D<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T, DESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "physBodyForce";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeGuoZhaoPhysBodyForce2D<T,DESCRIPTOR>(this->_sLattice.getBlockLattice(iC), this->_converter));
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeGuoZhaoPhysBodyForce2D<T, DESCRIPTOR>::operator()(T output[],
    const int input[])
{
  if (this->_sLattice.getLoadBalancer().rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF( this->_sLattice.getLoadBalancer().loc(input[0]) )(output, &input[1]);
  }
  else {
    return false;
  }
}

}  // end namespace olb

#endif
