/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#ifndef BLOCK_LATTICE_LOCAL_F_3D_HH
#define BLOCK_LATTICE_LOCAL_F_3D_HH


#include<cmath>
#include<math.h>

#include "blockLatticeLocalF3D.h"
#include "blockBaseF3D.h"
#include "core/blockLatticeStructure3D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
BlockLatticeFpop3D<T, DESCRIPTOR>::BlockLatticeFpop3D(BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, DESCRIPTOR::q)
{
  this->getName() = "fPop";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeFpop3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    output[iPop] =
      this->_blockLattice.get(input[0], input[1], input[2])[iPop]
      + descriptors::t<T,DESCRIPTOR>(iPop);
  }
  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticeDissipation3D<T, DESCRIPTOR>::BlockLatticeDissipation3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1), _converter(converter)
{
  this->getName() = "dissipation";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeDissipation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_blockLattice.get(input[0], input[1], input[2]).computeAllMomenta(rho,
      uTemp,
      pi);

  T PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                    + pi[5] * pi[5];
  }

  T nuLattice = _converter.getLatticeViscosity();
  T omega = 1. / _converter.getLatticeRelaxationTime();
  output[0] = PiNeqNormSqr * nuLattice
              * pow(omega * descriptors::invCs2<T,DESCRIPTOR>(), 2) / rho / 2.;

  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticePhysDissipation3D<T, DESCRIPTOR>::BlockLatticePhysDissipation3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
  int overlap,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _overlap(overlap),
    _converter(converter)
{
  this->getName() = "physDissipation";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysDissipation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_blockLattice.get(
    input[0]+_overlap, input[1]+_overlap, input[2]+_overlap
  ).computeAllMomenta(rho, uTemp, pi);

  T PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                    + pi[5] * pi[5];
  }

  T nuLattice = this->_converter.getLatticeViscosity();
  T omega = 1. / this->_converter.getLatticeRelaxationTime();
  T dt = this->_converter.getConversionFactorTime();
  output[0] = PiNeqNormSqr * nuLattice
              * pow(omega * descriptors::invCs2<T,DESCRIPTOR>() / rho, 2) / 2.
              * this->_converter.getPhysViscosity() / nuLattice / dt / dt;

  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticeEffevtiveDissipation3D<T, DESCRIPTOR>::BlockLatticeEffevtiveDissipation3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst,
  LESDynamics<T, DESCRIPTOR>& LESdynamics)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _converter(converter), _smagoConst(smagoConst), _LESdynamics(LESdynamics)
{
  this->getName() = "EffevtiveDissipation";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeEffevtiveDissipation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_blockLattice.get(input[0], input[1], input[2]).computeAllMomenta(rho,
      uTemp,
      pi);

  T PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                    + pi[5] * pi[5];
  }

  T omegaEff = _LESdynamics.getEffectiveOmega(this->_blockLattice.get(input[0], input[1], input[2]));
  T nuEff = ((1./omegaEff)-0.5)/descriptors::invCs2<T,DESCRIPTOR>();

  output[0] = PiNeqNormSqr * (nuEff)
              * pow(omegaEff * descriptors::invCs2<T,DESCRIPTOR>() / rho, 2)  / 2.;

  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticePhysEffevtiveDissipation3D<T, DESCRIPTOR>::BlockLatticePhysEffevtiveDissipation3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst,
  LESDynamics<T, DESCRIPTOR>& LESdynamics)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _converter(converter), _smagoConst(smagoConst), _LESdynamics(LESdynamics)
{
  this->getName() = "physEffevtiveDissipation";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysEffevtiveDissipation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_blockLattice.get(input[0], input[1], input[2]).computeAllMomenta(rho,
      uTemp,
      pi);

  T PiNeqNormSqr = pi[0] * pi[0] + 2. * pi[1] * pi[1] + pi[2] * pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2] * pi[2] + pi[3] * pi[3] + 2. * pi[4] * pi[4]
                    + pi[5] * pi[5];
  }

  T dt = 1./ _converter.getConversionFactorTime();
  T omegaEff = _LESdynamics.getEffectiveOmega(this->_blockLattice.get(input[0], input[1], input[2]));
  T nuEff = ((1./omegaEff)-0.5)/descriptors::invCs2<T,DESCRIPTOR>();  // BGK shear viscosity definition

  output[0] = PiNeqNormSqr * nuEff
              * pow(omegaEff * descriptors::invCs2<T,DESCRIPTOR>() / rho, 2) / 2.
              * _converter.getPhysViscosity() / _converter.getLatticeViscosity() / dt / dt;

  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticeDensity3D<T, DESCRIPTOR>::BlockLatticeDensity3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1)
{
  this->getName() = "density";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeDensity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = this->_blockLattice.get(input[0], input[1], input[2]).computeRho();
  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticeVelocity3D<T, DESCRIPTOR>::BlockLatticeVelocity3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3)
{
  this->getName() = "velocity";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeVelocity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho;
  this->_blockLattice.get(input[0], input[1], input[2]).computeRhoU(rho, output);
  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticeExternalVelocity3D<T, DESCRIPTOR>::BlockLatticeExternalVelocity3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3)
{
  this->getName() = "externalVelocity";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeExternalVelocity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T* ExtVel = this->_blockLattice.get(input[0], input[1], input[2]).template getFieldPointer<descriptors::VELOCITY>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    output[iVel] = ExtVel[iVel];
  }
  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticeFlux3D<T, DESCRIPTOR>::BlockLatticeFlux3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3)
{
  this->getName() = "flux";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeFlux3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  this->_blockLattice.get(input[0], input[1], input[2]).computeJ(output);
  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticeStrainRate3D<T, DESCRIPTOR>::BlockLatticeStrainRate3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 9)
{
  this->getName() = "strainRate";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeStrainRate3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_blockLattice.get(input[0], input[1], input[2]).computeAllMomenta(rho,
      uTemp,
      pi);

  T omega = 1. / this->_converter.getLatticeRelaxationTime();

  output[0] = -pi[0] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2.;
  output[1] = -pi[1] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2.;
  output[2] = -pi[2] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2.;
  output[3] = -pi[1] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2.;
  output[4] = -pi[3] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2.;
  output[5] = -pi[4] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2.;
  output[6] = -pi[2] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2.;
  output[7] = -pi[4] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2.;
  output[8] = -pi[5] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2.;

  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticePhysStrainRate3D<T, DESCRIPTOR>::BlockLatticePhysStrainRate3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
  int overlap,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 9),
    _overlap(overlap)
{
  this->getName() = "strainRate";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysStrainRate3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_blockLattice.get(
    input[0]+_overlap, input[1]+_overlap, input[2]+_overlap
  ).computeAllMomenta(rho, uTemp, pi);

  T omega = 1. / this->_converter.getLatticeRelaxationTime();
  T dt = this->_converter.getConversionFactorTime();

  output[0] = -pi[0] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[1] = -pi[1] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[2] = -pi[2] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[3] = -pi[1] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[4] = -pi[3] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[5] = -pi[4] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[6] = -pi[2] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[7] = -pi[4] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;
  output[8] = -pi[5] * omega * descriptors::invCs2<T,DESCRIPTOR>() / rho / 2. / dt;

  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticeGeometry3D<T, DESCRIPTOR>::BlockLatticeGeometry3D(BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
    BlockGeometryStructure3D<T>& blockGeometry, int material)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _blockGeometry(blockGeometry),
    _material(material)
{
  this->getName() = "geometry";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeGeometry3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = _blockGeometry.getMaterial(input[0], input[1], input[2]);

  if (_material != -1) {
    if ( util::nearZero(_material-output[0]) ) {
      output[0] = 1.;
      return true;
    }
    else {
      output[0] = 0.;
      return true;
    }
  }
  return true;
}


template<typename T, typename DESCRIPTOR>
BlockLatticeRank3D<T,DESCRIPTOR>::BlockLatticeRank3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1)
{
  this->getName() = "rank";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeRank3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = singleton::mpi().getRank() + 1;
  return true;
}


template<typename T, typename DESCRIPTOR>
BlockLatticeCuboid3D<T,DESCRIPTOR>::BlockLatticeCuboid3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, int iC)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1), _iC(iC)
{
  this->getName() = "cuboid";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeCuboid3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = _iC + 1;
  return true;
}



template<typename T, typename DESCRIPTOR>
BlockLatticePhysPressure3D<T, DESCRIPTOR>::BlockLatticePhysPressure3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
  int overlap,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 1),
    _overlap(overlap)
{
  this->getName() = "physPressure";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysPressure3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  // lattice pressure = c_s^2 ( rho -1 )
  T latticePressure = ( this->_blockLattice.get(input[0]+_overlap, input[1]+_overlap, input[2]+_overlap).computeRho() - 1.0) / descriptors::invCs2<T,DESCRIPTOR>();
  output[0] = this->_converter.getPhysPressure(latticePressure);

  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticePhysVelocity3D<T, DESCRIPTOR>::BlockLatticePhysVelocity3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
  int overlap,
  const UnitConverter<T,DESCRIPTOR>& converter,
  bool print)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 3),
    _overlap(overlap),
    _print(print)
{
  this->getName() = "physVelocity";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysVelocity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (_print) {
    std::cout << input[0] << " " << input[1] << " " << input[2] << " | "
              << singleton::mpi().getRank() << std::endl;
  }

  T rho;
  this->_blockLattice.get(
    input[0]+_overlap, input[1]+_overlap, input[2]+_overlap).computeRhoU(rho, output);
  output[0] = this->_converter.getPhysVelocity(output[0]);
  output[1] = this->_converter.getPhysVelocity(output[1]);
  output[2] = this->_converter.getPhysVelocity(output[2]);

  return true;
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysExternalPorosity3D<T,DESCRIPTOR>::BlockLatticePhysExternalPorosity3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
  int overlap,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,2),
    _overlap(overlap)
{
  this->getName() = "ExtPorosityField";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysExternalPorosity3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  this->_blockLattice.get(
    input[0]+_overlap, input[1]+_overlap, input[2]+_overlap
  ).template computeField<descriptors::POROSITY>(output);
  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticePhysExternalVelocity3D<T, DESCRIPTOR>::BlockLatticePhysExternalVelocity3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 3)
{
  this->getName() = "physVelExtField";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysExternalVelocity3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  this->_blockLattice.get(input[0], input[1], input[2]).template computeField<descriptors::VELOCITY>(output);
  output[0] = this->_converter.getPhysVelocity(output[0]);
  output[1] = this->_converter.getPhysVelocity(output[1]);
  output[2] = this->_converter.getPhysVelocity(output[2]);
  return true;
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysExternalParticleVelocity3D<T,DESCRIPTOR>::BlockLatticePhysExternalParticleVelocity3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,2)
{
  this->getName() = "ExtParticleVelocityField";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysExternalParticleVelocity3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  const T* velocity_numerator   = this->blockLattice.get(input).template getFieldPointer<descriptors::VELOCITY_NUMERATOR>();
  const T* velocity_denominator = this->blockLattice.get(input).template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>();

  if (velocity_denominator[0] > std::numeric_limits<T>::epsilon()) {
    output[0]=this->_converter.getPhysVelocity(velocity_numerator[0]/velocity_denominator[0]);
    output[1]=this->_converter.getPhysVelocity(velocity_numerator[1]/velocity_denominator[0]);
    output[2]=this->_converter.getPhysVelocity(velocity_numerator[2]/velocity_denominator[0]);
    return true;
  }
  output[0]=this->_converter.getPhysVelocity(velocity_numerator[0]);
  output[1]=this->_converter.getPhysVelocity(velocity_numerator[1]);
  output[2]=this->_converter.getPhysVelocity(velocity_numerator[2]);
  return true;
}


template<typename T, typename DESCRIPTOR>
BlockLatticePhysExternal3D<T, DESCRIPTOR>::BlockLatticePhysExternal3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
  T convFactorToPhysUnits, int offset, int size)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3),
    _convFactorToPhysUnits(convFactorToPhysUnits),
    _offset(offset), _size(size)
{
  this->getName() = "physExtField";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysExternal3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  const auto& cell = this->_blockLattice.get(input[0], input[1], input[2]);
  output[0] = _convFactorToPhysUnits*cell[_offset+0];
  output[1] = _convFactorToPhysUnits*cell[_offset+1];
  output[2] = _convFactorToPhysUnits*cell[_offset+2];
  return true;
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysBoundaryForce3D<T,DESCRIPTOR>::BlockLatticePhysBoundaryForce3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
  BlockIndicatorF3D<T>&                  indicatorF,
  const UnitConverter<T,DESCRIPTOR>&     converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice, converter, 3),
    _indicatorF(indicatorF),
    _blockGeometry(indicatorF.getBlockGeometryStructure())
{
  this->getName() = "physBoundaryForce";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysBoundaryForce3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = T();
  }

  if (_indicatorF(input)) {
    for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
      // Get direction
      const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
      // Get next cell located in the current direction
      // Check if the next cell is a fluid node
      if (_blockGeometry.get(input[0] + c[0], input[1] + c[1], input[2] + c[2]) == 1) {
        // Get f_q of next fluid cell where l = opposite(q)
        T f = this->_blockLattice.get(input[0] + c[0], input[1] + c[1], input[2] + c[2])[iPop];
        // Get f_l of the boundary cell
        // Add f_q and f_opp
        f += this->_blockLattice.get(input)[util::opposite<DESCRIPTOR>(iPop)];
        // Update force
        for (int i = 0; i < this->getTargetDim(); ++i) {
          output[i] -= c[i] * f;
        }
      }
    }
    for (int i = 0; i < this->getTargetDim(); ++i) {
      output[i] = this->_converter.getPhysForce(output[i]);
    }
  }
  return true;
}

template <typename T, typename DESCRIPTOR>
BlockLatticePSMPhysForce3D<T,DESCRIPTOR>::BlockLatticePSMPhysForce3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
  const UnitConverter<T,DESCRIPTOR>&     converter,
  int mode_)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice, converter, 3)
{
  this->getName() = "physPSMForce";
  mode = (Mode) mode_;
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePSMPhysForce3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = T();
  }

  T epsilon = 1. - *(this->_blockLattice.get(input).template getFieldPointer<descriptors::POROSITY>());

  //if ((epsilon > 1e-5 && epsilon < 1 - 1e-5)) {
  if ((epsilon > 1e-5)) {
    T rho, u[DESCRIPTOR::d], u_s[DESCRIPTOR::d];

    for (int i = 0; i < DESCRIPTOR::d; i++) {
      u_s[i] = (this->_blockLattice.get(input).template getFieldPointer<descriptors::VELOCITY_SOLID>())[i];
    }
    T paramA = this->_converter.getLatticeRelaxationTime() - 0.5;
    // speed up paramB
    T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);

    T omega_s;
    T omega = 1. / this->_converter.getLatticeRelaxationTime();

    this->_blockLattice.get(input).computeRhoU(rho, u);

    const T uSqr_s = util::normSqr<T,DESCRIPTOR::d>(u_s);
    T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      switch (mode) {
      case M2:
        omega_s = (lbDynamicsHelpers<T, DESCRIPTOR>::equilibrium(iPop, rho, u_s, uSqr_s)
                   - this->_blockLattice.get(input[0], input[1], input[2])[iPop])
                  + (1 - omega)
                  * (this->_blockLattice.get(input[0], input[1], input[2])[iPop]
                     - lbDynamicsHelpers<T, DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr));
        break;
      case M3:
        omega_s =
          (this->_blockLattice.get(input[0], input[1], input[2])[descriptors::opposite<DESCRIPTOR>(iPop)]
           - lbDynamicsHelpers<T, DESCRIPTOR>::equilibrium(
             descriptors::opposite<DESCRIPTOR>(iPop), rho, u_s, uSqr_s))
          - (this->_blockLattice.get(input[0], input[1], input[2])[iPop]
             - lbDynamicsHelpers<T, DESCRIPTOR>::equilibrium(iPop, rho, u_s, uSqr_s));
      }

      for (int i = 0; i < this->getTargetDim(); ++i) {
        output[i] -= descriptors::c<DESCRIPTOR>(iPop,i) * omega_s;
      }
    }

    for (int i = 0; i < this->getTargetDim(); ++i) {
      output[i] = this->_converter.getPhysForce(output[i] * paramB);
    }
  }
  return true;
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysWallShearStress3D<T,DESCRIPTOR>::BlockLatticePhysWallShearStress3D(
    BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
    BlockGeometryStructure3D<T>& blockGeometry,
    int overlap,
    int material,
    const UnitConverter<T,DESCRIPTOR>& converter,
    IndicatorF3D<T>& indicator)
    : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,1),
      _blockGeometry(blockGeometry),
      _overlap(overlap),
      _material(material)
      {
  this->getName() = "physWallShearStress";
  const T scaling = this->_converter.getConversionFactorLength() * 0.1;
  const T omega = 1. / this->_converter.getLatticeRelaxationTime();
  const T dt = this->_converter.getConversionFactorTime();
  _physFactor = -omega * descriptors::invCs2<T,DESCRIPTOR>() / dt * this->_converter.getPhysDensity() * this->_converter.getPhysViscosity();
  std::vector<int> discreteNormalOutwards(4, 0);

  for (int iX = 1 ; iX < _blockGeometry.getNx() - 1; iX++) {
    _discreteNormal.resize(_blockGeometry.getNx() - 2);
    _normal.resize(_blockGeometry.getNx() - 2);

    for (int iY = 1; iY < _blockGeometry.getNy() - 1; iY++) {
      _discreteNormal[iX-1].resize(_blockGeometry.getNy() - 2);
      _normal[iX-1].resize(_blockGeometry.getNy() - 2);

      for (int iZ = 1; iZ < _blockGeometry.getNz() - 1; iZ++) {
        _discreteNormal[iX-1][iY-1].resize(_blockGeometry.getNz() - 2);
        _normal[iX-1][iY-1].resize(_blockGeometry.getNz() - 2);

        if (_blockGeometry.get(iX, iY, iZ) == _material) {
          discreteNormalOutwards = _blockGeometry.getStatistics().getType(iX, iY, iZ);
          _discreteNormal[iX - 1][iY - 1][iZ - 1].resize(3);
          _normal[iX - 1][iY - 1][iZ - 1].resize(3);

          _discreteNormal[iX - 1][iY- 1][iZ- 1][0] = -discreteNormalOutwards[1];
          _discreteNormal[iX- 1][iY- 1][iZ- 1][1] = -discreteNormalOutwards[2];
          _discreteNormal[iX- 1][iY- 1][iZ- 1][2] = -discreteNormalOutwards[3];

          T physR[3];
          _blockGeometry.getPhysR(physR,iX, iY, iZ);
          Vector<T,3> origin(physR[0],physR[1],physR[2]);
          Vector<T,3> direction(-_discreteNormal[iX- 1][iY- 1][iZ- 1][0] * scaling,
                                -_discreteNormal[iX- 1][iY- 1][iZ- 1][1] * scaling,
                                -_discreteNormal[iX- 1][iY- 1][iZ- 1][2] * scaling);
          Vector<T,3> normal(0.,0.,0.);
          origin[0] = physR[0];
          origin[1] = physR[1];
          origin[2] = physR[2];

          indicator.normal(normal, origin, direction);
          normal.normalize();

          _normal[iX- 1][iY- 1][iZ- 1][0] = normal[0];
          _normal[iX- 1][iY- 1][iZ- 1][1] = normal[1];
          _normal[iX- 1][iY- 1][iZ- 1][2] = normal[2];
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysWallShearStress3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = T();

  if (input[0] + _overlap < 1 ||
      input[1] + _overlap < 1 ||
      input[2] + _overlap < 1 ||
      input[0] + _overlap >= _blockGeometry.getNx()-1 ||
      input[1] + _overlap >= _blockGeometry.getNy()-1 ||
      input[2] + _overlap >= _blockGeometry.getNz()-1 ){

#ifdef OLB_DEBUG
    std::cout << "Input address not mapped by _discreteNormal, overlap too small" << std::endl;
#endif
    return true;
  }

  if (_blockGeometry.get(input[0]+_overlap,input[1]+_overlap,input[2]+_overlap) == _material) {

    T traction[3];
    T stress[6];
    T rho = this->_blockLattice.get(input[0] + _overlap + _discreteNormal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][0],
                                    input[1] + _overlap + _discreteNormal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][1],
                                    input[2] + _overlap + _discreteNormal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][2]).computeRho();
     this->_blockLattice.get(input[0] + _overlap +   _discreteNormal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][0],
                      input[1] + _overlap +   _discreteNormal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][1],
                      input[2] + _overlap +   _discreteNormal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][2]).computeStress(stress);

    traction[0] = stress[0]*_physFactor/rho*_normal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][0] +
                  stress[1]*_physFactor/rho*_normal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][1] +
                  stress[2]*_physFactor/rho*_normal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][2];
    traction[1] = stress[1]*_physFactor/rho*_normal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][0] +
                  stress[3]*_physFactor/rho*_normal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][1] +
                  stress[4]*_physFactor/rho*_normal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][2];
    traction[2] = stress[2]*_physFactor/rho*_normal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][0] +
                  stress[4]*_physFactor/rho*_normal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][1] +
                  stress[5]*_physFactor/rho*_normal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][2];

    T traction_normal_SP;
    T tractionNormalComponent[3];
    // scalar product of traction and normal vector
    traction_normal_SP = traction[0] * _normal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][0] +
                         traction[1] * _normal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][1] +
                         traction[2] * _normal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][2];
    tractionNormalComponent[0] = traction_normal_SP * _normal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][0];
    tractionNormalComponent[1] = traction_normal_SP * _normal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][1];
    tractionNormalComponent[2] = traction_normal_SP * _normal[input[0]+_overlap-1][input[1]+_overlap-1][input[2]+_overlap-1][2];

    T WSS[3];
    WSS[0] = traction[0] - tractionNormalComponent[0];
    WSS[1] = traction[1] - tractionNormalComponent[1];
    WSS[2] = traction[2] - tractionNormalComponent[2];
    // magnitude of the wall shear stress vector
    output[0] = sqrt(WSS[0]*WSS[0] + WSS[1]*WSS[1] + WSS[2]*WSS[2]);

    return true;
  }
  else {
    return true;
  }
}


template <typename T, typename DESCRIPTOR>
BlockLatticePhysCorrBoundaryForce3D<T,DESCRIPTOR>::BlockLatticePhysCorrBoundaryForce3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
  BlockIndicatorF3D<T>&                  indicatorF,
  const UnitConverter<T,DESCRIPTOR>&     converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice, converter, 3),
    _indicatorF(indicatorF),
    _blockGeometry(indicatorF.getBlockGeometryStructure())
{
  this->getName() = "physCorrBoundaryForce";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysCorrBoundaryForce3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = T();
  }

  if (_indicatorF(input)) {
    for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
      // Get direction
      const Vector<int,3> c = descriptors::c<DESCRIPTOR>(iPop);
      // Get next cell located in the current direction
      // Check if the next cell is a fluid node
      if (_blockGeometry.get(input[0] + c[0], input[1] + c[1], input[2] + c[2]) == 1) {
        // Get f_q of next fluid cell where l = opposite(q)
        T f = this->_blockLattice.get(input[0] + c[0], input[1] + c[1], input[2] + c[2])[iPop];
        // Get f_l of the boundary cell
        // Add f_q and f_opp
        f += this->_blockLattice.get(input)[util::opposite<DESCRIPTOR>(iPop)];
        // Update force
        for (int i = 0; i < this->getTargetDim(); ++i) {
          output[i] -= c[i] * (f - 2. * descriptors::t<T,DESCRIPTOR>(iPop));
        }
      }
    }
    for (int i = 0; i < this->getTargetDim(); ++i) {
      output[i] = this->_converter.getPhysForce(output[i]);
    }
  }
  return true;
}

template<typename T, typename DESCRIPTOR, typename FIELD>
BlockLatticeField3D<T,DESCRIPTOR,FIELD>::BlockLatticeField3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, DESCRIPTOR::template size<FIELD>())
{
  this->getName() = "externalField";
}

template<typename T, typename DESCRIPTOR, typename FIELD>
bool BlockLatticeField3D<T,DESCRIPTOR,FIELD>::operator()(T output[], const int input[])
{
  this->_blockLattice.get(input[0], input[1], input[2]).template computeField<FIELD>(output);
  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticePorosity3D<T, DESCRIPTOR>::BlockLatticePorosity3D(BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1)
{
  this->getName() = "porosity";
}

// under construction
template<typename T, typename DESCRIPTOR>
bool BlockLatticePorosity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
#ifndef excludeDualDynamics
  this->_blockLattice.get(input[0], input[1], input[2]).template computeField<descriptors::POROSITY>(
    output);
#else
  this->_blockLattice.get(input[0], input[1], input[2]).template computeField<descriptors::POROSITY>(
    output);

#endif
  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticeVolumeFractionApproximation3D<T, DESCRIPTOR>::BlockLatticeVolumeFractionApproximation3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
    BlockGeometryStructure3D<T>& blockGeometry,
    IndicatorF3D<T>& indicator,
    int refinementLevel,
    const UnitConverter<T,DESCRIPTOR>& converter,
    bool insideOut)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _blockGeometry(blockGeometry), _indicator(indicator), _refinementLevel(refinementLevel), _converter(converter), _insideOut(insideOut),
    _physSubGridMinPhysRshift((1./ T(_refinementLevel * 2.) - 0.5) * _converter.getPhysDeltaX()),
    _physSubGridDeltaX(1. / T(_refinementLevel) * _converter.getPhysDeltaX()),
    _latticeSubGridVolume(1. / (_refinementLevel * _refinementLevel * _refinementLevel))
{
  this->getName() = "volumeFractionApproximation";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeVolumeFractionApproximation3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = 0.;
  T physR[3];
  bool inside[1];
  _blockGeometry.getPhysR(physR, input[0], input[1], input[2]);

  T subGridMinPhysR[3];
  subGridMinPhysR[0] = physR[0] + _physSubGridMinPhysRshift;
  subGridMinPhysR[1] = physR[1] + _physSubGridMinPhysRshift;
  subGridMinPhysR[2] = physR[2] + _physSubGridMinPhysRshift;

  for (int x = 0; x < _refinementLevel; x++) {
    for (int y = 0; y < _refinementLevel; y++) {
      for (int z = 0; z < _refinementLevel; z++) {
        physR[0] = subGridMinPhysR[0] + x * _physSubGridDeltaX;
        physR[1] = subGridMinPhysR[1] + y * _physSubGridDeltaX;
        physR[2] = subGridMinPhysR[2] + z * _physSubGridDeltaX;
        _indicator(inside, physR);
        if (!_insideOut) {
          if (inside[0]) {
            output[0] += _latticeSubGridVolume;
          }
        }
        else {
          if (!inside[0]) {
            output[0] += _latticeSubGridVolume;
          }
        }
      }
    }
  }
  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticePhysPermeability3D<T, DESCRIPTOR>::BlockLatticePhysPermeability3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 1)
{
  this->getName() = "permeability";
}

//TODO: consistency with 2D (181219)
//template<typename T, typename DESCRIPTOR>
//BlockLatticePhysPermeability3D<T, DESCRIPTOR>::BlockLatticePhysPermeability3D(
//  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
//  BlockGeometry3D<T>& blockGeometry, int material,
//  const UnitConverter<T>& converter)
//  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 1),
//    _blockGeometry(blockGeometry),
//    _material(material)
//{
//  this->getName() = "permeability";
//}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysPermeability3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T value;
  // get porosity
  this->_blockLattice.get(input[0], input[1], input[2]).template computeField<descriptors::POROSITY>(
    &value);
  // convert to physPermeability
//  output[0] = this->_converter.physPermeability(value);  // TODO converter MG
  // \todo Why do we need this???
  if (output[0] >= 42 && output[0] <= 42 && output[0] != 42) {
    output[0] = 999999;
  }
  if (std::isinf(output[0])) {
    output[0] = 1e6;
  }
  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticePhysCroppedPermeability3D<T, DESCRIPTOR>::BlockLatticePhysCroppedPermeability3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 1)
{
  this->getName() = "cropped_permeability";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysCroppedPermeability3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T value;
  // get porosity
  this->_blockLattice.get(input[0], input[1], input[2]).template computeField<descriptors::POROSITY>(
    &value);
  // convert to physPermeability
//  output[0] = _converter.physPermeability(value); // TODO converter MG
  // \todo Why do we need this???
  if (output[0] >= 42 && output[0] <= 42 && output[0] != 42) {
    output[0] = 1;
  }
  if (std::isinf(output[0])) {
    output[0] = 1;
  }
  if (output[0] > 1.) {
    output[0] = 1.;
  }
  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticePhysDarcyForce3D<T, DESCRIPTOR>::BlockLatticePhysDarcyForce3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, BlockGeometry3D<T>& blockGeometry,
  int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 3),
    _blockGeometry(blockGeometry),
    _material(material)
{
  this->getName() = "alphaU";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysDarcyForce3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  BlockLatticePhysPermeability3D<T, DESCRIPTOR> permeability(this->_blockLattice, this->_converter);
  BlockLatticeVelocity3D<T, DESCRIPTOR> velocity(this->_blockLattice);

  T nu = this->_converter.getPhysViscosity();
  permeability(output,input);
  T K = output[0];
  velocity(output,input);

  output[0] *= -nu / K;
  output[1] *= -nu / K;
  output[2] *= -nu / K;

  return true;
}


template<typename T, typename DESCRIPTOR>
BlockEuklidNorm3D<T, DESCRIPTOR>::BlockEuklidNorm3D(BlockF3D<T>& f)
  : BlockF3D<T>(f.getBlockStructure(), 1),
    _f(f)
{
  this->getName() = "EuklidNorm(" + f.getName() + ")";
}

template<typename T, typename DESCRIPTOR>
bool BlockEuklidNorm3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = T();
  T data[_f.getTargetDim()];
  _f(data,input);
  for (int i = 0; i < _f.getTargetDim(); ++i) {
    output[0] += data[i] * data[i];
  }
  output[0] = sqrt(output[0]);
  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticeInterpPhysVelocity3D<T, DESCRIPTOR>::BlockLatticeInterpPhysVelocity3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, UnitConverter<T,DESCRIPTOR> const& converter, Cuboid3D<T>* c, int overlap)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 3),
    _cuboid(c),
    _overlap(overlap)
{
  this->getName() = "BlockLatticeInterpVelocity3D";
}

template<typename T, typename DESCRIPTOR>
BlockLatticeInterpPhysVelocity3D<T, DESCRIPTOR>::BlockLatticeInterpPhysVelocity3D(
  const BlockLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& rhs) :
  BlockLatticePhysF3D<T, DESCRIPTOR>(rhs._blockLattice, rhs._converter, 3),
  _cuboid(rhs._cuboid)
{
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeInterpPhysVelocity3D<T, DESCRIPTOR>::operator()(T output[3], const T input[3])
{
  T u[3], rho, volume;
  T d[3], e[3];
  int latIntPos[3] = {0};
  T latPhysPos[3] = {T()};
  _cuboid->getFloorLatticeR(latIntPos, &input[0]);
  _cuboid->getPhysR(latPhysPos, latIntPos);

  T deltaRinv = 1. / _cuboid->getDeltaR();
  d[0] = (input[0] - latPhysPos[0]) * deltaRinv;
  d[1] = (input[1] - latPhysPos[1]) * deltaRinv;
  d[2] = (input[2] - latPhysPos[2]) * deltaRinv;

  e[0] = 1. - d[0];
  e[1] = 1. - d[1];
  e[2] = 1. - d[2];

  latIntPos[0]+=_overlap;
  latIntPos[1]+=_overlap;
  latIntPos[2]+=_overlap;

  this->_blockLattice.get(latIntPos[0], latIntPos[1],
                          latIntPos[2]).computeRhoU(rho, u);
  volume = e[0] * e[1] * e[2];
  output[0] = u[0] * volume;
  output[1] = u[1] * volume;
  output[2] = u[2] * volume;

  this->_blockLattice.get(latIntPos[0], latIntPos[1] + 1,
                          latIntPos[2]).computeRhoU(rho, u);
  volume = e[0] * d[1] * e[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0] + 1, latIntPos[1],
                          latIntPos[2]).computeRhoU(rho, u);
  volume = d[0] * e[1] * e[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0] + 1, latIntPos[1] + 1,
                          latIntPos[2]).computeRhoU(rho, u);
  volume = d[0] * d[1] * e[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0], latIntPos[1],
                          latIntPos[2] + 1).computeRhoU(rho,
                              u);
  volume = e[0] * e[1] * d[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0], latIntPos[1] + 1,
                          latIntPos[2] + 1).computeRhoU(rho,
                              u);
  volume = e[0] * d[1] * d[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0] + 1, latIntPos[1],
                          latIntPos[2] + 1).computeRhoU(rho,
                              u);
  volume = d[0] * e[1] * d[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  this->_blockLattice.get(latIntPos[0] + 1, latIntPos[1] + 1,
                          latIntPos[2] + 1).computeRhoU(rho,
                              u);
  volume = d[0] * d[1] * d[2];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;
  output[2] += u[2] * volume;

  output[0] = this->_converter.getPhysVelocity(output[0]);
  output[1] = this->_converter.getPhysVelocity(output[1]);
  output[2] = this->_converter.getPhysVelocity(output[2]);
}

template<typename T, typename DESCRIPTOR>
BlockLatticePorousMomentumLossForce3D<T, DESCRIPTOR>::BlockLatticePorousMomentumLossForce3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, BlockGeometryStructure3D<T>& blockGeometry,
  std::vector<SmoothIndicatorF3D<T,T,true>* >& indicator, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 7*indicator.size()), _blockGeometry(blockGeometry), _vectorOfIndicator(indicator)
{
  this->getName() = "physPorousMomentumLossForce";
}


template<typename T, typename DESCRIPTOR>
bool BlockLatticePorousMomentumLossForce3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  // iterate over all particles in _indicator
  for (size_t iInd=0; iInd!=_vectorOfIndicator.size(); iInd++) {

    int numVoxels = 0;
    int start[3] = {0}; 
    int end[3] = {0};   
    // check for intersection of cuboid and indicator
    Cuboid3D<T> tmpCuboid(_blockGeometry.getOrigin()[0], _blockGeometry.getOrigin()[1], _blockGeometry.getOrigin()[2], this->_converter.getPhysDeltaX(), _blockGeometry.getNx(), _blockGeometry.getNy(), _blockGeometry.getNz());
    T posXmin = _vectorOfIndicator[iInd]->getPos()[0] - _vectorOfIndicator[iInd]->getCircumRadius();
    T posXmax = _vectorOfIndicator[iInd]->getPos()[0] + _vectorOfIndicator[iInd]->getCircumRadius();
    T posYmin = _vectorOfIndicator[iInd]->getPos()[1] - _vectorOfIndicator[iInd]->getCircumRadius();
    T posYmax = _vectorOfIndicator[iInd]->getPos()[1] + _vectorOfIndicator[iInd]->getCircumRadius();
    T posZmin = _vectorOfIndicator[iInd]->getPos()[2] - _vectorOfIndicator[iInd]->getCircumRadius();
    T posZmax = _vectorOfIndicator[iInd]->getPos()[2] + _vectorOfIndicator[iInd]->getCircumRadius();
    if (tmpCuboid.checkInters(posXmin, posXmax, posYmin, posYmax, posZmin, posZmax, start[0], end[0], start[1], end[1], start[2], end[2])) {

      for(int k=0; k<3; k++) {
        start[k] -= 1;  
        if (start[k] < 0) start[k] = 0;
        end[k] += 2;    
        if (end[k] > _blockGeometry.getExtend()[k]) end[k] = _blockGeometry.getExtend()[k];
      }

      // iterate over cells in the constructed intersection box
      for (int iX = start[0]; iX < end[0]; iX++) {
        for (int iY = start[1]; iY < end[1]; iY++) {
          for (int iZ = start[2]; iZ < end[2]; iZ++) {
    
            // check if cell belongs to particle
            T inside[1] = {0.};
            T posIn[3] = {0.};
            _blockGeometry.getPhysR(posIn, iX, iY, iZ);
            (*(_vectorOfIndicator[iInd]))( inside, posIn);
            if ( !util::nearZero(inside[0]) && this->_blockGeometry.get(iX,iY,iZ)==1) {
              // compute momentum exchange force on particle
              T tmpForce[3] = {0.,0.,0.};
              tmpForce[0] += this->_blockLattice.get(iX, iY, iZ).template getFieldPointer<descriptors::VELOCITY_NUMERATOR>()[0];
              tmpForce[1] += this->_blockLattice.get(iX, iY, iZ).template getFieldPointer<descriptors::VELOCITY_NUMERATOR>()[1];
              tmpForce[2] += this->_blockLattice.get(iX, iY, iZ).template getFieldPointer<descriptors::VELOCITY_NUMERATOR>()[2];
              // reset external field for next timestep
              T reset_to_zero[3] = {0.,0.,0.};
              this->_blockLattice.get(iX, iY, iZ).template setField<descriptors::VELOCITY_NUMERATOR>(reset_to_zero);
              // convert force to SI units and compute torque
              numVoxels++;
              // division bei length of lattice cell necessary due to converter handling of force
              tmpForce[0] = this->_converter.getPhysForce(tmpForce[0]);
              tmpForce[1] = this->_converter.getPhysForce(tmpForce[1]);
              tmpForce[2] = this->_converter.getPhysForce(tmpForce[2]);
              output[0+iInd*7] += tmpForce[0];
              output[1+iInd*7] += tmpForce[1];
              output[2+iInd*7] += tmpForce[2];
              output[3+iInd*7] += (posIn[1]-_vectorOfIndicator[iInd]->getPos()[1])*tmpForce[2]
                                - (posIn[2]-_vectorOfIndicator[iInd]->getPos()[2])*tmpForce[1];
              output[4+iInd*7] += (posIn[2]-_vectorOfIndicator[iInd]->getPos()[2])*tmpForce[0]
                                - (posIn[0]-_vectorOfIndicator[iInd]->getPos()[0])*tmpForce[2];
              output[5+iInd*7] += (posIn[0]-_vectorOfIndicator[iInd]->getPos()[0])*tmpForce[1] 
                                - (posIn[1]-_vectorOfIndicator[iInd]->getPos()[1])*tmpForce[0];
            }
          }
        }
      }

    }
    output[6+iInd*7] = numVoxels;

  }
  return true;
}


template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
BlockLatticePhysTemperature3D<T,DESCRIPTOR,TDESCRIPTOR>::BlockLatticePhysTemperature3D
(BlockLatticeStructure3D<T,TDESCRIPTOR>& blockLattice, ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter)
  : BlockLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR>(blockLattice,converter,1)
{
  this->getName() = "physTemperature";
}


template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
bool BlockLatticePhysTemperature3D<T,DESCRIPTOR,TDESCRIPTOR>::operator() (T output[], const int input[])
{
  T latticeTemperature = this->_blockLattice.get( input[0], input[1], input[2]).computeRho();
  output[0] = this->_converter.getPhysTemperature(latticeTemperature);

  return true;
}


template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
BlockLatticePhysHeatFlux3D<T,DESCRIPTOR,TDESCRIPTOR>::BlockLatticePhysHeatFlux3D
(BlockLatticeStructure3D<T,TDESCRIPTOR>& blockLattice, ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter)
  : BlockLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR>(blockLattice,converter,3),
    _temp(converter.getLatticeSpecificHeatCapacity(converter.getPhysSpecificHeatCapacity())*(converter.getLatticeThermalRelaxationTime() - 0.5) / converter.getLatticeThermalRelaxationTime())
{
  this->getName() = "physHeatFlux";
}

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
bool BlockLatticePhysHeatFlux3D<T,DESCRIPTOR,TDESCRIPTOR>::operator() (T output[], const int input[])
{
  T temperature, extVel[3], j[3];
  this->_blockLattice.get( input[0], input[1], input[2] ).computeRhoU(temperature,extVel);
  this->_blockLattice.get( input[0], input[1], input[2] ).computeJ(j);
  output[0]= this->_converter.getPhysHeatFlux((j[0] - temperature * extVel[0])*_temp);
  output[1]= this->_converter.getPhysHeatFlux((j[1] - temperature * extVel[1])*_temp);
  output[2]= this->_converter.getPhysHeatFlux((j[2] - temperature * extVel[2])*_temp);

  return true;
}


template<typename T, typename DESCRIPTOR>
BlockLatticePhysBoundaryDistance3D<T, DESCRIPTOR>::BlockLatticePhysBoundaryDistance3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, BlockGeometryStructure3D<T>& blockGeometry, XMLreader const& xmlReader)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1), _blockGeometry(blockGeometry)
{
  this->getName() = "physBoundaryDistance";

  for (XMLreader* child : xmlReader) {
    // clout << "iterator to xml-child: " << child->getName() << std::endl;
    _tmpIndicator = createIndicatorSphere3D<T>(*child);
    _indicatorList.push_back(_tmpIndicator);
  }
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysBoundaryDistance3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T minDistance = std::numeric_limits<T>::max();
  T origin[3];
  _blockGeometry.getPhysR(origin, input);
  for (auto &indicator : _indicatorList) {
    bool inside[1] = {false};
    (*indicator)(inside, origin);
    if (inside[0]) {
      minDistance = -1;
      break;
    }
    T distance = 0;
    indicator->distance(distance, origin);
    // clout << "sphere distance = " << distance << endl;
    if ( distance < minDistance ) {
      minDistance = distance;
    }
  }
  // clout << "min distance = " << minDistance << endl;

  output[0] = minDistance;
  return true;
}


template<typename T, typename DESCRIPTOR>
BlockLatticePhysPoreSizeDistribution3D<T, DESCRIPTOR>::BlockLatticePhysPoreSizeDistribution3D(
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice, BlockGeometryStructure3D<T>& blockGeometry, int material, XMLreader const& xmlReader)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1), _blockGeometry(blockGeometry), _material(material),
    _distanceFunctor(blockLattice, blockGeometry, xmlReader), _distanceCache(_distanceFunctor)
{
  this->getName() = "physPoreSizeDistribution";

  for (XMLreader* child : xmlReader) {
    // clout << "iterator to xml-child: " << child->getName() << std::endl;
    _tmpIndicator = createIndicatorSphere3D<T>(*child);
    _indicatorList.push_back(_tmpIndicator);
  }

}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysPoreSizeDistribution3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (this->_blockGeometry.get(input[0],input[1],input[2]) == _material) {
    T localDistance[1] = {0.};
    _distanceCache(localDistance, input);
    // cout << localDistance[0] << endl;

    // filter by local maximum (compare to 26 neighbours)
    for (int iPop = 1; iPop < 27; iPop++) {
      int neighbourInput[3] = {0,0,0};
      for (int iDim = 0; iDim < 3; iDim++) {
        neighbourInput[iDim] = input[iDim] + descriptors::c<descriptors::D3Q27<>>(iPop,iDim);
        // cout << neighbourInput[iDim] << " ";
      }
      // cout << iPop << ": ";

      T neighbourDistance[1] = {0.};
      if ( _distanceCache(neighbourDistance, neighbourInput) ) {
        if ( neighbourDistance[0] > localDistance[0] ) {
          // cout << "neighbour larger ";
          // cout << iPop << std::endl;
          output[0] = -1;
          return true;
        }
      }
      else {
        // cout << "distance not calculated ";
        // cout << iPop << std::endl;
        output[0] = -1;
        return true;
      }
    }

    output[0] = localDistance[0];
    return true;
  }
  else {
    output[0] = -1;
    return false;
  }
}


template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
BlockLatticePhysTauFromBoundaryDistance3D<T,DESCRIPTOR,TDESCRIPTOR>::BlockLatticePhysTauFromBoundaryDistance3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockGeometryStructure3D<T>& blockGeometry, XMLreader const& xmlReader, ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter, const T p, const T T_avg, const T c_p, const T beta, const T lambda_0, const T sigma, const T p_0, const T n_0)
  : BlockLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR>(blockLattice,converter,1), _blockGeometry(blockGeometry), _distanceFunctor(blockLattice, blockGeometry, xmlReader), _tmp1(lambda_0 * 287.058 * T_avg / p / c_p), _tmp2(2. * beta / (sqrt(2.) * M_PI * sigma * sigma * p * n_0 / p_0))
{
  this->getName() = "physTauFromBoundaryDistance";
}


template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
bool BlockLatticePhysTauFromBoundaryDistance3D<T,DESCRIPTOR,TDESCRIPTOR>::operator() (T output[], const int input[])
{
  T L[1] = {0.};
  _distanceFunctor(L, input);
  if ( L[0] < this->_converter.getPhysDeltaX() ) {
    L[0] = this->_converter.getPhysDeltaX();
  }

  const T alpha = _tmp1 / ( 1. + _tmp2 / L[0] );

  output[0] = alpha / this->_converter.getConversionFactorViscosity() * descriptors::invCs2<T,TDESCRIPTOR>() + 0.5;

  // std::cout << L[0] << " " << alpha << " " << output[0] << std::endl;

  return true;
}


template<typename T, typename DESCRIPTOR, bool HLBM>
BlockLatticeIndicatorSmoothIndicatorIntersection3D<T, DESCRIPTOR, HLBM>::BlockLatticeIndicatorSmoothIndicatorIntersection3D (
  BlockLatticeStructure3D<T, DESCRIPTOR>& blockLattice,
  BlockGeometryStructure3D<T>& blockGeometry,
  IndicatorF3D<T>& normalInd,
  SmoothIndicatorF3D<T,T,HLBM>& smoothInd )
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1),
    _blockGeometry(blockGeometry), _normalInd(normalInd), _smoothInd(smoothInd)
{
  this->getName() = "Indicator-SmoothIndicator Intersection";
}

template<typename T, typename DESCRIPTOR, bool HLBM>
bool BlockLatticeIndicatorSmoothIndicatorIntersection3D<T, DESCRIPTOR, HLBM>::operator()(T output[], const int input[])
{
  output[0] = 0.;
  int start[3] = {0};
  int end[3] = {0};
  // check for intersection of cuboid and smoothIndicator
  Cuboid3D<T> tmpCuboid(_blockGeometry.getOrigin()[0], _blockGeometry.getOrigin()[1], _blockGeometry.getOrigin()[2], _blockGeometry.getDeltaR(), _blockGeometry.getNx(), _blockGeometry.getNy(), _blockGeometry.getNz());
    T posXmin = _smoothInd.getPos()[0] - _smoothInd.getCircumRadius();
    T posXmax = _smoothInd.getPos()[0] + _smoothInd.getCircumRadius();
    T posYmin = _smoothInd.getPos()[1] - _smoothInd.getCircumRadius();
    T posYmax = _smoothInd.getPos()[1] + _smoothInd.getCircumRadius();
    T posZmin = _smoothInd.getPos()[2] - _smoothInd.getCircumRadius();
    T posZmax = _smoothInd.getPos()[2] + _smoothInd.getCircumRadius();
    if (tmpCuboid.checkInters(posXmin, posXmax, posYmin, posYmax, posZmin, posZmax, start[0], 
                              end[0], start[1], end[1], start[2], end[2])) {
  
      for(int k=0; k<3; k++) {
        start[k] -= 1;
        if (start[k] < 0) start[k] = 0;
        end[k] += 2;
        if (end[k] > _blockGeometry.getExtend()[k]) end[k] = _blockGeometry.getExtend()[k];
      }

    // iterate over cells in the constructed intersection box
    for (int iX = start[0]; iX < end[0]; iX++) {
      for (int iY = start[1]; iY < end[1]; iY++) {
        for (int iZ = start[2]; iZ < end[2]; iZ++) {

          // check if cell belongs to particle
          T insideT[1] = {0.};
          T posIn[3] = {0.};
          _blockGeometry.getPhysR(posIn, iX, iY, iZ);
          _smoothInd( insideT, posIn);
          if ( !util::nearZero(insideT[0]) && this->_blockGeometry.get(iX,iY,iZ)==1) {
            // Return true if at least one cell is found to be inside both A and B
            bool insideBool[1] = {false};
            _normalInd(insideBool, posIn);
            if (insideBool[0]) {
              output[0] = 1.;
              return true;
            }
          }
        }
      }
    }
  }

  return true;
}

/////////// BlockLatticeGuoZhaoEpsilon3D /////////////////////////////////////////////

template <typename T, typename DESCRIPTOR>
BlockLatticeGuoZhaoEpsilon3D<T,DESCRIPTOR>::BlockLatticeGuoZhaoEpsilon3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 1)
{
  this->getName() = "epsilon";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeGuoZhaoEpsilon3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = *this->_blockLattice.get( input[0], input[1], input[2] ).template getFieldPointer<descriptors::EPSILON>();
  return true;
}

/////////// BlockLatticeGuoZhaoPhysK3D /////////////////////////////////////////////

template <typename T, typename DESCRIPTOR>
BlockLatticeGuoZhaoPhysK3D<T,DESCRIPTOR>::BlockLatticeGuoZhaoPhysK3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice, converter, 1)
{
  this->getName() = "physK";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeGuoZhaoPhysK3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = *this->_blockLattice.get( input[0], input[1], input[2] )
              .template getFieldPointer<descriptors::K>()
              * this->_converter.getConversionFactorLength()
              * this->_converter.getConversionFactorLength();
  return true;
}

/////////// BlockLatticeGuoZhaoPhysBodyForce3D /////////////////////////////////////////////

template <typename T, typename DESCRIPTOR>
BlockLatticeGuoZhaoPhysBodyForce3D<T,DESCRIPTOR>::BlockLatticeGuoZhaoPhysBodyForce3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice, converter, 3)
{
  this->getName() = "physBodyForce";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeGuoZhaoPhysBodyForce3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = this->_converter.getPhysForce(
                *this->_blockLattice.get( input[0], input[1], input[2] )
                .template getFieldPointer<descriptors::BODY_FORCE>() );
  output[1] = this->_converter.getPhysForce(
                *this->_blockLattice.get( input[0], input[1], input[2] )
                .template getFieldPointer<descriptors::BODY_FORCE>(1) );
  output[2] = this->_converter.getPhysForce(
                *this->_blockLattice.get( input[0], input[1], input[2] )
                .template getFieldPointer<descriptors::BODY_FORCE>(2) );
  return true;
}

/////////// BlockLatticeTimeStepScale3D /////////////////////////////////////////////

template <typename T, typename DESCRIPTOR>
BlockLatticeTimeStepScale3D<T,DESCRIPTOR>::BlockLatticeTimeStepScale3D
(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, T oldTau, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, DESCRIPTOR::q), _tau_old(oldTau), _converter(converter)
{
  this->getName() = "latticeTimeStepScale";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeTimeStepScale3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{

  Cell<T,DESCRIPTOR>& cell = this->_blockLattice.get(input[0], input[1], input[2]);
  Dynamics<T,DESCRIPTOR>* dynamics = this->_blockLattice.getDynamics(input[0], input[1], input[2]);
  T rho_old, rho_new, u_old[DESCRIPTOR::d], u_new[DESCRIPTOR::d];;
  cell.computeRhoU(rho_old,u_old);
  T fNeq[DESCRIPTOR::q];

  T uSqr_old = util::normSqr<T,DESCRIPTOR::d>(u_old);

  T tau_new = _converter.getLatticeRelaxationTime();

  T physDeltaT_new = _converter.getPhysDeltaT();
  T physDeltaT_old = physDeltaT_new * (_tau_old-0.5) / (tau_new-0.5);
  T neqScaling = ((physDeltaT_new*tau_new)/(physDeltaT_old*_tau_old));

  for (int iDim=0; iDim < DESCRIPTOR::d; iDim++) {
    u_new[iDim] = u_old[iDim] * physDeltaT_new / physDeltaT_old ;
  }

  rho_new = (rho_old-1.0) * (physDeltaT_new / physDeltaT_old) * (physDeltaT_new / physDeltaT_old) + 1.0;

  T uSqr_new = util::normSqr<T,DESCRIPTOR::d>(u_new);


  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    fNeq[iPop] = cell[iPop] - dynamics -> computeEquilibrium(iPop,rho_old,u_old,uSqr_old);
//    output[iPop] = cell[iPop];
    output[iPop] = (dynamics -> computeEquilibrium(iPop,rho_old,u_old,uSqr_old)+descriptors::t<T,DESCRIPTOR>(iPop)// * ((dynamics -> computeEquilibrium(iPop,rho_old,u_new,uSqr_new)+descriptors::t<T,DESCRIPTOR>(iPop))/(dynamics -> computeEquilibrium(iPop,rho_old,u_old,uSqr_old)+descriptors::t<T,DESCRIPTOR>(iPop)))
        + fNeq[iPop]*neqScaling)
        * ((dynamics -> computeEquilibrium(iPop,rho_new,u_new,uSqr_new)+descriptors::t<T,DESCRIPTOR>(iPop))/(dynamics -> computeEquilibrium(iPop,rho_old,u_old,uSqr_old)+descriptors::t<T,DESCRIPTOR>(iPop))) - descriptors::t<T,DESCRIPTOR>(iPop);
  }
  return true;
}



}  // end namespace olb

#endif
