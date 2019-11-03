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

#ifndef BLOCK_LATTICE_LOCAL_F_2D_HH
#define BLOCK_LATTICE_LOCAL_F_2D_HH

#include<vector>
#include<cmath>

#include "blockLatticeLocalF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "core/blockLattice2D.h"
#include "communication/mpiManager.h"
#include "core/blockLatticeStructure2D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity

namespace olb {


template <typename T, typename DESCRIPTOR>
BlockLatticeDissipation2D<T,DESCRIPTOR>::BlockLatticeDissipation2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1), _converter(converter)
{
  this->getName() = "dissipation";
}

// todo: get functor working
template <typename T, typename DESCRIPTOR>
bool BlockLatticeDissipation2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  //  int lociz= input[3];
  //  if ( this->_blockLattice.get_load().rank(globIC) == singleton::mpi().getRank() ) {
  //    // local coords are given, fetch local cell and compute value(s)
  //    T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  //    int overlap = this->_blockLattice.getOverlap();
  //    this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap).computeAllMomenta(rho, uTemp, pi);

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

template <typename T, typename DESCRIPTOR>
BlockLatticePhysDissipation2D<T,DESCRIPTOR>::BlockLatticePhysDissipation2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1)
{
  this->getName() = "physDissipation";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysDissipation2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_blockLattice.get( input[0], input[1] ).computeAllMomenta(rho, uTemp, pi);

  T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] + pi[5]*pi[5];
  }

  T nuLattice = this->_converter.getLatticeViscosity();
  T omega = 1. / this->_converter.getLatticeRelaxationTime();
  T dt = this->_converter.getConversionFactorTime();
  output[0] = PiNeqNormSqr*nuLattice*pow(omega*descriptors::invCs2<T,DESCRIPTOR>()/rho,2)/2.*this->_converter.getPhysViscosity()/this->_converter.getLatticeViscosity()/dt/dt;

  return true;
}

template <typename T, typename DESCRIPTOR>
BlockLatticeDensity2D<T,DESCRIPTOR>::BlockLatticeDensity2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1)
{
  this->getName() = "density";
}


template <typename T, typename DESCRIPTOR>
bool BlockLatticeDensity2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = this->_blockLattice.get( input[0], input[1] ).computeRho();
  return true;
}


template <typename T, typename DESCRIPTOR>
BlockLatticeVelocity2D<T,DESCRIPTOR>::BlockLatticeVelocity2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,2)
{
  this->getName() = "velocity";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeVelocity2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T rho;
  this->_blockLattice.get(input[0], input[1]).computeRhoU(rho,output);
  return true;
}

/*template <typename T, typename DESCRIPTOR>
BlockLatticeStrainRate2D<T,DESCRIPTOR>::BlockLatticeStrainRate2D
  (BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,4)
{ this->getName() = "strainRate"; }

template <typename T, typename DESCRIPTOR>
std::vector<T> BlockLatticeStrainRate2D<T,DESCRIPTOR>::operator() (std::vector<int> input) {

  T strainRate[4];
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_blockLattice.get( input[0], input[1] ).computeAllMomenta(rho, uTemp, pi);

  T omega = this->_converter.getOmega();

  strainRate[0] = -pi[0]*omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2.;
  strainRate[1] = -pi[1]*omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2.;
  strainRate[2] = -pi[1]*omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2.;
  strainRate[3] = -pi[2]*omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2.;

  //cout << pi[0] << " " << pi[1] << " " << pi[2] << " " << descriptors::invCs2<T,DESCRIPTOR>() << endl;

  std::vector<T> output(strainRate, strainRate+4); // first adress, last adress
  return output;
}*/

template <typename T, typename DESCRIPTOR>
BlockLatticePhysStrainRate2D<T,DESCRIPTOR>::BlockLatticePhysStrainRate2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,4)
{
  this->getName() = "strainRate";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysStrainRate2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_blockLattice.get( input[0], input[1] ).computeAllMomenta(rho, uTemp, pi);

  T omega = 1. / this->_converter.getLatticeRelaxationTime();
  T dt = this->_converter.getConversionFactorTime();

  output[0] = -pi[0]*omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2./dt;
  output[1] = -pi[1]*omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2./dt;
  output[2] = -pi[1]*omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2./dt;
  output[3] = -pi[2]*omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2./dt;

  return true;
}

template <typename T, typename DESCRIPTOR>
BlockLatticeGeometry2D<T,DESCRIPTOR>::BlockLatticeGeometry2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, BlockGeometryStructure2D<T>& blockGeometry, int material)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "geometry";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeGeometry2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  int materialTmp = _blockGeometry.getMaterial( input[0], input[1] );

  if (_material != -1) {
    if (_material == materialTmp) {
      output[0] = T(1);
      return true;
    }
    else {
      output[0] = T();
      return true;
    }
  }
  output[0]=T(materialTmp);
  return false;
}


template <typename T, typename DESCRIPTOR>
BlockLatticeRank2D<T,DESCRIPTOR>::BlockLatticeRank2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1)
{
  this->getName() = "rank";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeRank2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = singleton::mpi().getRank() + 1;
  return false;
}


template <typename T, typename DESCRIPTOR>
BlockLatticeCuboid2D<T,DESCRIPTOR>::BlockLatticeCuboid2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const int iC)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1), _iC(iC)
{
  this->getName() = "cuboid";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeCuboid2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = _iC + 1;
  return false;
}


template <typename T, typename DESCRIPTOR>
BlockLatticePhysPressure2D<T,DESCRIPTOR>::BlockLatticePhysPressure2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1)
{
  this->getName() = "physPressure";
}


template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysPressure2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  // lattice pressure = c_s^2 ( rho -1 )
  T latticePressure = ( this->_blockLattice.get( input[0], input[1] ).computeRho() - 1.0 ) / descriptors::invCs2<T,DESCRIPTOR>();
  output[0] = this->_converter.getPhysPressure(latticePressure);

  return true;
}


template <typename T, typename DESCRIPTOR>
BlockLatticePhysVelocity2D<T,DESCRIPTOR>::BlockLatticePhysVelocity2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2)
{
  this->getName() = "physVelocity";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysVelocity2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T rho;
  this->_blockLattice.get( input[0], input[1] ).computeRhoU(rho,output);
  output[0]=this->_converter.getPhysVelocity(output[0]);
  output[1]=this->_converter.getPhysVelocity(output[1]);
  return true;
}


template<typename T, typename DESCRIPTOR, typename FIELD>
BlockLatticeField2D<T,DESCRIPTOR,FIELD>::BlockLatticeField2D(
  BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF2D<T, DESCRIPTOR>(blockLattice, DESCRIPTOR::template size<FIELD>())
{
  this->getName() = "extField";
}

template<typename T, typename DESCRIPTOR, typename FIELD>
bool BlockLatticeField2D<T,DESCRIPTOR,FIELD>::operator()(
  T output[], const int input[])
{
  this->_blockLattice.get(input[0], input[1]).template computeField<FIELD>(output);
  return true;
}


template <typename T, typename DESCRIPTOR>
BlockLatticePhysExternalPorosity2D<T,DESCRIPTOR>::BlockLatticePhysExternalPorosity2D(
  BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
  int overlap,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice, converter, 2),
    _overlap(overlap)
{
  this->getName() = "ExtPorosityField";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysExternalPorosity2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  this->_blockLattice.get(
    input[0]+_overlap, input[1]+_overlap
  ).template computeField<descriptors::POROSITY>(output);
  return true;
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysExternalVelocity2D<T,DESCRIPTOR>::BlockLatticePhysExternalVelocity2D(
  BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
  int overlap,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice, converter, 2),
    _overlap(overlap)
{
  this->getName() = "ExtVelocityField";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysExternalVelocity2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  this->_blockLattice.get(
    input[0]+_overlap, input[1]+_overlap
  ).template computeField<descriptors::VELOCITY>(output);
  output[0] = this->_converter.getPhysVelocity(output[0]);
  output[1] = this->_converter.getPhysVelocity(output[1]);
  return true;
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysExternalParticleVelocity2D<T,DESCRIPTOR>::BlockLatticePhysExternalParticleVelocity2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2)
{
  this->getName() = "ExtParticleVelocityField";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysExternalParticleVelocity2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  const T* velocity_numerator   = this->blockLattice.get(input).template getFieldPointer<descriptors::VELOCITY_NUMERATOR>();
  const T* velocity_denominator = this->blockLattice.get(input).template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>();

  if (velocity_denominator[0] > std::numeric_limits<T>::epsilon()) {
    output[0]=this->_converter.getPhysVelocity(velocity_numerator[0]/velocity_denominator[0]);
    output[1]=this->_converter.getPhysVelocity(velocity_numerator[1]/velocity_denominator[0]);
    return true;
  }
  output[0]=this->_converter.getPhysVelocity(velocity_numerator[0]);
  output[1]=this->_converter.getPhysVelocity(velocity_numerator[1]);
  return true;
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysWallShearStress2D<T,DESCRIPTOR>::BlockLatticePhysWallShearStress2D(
    BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
    BlockGeometryStructure2D<T>& blockGeometry,
    int overlap,
    int material,
    const UnitConverter<T,DESCRIPTOR>& converter,
    IndicatorF2D<T>& indicator)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1),
    _blockGeometry(blockGeometry), _overlap(overlap), _material(material)
{
  this->getName() = "physWallShearStress";
  const T scaling = this->_converter.getConversionFactorLength() * 0.1;
  const T omega = 1. / this->_converter.getLatticeRelaxationTime();
  const T dt = this->_converter.getConversionFactorTime();
  _physFactor = -omega * descriptors::invCs2<T,DESCRIPTOR>() / dt * this->_converter.getPhysDensity() * this->_converter.getPhysViscosity();
  std::vector<int> discreteNormalOutwards(3, 0);

  for (int iX = 1; iX < _blockGeometry.getNx() - 1; iX++) {
    _discreteNormal.resize(_blockGeometry.getNx() - 2);
    _normal.resize(_blockGeometry.getNx() - 2);

    for (int iY = 1; iY < _blockGeometry.getNy() - 1; iY++) {
      _discreteNormal[iX-1].resize(_blockGeometry.getNy() - 2);
      _normal[iX-1].resize(_blockGeometry.getNy() - 2);

      if (_blockGeometry.get(iX, iY) == _material) {
        discreteNormalOutwards = _blockGeometry.getStatistics().getType(iX, iY);
        _discreteNormal[iX-1][iY-1].resize(2);
        _normal[iX-1][iY-1].resize(2);

        _discreteNormal[iX-1][iY-1][0] = -discreteNormalOutwards[1];
        _discreteNormal[iX-1][iY-1][1] = -discreteNormalOutwards[2];

        T physR[2];
        _blockGeometry.getPhysR(physR, iX, iY);
        Vector<T,2> origin(physR[0],physR[1]);
        Vector<T,2> direction(-_discreteNormal[iX-1][iY-1][0] * scaling,
                              -_discreteNormal[iX-1][iY-1][1] * scaling);
        Vector<T,2> normal(0.,0.);
        origin[0] = physR[0];
        origin[1] = physR[1];

        indicator.normal(normal, origin, direction);
        normal.normalize();

        _normal[iX-1][iY-1][0] = normal[0];
        _normal[iX-1][iY-1][1] = normal[0];
      }
    }
  }
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysWallShearStress2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = T();

  if (input[0] + _overlap < 1 ||
      input[1] + _overlap < 1 ||
      input[0] + _overlap >= _blockGeometry.getNx()-1 ||
      input[1] + _overlap >= _blockGeometry.getNy()-1){
#ifdef OLB_DEBUG
    std::cout << "Input address not mapped by _discreteNormal, overlap too small" << std::endl;
#endif
    return true;
  }

  if (this->_blockGeometry.get(input[0]+_overlap,input[1]+_overlap)==_material) {

    T traction[2];
    T stress[3];
    T rho = this->_blockLattice.get(input[0] + _overlap + _discreteNormal[input[0]+_overlap-1][input[1]+_overlap-1][0],
                                    input[1] + _overlap + _discreteNormal[input[0]+_overlap-1][input[1]+_overlap-1][1]).computeRho();

    this->_blockLattice.get(input[0] + _overlap +   _discreteNormal[input[0]+_overlap-1][input[1]+_overlap-1][0],
                            input[1] + _overlap +   _discreteNormal[input[0]+_overlap-1][input[1]+_overlap-1][1]).computeStress(stress);

    traction[0] = stress[0]*_physFactor/rho*_normal[input[0]+_overlap-1][input[1]+_overlap-1][0] +
                  stress[1]*_physFactor/rho*_normal[input[0]+_overlap-1][input[1]+_overlap-1][1];
    traction[1] = stress[1]*_physFactor/rho*_normal[input[0]+_overlap-1][input[1]+_overlap-1][0] +
                  stress[2]*_physFactor/rho*_normal[input[0]+_overlap-1][input[1]+_overlap-1][1];

    T traction_normal_SP;
    T tractionNormalComponent[2];
    // scalar product of traction and normal vector
    traction_normal_SP = traction[0] * _normal[input[0]+_overlap-1][input[1]+_overlap-1][0] +
                           traction[1] * _normal[input[0]+_overlap-1][input[1]+_overlap-1][1];
    tractionNormalComponent[0] = traction_normal_SP * _normal[input[0]+_overlap-1][input[1]+_overlap-1][0];
    tractionNormalComponent[1] = traction_normal_SP * _normal[input[0]+_overlap-1][input[1]+_overlap-1][1];

    T WSS[2];
    WSS[0] = traction[0] - tractionNormalComponent[0];
    WSS[1] = traction[1] - tractionNormalComponent[1];
    // magnitude of the wall shear stress vector
    output[0] = sqrt(WSS[0]*WSS[0] + WSS[1]*WSS[1]);
    return true;
  }
  else {
    return true;
  }
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysBoundaryForce2D<T,DESCRIPTOR>::BlockLatticePhysBoundaryForce2D(
  BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
  BlockIndicatorF2D<T>&                  indicatorF,
  const UnitConverter<T,DESCRIPTOR>&     converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice, converter, 2),
    _indicatorF(indicatorF),
    _blockGeometry(indicatorF.getBlockGeometryStructure())
{
  this->getName() = "physBoundaryForce";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysBoundaryForce2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = T();
  }

  if (_indicatorF(input)) {
    for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
      // Get direction
      // Get next cell located in the current direction
      // Check if the next cell is a fluid node
      if (_blockGeometry.get(input[0] + descriptors::c<DESCRIPTOR >(iPop,0), input[1] + descriptors::c<DESCRIPTOR >(iPop,1)) == 1) {
        // Get f_q of next fluid cell where l = opposite(q)
        T f = this->_blockLattice.get(input[0] + descriptors::c<DESCRIPTOR >(iPop,0), input[1] + descriptors::c<DESCRIPTOR >(iPop,1))[iPop];
        // Get f_l of the boundary cell
        // Add f_q and f_opp
        f += this->_blockLattice.get(input[0], input[1])[util::opposite<DESCRIPTOR >(iPop)];
        // Update force
        for (int i = 0; i < this->getTargetDim(); ++i) {
          output[i] -= descriptors::c<DESCRIPTOR >(iPop,i) * f;
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
BlockLatticePSMPhysForce2D<T,DESCRIPTOR>::BlockLatticePSMPhysForce2D(
  BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
  const UnitConverter<T,DESCRIPTOR>&     converter,
  int mode_)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice, converter, 2)
{
  this->getName() = "physPSMForce";
  mode = (Mode) mode_;
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePSMPhysForce2D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = T();
  }

  T epsilon = 1. - *(this->_blockLattice.get(input[0], input[1]).template getFieldPointer<descriptors::POROSITY>());

  //if ((epsilon > 1e-5 && epsilon < 1 - 1e-5)) {
  if ((epsilon > 1e-5)) {
    T rho, u[DESCRIPTOR::d], u_s[DESCRIPTOR::d];

    for (int i = 0; i < DESCRIPTOR::d; i++) {
      u_s[i] = (this->_blockLattice.get(input[0], input[1]).template getFieldPointer<descriptors::VELOCITY_SOLID>())[i];
    }
    T paramA = this->_converter.getLatticeRelaxationTime() - 0.5;
    // speed up paramB
    T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);

    T omega_s;
    T omega = 1. / this->_converter.getLatticeRelaxationTime();

    this->_blockLattice.get(input[0], input[1]).computeRhoU(rho, u);

    const T uSqr_s = util::normSqr<T,DESCRIPTOR::d>(u_s);
    T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      switch (mode) {
        case M2:
          omega_s = (lbDynamicsHelpers<T, DESCRIPTOR>::equilibrium(iPop, rho, u_s, uSqr_s)
              - this->_blockLattice.get(input[0], input[1])[iPop])
              + (1 - omega)
                  * (this->_blockLattice.get(input[0], input[1])[iPop]
                      - lbDynamicsHelpers<T, DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr));
          break;
        case M3:
          omega_s =
              (this->_blockLattice.get(input[0], input[1])[descriptors::opposite<DESCRIPTOR>(iPop)]
                  - lbDynamicsHelpers<T, DESCRIPTOR>::equilibrium(
                      descriptors::opposite<DESCRIPTOR>(iPop), rho, u_s, uSqr_s))
                  - (this->_blockLattice.get(input[0], input[1])[iPop]
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
BlockLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR>::BlockLatticePhysCorrBoundaryForce2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,BlockGeometry2D<T>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2),
    _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "physCorrBoundaryForce";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysCorrBoundaryForce2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  //  int lociz= input[3];

  //  std::vector<T> force(3, T());
  //  if ( this->_blockLattice.get_load().rank(globIC) == singleton::mpi().getRank() ) {
  //    int globX = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosX() + locix;
  //    int globY = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosY() + lociy;
  //    int globZ = (int)this->_blockLattice.get_cGeometry().get(globIC).get_globPosZ() + lociz;
  //    if(BlockGeometry.getMaterial(globX, globY, globZ) == material) {
  //      for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
  //        // Get direction
  //        const int* c = DESCRIPTOR::c[iPop];
  //        // Get next cell located in the current direction
  //        // Check if the next cell is a fluid node
  //        if (BlockGeometry.getMaterial(globX + c[0], globY + c[1], globZ + c[2]) == 1) {
  //          int overlap = this->_blockLattice.getOverlap();
  //          // Get f_q of next fluid cell where l = opposite(q)
  //          T f = this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap + c[0], lociy+overlap + c[1], lociz+overlap + c[2])[iPop];
  //          // Get f_l of the boundary cell
  //          // Add f_q and f_opp
  //          f += this->_blockLattice.getExtendedBlockLattice(this->_blockLattice.get_load().loc(globIC)).get(locix+overlap, lociy+overlap, lociz+overlap)[util::opposite<DESCRIPTOR >(iPop)];
  //          // Update force
  //          force[0] -= c[0]*(f-2.*descriptors::t<T,DESCRIPTOR>(iPop));
  //          force[1] -= c[1]*(f-2.*descriptors::t<T,DESCRIPTOR>(iPop));
  //          force[2] -= c[2]*(f-2.*descriptors::t<T,DESCRIPTOR>(iPop));
  //        }
  //      }
  //      force[0]=this->_converter.physForce(force[0]);
  //      force[1]=this->_converter.physForce(force[1]);
  //      force[2]=this->_converter.physForce(force[2]);
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


template <typename T, typename DESCRIPTOR>
BlockLatticePorosity2D<T,DESCRIPTOR>::BlockLatticePorosity2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, BlockGeometryStructure2D<T>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "porosity";
}

// under construction
template <typename T, typename DESCRIPTOR>
bool BlockLatticePorosity2D<T,DESCRIPTOR>::operator()(T output[], const int input[])
{
  return false;
}

template<typename T, typename DESCRIPTOR>
BlockLatticeVolumeFractionApproximation2D<T, DESCRIPTOR>::BlockLatticeVolumeFractionApproximation2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
    BlockGeometryStructure2D<T>& blockGeometry,
    IndicatorF2D<T>& indicator,
    int refinementLevel,
    const UnitConverter<T,DESCRIPTOR>& converter,
    bool insideOut)
  : BlockLatticeF2D<T, DESCRIPTOR>(blockLattice, 1),
    _blockGeometry(blockGeometry), _indicator(indicator), _refinementLevel(refinementLevel), _converter(converter), _insideOut(insideOut),
    _physSubGridMinPhysRshift((1./ T(_refinementLevel * 2.) - 0.5) * _converter.getPhysDeltaX()),
    _physSubGridDeltaX(1. / T(_refinementLevel) * _converter.getPhysDeltaX()),
    _latticeSubGridVolume(1. / (_refinementLevel * _refinementLevel))
{
  this->getName() = "volumeFractionApproximation";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeVolumeFractionApproximation2D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = 0.;
  T physR[2];
  bool inside[1];
  _blockGeometry.getPhysR(physR, input[0], input[1]);

  T subGridMinPhysR[2];
  subGridMinPhysR[0] = physR[0] + _physSubGridMinPhysRshift;
  subGridMinPhysR[1] = physR[1] + _physSubGridMinPhysRshift;

  for (int x = 0; x < _refinementLevel; x++) {
    for (int y = 0; y < _refinementLevel; y++) {
      physR[0] = subGridMinPhysR[0] + x * _physSubGridDeltaX;
      physR[1] = subGridMinPhysR[1] + y * _physSubGridDeltaX;
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
  return true;
}



template<typename T, typename DESCRIPTOR>
BlockLatticeVolumeFractionPolygonApproximation2D<T, DESCRIPTOR>::BlockLatticeVolumeFractionPolygonApproximation2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
    BlockGeometryStructure2D<T>& blockGeometry,
    IndicatorF2D<T>& indicator,
    const UnitConverter<T,DESCRIPTOR>& converter,
    bool insideOut)
  : BlockLatticeF2D<T, DESCRIPTOR>(blockLattice, 1),
    _blockGeometry(blockGeometry), _indicator(indicator), _converter(converter), _insideOut(insideOut)
{
  this->getName() = "volumeFractionPolygonApproximation";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeVolumeFractionPolygonApproximation2D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  output[0] = 0.;
  T physR[2];

  bool cornerXMYM_inside[1];
  bool cornerXMYP_inside[1];
  bool cornerXPYM_inside[1];
  bool cornerXPYP_inside[1];
  _blockGeometry.getPhysR(physR, input[0], input[1]);

  T cornerXMYM[2];
  T cornerXMYP[2];
  T cornerXPYM[2];
  T cornerXPYP[2];

  cornerXMYM[0] = physR[0] - 0.5*_converter.getPhysDeltaX();
  cornerXMYM[1] = physR[1] - 0.5*_converter.getPhysDeltaX();

  cornerXMYP[0] = physR[0] - 0.5*_converter.getPhysDeltaX();
  cornerXMYP[1] = physR[1] + 0.5*_converter.getPhysDeltaX();

  cornerXPYM[0] = physR[0] + 0.5*_converter.getPhysDeltaX();
  cornerXPYM[1] = physR[1] - 0.5*_converter.getPhysDeltaX();

  cornerXPYP[0] = physR[0] + 0.5*_converter.getPhysDeltaX();
  cornerXPYP[1] = physR[1] + 0.5*_converter.getPhysDeltaX();

  _indicator(cornerXMYM_inside, cornerXMYM);
  _indicator(cornerXMYP_inside, cornerXMYP);
  _indicator(cornerXPYM_inside, cornerXPYM);
  _indicator(cornerXPYP_inside, cornerXPYP);

  if ((cornerXMYM_inside[0] && cornerXMYP_inside[0] && cornerXPYM_inside[0] && cornerXPYP_inside[0]) ||
      (!cornerXMYM_inside[0] && !cornerXMYP_inside[0] && !cornerXPYM_inside[0] && !cornerXPYP_inside[0]) ) {
    if (!_insideOut) {
      if (cornerXMYM_inside[0]) {
        output[0] = 1.0;
      }
    }
    else {
      if (!cornerXMYM_inside[0]) {
        output[0] = 1.0;
      }
    }
  } else {
    Vector<T,2> cornerXMYM_vec(physR[0] - 0.5*_converter.getPhysDeltaX(), physR[1] - 0.5*_converter.getPhysDeltaX());
    Vector<T,2> cornerXPYP_vec(physR[0] + 0.5*_converter.getPhysDeltaX(), physR[1] + 0.5*_converter.getPhysDeltaX());

    Vector<T,2> directionXP(_converter.getPhysDeltaX(), 0);
    Vector<T,2> directionXM(-_converter.getPhysDeltaX(), 0);
    Vector<T,2> directionYP(0, _converter.getPhysDeltaX());
    Vector<T,2> directionYM(0, -_converter.getPhysDeltaX());

    T distanceXP = 1.01 * _converter.getPhysDeltaX();
    T distanceXM = 1.01 * _converter.getPhysDeltaX();
    T distanceYP = 1.01 * _converter.getPhysDeltaX();
    T distanceYM = 1.01 * _converter.getPhysDeltaX();

    if ((cornerXMYM_inside[0] && !cornerXMYP_inside[0]) ||
        (!cornerXMYM_inside[0] && cornerXMYP_inside[0]) ) {
      _indicator.distance(distanceYP, cornerXMYM, directionYP);
    }

    if ((cornerXMYM_inside[0] && !cornerXPYM_inside[0]) ||
        (!cornerXMYM_inside[0] && cornerXPYM_inside[0]) ) {
      _indicator.distance(distanceXP, cornerXMYM, directionXP);
    }

    if ((cornerXPYP_inside[0] && !cornerXMYP_inside[0]) ||
        (!cornerXPYP_inside[0] && cornerXMYP_inside[0]) ) {
      _indicator.distance(distanceXM, cornerXPYP, directionXM);
    }

    if ((cornerXPYP_inside[0] && !cornerXPYM_inside[0]) ||
        (!cornerXPYP_inside[0] && cornerXPYM_inside[0]) ) {
      _indicator.distance(distanceYM, cornerXPYP, directionYM);
    }

    T volumeFraction = 0.0;

    if (distanceXP < _converter.getPhysDeltaX() && distanceXM < _converter.getPhysDeltaX()) {
      volumeFraction = distanceXP * _converter.getPhysDeltaX();
      volumeFraction += 0.5 * (_converter.getPhysDeltaX() - distanceXM - distanceXP) * _converter.getPhysDeltaX();
      volumeFraction /= _converter.getPhysDeltaX() * _converter.getPhysDeltaX();
      if (!cornerXMYM_inside[0]) {
        volumeFraction = 1.0 - volumeFraction;
      }
    }

    if (distanceYP < _converter.getPhysDeltaX() && distanceYM < _converter.getPhysDeltaX()) {
      volumeFraction = distanceYP * _converter.getPhysDeltaX();
      volumeFraction += 0.5 * (_converter.getPhysDeltaX() - distanceYM - distanceYP) * _converter.getPhysDeltaX();
      volumeFraction /= _converter.getPhysDeltaX() * _converter.getPhysDeltaX();
      if (!cornerXMYM_inside[0]) {
        volumeFraction = 1.0 - volumeFraction;
      }
    }

    if (distanceXP < _converter.getPhysDeltaX() && distanceYP < _converter.getPhysDeltaX()) {
      volumeFraction = 0.5 * distanceXP * distanceYP;
      volumeFraction /= _converter.getPhysDeltaX() * _converter.getPhysDeltaX();
      if (!cornerXMYM_inside[0]) {
        volumeFraction = 1.0 - volumeFraction;
      }
    }

    if (distanceXM < _converter.getPhysDeltaX() && distanceYM < _converter.getPhysDeltaX()) {
      volumeFraction = 0.5 * distanceXM * distanceYM;
      volumeFraction /= _converter.getPhysDeltaX() * _converter.getPhysDeltaX();
      if (!cornerXPYP_inside[0]) {
        volumeFraction = 1.0 - volumeFraction;
      }
    }

    if (distanceXM < _converter.getPhysDeltaX() && distanceYP < _converter.getPhysDeltaX()) {
      volumeFraction = 0.5 * (_converter.getPhysDeltaX() - distanceXM) * (_converter.getPhysDeltaX() - distanceYP);
      volumeFraction /= _converter.getPhysDeltaX() * _converter.getPhysDeltaX();
      if (!cornerXMYP_inside[0]) {
        volumeFraction = 1.0 - volumeFraction;
      }
    }

    if (distanceYM < _converter.getPhysDeltaX() && distanceXP < _converter.getPhysDeltaX()) {
      volumeFraction = 0.5 * (_converter.getPhysDeltaX() - distanceXP) * (_converter.getPhysDeltaX() - distanceYM);
      volumeFraction /= _converter.getPhysDeltaX() * _converter.getPhysDeltaX();
      if (!cornerXPYM_inside[0]) {
        volumeFraction = 1.0 - volumeFraction;
      }
    }

    if (!_insideOut) {
        output[0] = volumeFraction;
    }
    else {
        output[0] = 1.0 - volumeFraction;
    }

  }
  return true;
}



template <typename T, typename DESCRIPTOR>
BlockLatticePhysPermeability2D<T,DESCRIPTOR>::BlockLatticePhysPermeability2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, BlockGeometryStructure2D<T>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "permeability";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysPermeability2D<T,DESCRIPTOR>::operator()(T output[], const int input[])
{
  return false;
}


template <typename T, typename DESCRIPTOR>
BlockLatticePhysDarcyForce2D<T,DESCRIPTOR>::BlockLatticePhysDarcyForce2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, BlockGeometry2D<T>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2),
    _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "alphaU";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysDarcyForce2D<T,DESCRIPTOR>::operator()(T output[], const int input[])
{
  BlockLatticePhysPermeability2D<T,DESCRIPTOR> permeability(this->_blockLattice,this->_blockGeometry,this->_material,this->_converter);
  BlockLatticeVelocity2D<T,DESCRIPTOR> velocity(this->_blockLattice);

  T nu = this->_converter.getPhysViscosity();
  permeability(output,input);
  T K = output[0];
  velocity(output,input);

  output[0] *= -nu/K;
  output[1] *= -nu/K;

  return true;
}


template <typename T, typename DESCRIPTOR>
BlockLatticeAverage2D<T,DESCRIPTOR>::BlockLatticeAverage2D
(BlockLatticeF2D<T,DESCRIPTOR>& f, BlockGeometry2D<T>& blockGeometry,
 int material, T radius)
  : BlockLatticeF2D<T,DESCRIPTOR>(f.getBlockLattice(), f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material), _radius(radius)
{
  this->getName() = "Average("+f.getName()+")";
}


template <typename T, typename DESCRIPTOR>
bool BlockLatticeAverage2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //  CuboidGeometry2D<T>& cGeometry = f.getBlockLattice2D().get_cGeometry();
  //  loadBalancer& load = f.getBlockLattice2D().get_load();

  //  //create boolean indicator functor isInSphere
  //  std::vector<T> center(3,0);
  //  center[0] = (int)cGeometry.get(load.glob(input[0])).get_globPosX() + input[1];
  //  center[1] = (int)cGeometry.get(load.glob(input[0])).get_globPosY() + input[2];
  //  center[2] = (int)cGeometry.get(load.glob(input[0])).get_globPosZ() + input[3];
  //  SphereAnalyticalF2D<bool,T> isInSphere(center,radius);

  // iterate over all cuboids & points and test for material && isInSphere
  //  std::vector<T> tmp( this->_n, T() );
  //  int numVoxels(0);
  //  if (this->_blockGeometry.getMaterial(center[0],center[1],center[2]) == material) {
  //    for (int iC=0; iC<load.size(); iC++) {
  //      int nX = cGeometry.get(load.glob(iC)).getNx();
  //      int nY = cGeometry.get(load.glob(iC)).getNy();
  //      int nZ = cGeometry.get(load.glob(iC)).getNz();
  //      for (int iX=0; iX<nX; ++iX) {
  //        for (int iY=0; iY<nY; ++iY) {
  //          for (int iZ=0; iZ<nZ; iZ++) {
  //            std::vector<T> glob(3,0);
  //            glob[0] = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
  //            glob[1] = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
  //            glob[2] = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
  //            if (this->_blockGeometry.getMaterial(glob[0],glob[1],glob[2]) == material
  //                && isInSphere(glob)[0]==true) {
  //              for (unsigned iD=0; iD<f(load.glob(0),0,0,0).size(); iD++) {
  //                tmp[iD]+=f(load.glob(iC),iX,iY,iZ)[iD];
  //              }
  //              numVoxels++;
  //            }
  //          }
  //        }
  //      }
  //    }

  //#ifdef PARALLEL_MODE_MPI
  //    singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
  //#endif
  //    for (int iD=0; iD<f.getTargetDim(); iD++) {
  //#ifdef PARALLEL_MODE_MPI
  //      singleton::mpi().reduceAndBcast(tmp[iD], MPI_SUM);
  //#endif
  //      if (numVoxels>0) {
  //        tmp[iD] /= numVoxels;
  //      }
  //    }
  //  }
  //  return tmp;

  return false;
}


template <typename T, typename DESCRIPTOR>
BlockEuklidNorm2D<T,DESCRIPTOR>::BlockEuklidNorm2D(BlockF2D<T>& f)
  : BlockF2D<T>(f.getBlockStructure(),1), _f(f)
{
  this->getName() = "l2("+f.getName()+")";
}

template <typename T, typename DESCRIPTOR>
bool BlockEuklidNorm2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = T();  // flash output, since this methods adds values.
  T data[_f.getTargetDim()];
  _f(data,input);
  for ( int i = 0; i < _f.getTargetDim(); ++i) {
    output[0] += data[i]*data[i];
  }
  output[0] = sqrt(output[0]);
  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticePorousMomentumLossForce2D<T, DESCRIPTOR>::BlockLatticePorousMomentumLossForce2D(
  BlockLatticeStructure2D<T, DESCRIPTOR>& blockLattice, BlockGeometryStructure2D<T>& blockGeometry,
  std::vector<SmoothIndicatorF2D<T,T,true>* >& indicator, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T, DESCRIPTOR>(blockLattice, converter, 4*indicator.size()), _blockGeometry(blockGeometry), _vectorOfIndicator(indicator)
{
  this->getName() = "physPorousMomentumLossForce";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePorousMomentumLossForce2D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  // iterate over all particles in _indicator
  for (size_t iInd=0; iInd!=_vectorOfIndicator.size(); iInd++) {

    int numVoxels = 0;
    int start[2] = {0};
    int end[2] = {0};
    // check for intersection of cuboid and indicator
    Cuboid2D<T> tmpCuboid(_blockGeometry.getOrigin()[0], _blockGeometry.getOrigin()[1], this->_converter.getPhysDeltaX(), _blockGeometry.getNx(), _blockGeometry.getNy());
    T posXmin = _vectorOfIndicator[iInd]->getPos()[0] - _vectorOfIndicator[iInd]->getCircumRadius();
    T posXmax = _vectorOfIndicator[iInd]->getPos()[0] + _vectorOfIndicator[iInd]->getCircumRadius();
    T posYmin = _vectorOfIndicator[iInd]->getPos()[1] - _vectorOfIndicator[iInd]->getCircumRadius();
    T posYmax = _vectorOfIndicator[iInd]->getPos()[1] + _vectorOfIndicator[iInd]->getCircumRadius();
    if (tmpCuboid.checkInters(posXmin, posXmax, posYmin, posYmax, start[0], end[0], start[1], end[1])) {

      T invDeltaX = 1./this->_converter.getPhysDeltaX();
      for(int k=0; k<2; k++) {
        start[k] -= 1;
        if (start[k] < 0) start[k] = 0;
        end[k] += 2;
        if (end[k] > _blockGeometry.getExtend()[k]) end[k] = _blockGeometry.getExtend()[k];
      }
      // iterate over cells in the constructed intersection box
      for (int iX = start[0]; iX < end[0]; iX++) {
        for (int iY = start[1]; iY < end[1]; iY++) {
          // check if cell belongs to particle
          T inside[1] = {0.};
          T posIn[2] = {0.};
          _blockGeometry.getPhysR(posIn, iX, iY);
          (*(_vectorOfIndicator[iInd]))( inside, posIn);
          if ( !util::nearZero(inside[0]) && this->_blockGeometry.get(iX,iY)==1) {
            // compute momentum exchange force on particle
            T tmpForce[2] = {0.,0.};
            tmpForce[0] += this->_blockLattice.get(iX, iY).template getFieldPointer<descriptors::VELOCITY_NUMERATOR>()[0];
            tmpForce[1] += this->_blockLattice.get(iX, iY).template getFieldPointer<descriptors::VELOCITY_NUMERATOR>()[1];
            // reset external field for next timestep
            T reset_to_zero[2] = {0.,0.};
            this->_blockLattice.get(iX, iY).template setField<descriptors::VELOCITY_NUMERATOR>(reset_to_zero);
            // convert force to SI units and compute torque
            numVoxels++;
            // division bei length of lattice cell necessary due to converter handling of force
            tmpForce[0] = this->_converter.getPhysForce(tmpForce[0])*invDeltaX;
            tmpForce[1] = this->_converter.getPhysForce(tmpForce[1])*invDeltaX;
            output[0+iInd*4] += tmpForce[0];
            output[1+iInd*4] += tmpForce[1];
            output[2+iInd*4] += (posIn[0]-_vectorOfIndicator[iInd]->getPos()[0])*tmpForce[1] - (posIn[1]-_vectorOfIndicator[iInd]->getPos()[1])*tmpForce[0];
          }
        }
      }

    }
    output[3+iInd*4] = numVoxels;

  }
  return true;
}

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
BlockLatticePhysTemperature2D<T,DESCRIPTOR,TDESCRIPTOR>::BlockLatticePhysTemperature2D
(BlockLatticeStructure2D<T,TDESCRIPTOR>& blockLattice, ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter)
  : BlockLatticeThermalPhysF2D<T,DESCRIPTOR,TDESCRIPTOR>(blockLattice,converter,1)
{
  this->getName() = "physTemperature";
}


template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
bool BlockLatticePhysTemperature2D<T,DESCRIPTOR,TDESCRIPTOR>::operator() (T output[], const int input[])
{
  T latticeTemperature = this->_blockLattice.get( input[0], input[1] ).computeRho();
  output[0] = this->_converter.getPhysTemperature(latticeTemperature);

  return true;
}

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
BlockLatticePhysHeatFlux2D<T,DESCRIPTOR,TDESCRIPTOR>::BlockLatticePhysHeatFlux2D
(BlockLatticeStructure2D<T,TDESCRIPTOR>& blockLattice, ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter)
  : BlockLatticeThermalPhysF2D<T,DESCRIPTOR,TDESCRIPTOR>(blockLattice,converter,2),
    _temp(converter.getLatticeSpecificHeatCapacity(converter.getPhysSpecificHeatCapacity())*(converter.getLatticeThermalRelaxationTime() - 0.5) / converter.getLatticeThermalRelaxationTime())
{
  this->getName() = "physHeatFlux";

  // std::cout << _temp << " " << converter.getConversionFactorHeatFlux() << " " << 0.426 / _temp / converter.getConversionFactorHeatFlux() << std::endl;
}

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
bool BlockLatticePhysHeatFlux2D<T,DESCRIPTOR,TDESCRIPTOR>::operator() (T output[], const int input[])
{
  T temperature, extVel[2], j[2];
  this->_blockLattice.get( input[0], input[1] ).computeRhoU(temperature,extVel);
  this->_blockLattice.get( input[0], input[1] ).computeJ(j);

  output[0] = this->_converter.getPhysHeatFlux((j[0] - temperature * extVel[0])*_temp);
  output[1] = this->_converter.getPhysHeatFlux((j[1] - temperature * extVel[1])*_temp);

  return true;
}

template<typename T, typename DESCRIPTOR, bool HLBM>
BlockLatticeIndicatorSmoothIndicatorIntersection2D<T,DESCRIPTOR,HLBM>::BlockLatticeIndicatorSmoothIndicatorIntersection2D (
  BlockLatticeStructure2D<T, DESCRIPTOR>& blockLattice,
  BlockGeometryStructure2D<T>& blockGeometry,
  IndicatorF2D<T>& normalInd, SmoothIndicatorF2D<T,T,HLBM>& smoothInd )
  : BlockLatticeF2D<T, DESCRIPTOR>(blockLattice, 1),
    _blockGeometry(blockGeometry), _normalInd(normalInd), _smoothInd(smoothInd)
{
  this->getName() = "Indicator-SmoothIndicator Intersection";
}

template<typename T, typename DESCRIPTOR, bool HLBM>
bool BlockLatticeIndicatorSmoothIndicatorIntersection2D<T, DESCRIPTOR,HLBM>::operator()(T output[], const int input[])
{
  output[0] = 0.;
  int start[2] = {0};
  int end[2] = {0};
  // check for intersection of cuboid and smoothIndicator
  Cuboid2D<T> tmpCuboid(_blockGeometry.getOrigin()[0], _blockGeometry.getOrigin()[1], _blockGeometry.getDeltaR(), _blockGeometry.getNx(), _blockGeometry.getNy());
  T posXmin = _smoothInd.getPos()[0] - _smoothInd.getCircumRadius();
  T posXmax = _smoothInd.getPos()[0] + _smoothInd.getCircumRadius();
  T posYmin = _smoothInd.getPos()[1] - _smoothInd.getCircumRadius();
  T posYmax = _smoothInd.getPos()[1] + _smoothInd.getCircumRadius();
  if (tmpCuboid.checkInters(posXmin, posXmax, posYmin, posYmax, start[0], end[0], start[1], end[1])) {

    for(int k=0; k<2; k++) {
      start[k] -= 1;
      if (start[k] < 0) start[k] = 0;
      end[k] += 2;
      if (end[k] > _blockGeometry.getExtend()[k]) end[k] = _blockGeometry.getExtend()[k];
    }

    // iterate over cells in the constructed intersection box
    for (int iX = start[0]; iX < end[0]; iX++) {
      for (int iY = start[1]; iY < end[1]; iY++) {

        // check if cell belongs to particle
        T insideT[1] = {0.};
        T posIn[2] = {0.};
        _blockGeometry.getPhysR(posIn, iX, iY);
        _smoothInd( insideT, posIn);
        if ( !util::nearZero(insideT[0]) && this->_blockGeometry.get(iX,iY)==1) {
          // Return 1 if at least one cell is found to be inside both A and B
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

  return true;
}

/////////// BlockLatticeGuoZhaoEpsilon2D /////////////////////////////////////////////

template <typename T, typename DESCRIPTOR>
BlockLatticeGuoZhaoEpsilon2D<T,DESCRIPTOR>::BlockLatticeGuoZhaoEpsilon2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice, 1)
{
  this->getName() = "epsilon";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeGuoZhaoEpsilon2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = *this->_blockLattice.get( input[0], input[1] ).template getFieldPointer<descriptors::EPSILON>();
  return true;
}

/////////// BlockLatticeGuoZhaoPhysK2D /////////////////////////////////////////////

template <typename T, typename DESCRIPTOR>
BlockLatticeGuoZhaoPhysK2D<T,DESCRIPTOR>::BlockLatticeGuoZhaoPhysK2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice, converter, 1)
{
  this->getName() = "physK";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeGuoZhaoPhysK2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = *this->_blockLattice.get( input[0], input[1] )
              .template getFieldPointer<descriptors::K>()
              * this->_converter.getConversionFactorLength()
              * this->_converter.getConversionFactorLength();
  return true;
}

/////////// BlockLatticeGuoZhaoPhysBodyForce2D /////////////////////////////////////////////

template <typename T, typename DESCRIPTOR>
BlockLatticeGuoZhaoPhysBodyForce2D<T,DESCRIPTOR>::BlockLatticeGuoZhaoPhysBodyForce2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice, converter, 2)
{
  this->getName() = "physBodyForce";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeGuoZhaoPhysBodyForce2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  output[0] = this->_converter.getPhysForce( *this->_blockLattice.get( input[0], input[1] )
              .template getFieldPointer<descriptors::BODY_FORCE>() );
  output[1] = this->_converter.getPhysForce( *this->_blockLattice.get( input[0], input[1] )
              .template getFieldPointer<descriptors::BODY_FORCE>(1) );
  return true;
}

} // end namespace olb

#endif

