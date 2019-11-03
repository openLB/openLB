/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2019 Davide Dapelo
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

/* Models for Lagrangian back-coupling methods -- generic implementation.
 */

#ifndef LB_BACK_COUPLING_MODELS_HH
#define LB_BACK_COUPLING_MODELS_HH

namespace olb {


////////////////////// Class BaseBackCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
BaseBackCouplingModel<T,Lattice,Particle>::BaseBackCouplingModel (
                 UnitConverter<T, Lattice>& converter,
                         SuperLattice3D<T, Lattice>& sLattice,
                         SuperGeometry3D<T>& sGeometry )
         : _converter(converter),
           _sGeometry(sGeometry),
           _sLattice(sLattice)
{
  _zeroAnalytical = std::make_shared<AnalyticalConst3D<T, T> > (T());
  _zeroField = std::make_shared<AnalyticalComposed3D<T, T> > (*_zeroAnalytical, *_zeroAnalytical, *_zeroAnalytical);
}

template<typename T, typename Lattice, template<typename V> class Particle>
void BaseBackCouplingModel<T,Lattice,Particle>::resetExternalField(int material)
{
  // resets external field
  this->_sLattice.template defineField<descriptors::FORCE>(this->_sGeometry, material, *_zeroField);

  // NECESSARY to communicate values before using them in operator()
  this->_sLattice.communicate();
}


////////////////////// Class CubicDeltaBackCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
CubicDeltaBackCouplingModel<T,Lattice,Particle>::CubicDeltaBackCouplingModel (
                 UnitConverter<T, Lattice>& converter,
                         SuperLattice3D<T, Lattice>& sLattice,
                         SuperGeometry3D<T>& sGeometry )
         : BaseBackCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry)
{
  _cubicDeltaFunctional = std::make_shared<SuperLatticeSmoothDiracDelta3D<T, Lattice> > (
      this->_sLattice, this->_converter, this->_sGeometry );
}

template<typename T, typename Lattice, template<typename V> class Particle>
bool CubicDeltaBackCouplingModel<T,Lattice,Particle>::operator() (Particle<T>* p, int globic, int material, int subCycles)
{
  int locIC = this->_sLattice.getLoadBalancer().loc(globic);

  // reading the force from the value stored inside the particle
  std::vector<T> physForceP = p->getStoreForce(); // physical force acting on the particle

  T latticeForceP[3] = {T(), T(), T()}; // dimensionless force acting on the particle
  latticeForceP[0] = physForceP[0] / this->_converter.getConversionFactorForce();
  latticeForceP[1] = physForceP[1] / this->_converter.getConversionFactorForce();
  latticeForceP[2] = physForceP[2] / this->_converter.getConversionFactorForce();

  T physPosP[3] = {T(), T(), T()}; // particle's physical position
  physPosP[0] = (p->getPos()[0]);
  physPosP[1] = (p->getPos()[1]);
  physPosP[2] = (p->getPos()[2]);

  // particle's dimensionless position, rounded at neighbouring voxel
  int latticeRoundedPosP[3] = {0, 0, 0};
  this->_sLattice.getCuboidGeometry().get(globic).getLatticeR (
           latticeRoundedPosP, physPosP );

  // smooth Dirac delta
  this->_cubicDeltaFunctional->operator() (_delta, physPosP, globic);

  T tempDelta = T();
  T F[3] = {T(), T(), T()}; // dimensionless smoothed force

  for (int i = -_range; i <= _range; ++i) {
    for (int j = -_range; j <= _range; ++j) {
      for (int k = -_range; k <= _range; ++k) {
        if (this->_sGeometry.getBlockGeometry(locIC).getMaterial(
              latticeRoundedPosP[0] + i, latticeRoundedPosP[1] + j,
              latticeRoundedPosP[2] + k) == material) {

          tempDelta = _delta[i + _range][j + _range][k + _range];

          F[0] = -latticeForceP[0] * tempDelta / (T)(subCycles);
          F[1] = -latticeForceP[1] * tempDelta / (T)(subCycles);
          F[2] = -latticeForceP[2] * tempDelta / (T)(subCycles);

          this->_sLattice.getBlockLattice(locIC).get (
                   latticeRoundedPosP[0] + i,
                   latticeRoundedPosP[1] + j,
                   latticeRoundedPosP[2] + k ).template addField<descriptors::FORCE>( F );
        }
      }
    }
  }
  return true;
}


////////////////////// Class LocalDeltaBackCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
LocalBackCouplingModel<T,Lattice,Particle>::LocalBackCouplingModel (
                 UnitConverter<T, Lattice>& converter,
                         SuperLattice3D<T, Lattice>& sLattice,
                         SuperGeometry3D<T>& sGeometry )
         : BaseBackCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry)
{}

template<typename T, typename Lattice, template<typename V> class Particle>
bool LocalBackCouplingModel<T,Lattice,Particle>::operator() (Particle<T>* p, int globic, int material, int subCycles)
{
  int locIC = this->_sLattice.getLoadBalancer().loc(globic);

  // reading the force from the value stored inside the particle
  std::vector<T> physForceP = p->getStoreForce(); // physical force acting on the particle

  T latticeForceP[3] = {T(), T(), T()}; // dimensionless force acting on the particle
  latticeForceP[0] = physForceP[0] / this->_converter.getConversionFactorForce();
  latticeForceP[1] = physForceP[1] / this->_converter.getConversionFactorForce();
  latticeForceP[2] = physForceP[2] / this->_converter.getConversionFactorForce();

  T physPosP[3] = {T(), T(), T()}; // particle's physical position
  physPosP[0] = (p->getPos()[0]);
  physPosP[1] = (p->getPos()[1]);
  physPosP[2] = (p->getPos()[2]);

  // particle's dimensionless position, rounded at neighbouring voxel
  int latticeRoundedPosP[3] = {0, 0, 0};
  this->_sLattice.getCuboidGeometry().get(globic).getLatticeR (
           latticeRoundedPosP, physPosP );

  if (this->_sGeometry.getBlockGeometry(locIC).getMaterial(
        latticeRoundedPosP[0], latticeRoundedPosP[1],
        latticeRoundedPosP[2]) == material) {

    T F[3] = {T(), T(), T()}; // dimensionless smoothed force
    F[0] = -latticeForceP[0] / (T)(subCycles);
    F[1] = -latticeForceP[1] / (T)(subCycles);
    F[2] = -latticeForceP[2] / (T)(subCycles);

    this->_sLattice.getBlockLattice(locIC).get (
             latticeRoundedPosP[0],
             latticeRoundedPosP[1],
             latticeRoundedPosP[2] ).template addField<descriptors::FORCE>( F );
  }

  return true;
}


////////////////////// Class NonLocalBaseBackCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
NonLocalBaseBackCouplingModel<T,Lattice,Particle>::NonLocalBaseBackCouplingModel (
                         UnitConverter<T, Lattice>& converter,
                         SuperLattice3D<T, Lattice>& sLattice,
                         SuperGeometry3D<T>& sGeometry,
                         SmoothingFunctional<T, Lattice>& smoothingFunctional )
         : BaseBackCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry),
           _smoothingFunctional(smoothingFunctional)
{}

template<typename T, typename Lattice, template<typename V> class Particle>
bool NonLocalBaseBackCouplingModel<T,Lattice,Particle>::operator() (Particle<T>* p, int globic, int material, int subCycles)
{
  int locIC = this->_sLattice.getLoadBalancer().loc(globic);

  // reading the force from the value stored inside the particle
  std::vector<T> physForceP = p->getStoreForce(); // physical force acting on the particle
  T latticeForceP[3] = {T(), T(), T()}; // dimensionless force acting on the particle
  latticeForceP[0] = physForceP[0] / this->_converter.getConversionFactorForce();
  latticeForceP[1] = physForceP[1] / this->_converter.getConversionFactorForce();
  latticeForceP[2] = physForceP[2] / this->_converter.getConversionFactorForce();

  // Updating force through voxels within kernel smoothing length from the bubble's position
  for (int i=0; i<this->_smoothingFunctional.getSize(); i++) {

    // Position of the iterated voxel
    int iLatticePosF[3] = {0, 0, 0};
    this->_smoothingFunctional.getLatticePos(iLatticePosF, i);

    // Updating iterated voxel
    if (this->_sGeometry.getBlockGeometry(locIC).getMaterial(
          iLatticePosF[0], iLatticePosF[1],
          iLatticePosF[2]) == material) {

      // Weighted force acting on the iterated voxel
      T F[3] = {T(), T(), T()}; // dimensionless smoothed force
      F[0] = -latticeForceP[0] * this->_smoothingFunctional.getWeight(i) / (T)(subCycles);
      F[1] = -latticeForceP[1] * this->_smoothingFunctional.getWeight(i) / (T)(subCycles);
      F[2] = -latticeForceP[2] * this->_smoothingFunctional.getWeight(i) / (T)(subCycles);

      this->_sLattice.getBlockLattice(locIC).get (
              iLatticePosF[0], iLatticePosF[1], iLatticePosF[2] ).template addField<descriptors::FORCE>( F );
    }
  }
  return true;
}

}

#endif
