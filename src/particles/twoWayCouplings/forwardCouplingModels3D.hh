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

/* Models for Lagrangian forward-coupling methods -- generic implementation.
 */

#ifndef LB_FORWARD_COUPLING_MODELS_HH
#define LB_FORWARD_COUPLING_MODELS_HH

namespace olb {


////////////////////// Class LocalBaseCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
LocalBaseForwardCouplingModel<T,Lattice,Particle>::LocalBaseForwardCouplingModel (
                     UnitConverter<T, Lattice>& converter,
                     SuperLattice3D<T, Lattice>& sLattice,
                     SuperGeometry3D<T>& sGeometry,
                     DragModel<T, Particle>& dragModel )
         : _sGeometry(sGeometry),
           _converter(converter),
           _sLattice(sLattice),
           _dragModel(dragModel)
{
  this->_interpLatticeDensity = std::make_shared<SuperLatticeInterpDensity3Degree3D<T, Lattice> > (
                                this->_sLattice, _sGeometry, this->_converter );
  this->_interpLatticeVelocity = std::make_shared<SuperLatticeInterpPhysVelocity3D<T, Lattice> > (
                                 this->_sLattice, this->_converter );
}

template<typename T, typename Lattice, template<typename V> class Particle>
bool LocalBaseForwardCouplingModel<T,Lattice,Particle>::operator() (Particle<T>* p, int globic)
{
  /// Getting the particle and its containing cell's position
  T physPosP[3] = {T(), T(), T()}; // particle's physical position
  physPosP[0] = (p->getPos()[0]);
  physPosP[1] = (p->getPos()[1]);
  physPosP[2] = (p->getPos()[2]);
  // particle's dimensionless position, rounded at neighbouring voxel
  int latticeRoundedPosP[3] = {0, 0, 0};
  this->_sLattice.getCuboidGeometry().get(globic).getLatticeR (
                  latticeRoundedPosP, physPosP );
  // { globic, latticeRoundedP[0, ..., 2] }
  int globicFull[4] = { globic,
                        latticeRoundedPosP[0],
                        latticeRoundedPosP[1],
                        latticeRoundedPosP[2] };

  // Particle's velocity
  T physVelP[3] = {T(), T(), T()}; // Physical
  physVelP[0] = (p->getVel()[0]);
  physVelP[1] = (p->getVel()[1]);
  physVelP[2] = (p->getVel()[2]);
  // Lattice
  T latticeVelP[3] = {T(), T(), T()}; // particle's dimensionless velocity
  latticeVelP[0] = this->_converter.getLatticeVelocity(physVelP[0]);
  latticeVelP[1] = this->_converter.getLatticeVelocity(physVelP[1]);
  latticeVelP[2] = this->_converter.getLatticeVelocity(physVelP[2]);

  // Lattice's velocity at particle's location
  T physVelF[3] = {T(), T(), T()}; // Physical
  this->_interpLatticeVelocity->operator() (physVelF, physPosP, globic);
  // Lattice
  T latticeVelF[3] = {T(), T(), T()}; // Lattice's dimensionless velocity at particle's location
  latticeVelF[0] = this->_converter.getLatticeVelocity(physVelF[0]);
  latticeVelF[1] = this->_converter.getLatticeVelocity(physVelF[1]);
  latticeVelF[2] = this->_converter.getLatticeVelocity(physVelF[2]);

  // Computing fluid-particle momentum transfer
  T gF[3] = {T(), T(), T()}; // force density gF
  this->_momentumExchange->operator() ( gF, latticeVelF, latticeVelP,
                           physPosP, latticeRoundedPosP, globic);

  // Computing drag coefficient
  T Cd = this->_dragModel(p, latticeVelF, latticeVelP, globicFull);

  /// Computing drag force in dimensionless units
  T latticePRad = p->getRad() / _converter.getConversionFactorLength();
  T latticeForceP[3] = {T(), T(), T()}; // dimensionless force acting on the particle
  latticeForceP[0] = .5 * Cd * M_PI*pow(latticePRad,2) * gF[0] * (latticeVelF[0] - latticeVelP[0]);
  latticeForceP[1] = .5 * Cd * M_PI*pow(latticePRad,2) * gF[1] * (latticeVelF[1] - latticeVelP[1]);
  latticeForceP[2] = .5 * Cd * M_PI*pow(latticePRad,2) * gF[2] * (latticeVelF[2] - latticeVelP[2]);

  /// Computing physical drag force
  std::vector<T> physForceP(3, T()); // physical force acting on the particle
  physForceP[0] = latticeForceP[0] * this->_converter.getConversionFactorForce();
  physForceP[1] = latticeForceP[1] * this->_converter.getConversionFactorForce();
  physForceP[2] = latticeForceP[2] * this->_converter.getConversionFactorForce();

  /// Updating the particle
  p->setStoreForce(physForceP);
  return true;
}


////////////////////// Class NaiveForwardCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
NaiveForwardCouplingModel<T,Lattice,Particle>::NaiveForwardCouplingModel (
                     UnitConverter<T, Lattice>& converter,
                     SuperLattice3D<T, Lattice>& sLattice,
                     SuperGeometry3D<T>& sGeometry,
                     DragModel<T, Particle>& dragModel )
         : LocalBaseForwardCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry, dragModel)
{
  this->_momentumExchange = std::make_shared<NaiveMomentumExchange<T, Lattice> > (
                            this->_converter, this->_sLattice, this->_interpLatticeDensity );
}


////////////////////// Class LaddForwardCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
LaddForwardCouplingModel<T,Lattice,Particle>::LaddForwardCouplingModel (
                     UnitConverter<T, Lattice>& converter,
                     SuperLattice3D<T, Lattice>& sLattice,
                     SuperGeometry3D<T>& sGeometry,
                     DragModel<T, Particle>& dragModel )
         : LocalBaseForwardCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry, dragModel)
{
  this->_momentumExchange = std::make_shared<LaddMomentumExchange<T, Lattice> > (
                            this->_converter, this->_sLattice,
                            this->_interpLatticeDensity, this->_interpLatticeVelocity );
}


////////////////////// Class NonLocalBaseCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
NonLocalBaseForwardCouplingModel<T,Lattice,Particle>::NonLocalBaseForwardCouplingModel (
                     UnitConverter<T, Lattice>& converter,
                     SuperLattice3D<T, Lattice>& sLattice,
                     SuperGeometry3D<T>& sGeometry,
                     DragModel<T, Particle>& dragModel,
                     SmoothingFunctional<T, Lattice>& smoothingFunctional )
         : _sGeometry(sGeometry),
           _converter(converter),
           _sLattice(sLattice),
           _dragModel(dragModel),
           _smoothingFunctional(smoothingFunctional)
{
  this->_interpLatticeDensity = std::make_shared<SuperLatticeInterpDensity3Degree3D<T, Lattice> > (
                                this->_sLattice, _sGeometry, this->_converter );
  this->_interpLatticeVelocity = std::make_shared<SuperLatticeInterpPhysVelocity3D<T, Lattice> > (
                                 this->_sLattice, this->_converter );
}

template<typename T, typename Lattice, template<typename V> class Particle>
bool NonLocalBaseForwardCouplingModel<T,Lattice,Particle>::operator() (Particle<T>* p, int globic)
{
  /// Getting the particle and its containing cell's position
  T physPosP[3] = {T(), T(), T()}; // particle's physical position
  physPosP[0] = (p->getPos()[0]);
  physPosP[1] = (p->getPos()[1]);
  physPosP[2] = (p->getPos()[2]);
  //std::cout <<  "globic=" << globic << "   physPosP=(" << physPosP[0] << ", " << physPosP[1] << ", " << physPosP[2] << ")" << std::endl; // <---
  // particle's dimensionless position, rounded at neighbouring voxel
  int latticeRoundedPosP[3] = {0, 0, 0};
  this->_sLattice.getCuboidGeometry().get(globic).getLatticeR (
           latticeRoundedPosP, physPosP );
  // { globic, latticeRoundedP[0, ..., 2] }
  int globicFull[4] = { globic,
                        latticeRoundedPosP[0],
                        latticeRoundedPosP[1],
                        latticeRoundedPosP[2] };

  // Particle's velocity
  T physVelP[3] = {T(), T(), T()}; // Physical
  physVelP[0] = (p->getVel()[0]);
  physVelP[1] = (p->getVel()[1]);
  physVelP[2] = (p->getVel()[2]);
  // Lattice
  T latticeVelP[3] = {T(), T(), T()}; // particle's dimensionless velocity
  latticeVelP[0] = this->_converter.getLatticeVelocity(physVelP[0]);
  latticeVelP[1] = this->_converter.getLatticeVelocity(physVelP[1]);
  latticeVelP[2] = this->_converter.getLatticeVelocity(physVelP[2]);
    //std::cout <<  "globic=" << globic << "   latticeVelP=(" << latticeVelP[0] << ", " << latticeVelP[1] << ", " << latticeVelP[2] << ")" << std::endl; // <---

  // Update the smoothing functional
  if ( ! this->_smoothingFunctional.update(physPosP, globic)) {
    std::cout << "ERROR: no lattice point enclosed in particle's kernel length!" << std::endl;
    return false;
  }

  /// Computing dimensionless drag force through smoothed average within the kernel smoothing length
  T latticeForceP[3] = {T(), T(), T()}; // dimensionless force acting on the particle
  for (int i=0; i<this->_smoothingFunctional.getSize(); i++) {

    // Position of the iterated voxel
    int iLatticePosF[3] = {0, 0, 0};
    this->_smoothingFunctional.getLatticePos(iLatticePosF, i);
    // Physical
    T iPhysPosF[3] = {T(), T(), T()};
    this->_sLattice.getCuboidGeometry().get(globic).getPhysR (
           iPhysPosF, iLatticePosF );

    // Fluid velocity at the iterated voxel
    T iPhysVelF[3] = {T(), T(), T()}; // Physical
    this->_interpLatticeVelocity->operator() (iPhysVelF, iPhysPosF, globic);
    //std::cout <<  "globic=" << globic << "   iPhysVelF=(" << iPhysVelF[0] << ", " << iPhysVelF[1] << ", " << iPhysVelF[2] << ")" << std::endl; // <---
    // Lattice
    T iLatticeVelF[3] = {T(), T(), T()}; // Lattice's dimensionless velocity at particle's location
    iLatticeVelF[0] = this->_converter.getLatticeVelocity(iPhysVelF[0]);
    iLatticeVelF[1] = this->_converter.getLatticeVelocity(iPhysVelF[1]);
    iLatticeVelF[2] = this->_converter.getLatticeVelocity(iPhysVelF[2]);

    // Computing fluid-particle momentum transfer
    T gF[3] = {T(), T(), T()}; // force density gF
    this->_momentumExchange->operator() ( gF, iLatticeVelF, latticeVelP,
                      physPosP, iLatticePosF, globic);
    //std::cout <<  "globic=" << globic << "   gF=(" << gF[0] << ", " << gF[1] << ", " << gF[2] << ")" << std::endl; // <---

    // Computing drag coefficient
    T Cd = this->_dragModel(p, iLatticeVelF, latticeVelP, globicFull);
    //std::cout <<  "globic=" << globic << "   Cd=" << Cd << std::endl; //<---

    /// Computing drag force in dimensionless units
    T latticePRad = p->getRad() / _converter.getConversionFactorLength();
    latticeForceP[0] += .5 * Cd * M_PI*pow(latticePRad,2) * gF[0] * (iLatticeVelF[0] - latticeVelP[0])
                      * this->_smoothingFunctional.getWeight(i);
    latticeForceP[1] += .5 * Cd * M_PI*pow(latticePRad,2) * gF[1] * (iLatticeVelF[1] - latticeVelP[1])
                      * this->_smoothingFunctional.getWeight(i);
    latticeForceP[2] += .5 * Cd * M_PI*pow(latticePRad,2) * gF[2] * (iLatticeVelF[2] - latticeVelP[2])
                      * this->_smoothingFunctional.getWeight(i);
  }

  /// Computing physical drag force
  std::vector<T> physForceP(3, T()); // physical force acting on the particle
  physForceP[0] = latticeForceP[0] * this->_converter.getConversionFactorForce();
  physForceP[1] = latticeForceP[1] * this->_converter.getConversionFactorForce();
  physForceP[2] = latticeForceP[2] * this->_converter.getConversionFactorForce();

  /// Updating the particle
  p->setStoreForce(physForceP);
  //std::cout <<  "globic=" << globic << "   physForceP=(" << physForceP[0] << ", " << physForceP[1] << ", " << physForceP[2] << ")" << std::endl; // <---
  return true;
}


////////////////////// Class NaiveNonLocalForwardCouplingModel ////////////////////////

template<typename T, typename Lattice, template<typename V> class Particle>
NaiveNonLocalForwardCouplingModel<T,Lattice,Particle>::NaiveNonLocalForwardCouplingModel (
                     UnitConverter<T, Lattice>& converter,
                     SuperLattice3D<T, Lattice>& sLattice,
                     SuperGeometry3D<T>& sGeometry,
                     DragModel<T, Particle>& dragModel,
                     SmoothingFunctional<T, Lattice>& smoothingFunctional )
         : NonLocalBaseForwardCouplingModel<T,Lattice,Particle>(converter, sLattice, sGeometry, dragModel, smoothingFunctional)
{
  this->_momentumExchange = std::make_shared<NaiveMomentumExchange<T, Lattice> > (
                            this->_converter, this->_sLattice, this->_interpLatticeDensity );
}


}

#endif
