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

/* Helper functionals for Lagrangian two-way coupling methods -- header file.
 */

#ifndef LB_TWO_WAY_HELPER_FUNCTIONALS_H
#define LB_TWO_WAY_HELPER_FUNCTIONALS_H

#include "functors/lattice/reductionF3D.h"

namespace olb {

/** Data structure for smoothing functionals.
  * Stores the lattice position of a cell within smoothing kernel length
  * and the related multiplicative weight.
  */
template<typename T>
struct LatticePosAndWeight {
  int globic = 0;
  int latticePos[3] = {0, 0, 0};
  T weight = T();
};

/** Abstract class for particle Reynolds number computation within drag model.
  * Its raison d'etre consists of not being templetized in Lattice.
  */
template<typename T, template<typename V> class Particle>
class ParticleReynoldsNumber {
public:
  /// Returns the particle Reynolds number. globicFull = { globic, latticeRoundedP[0, ..., 2] }
  virtual T operator() (Particle<T>* p, T magU, int globicFull[])=0;
  /// Destructor
  virtual ~ParticleReynoldsNumber() {};
protected:
  T _RePmin = 0.01;
};

/** Abstract class for particle Reynolds number computation within drag model.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class ParticleReynoldsNumberBase : public ParticleReynoldsNumber<T, Particle> {
public:
  /// Destructor
  virtual ~ParticleReynoldsNumberBase() {};
protected:
  /// Constructor
  ParticleReynoldsNumberBase(UnitConverter<T, Lattice>& converter);
  UnitConverter<T, Lattice>& _converter; // reference to a UnitConverter
};

/** Class class for Newtonian particle Reynolds number computation within drag model.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class NewtonianParticleReynoldsNumber: public ParticleReynoldsNumberBase<T,Lattice,Particle> {
public:
  /// Constructor
  NewtonianParticleReynoldsNumber(UnitConverter<T, Lattice>& converter);
  /// Destructor
  ~NewtonianParticleReynoldsNumber() {};
  /// Returns the particle Reynolds number. globicFull = { globic, latticeRoundedP[0, ..., 2] }
  virtual T operator() (Particle<T>* p, T magU, int globicFull[]) override;
};

/** Class class for power-law particle Reynolds number computation within drag model.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class PowerLawParticleReynoldsNumber: public ParticleReynoldsNumberBase<T,Lattice,Particle> {
public:
  /// Constructor
  PowerLawParticleReynoldsNumber (
        UnitConverter<T, Lattice>& converter, SuperLattice3D<T, Lattice>& sLattice );
  /// Destructor
  ~PowerLawParticleReynoldsNumber() {};
  /// Returns the particle Reynolds number. globicFull = { globic, latticeRoundedP[0, ..., 2] }
  virtual T operator() (Particle<T>* p, T magU, int globicFull[]) override;
protected:
  SuperLattice3D<T, Lattice>& _sLattice; // reference to a lattice
};

/** Abstact class for all the local forward-coupling models,
  * viz., momentum coupling from fluid to particle.
  * Input parameters in attice units.
  */
template<typename T, typename Lattice>
class TwoWayHelperFunctional {
public:
  /// Computes the momentum transfer from fluid to particle.
  virtual bool operator() ( T gF[], T latticeVelF[], T latticeVelP[],
                            T physPosP[], int latticeRoundedP[],
                            int globic )=0;
  virtual ~TwoWayHelperFunctional();
protected:
  /// Constructor
  TwoWayHelperFunctional ( UnitConverter<T, Lattice>& converter,
                           SuperLattice3D<T, Lattice>& sLattice );
  UnitConverter<T, Lattice>& _converter; // reference to a UnitConverter
  SuperLattice3D<T, Lattice>& _sLattice; // reference to a lattice
  std::shared_ptr<SuperLatticeInterpDensity3Degree3D<T, Lattice> > _interpLatticeDensity;
  std::shared_ptr<SuperLatticeInterpPhysVelocity3D<T, Lattice> > _interpLatticeVelocity ;
};

/** Naive way
  */
template<typename T, typename Lattice>
class NaiveMomentumExchange : public TwoWayHelperFunctional<T, Lattice> {
public:
  /// Constructor
  NaiveMomentumExchange ( UnitConverter<T, Lattice>& converter,
                          SuperLattice3D<T, Lattice>& sLattice,
                          std::shared_ptr<SuperLatticeInterpDensity3Degree3D<T, Lattice> > interpLatticeDensity );
  /// Computes the momentum transfer from fluid to particle.
  virtual bool operator() ( T gF[], T latticeVelF[], T latticeVelP[],
                            T physPosP[], int latticeRoundedP[],
                            int globic ) override;
};

/** Using Ladd mechanism
  */
template<typename T, typename Lattice>
class LaddMomentumExchange : public TwoWayHelperFunctional<T, Lattice> {
public:
  /// Constructor
  LaddMomentumExchange ( UnitConverter<T, Lattice>& converter,
                         SuperLattice3D<T, Lattice>& sLattice,
                         std::shared_ptr<SuperLatticeInterpDensity3Degree3D<T, Lattice> > interpLatticeDensity,
                         std::shared_ptr<SuperLatticeInterpPhysVelocity3D<T, Lattice> > interpLatticeVelocity );
  /// Computes the momentum transfer from fluid to particle.
  virtual bool operator() ( T gF[], T latticeVelF[], T latticeVelP[],
                            T physPosP[], int latticeRoundedP[],
                            int globic ) override;
};

/** Abstact class for all the smoothing functionals.
  */
template<typename T, typename Lattice>
class SmoothingFunctional {
public:
  // Returns the size of _latticePosAndWeight
  int getSize();
  // Returns the lattice position of the i-th element of _latticePosAndWeight
  void getLatticePos(int latticePos[], int i);
  // Returns the globic of the i-th element of _latticePosAndWeight
  int getGlobic(int i);
  // Returns the weight relative to the i-th element of _latticePosAndWeight
  T getWeight(int i);
  // Rebuilds _latticePosAndWeight with the new cells within _kernelLength from the bubble's position
  bool update(T physPosP[], int globic);
protected:
  /// Constructor
  SmoothingFunctional(T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice3D<T, Lattice>& sLattice);
  /// The actual smoothing function
  virtual T smoothingFunction(T delta)=0;
  /// Returns the weight for smoothing.
  virtual T compute(T physPosP[], T physPosL[])=0;
  T _kernelLength; // Kernel's smoothing length.
  UnitConverter<T, Lattice>& _converter;  // reference to a UnitConverter
  SuperLattice3D<T, Lattice>& _sLattice; // reference to a lattice
  // positions and weights of the cells within _kernelLength from bubble's position
  std::deque<LatticePosAndWeight<T> > _latticePosAndWeight;
};

/** Abstact class for all the linear-averaging smoothing functionals.
  */
template<typename T, typename Lattice>
class LinearAveragingSmoothingFunctional : public SmoothingFunctional<T, Lattice> {
protected:
  /// Constructor
  LinearAveragingSmoothingFunctional(T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice3D<T, Lattice>& sLattice);
  /// Returns the weight for smoothing.
  virtual T compute(T physPosP[], T physPosL[]) override;
};

/** Abstact class for all the volume-averaging smoothing functionals.
  */
template<typename T, typename Lattice>
class VolumeAveragingSmoothingFunctional : public SmoothingFunctional<T, Lattice> {
protected:
  /// Constructor
  VolumeAveragingSmoothingFunctional(T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice3D<T, Lattice>& sLattice);
  /// Returns the weight for smoothing.
  virtual T compute(T physPosP[], T physPosL[]) override;
};

/** Smoothing functional as in Deen et al (2004), Chem. Eng. Sci 59.
  */
template<typename T, typename Lattice>
class DeenSmoothingFunctional : public LinearAveragingSmoothingFunctional<T, Lattice> {
public:
  /// Constructor
  DeenSmoothingFunctional(T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice3D<T, Lattice>& sLattice);
protected:
  /// The actual smoothing function
  virtual T smoothingFunction(T delta) override;
};

/** Stepwise smoothing functional.
  */
template<typename T, typename Lattice>
class StepSmoothingFunctional : public VolumeAveragingSmoothingFunctional<T, Lattice> {
public:
  /// Constructor
  StepSmoothingFunctional(T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice3D<T, Lattice>& sLattice);
protected:
  /// The actual smoothing function
  virtual T smoothingFunction(T delta) override;
};

}

#endif
