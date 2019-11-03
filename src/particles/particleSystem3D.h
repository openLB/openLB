/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Thomas Henn
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

#ifndef PARTICLE_SYSTEM_3D_H
#define PARTICLE_SYSTEM_3D_H

// Enables particle collision functions
// CollisionModels: Can be used as an alternative to a mechanic contact force
// #define CollisionModels
// CollisionModelsCombindedWithMechContactForce: Can be used in combination with
// a mechanic contact force
// #define CollisionModelsCombindedWithMechContactForce

#include <set>
#include <vector>
#include <list>
#include <deque>
#include "contactDetection/contactDetection.h"
#include "geometry/superGeometry3D.h"
#include "core/blockLatticeStructure3D.h"
#include "forces/force3D.h"
#include "functors/analytical/analyticalF.h"
#include "boundaries/boundary3D.h"
#include "functors/lattice/latticeFrameChangeF3D.h"
#include "utilities/vectorHelpers.h"
#include "twoWayCouplings/twoWayCouplings3D.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
class SuperLatticeInterpPhysVelocity3D;

template<typename T, template<typename U> class PARTICLETYPE>
class Force3D;
template<typename T, template<typename U> class PARTICLETYPE>
class Boundary3D;
template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D;
template<typename T, template<typename U> class PARTICLETYPE>
class SuperParticleSystem3D;
template<typename T, template<typename U> class PARTICLETYPE>
class SuperParticleSysVtuWriter;
template<typename T>
class SuperParticleSysVtuWriterMag;

template<typename T, template<typename U> class PARTICLETYPE>
class SimulateParticles {
public:
  SimulateParticles(ParticleSystem3D<T, PARTICLETYPE>* ps)
    : _pSys(ps)
  {
  }
  inline void simulate(T dT, bool scale = false)
  {
    _pSys->computeForce();
    _pSys->explicitEuler(dT, scale);
    //_pSys->rungeKutta4(dT);
  }
  inline void simulate(T dT, std::set<int> sActivityOfParticle, bool scale = false) {simulate(dT, scale);}
  inline void simulateWithTwoWayCoupling_Mathias ( T dT,
                                    ForwardCouplingModel<T,PARTICLETYPE>& forwardCoupling,
                                    BackCouplingModel<T,PARTICLETYPE>& backCoupling,
                                    int material, int subSteps = 1, bool scale = false )
  {
    for (int iSubStep=1; iSubStep<=subSteps; iSubStep++) {
      if (! _pSys->executeForwardCoupling(forwardCoupling) ) {
        std::cout << " on substep " << iSubStep << std::endl;
        singleton::exit(1);
      }
      _pSys->computeForce();
      _pSys->explicitEuler(dT/(T)(subSteps), scale);
      //_pSys->rungeKutta4(dT/(T)(subSteps));
    }
    _pSys->executeBackwardCoupling(backCoupling, material);
  }
  inline void simulateWithTwoWayCoupling_Davide ( T dT,
                                    ForwardCouplingModel<T,PARTICLETYPE>& forwardCoupling,
                                    BackCouplingModel<T,PARTICLETYPE>& backCoupling,
                                    int material, int subSteps = 1, bool scale = false )
  {
    for (int iSubStep=1; iSubStep<=subSteps; iSubStep++) {
      if (! _pSys->executeForwardCoupling(forwardCoupling) ) {
        std::cout << " on substep " << iSubStep << std::endl;
        singleton::exit(1);
      }
      _pSys->executeBackwardCoupling(backCoupling, material, subSteps);
      _pSys->computeForce();
      _pSys->explicitEuler(dT/(T)(subSteps), scale);
      //_pSys->rungeKutta4(dT/(T)(subSteps));
    }
  }

private:
  ParticleSystem3D<T, PARTICLETYPE>* _pSys;

};


template<typename T>
class SimulateParticles<T, RotatingParticle3D> {
public:
  SimulateParticles(ParticleSystem3D<T, RotatingParticle3D>* ps)
    : _pSys(ps)
  {
  }
  inline void simulate(T dT, bool scale = false)
  {
    _pSys->computeForce();
    _pSys->explicitEuler(dT, scale);
    _pSys->integrateTorque(dT);
  }
private:
  ParticleSystem3D<T, RotatingParticle3D>* _pSys;
};

template<typename T>
class SimulateParticles<T, MagneticParticle3D> {
public:
  SimulateParticles(ParticleSystem3D<T, MagneticParticle3D>* ps)
    : _pSys(ps)
  {
  }

  inline void simulate(T dT, bool scale = false)
  {
    _pSys->resetMag();
    _pSys->computeForce();
    _pSys->explicitEuler(dT, scale);
    _pSys->integrateTorqueMag(dT);

#ifdef CollisionModels
    _pSys->partialElasticImpact(0.67);
#endif
  }
  inline void simulateWithTwoWayCoupling_Mathias ( T dT,
                                    ForwardCouplingModel<T,MagneticParticle3D>& forwardCoupling,
                                    BackCouplingModel<T,MagneticParticle3D>& backCoupling,
                                    int material, int subSteps = 1, bool scale = false )
  {
    for (int iSubStep=1; iSubStep<=subSteps; iSubStep++) {
      if (! _pSys->executeForwardCoupling(forwardCoupling) ) {
        std::cout << " on substep " << iSubStep << std::endl;
        singleton::exit(1);
      }
      _pSys->computeForce();
      _pSys->explicitEuler(dT/(T)(subSteps), scale);
      //_pSys->rungeKutta4(dT/(T)(subSteps));
    }
    _pSys->executeBackwardCoupling(backCoupling, material);
  }
  inline void simulateWithTwoWayCoupling_Davide ( T dT,
                                    ForwardCouplingModel<T,MagneticParticle3D>& forwardCoupling,
                                    BackCouplingModel<T,MagneticParticle3D>& backCoupling,
                                    int material, int subSteps = 1, bool scale = false )
  {
    for (int iSubStep=1; iSubStep<=subSteps; iSubStep++) {
      if (! _pSys->executeForwardCoupling(forwardCoupling) ) {
        std::cout << " on substep " << iSubStep << std::endl;
        singleton::exit(1);
      }
      _pSys->computeForce();
      _pSys->explicitEuler(dT/(T)(subSteps), scale);
      //_pSys->rungeKutta4(dT/(T)(subSteps));
      _pSys->executeBackwardCoupling(backCoupling, material, subSteps);
    }
  }

//   multiple collision models
  inline void simulate(T dT, std::set<int> sActivityOfParticle , bool scale = false)
  {

    _pSys->resetMag(sActivityOfParticle);
    _pSys->computeForce(sActivityOfParticle);
    _pSys->explicitEuler(dT, sActivityOfParticle, scale);
    _pSys->integrateTorqueMag(dT, sActivityOfParticle);

#ifdef CollisionModels
    _pSys->partialElasticImpact(0.67);
#endif

#ifdef CollisionModelsCombindedWithMechContactForce
    _pSys->partialElasticImpactForCombinationWithMechContactForce(0.67);
#endif
  }

private:
  ParticleSystem3D<T, MagneticParticle3D>* _pSys;

};


template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D {
public:

  /// Default constructor for ParticleSystem
  ParticleSystem3D() = default;
  /// Constructor for ParticleSystem
  ParticleSystem3D(int iGeometry, SuperGeometry3D<T>& superGeometry);
  /// Copy constructor for ParticleSystem
  ParticleSystem3D(const ParticleSystem3D<T, PARTICLETYPE>& pS);
  /// Move constructor for ParticleSystem
  ParticleSystem3D(ParticleSystem3D<T, PARTICLETYPE> && pS);
  /// Destructor for ParticleSystem
  virtual ~ParticleSystem3D() {}

  virtual void simulate(T deltatime, bool scale = false);
  virtual void simulateWithTwoWayCoupling_Mathias ( T dT,
                                    ForwardCouplingModel<T,PARTICLETYPE>& forwardCoupling,
                                    BackCouplingModel<T,PARTICLETYPE>& backCoupling,
                                    int material, int subSteps, bool scale );
  virtual void simulateWithTwoWayCoupling_Davide ( T dT,
                                    ForwardCouplingModel<T,PARTICLETYPE>& forwardCoupling,
                                    BackCouplingModel<T,PARTICLETYPE>& backCoupling,
                                    int material, int subSteps, bool scale );
  // multiple collision models
  virtual void simulate(T deltatime, std::set<int> sActivityOfParticle, bool scale = false);

  void printDeep(std::string message="");

  /// Get number of particles in ParticleSystem
  int size();
  /// Get number of particles including shadow particles in ParticleSystem
  int sizeInclShadow() const;
  /// Get number of active particles in ParticleSystem
  int numOfActiveParticles();
  /// Get number of linked forces in ParticleSystem
  int numOfForces();
  /// Get number of particles in vincinity of material number mat
  int countMaterial(int mat = 1);

  /// Add a particle to ParticleSystem
  void addParticle(PARTICLETYPE<T>& p);
  /// Removes all particles from system
  void clearParticles();
  /// Add a force to ParticleSystem
  void addForce(std::shared_ptr<Force3D<T, PARTICLETYPE> > pF);
  /// Add a boundary to ParticleSystem
  void addBoundary(std::shared_ptr<Boundary3D<T, PARTICLETYPE> > pF);

  /// Get reference to a particle in the ParticleSystem
  /// runs over all particles incl. shadow particles
  PARTICLETYPE<T>& operator[](const int i);
  const PARTICLETYPE<T>& operator[](const int i) const;

  /// Set velocity of all particles to fluid velocity
  template<typename DESCRIPTOR>
  void setVelToFluidVel(SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR> &);

  /// Set particle velocity to analytical velocity (e.g. as inital condition
  void setVelToAnalyticalVel(AnalyticalConst3D<T, T>&);

  /// Set global coordinates and extends of Particlesystem (SI units)
  void setPosExt(std::vector<T> physPos, std::vector<T> physExtend);

  /// Get global coordinates and extends of Particlesystem (SI units)
  const std::vector<T>& getPhysPos();
  const std::vector<T>& getPhysExtend();

  /// Save particle positions to file
  void saveToFile(std::string name);

  /// Compute all forces on particles
  void computeForce();
  // multiple collision models
  void computeForce(std::set<int> sActivityOfParticle) {computeForce();};

  /// Stores old particle positions - is used in ..._ActExt
  void setStoredValues();

  /// Sorts the vector of neighbour Particles by increasing distance
  struct getMinDistPart {
    bool operator() (std::pair<size_t, T> i, std::pair<size_t, T> j) { return (i.second < j.second);}
  } getMinDistPartObj;
  void getMinDistParticle (std::vector<std::pair<size_t, T>> ret_matches);

  /// Compute boundary contact
  void computeBoundary();

  /// Set boundary detection algorithm (for future features)
  void setContactDetection(ContactDetection<T, PARTICLETYPE>& contactDetection);
  ContactDetection<T, PARTICLETYPE>* getContactDetection();

  /// Particle-Fluid interaction for subgrid scale particles
  //  template<template<typename V> class DESCRIPTOR>
  //  void particleOnFluid(BlockLatticeStructure3D<T, DESCRIPTOR>& bLattice,
  //                       Cuboid3D<T>& cuboid, int overlap, T eps,
  //                       BlockGeometryStructure3D<T>& bGeometry);
  //  template<template<typename V> class DESCRIPTOR>
  //  void resetFluid(BlockLatticeStructure3D<T, DESCRIPTOR>& bLattice,
  //                  Cuboid3D<T>& cuboid, int overlap);

  friend class SuperParticleSystem3D<T, PARTICLETYPE>;
  friend class SuperParticleSysVtuWriter<T, PARTICLETYPE>;
  friend class SuperParticleSysVtuWriterMag<T>;
  friend class SimulateParticles<T, PARTICLETYPE>;

  //std::map<T, int> radiusDistribution();

  /// Integration method: explicit Euler
  /// if scale = true, velocity is scaled to maximal velocity
  /// maximal velocity = _superGeometry.getCuboidGeometry().getMaxDeltaR()/dT
  void explicitEuler(T dT, bool scale = false);
  // multiple collision models
  void explicitEuler(T dT, std::set<int> sActivityOfParticle, bool scale = false) {explicitEuler(dT, scale);};
  bool executeForwardCoupling(ForwardCouplingModel<T,PARTICLETYPE>& forwardCoupling);
  bool executeBackwardCoupling(BackCouplingModel<T,PARTICLETYPE>& backwardCoupling, int material, int subSteps=1);

  ContactDetection<T, PARTICLETYPE>* getDetection()
  {
    return _contactDetection;
  }

  int getIGeometry()
  {
    return _iGeometry;
  }
  /// returns deque of pointer to particles (not shadow particles)
  /// contained in a particleSystem3D
  std::deque<PARTICLETYPE<T>*> getParticlesPointer();
  /// returns deque of pointer to all particles (incl. shadow particles)
  /// contained in a particleSystem3D
  std::deque<PARTICLETYPE<T>*> getAllParticlesPointer();
  /// returns deque of pointer to all shadow particles
  /// contained in a particleSystem3D
  std::deque<PARTICLETYPE<T>*> getShadowParticlesPointer();

  /// returns deque of particles (no shadow particles)
  /// contained in a particleSystem3D
  std::deque<PARTICLETYPE<T>> getParticles()
  {
    return _particles;
  }

  /// returns shared pointer of forces
  std::list<std::shared_ptr<Force3D<T, PARTICLETYPE> > > getForcesPointer();

  /// Deque of Lists of agglomerated particles
  std::deque<std::list<PARTICLETYPE<T>*>> _Agglomerates;

protected:
  void integrateTorque(T dT);
  void integrateTorqueMag(T dT) {};
  // multiple collision models
  void integrateTorqueMag(T dT, std::set<int> sActivityOfParticle) {};
  void resetMag() {};
  // multiple collision models
  void resetMag(std::set<int> sActivityOfParticle) {};

  /// Collision models: Todo: enable for paralle mode
  /// Resets existing particle overlaps in the event of a collision
  void setOverlapZero() {};
  /// For the combined use of setOverlapZero() and a mechanic contact force
  void setOverlapZeroForCombinationWithMechContactForce() {};
  /// Resets existing particle overlaps in the event of a collision
  /// and applies the physics of an partial elastic impact
  void partialElasticImpact(T restitutionCoeff) {};
  /// Applies the physics of an partial elastic impact while multiple
  /// particle overlapping only to the particles with the least separation distance
  void partialElasticImpactV2(T restitutionCoeff) {};
  /// For the combined use of partialElasticImpact() and a mechanic contact force
  void partialElasticImpactForCombinationWithMechContactForce(T restitutionCoeff) {};
  
  /// Detects and manages particle agglomerates
  void findAgglomerates() {};
  /// Adds new generated particles to the list of non agglomerated Particles
  void initAggloParticles() {};

  void addShadowParticle(PARTICLETYPE<T>& p);

  mutable OstreamManager clout;
  int _iGeometry = -1;
  SuperGeometry3D<T>& _superGeometry;
  ContactDetection<T, PARTICLETYPE>* _contactDetection;
  SimulateParticles<T, PARTICLETYPE> _sim;

  std::deque<PARTICLETYPE<T> > _particles;

  std::deque<PARTICLETYPE<T> > _shadowParticles;
  std::list<std::shared_ptr<Force3D<T, PARTICLETYPE> > > _forces;
  std::list<std::shared_ptr<Boundary3D<T, PARTICLETYPE> > > _boundaries;

  std::vector<T> _physPos;
  std::vector<T> _physExtend;


  /// Integration methods, each need a special template particle
  void velocityVerlet1(T dT);
  void velocityVerlet2(T dT);
//  void implicitEuler(T dT, AnalyticalF3D<T,T>& getvel);
//  void adamBashforth4(T dT);
//  void predictorCorrector1(T dT);
//  void predictorCorrector2(T dT);
  void rungeKutta4_1(T dt);
  void rungeKutta4_2(T dt);
  void rungeKutta4_3(T dt);
  void rungeKutta4_4(T dt);
  void rungeKutta4(T dT);
  void updateParticleDistribution();
};

// Magnetic particle type
template<>
void ParticleSystem3D<double, MagneticParticle3D>::integrateTorqueMag(double dT);
template<>
void ParticleSystem3D<double, MagneticParticle3D>::integrateTorqueMag(double dT, std::set<int> sActivityOfParticle);
template<>
void ParticleSystem3D<double, MagneticParticle3D>::computeForce();
template<>
void ParticleSystem3D<double, MagneticParticle3D>::computeForce(std::set<int> sActivityOfParticle);
template<>
void ParticleSystem3D<double, MagneticParticle3D>::resetMag();
template<>
void ParticleSystem3D<double, MagneticParticle3D>::resetMag(std::set<int> sActivityOfParticle);
template<>
void ParticleSystem3D<double, MagneticParticle3D>::explicitEuler(double dT, bool scale);
template<>
void ParticleSystem3D<double, MagneticParticle3D>::explicitEuler(double dT, std::set<int> sActivityOfParticle, bool scale);

template<>
void ParticleSystem3D<double, MagneticParticle3D>::setOverlapZero();
template<>
void ParticleSystem3D<double, MagneticParticle3D>::setOverlapZeroForCombinationWithMechContactForce();
template<>
void ParticleSystem3D<double, MagneticParticle3D>::partialElasticImpact(double restitutionCoeff);
template<>
void ParticleSystem3D<double, MagneticParticle3D>::partialElasticImpactV2(double restitutionCoeff);
template<>
void ParticleSystem3D<double, MagneticParticle3D>::partialElasticImpactForCombinationWithMechContactForce(double restitutionCoeff);
template<>
void ParticleSystem3D<double, MagneticParticle3D>::findAgglomerates();
template<>
void ParticleSystem3D<double, MagneticParticle3D>::initAggloParticles();

}  //namespace olb
#endif /* PARTICLE_SYSTEM_3D_H */
