/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Mathias J. Krause
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

#ifndef PARTICLE_3D_H
#define PARTICLE_3D_H

#include <set>
#include <vector>
#include <list>
#include <deque>
#include <string>
#include <iostream>

namespace olb {

template<typename T>
class Particle3D {
public:
  Particle3D();
  Particle3D(std::vector<T> pos, T mas = 1., T rad = 1., int id = 0);
  Particle3D(std::vector<T> pos, std::vector<T> vel,
             T mas = 1., T rad = 1., int id = 0);
  Particle3D(const Particle3D<T>& p);

  inline void setPos(std::vector<T> pos) {
    _pos = pos;
  }
  inline void setStoredPos(std::vector<T> pos) {
    _storePos = pos;
  }
  inline std::vector<T>& getStoredPos() {
    return _storePos;
  }
  inline std::vector<T>& getPos() {
    return _pos;
  }
  inline const std::vector<T>& getPos() const {
    return _pos;
  }
  inline void setVel(std::vector<T> vel) {
    _vel = vel;
  }
  inline void setStoredVel(std::vector<T> vel) {
    _storeVel = vel;
  }
  inline std::vector<T>& getStoredVel() {
    return _storeVel;
  }
  inline int getID() {
    return _id;
  }
  inline void setID(int id) {
    _id = id;
  }
  inline std::vector<T>& getVel() {
    return _vel;
  }
  inline const std::vector<T>& getVel() const {
    return _vel;
  }

  inline void addForce(std::vector<T>& frc);

  inline void setForce(std::vector<T>& frc);
  inline void resetForce();
  inline std::vector<T>& getForce() {
    return _force;
  }
  inline const std::vector<T>& getForce() const {
    return _force;
  }

  inline void setStoreForce(std::vector<T>& storeForce);
  inline void resetStoreForce();
  inline std::vector<T>& getStoreForce()
  {
    return _storeForce;
  }
  inline const std::vector<T>& getStoreForce() const {
    return _storeForce;
  }

  // RK4
//  inline void setPos(std::vector<T> pos, int i) {
//    _positions[i] = pos;
//  }
//  inline std::vector<T>& getPos(int i) {
//    return _positions[i];
//  }
//  inline const std::vector<T>& getPos(int i) const {
//    return _positions[i];
//  }
//  inline void setVel(std::vector<T> vel, int i) {
//    _velocities[i] = vel;
//  }
//  inline std::vector<T>& getVel(int i) {
//    return _velocities[i];
//  }
//  inline const std::vector<T>& getVel(int i) const {
//    return _velocities[i];
//  }
//  inline void setForce(std::vector<T> force, int i) {
//    _forces[i] = force;
//  }
//  inline std::vector<T>& getForce(int i) {
//    return _forces[i];
//  }
//  inline const std::vector<T>& getForce(int i) const {
//    return _forces[i];
//  }

// write particle properties into serial
  void serialize(T serial[]);
  // write data into particle properties
  void unserialize(T*);
  void print();
  void printDeep(std::string message);

  inline const T& getMass() {
    return _mas;
  }

  inline const T& getInvMass() {
    return _invMas;
  }

  inline const T& getMass() const {
    return _mas;
  }
  inline void setMass(T m) {
    _mas = m;
    _invMas = 1. / _mas;
  }
  inline const T& getRad() {
    return _rad;
  }
  inline const T& getRad() const {
    return _rad;
  }
  inline void setRad(T r) {
    _rad = r;
  }
  inline const int& getCuboid() {
    return _cuboid;
  }
  inline void setCuboid(int c) {
    _cuboid = c;
  }
  inline const bool& getActive() {
    return _active;
  }
  inline const bool& getActive() const {
    return _active;
  }
  inline void setActive(bool act) {
    _active = act;
    if (!act) {
      _vel[0] = 0;
      _vel[1] = 0;
      _vel[2] = 0;
    }
  }

  // static const int serialPartSize = 14; // pos, vel, force, mas, rad, cuboid, id, active
  static const int serialPartSize = 17; // pos, vel, force, mas, rad, cuboid, id, active, storeForce
  // static const int serialPartSize = 23; // pos, vel, force, mas, rad, cuboid, id, active, storeForce, storePos, storeVel

  std::vector<std::pair<size_t, T> > _verletList;

protected:
  std::vector<T> _pos;
  std::vector<T> _vel;
  std::vector<T> _force;
  T _invMas;
  T _mas;
  T _rad;
  ///globIC
  int _cuboid;
  int _id;
  bool _active;
  std::vector<T> _storePos;
  std::vector<T> _storeVel;
  std::vector<T> _storeForce;
  // RK4
//  std::vector<std::vector<T> > _positions;
//  std::vector<std::vector<T> > _velocities;
//  std::vector<std::vector<T> > _forces;

};

/*
 * Electric Particles
 */
template<typename T>
class ElParticle3D: public Particle3D<T> {
public:
  ElParticle3D();
  ElParticle3D(std::vector<T> pos, T mas = 1., T rad = 1., T charge = 1.);
  ElParticle3D(std::vector<T> pos, std::vector<T> vel,
               T mas = 1., T rad = 1., T charge = 1.);
  ElParticle3D(const ElParticle3D<T>& p);
  virtual ~ElParticle3D() {
  }
  ;
  virtual void serialize(T serial[]);
  virtual void unserialize(T*);
  static const int serialPartSize = 10;
  T _charge;
};

/*
 * Particles for Agglomeration
 */

template<typename T>
class AggParticle3D: public Particle3D<T> {
public:
  AggParticle3D();
  AggParticle3D(std::vector<T> pos, T mas = 1., T rad = 1.);
  AggParticle3D(std::vector<T> pos, std::vector<T> vel, T mas = 1., T rad = 1.);
  AggParticle3D(const Particle3D<T>& p);

  inline void setMass(T mas) {
    this->_mas = mas;
  }
  inline void setRad(T rad) {
    this->_rad = rad;
  }
  inline const bool& getAggl() {
    return _aggl;
  }
  inline const bool& getAggl() const {
    return _aggl;
  }
  inline void setAggl(bool aggl) {
    _aggl = aggl;
    if (aggl) {
      this->setActive(false);
    }
  }

  static const int serialPartSize = 14;
  void serialize(T serial[]);
  void unserialize(T*);

private:
  bool _aggl;
};

/*
 * Rotating Particles
 */
template<typename T>
class RotatingParticle3D: public Particle3D<T> {
public:
  RotatingParticle3D();
  RotatingParticle3D(std::vector<T> pos, T mas = 1., T rad = 1.);
  RotatingParticle3D(std::vector<T> pos, std::vector<T> vel, T mas = 1., T rad = 1.);
  RotatingParticle3D(const RotatingParticle3D<T>& p);

  static const int serialPartSize = 19;
  void serialize(T serial[]);
  void unserialize(T*);

  inline std::vector<T>& getAVel() {
    return _aVel;
  }
  inline const std::vector<T>& getAVel() const {
    return _aVel;
  }
  inline std::vector<T>& getTorque() {
    return _torque;
  }
  inline const std::vector<T>& getTorque() const {
    return _torque;
  }
private:
  std::vector<T> _aVel;
  std::vector<T> _torque;
};

// NEW
/*
 * Particles for magnetic force
 */

template<typename T>
class MagneticParticle3D: public Particle3D<T> {
private:
  std::vector<T> _dMoment;
  std::vector<T> _aVel;
  std::vector<T> _torque;
  T _magnetisation;

  // _sActivity values range from 1 to 3
  // Gives information about the poisiton of the particle in the simulation area
  // One collision model:
  // 1: Particle is not deposited
  // 2: Not in use
  // 3: Particle is deposited on mag. wire.
  // Multiple collision models:
  // 1: Particle is in sufficiently large distance to the inhom. mag. field of the wire
  //    Particle collsions are managed by a collision model function
  // 2: Particle is within the inhom. mag. field of the wire
  //    Particle collsions are managed by a mechanic contact force
  //    The sActivity change (1 --> 2) is done in explicitEuler()
  // 3: Particle is deposited on mag. wire.
  //    Particle collsions are managed by a mechanic contact force
  int _sActivity;

  // Assigns an agglomerate to the particle
  typename std::deque<std::list<MagneticParticle3D<T>*>>::iterator _aggloItr;
  // Assigns the particle position in an Agglomerate to a Particle
  // typename std::list<MagneticParticle3D<T>*>::iterator _aggloItrPos;

public:

  MagneticParticle3D();
  MagneticParticle3D(std::vector<T> pos, T mas = 1., T rad = 1., int id = 0);
  MagneticParticle3D(std::vector<T> pos, std::vector<T> vel, T mas = 1., T rad =
                       1., int id = 0);
  MagneticParticle3D(const MagneticParticle3D<T>& p);
  MagneticParticle3D(std::vector<T> pos, std::vector<T> vel, T mas, T rad, int id,
                     std::vector<T> dMoment, std::vector<T> aVel, std::vector<T> torque, T magnetisation);
  MagneticParticle3D(std::vector<T> pos, std::vector<T> vel, T mas, T rad, int id,
                     std::vector<T> dMoment, std::vector<T> aVel, std::vector<T> torque, T magnetisation, int sActivity);

  static const int serialPartSize = 28;
  void serialize(T serial[]);
  void unserialize(T*);

  // Set torque zero
  inline void resetTorque();

  // Set and get orientation of magnetic dipolemoment
  inline void setMoment(std::vector<T> moment);
  inline std::vector<T>& getMoment() {
    return _dMoment;
  }
  inline const std::vector<T>& getMoment() const {
    return _dMoment;
  }

  // Set and get angular velocity
  inline void setAVel(std::vector<T> aVel);
  inline std::vector<T>& getAVel() {
    return _aVel;
  }
  inline const std::vector<T>& getAVel() const {
    return _aVel;
  }

  // Set and get torque
  inline void setTorque(std::vector<T> torque);
  inline std::vector<T>& getTorque() {
    return _torque;
  }
  inline const std::vector<T>& getTorque() const {
    return _torque;
  }

  // Set and get magnetisation
  inline void setMagnetisation(T magnetisation);
  inline T& getMagnetisation() {
    return _magnetisation;

  }
  inline const T& getMagnetisation() const {
    return _magnetisation;
  }

  // Set and get sActivity
  inline void setSActivity(int sActivity);
  inline int& getSActivity() {
    return _sActivity;
  }

  // Set and get aggloItr
  inline void setAggloItr(typename std::deque<std::list<MagneticParticle3D<T>*>>::iterator aggloItr);
  inline typename std::deque<std::list<MagneticParticle3D<T>*>>::iterator& getAggloItr() {
    return _aggloItr;
  }

};
}
#endif /* PARTICLE_3D_H */

