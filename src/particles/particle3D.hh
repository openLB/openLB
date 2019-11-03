/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn
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

#ifndef PARTICLE_3D_HH
#define PARTICLE_3D_HH

#include <string>
#include <iostream>
#include <set>
#include <vector>
#include <list>
#include <deque>

#include "particle3D.h"

namespace olb {

template<typename T>
Particle3D<T>::Particle3D()
  : _pos(3, 0.),
    _vel(3, 0.),
    _force(3, 0.),
    _mas(0),
    _rad(0),
    _cuboid(0),
    _id(0),
    _active(false),
    _storeForce(3, 0.)
{
}

template<typename T>
Particle3D<T>::Particle3D(std::vector<T> pos, T mas, T rad, int id)
  : _pos(pos),
    _vel(3, 0.),
    _force(3, 0.),
    _mas(mas),
    _rad(rad),
    _cuboid(0),
    _id(id),
    _active(true),
    _storeForce(3, 0.)
{
  _invMas = 1. / _mas;
  // RK4
//  _positions = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
//  _velocities = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
//  _forces = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
}

template<typename T>
Particle3D<T>::Particle3D(const Particle3D<T>& p)
  : _pos(p._pos),
    _vel(p._vel),
    _force(p._force),
    _mas(p._mas),
    _rad(p._rad),
    _cuboid(p._cuboid),
    _id(p._id),
    _active(p._active),
    _storeForce(p._storeForce)
{
  _invMas = 1. / _mas;
  // RK4
//  _positions = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
//  _velocities = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
//  _forces = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
}

template<typename T>
Particle3D<T>::Particle3D(std::vector<T> pos, std::vector<T> vel, T mas, T rad, int id)
  : _pos(pos),
    _vel(vel),
    _force(12, 0.),
    _mas(mas),
    _rad(rad),
    _cuboid(0),
    _id(id),
    _active(true),
    _storeForce(3, 0.)
{
  _vel.resize(12, 0.);
  // RK4
//  _positions = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
//  _velocities = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
//  _forces = std::vector<std::vector<T> > (4, std::vector<T> (3, T() ));
}

template<typename T>
inline void Particle3D<T>::addForce(std::vector<T>& force)
{
  for (int i = 0; i < 3; i++) {
    _force[i] += force[i];
  }
}
// set and get force
template<typename T>
inline void Particle3D<T>::setForce(std::vector<T>& force)
{
  _force = force;
}
template<typename T>
inline void Particle3D<T>::resetForce()
{
  for (int i = 0; i < 3; i++) {
    _force[i] = 0.;
  }
}

// set and get storedForce
template<typename T>
inline void Particle3D<T>::setStoreForce(std::vector<T>& storeForce)
{
  for (int i = 0; i < 3; i++) {
    _storeForce[i] = storeForce[i];
  }
}

template<typename T>
inline void Particle3D<T>::resetStoreForce()
{
  for (int i = 0; i < 3; i++) {
    _storeForce[i] = T(0);
  }
}

template<typename T>
void Particle3D<T>::serialize(T serial[])
{
  for (int i = 0; i < 3; i++) {
    serial[i] = _pos[i];
    serial[i + 3] = _vel[i];
    serial[i + 6] = _force[i];
  }
  serial[9] = _mas;
  serial[10] = _rad;
  serial[11] = _cuboid;
  serial[12] = _active;
  serial[13] = _id;

  for (int i = 0; i < 3; i++) {
    serial[i + 14] = _storeForce[i];
  }

  // for (int i = 0; i < 17; i++) {
  //   cout << "serialize " << i << ": " << serial[i]  << " tn: " << typeid(serial[i]).name() << endl;
  // }
}

template<typename T>
void Particle3D<T>::unserialize(T* data)
{
  for (int i = 0; i < 3; i++) {
    _pos[i] = data[i];
    _vel[i] = data[i + 3];
    _force[i] = data[i + 6];
  }
  _mas = data[9];
  _rad = data[10];
  _cuboid = int(data[11]);
  _active = data[12];
  _invMas = 1. / _mas;
  _id = data[13];

  for (int i = 0; i < 3; i++) {
    _storeForce[i] = data[i + 14];
  }

  // for (int i = 0; i < 17; i++) {
  //   cout << "unserialize " << i << ": " << data[i] << " tn: " << typeid(data[i]).name() << endl;
  // }
}

template<typename T>
void Particle3D<T>::print()
{
  std::cout << "Pos " << _pos[0] << " " << _pos[1] << " " << _pos[2] << " " << "Vel "
            << _vel[0] << " " << _vel[1] << " " << _vel[2] << " " << _cuboid
            << std::endl;
}

template<typename T>
void Particle3D<T>::printDeep(std::string message)
{
    std::cout << message <<  " Particle ID " << this->getID()
              << " rad "         << this->getRad()
              << " mass "        << this->getMass()
              << " invMass "     << this->getInvMass()
              << " force "       << this->getForce()[0] << " " << this->getForce()[1] << " " << this->getForce()[2]
              << " storeForce "  << this->getStoreForce()[0] << " " << this->getStoreForce()[1] << " " << this->getStoreForce()[2]
              << " ";
        this->print();
}

template<typename T>
ElParticle3D<T>::ElParticle3D()
  : Particle3D<T>(),
    _charge(1.)
{
}

template<typename T>
ElParticle3D<T>::ElParticle3D(std::vector<T> pos, T mas, T rad, T charge)
  : Particle3D<T>(pos, mas, rad),
    _charge(charge)
{
}

template<typename T>
ElParticle3D<T>::ElParticle3D(std::vector<T> pos, std::vector<T> vel, T mas,
                              T rad, T charge)
  : Particle3D<T>(pos, vel, mas, rad),
    _charge(charge)
{
}

template<typename T>
ElParticle3D<T>::ElParticle3D(const ElParticle3D<T>& p)
  : Particle3D<T>(p),
    _charge(p._charge)
{
}

template<typename T>
void ElParticle3D<T>::serialize(T serial[])
{
  serial[0] = this->_pos[0];
  serial[1] = this->_pos[1];
  serial[2] = this->_pos[2];
  serial[3] = this->_vel[0];
  serial[4] = this->_vel[1];
  serial[5] = this->_vel[2];
  serial[6] = this->_rad;
  serial[7] = this->_mas;
  serial[8] = (double) this->_active;
  serial[9] = _charge;
}

template<typename T>
void ElParticle3D<T>::unserialize(T* data)
{
  this->_pos[0] = data[0];
  this->_pos[1] = data[1];
  this->_pos[2] = data[2];
  this->_vel[0] = data[3];
  this->_vel[1] = data[4];
  this->_vel[2] = data[5];
  this->_rad = data[6];
  this->_mas = data[7];
  this->_active = (bool) data[8];
  _charge = data[9];
}

template<typename T>
AggParticle3D<T>::AggParticle3D()
  : Particle3D<T>::Particle3D()
{
  _aggl = false;
}

template<typename T>
AggParticle3D<T>::AggParticle3D(std::vector<T> pos, T mas, T rad)
  : Particle3D<T>::Particle3D(pos, mas, rad)
{
  _aggl = false;
}

template<typename T>
AggParticle3D<T>::AggParticle3D(const Particle3D<T>& p)
  : Particle3D<T>::Particle3D(p)
{
  _aggl = false;
}

template<typename T>
AggParticle3D<T>::AggParticle3D(std::vector<T> pos, std::vector<T> vel, T mas,
                                T rad)
  : Particle3D<T>::Particle3D(pos, vel, mas, rad)
{
  _aggl = false;
}

template<typename T>
void AggParticle3D<T>::serialize(T serial[])
{
  for (int i = 0; i < 3; i++) {
    serial[i] = this->_pos[i];
    serial[i + 3] = this->_vel[i];
    serial[i + 6] = this->_force[i];
  }
  serial[9] = this->_mas;
  serial[10] = this->_rad;
  serial[11] = this->_cuboid;
  serial[12] = (double) this->_active;
  serial[13] = (double) _aggl;
}

template<typename T>
void AggParticle3D<T>::unserialize(T* data)
{
  for (int i = 0; i < 3; i++) {
    this->_pos[i] = data[i];
    this->_vel[i] = data[i + 3];
    this->_force[i] = data[i + 6];
  }
  this->_mas = data[9];
  this->_rad = data[10];
  this->_cuboid = int(data[11]);
  this->_active = (bool) data[12];
  _aggl = (bool) data[13];
}

template<typename T>
RotatingParticle3D<T>::RotatingParticle3D()
  : Particle3D<T>::Particle3D(), _aVel(3, T()), _torque(3, T())
{
}

template<typename T>
RotatingParticle3D<T>::RotatingParticle3D(std::vector<T> pos, T mas, T rad)
  : Particle3D<T>::Particle3D(pos, mas, rad), _aVel(3, T()), _torque(3, T())
{
}

template<typename T>
RotatingParticle3D<T>::RotatingParticle3D(const RotatingParticle3D<T>& p)
  : Particle3D<T>::Particle3D(p), _aVel(p.getAVel()), _torque(p._torque)
{
}

template<typename T>
RotatingParticle3D<T>::RotatingParticle3D(std::vector<T> pos, std::vector<T> vel, T mas,
    T rad)
  : Particle3D<T>::Particle3D(pos, vel, mas, rad), _aVel(3, T()), _torque(3, T())
{
}


template<typename T>
void RotatingParticle3D<T>::serialize(T serial[])
{
  for (int i = 0; i < 3; i++) {
    serial[i] = this->_pos[i];
    serial[i + 3] = this->_vel[i];
    serial[i + 6] = this->_force[i];
  }
  serial[9] = this->_mas;
  serial[10] = this->_rad;
  serial[11] = this->_cuboid;
  serial[12] = (double) this->_active;
  serial[13] = (double) _aVel[0];
  serial[14] = (double) _aVel[1];
  serial[15] = (double) _aVel[2];
  serial[16] = (double) _torque[0];
  serial[17] = (double) _torque[1];
  serial[18] = (double) _torque[2];
}

template<typename T>
void RotatingParticle3D<T>::unserialize(T* data)
{
  for (int i = 0; i < 3; i++) {
    this->_pos[i] = data[i];
    this->_vel[i] = data[i + 3];
    this->_force[i] = data[i + 6];
  }
  this->_mas = data[9];
  this->_rad = data[10];
  this->_cuboid = int(data[11]);
  this->_active = (bool) data[12];
  _aVel[0] = (bool) data[13];
  _aVel[1] = (bool) data[14];
  _aVel[2] = (bool) data[15];
  _torque[0] = (bool) data[16];
  _torque[1] = (bool) data[17];
  _torque[2] = (bool) data[18];
}


template<typename T>
MagneticParticle3D<T>::MagneticParticle3D()
  : Particle3D<T>::Particle3D(), _dMoment(3, T()), _aVel(3, T()), _torque(3, T()), _magnetisation(T())
{

}

template<typename T>
MagneticParticle3D<T>::MagneticParticle3D(const MagneticParticle3D<T>& p)
  : Particle3D<T>::Particle3D(p), _dMoment(p._dMoment), _aVel(p._aVel), _torque(p._torque), _magnetisation(p._magnetisation), _sActivity(p._sActivity) 
{
}

template<typename T>
MagneticParticle3D<T>::MagneticParticle3D(std::vector<T> pos, std::vector<T> vel, T mas,
    T rad, int id)
  : Particle3D<T>::Particle3D(pos, vel, mas, rad), _dMoment(3, T()), _aVel(3, T()), _torque(3, T()), _magnetisation(T())
{
}

template<typename T>
MagneticParticle3D<T>::MagneticParticle3D(std::vector<T> pos, std::vector<T> vel, T mas, T rad, int id,
    std::vector<T> dMoment, std::vector<T> aVel, std::vector<T> torque, T magnetisation)
  : Particle3D<T>::Particle3D(pos, vel, mas, rad, id),
    _dMoment(dMoment), _aVel(aVel), _torque(torque), _magnetisation(magnetisation)
{
}

template<typename T>
MagneticParticle3D<T>::MagneticParticle3D(std::vector<T> pos, std::vector<T> vel, T mas, T rad, int id,
    std::vector<T> dMoment, std::vector<T> aVel, std::vector<T> torque, T magnetisation, int sActivity)
  : Particle3D<T>::Particle3D(pos, vel, mas, rad, id),
    _dMoment(dMoment), _aVel(aVel), _torque(torque), _magnetisation(magnetisation), _sActivity(sActivity)
{
}

template<typename T>
inline void MagneticParticle3D<T>::resetTorque()
{
  for (int i = 0; i < 3; i++) {
    _torque[i] = 0.;
  }
}

template<typename T>
inline void MagneticParticle3D<T>::setMoment(std::vector<T> moment)
{
  _dMoment = moment;
//  std::cout<< "Setting moment: "<< _dMoment[0] << " " << _dMoment[1] << " " <<_dMoment[2] << std::endl;
}

template<typename T>
inline void MagneticParticle3D<T>::setAVel(std::vector<T> aVel)
{
  _aVel = aVel;
}

template<typename T>
inline void MagneticParticle3D<T>::setTorque(std::vector<T> torque)
{
  _torque = torque;
  //std::cout<< "Setting torque: "<< _torque[0] << " " << _torque[1] << " " <<_torque[2] << std::endl;
}

template<typename T>
inline void MagneticParticle3D<T>::setMagnetisation(T magnetisation)
{
  _magnetisation = magnetisation;
  //std::cout<< "Setting magnetisation: "<< _magnetisation << std::endl;
}

template<typename T>
inline void MagneticParticle3D<T>::setSActivity(int sActivity)
{
  _sActivity = sActivity;
}

template<typename T>
inline void MagneticParticle3D<T>::setAggloItr(typename std::deque<std::list<MagneticParticle3D<T>*>>::iterator aggloItr)
{
  _aggloItr = aggloItr;
}



template<typename T>
void MagneticParticle3D<T>::serialize(T serial[])
{

  for (int i = 0; i < 3; i++) {
    serial[i] = this->_pos[i];
    serial[i + 3] = this->_vel[i];
    serial[i + 6] = this->_force[i];
  }
  serial[9] = this->_mas;
  serial[10] = this->_rad;
  serial[11] = (double) this->_cuboid;
  serial[12] = (double) this->_active;
  serial[13] = (double) this->_id;
  for (int i = 0; i < 3; i++) {
    serial[i + 14] = this->_storeForce[i];
    serial[i + 17] = _dMoment[i];
    serial[i + 20] = _aVel[i];
    serial[i + 23] = _torque[i];
  }
  serial[26] = _magnetisation;
  serial[27] = (double) _sActivity;

}

template<typename T>
void MagneticParticle3D<T>::unserialize(T* data)
{

  for (int i = 0; i < 3; i++) {
    this->_pos[i] = data[i];
    this->_vel[i] = data[i + 3];
    this->_force[i] = data[i + 6];
  }
  this->_mas = data[9];
  this->_rad = data[10];
  this->_cuboid = (int) data[11];
  this->_active = (bool) data[12];
  this->_id = (int) data[13];
  for (int i = 0; i < 3; i++) {
    this->_storeForce[i] = data[i + 14];
    _dMoment[i] = data[i + 17];
    _aVel[i] = data[i + 20];
    _torque[i] = data[i + 23];
  }
  _magnetisation = data[26];
  _sActivity = (int) data[27];

}

}

#endif /* PARTICLE_3D_HH */
