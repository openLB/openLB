/*
 *  Copyright (C) 2015 Marie-Luise Maier, Mathias J. Krause, Sascha Janz
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

/** Alberto Di Renzo, Francesco Paolo Di Maio:
 * "Comparison of contact-force models for the simulation of collisions in
 * DEM-based granular ow codes",
 * Chemical Engineering Science 59 (2004) 525 - 541
 */

#ifndef MagneticForceForMagP3D_HH
#define MagneticForceForMagP3D_HH

#include <cmath>
#include <vector>
#include "magneticForceForMagP3D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
MagneticForceForMagP3D<T, PARTICLETYPE, DESCRIPTOR>::MagneticForceForMagP3D(
  AnalyticalF3D<T, T>& getMagForce, AnalyticalF3D<T, T>& getMagField, T scale) :
  Force3D<T, PARTICLETYPE>(), _getMagForce(getMagForce), _getMagField(getMagField), _scale(scale)
{ }

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
void MagneticForceForMagP3D<T, PARTICLETYPE, DESCRIPTOR>::applyForce(
  typename std::deque<PARTICLETYPE<T> >::iterator p, int pInt,
  ParticleSystem3D<T, PARTICLETYPE>& pSys)
{

  T m_p = p->getMagnetisation();
  T mu_0 = 4 * 3.14159265e-7;
  T mu_i = 4. / 3.*M_PI * pow(p->getRad(), 3) * m_p; // norm mag. dipole moment

  T pos[3] = { T(), T(), T() };
  pos[0] = p->getPos()[0];
  pos[1] = p->getPos()[1];
  pos[2] = p->getPos()[2];

  T forceHelp[3] = { T(), T(), T() };
  _getMagForce(forceHelp, pos);
  T fieldHelp[3] = { T(), T(), T() };
  _getMagField(fieldHelp, pos);

  Vector<T, 3> dMom(p->getMoment()); // orientation vector mag. dipole moment
  dMom *= mu_i; // vector mag. dipole moment

  Vector<T, 3> fieldVec(fieldHelp); // H-Field
  Vector<T, 3> trq = crossProduct3D(dMom, fieldVec); // T = mu_0 mu x H
  trq *= mu_0;

  p->getTorque()[0] += trq[0];
  p->getTorque()[1] += trq[1];
  p->getTorque()[2] += trq[2];
  p->getForce()[0] += forceHelp[0] * _scale;
  p->getForce()[1] += forceHelp[1] * _scale;
  p->getForce()[2] += forceHelp[2] * _scale;

}

}
#endif
