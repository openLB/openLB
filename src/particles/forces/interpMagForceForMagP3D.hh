/*  This file is part of the OpenLB library
 *
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

/// Magnetic field that creates magnetization in wire has to be orthogonal to the wire.
/// to calculate the magnetic force on particles around a cylinder
/// (J. Lindner et al. / Computers and Chemical Engineering 54 (2013) 111-121)
#ifndef InterpMagForceForMagP3D_HH
#define InterpMagForceForMagP3D_HH

#include <cmath>
#include <vector>
#include "interpMagForceForMagP3D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
InterpMagForceForMagP3D<T, PARTICLETYPE, DESCRIPTOR>::InterpMagForceForMagP3D(T scale, T scaleT)
  : Force3D<T, PARTICLETYPE>(),
    _scale(scale),
    _scaleT(scaleT)
{
//  this->_name = "InterpMagForceForMagP3D";
}

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
void InterpMagForceForMagP3D<T, PARTICLETYPE, DESCRIPTOR>::applyForce(
  typename std::deque<PARTICLETYPE<T> >::iterator p, int pInt,
  ParticleSystem3D<T, PARTICLETYPE>& pSys)
{

  std::vector < std::pair<size_t, T> > ret_matches;
  pSys.getDetection()->getMatches(pInt, ret_matches);
  // std::random_shuffle ( ret_matches.begin(), ret_matches.end() );

  /*const*/ PARTICLETYPE<T>* p2 = NULL;
  typename std::vector<std::pair<size_t, T> >::iterator it =
    ret_matches.begin();

  Vector<T, 3> dMom_1 = Vector<T, 3>((p->getMoment()));
  if (dMom_1.norm() > std::numeric_limits<T>::epsilon()) {
    T m_p1 = p->getMagnetisation();

    for (; it != ret_matches.end(); it++) {
      if (it->second >= std::numeric_limits < T > ::epsilon()) {

        p2 = &pSys[it->first];

        // No magnetic force between particles which are in physcal contact
        // if ((p2->getRad() + p->getRad()) <= std::sqrt(it->second)) {

        // get neighbour particle moment
        Vector<T, 3> dMom_2((p2->getMoment()));
        if (dMom_2.norm() > std::numeric_limits<T>::epsilon()) {
          T m_p2 = p2->getMagnetisation();

          // given moment magnitudes as scale factors
          T m_i_scaleFactor = dMom_1.norm();
          T m_j_scaleFactor = dMom_2.norm();

          // normalised moment directions
          Vector<T, 3> n_i(dMom_1);
          n_i.normalize();
          Vector<T, 3> n_j(dMom_2);
          n_j.normalize();

          // constants
          T mu_0 = 4 * M_PI * 1.e-7; // magnetic constant
          T mu_i = 4. / 3 * M_PI * pow(p->getRad(), 3) * m_p1 * m_i_scaleFactor ; // magnetic moment of particle i
          T mu_j = 4. / 3 * M_PI * pow(p2->getRad(), 3) * m_p2 * m_j_scaleFactor ; // of particle j

          //vector from particle1 to particle2
          Vector<T, 3> r_ij;
          //S Das wäre der Richtungsvektor von p2 zu p und nicht von p zu p2, aber die Rechung unten berücksichtigt das wohl, da richtige Ergebnisse
          r_ij[0] = p->getPos()[0] - p2->getPos()[0];
          r_ij[1] = p->getPos()[1] - p2->getPos()[1];
          r_ij[2] = p->getPos()[2] - p2->getPos()[2];

          // distance from particle1 to particle2
          T r = r_ij.norm();

          // normalised direction vector
          Vector<T, 3> t_ij(r_ij);
          t_ij.normalize();

          // FORCE of Dipole 1 on Dipole 2
          Vector<T, 3> mi = {n_i};
          mi *= mu_i;
          Vector<T, 3> mj = {n_j};
          mj *= mu_j;
          Vector<T, 3> rn = {r_ij};
          rn *= -1. ;
          normalize(rn);
          T scalar_termF = (3 * mu_0 ) / (4 * M_PI * std::pow(r, 4));
          Vector<T, 3> force = mj * (mi * rn) + mi * (mj * rn) + rn * (mi * mj) - 5 * rn * (mi * rn) * (mj * rn) ;
          force *= scalar_termF;

          // TORQUE of Dipole 1 on Dipole 2
          // T_ij = -[(mu_0*mu_i*mu_j)/(4*M_PI*r^3)][n_i x n_j - 3(n_j * t_ij)n_i x t_ij]
          T scalar_termT = - (mu_0 * mu_i * mu_j) / (4 * M_PI * std::pow(r, 3));
          Vector<T, 3> torque = crossProduct3D(n_i, n_j);
          torque -= crossProduct3D( 3 * (n_j * t_ij) * n_i, t_ij);
          torque *= scalar_termT;

          // Force an torque for overlapping particles 
          if ((p2->getRad() + p->getRad()) >= std::sqrt(it->second)) {

            normalize(force);
            torque *= 0;
            force *= mu_0 * pow((mu_i + mu_j) / 2., 2.) / (4 * M_PI * pow(r / 2, 4.)) ;
          }

          p->getTorque()[0] += 0.5 * torque[0] * _scaleT;
          p->getTorque()[1] += 0.5 * torque[1] * _scaleT;
          p->getTorque()[2] += 0.5 * torque[2] * _scaleT;
          p->getForce()[0] -= 0.5 * force[0] * _scale ;
          p->getForce()[1] -= 0.5 * force[1] * _scale ;
          p->getForce()[2] -= 0.5 * force[2] * _scale ;
          p2->getTorque()[0] -= 0.5 * torque[0] * _scaleT;
          p2->getTorque()[1] -= 0.5 * torque[1] * _scaleT;
          p2->getTorque()[2] -= 0.5 * torque[2] * _scaleT;
          p2->getForce()[0] += 0.5 * force[0] * _scale ;
          p2->getForce()[1] += 0.5 * force[1] * _scale ;
          p2->getForce()[2] += 0.5 * force[2] * _scale ;

        }
        // }
      }
    }
  }
}

}

#endif // INTERPMAGFORCEFORMAGP3D_H_
