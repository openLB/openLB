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
 *
 * validation paper for contact:
 * H. Kruggel-Emden, E. Simsek, S. Rickelt, S. Wirtz, V. Scherer: Review and
 * extension of normal force models for the Discrete Element Method, Powder
 * Technology 171 (2007) 157-173
 */

#ifndef MagneticForceForMagP3D_H
#define MagneticForceForMagP3D_H

#include <cmath>
#include "functors/analytical/analyticalBaseF.h"
#include "particles/particleSystem3D.h"
#include "force3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D;

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
class MagneticForceForMagP3D: public Force3D<T, PARTICLETYPE> {

public:

  MagneticForceForMagP3D(AnalyticalF3D<T, T>& getMagForce, AnalyticalF3D<T, T>& getMagField, T scale = T(1.));
  ~MagneticForceForMagP3D() override {};
  void applyForce(typename std::deque<PARTICLETYPE<T> >::iterator p, int pInt,
                  ParticleSystem3D<T, PARTICLETYPE>& psSys) override;
private:
  AnalyticalF3D<T, T>& _getMagForce;
  AnalyticalF3D<T, T>& _getMagField;
  T _scale;
};

}

#endif
