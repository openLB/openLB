/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Marie-Luise Maier, Mathias J. Krause
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
#ifndef InterpMagForceForMagP3D_H
#define InterpMagForceForMagP3D_H

#include <cmath>
#include "functors/lattice/superLatticeLocalF3D.h"
#include "particles/particleSystem3D.h"
#include "force3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D;

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
class InterpMagForceForMagP3D : public Force3D<T, PARTICLETYPE> {
  
public:
  InterpMagForceForMagP3D(T scale = T(1), T scaleT = T(1));
  InterpMagForceForMagP3D(InterpMagForceForMagP3D<T, PARTICLETYPE, DESCRIPTOR>& f);

  ~InterpMagForceForMagP3D() override {};
  void applyForce(typename std::deque<PARTICLETYPE<T> >::iterator p,
                  int pInt, ParticleSystem3D<T, PARTICLETYPE>& psSys) override ;

private:
  T _scale;
  T _scaleT;
  T _sRad2;
};

}

#endif // INTERPMAGFORCEFORMAGP3D_H_
