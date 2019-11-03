/*
 *  Copyright (C) 2018 Marie-Luise Maier
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

#ifndef MagneticForceFromHField3D_H
#define MagneticForceFromHField3D_H

#include <cmath>
#include "functors/lattice/superLatticeLocalF3D.h"
#include "functors/analytical/interpolationF3D.h"
#include "particles/particleSystem3D.h"
#include "force3D.h"

namespace olb {

/*
template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D;

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
class MagneticForceFromHField3D: public Force3D<T, PARTICLETYPE> {

public:
  MagneticForceFromHField3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
      SuperLatticeField3D<T, DESCRIPTOR>& sLatticeExField, T Mp, T scale=T(1.) );
  ~MagneticForceFromHField3D() override {};
  void applyForce(typename std::deque<PARTICLETYPE<T> >::iterator p, int pInt,
                  ParticleSystem3D<T, PARTICLETYPE>& psSys) override;
private:
  SuperLattice3D<T, DESCRIPTOR>& _sLattice;
  SuperLatticeField3D<T, DESCRIPTOR>& _sLatticeHField;
  T _Mp;
  T _scale;
};
*/

}

#endif
