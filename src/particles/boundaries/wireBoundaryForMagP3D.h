/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Thomas Henn, Mathias J. Krause
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

#ifndef WireBoundaryForMagP3D_H
#define WireBoundaryForMagP3D_H

#include <set>

#include "geometry/superGeometry3D.h"
#include "particles/particleSystem3D.h"
#include "boundary3D.h"


namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D;

/*
 * Particle boundary based on the fluids material number.
 **/

template<typename T, template<typename U> class PARTICLETYPE>
class WireBoundaryForMagP3D : public Boundary3D<T, PARTICLETYPE> {

public:
  /// Copy constructor
  WireBoundaryForMagP3D(WireBoundaryForMagP3D<T, PARTICLETYPE>& f);
  /// Constructor with set of material numbers
  WireBoundaryForMagP3D(SuperGeometry3D<T>& sg, std::set<int> materials);

  ~WireBoundaryForMagP3D() override {};
  
  /// Add a single material number
  void addMaterial(int mat)
  {
    _materials.insert(mat);
  }
  /// Add several material numbers
  void addMaterial(std::vector<int> mats)
  {
    for (unsigned i = 0; i < mats.size(); ++i) {
      _materials.insert(mats[i]);
    }
  }
  /// Apply the boundary condition
  void applyBoundary(
    typename std::deque<PARTICLETYPE<T> >::iterator& p,
    ParticleSystem3D<T, PARTICLETYPE>& psSys) override;

private:
  SuperGeometry3D<T>& _sg;
  std::set<int> _materials;
  std::set<int>::iterator _matIter;
  //  int x0, y0, z0;
};

}

#endif /* WireBoundaryForMagP3D */
