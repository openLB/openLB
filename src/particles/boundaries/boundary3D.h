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


#ifndef BOUNDARY3D_H
#define BOUNDARY3D_H

#include "particles/particleSystem3D.h"
#include "particles/particleSystem3D.h"
#include "io/ostreamManager.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D;


/**
 * Prototype for all particle boundaries
 */

template<typename T, template<typename U> class PARTICLETYPE>
class Boundary3D {
public:
  Boundary3D();
  Boundary3D(Boundary3D<T, PARTICLETYPE>&);
  Boundary3D(const Boundary3D<T, PARTICLETYPE>&);

  virtual ~Boundary3D() {};
  virtual void applyBoundary(typename std::deque<PARTICLETYPE<T> >::iterator& p, ParticleSystem3D<T, PARTICLETYPE>& psSys)=0;

protected:
  mutable OstreamManager clout;
};

}
#endif /* BOUNDARY3D_H */
