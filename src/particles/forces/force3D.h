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


#ifndef FORCE_3D_H
#define FORCE_3D_H

#include <deque>
#include "io/ostreamManager.h"
#include "particles/particleSystem3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D;


/**
 * Prototype for all particle forces
 */

template<typename T, template<typename U> class PARTICLETYPE>
class Force3D {
public:
  Force3D();
  Force3D(Force3D<T, PARTICLETYPE>&);
  Force3D(const Force3D<T, PARTICLETYPE>&);

  virtual ~Force3D() {};
  virtual void applyForce(typename std::deque<PARTICLETYPE<T> >::iterator p, int pInt, ParticleSystem3D<T, PARTICLETYPE>& psSys)=0;

protected:
  mutable OstreamManager clout;

};

}
#endif /* FORCE3D_H */
