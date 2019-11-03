/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Mathias J. Krause, Marie-Luise Maier
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

#include "dynamics/latticeDescriptors.h"
 
#include "functors/lattice/superLatticeLocalF3D.h"
#include "functors/lattice/superLatticeLocalF3D.hh"
#include "particle3D.h"
#include "particle3D.hh"
#include "superParticleSystem3D.h"
#include "superParticleSystem3D.hh"

namespace olb {

template class SuperParticleSystem3D<double, Particle3D>;

template class SuperParticleSystem3D<double, MagneticParticle3D>;

template<>
template<>
void SuperParticleSystem3D<double, Particle3D>::
setVelToFluidVel<descriptors::D3Q19<>>(
  SuperLatticeInterpPhysVelocity3D<double, descriptors::D3Q19<>>& fVel)
{
  for (auto pS : _pSystems) {
    pS->setVelToFluidVel(fVel);
  }
};

template<>
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::
setVelToFluidVel<descriptors::D3Q19<>>(
  SuperLatticeInterpPhysVelocity3D<double, descriptors::D3Q19<>>& fVel)
{
  for (auto pS : _pSystems) {
    pS->setVelToFluidVel(fVel);
  }
};

#ifndef OLB_PRECOMPILED
template<>
template<>
void SuperParticleSystem3D<double, Particle3D>::
setVelToFluidVel<descriptors::D3Q19<descriptors::FORCE>>(
  SuperLatticeInterpPhysVelocity3D<double, descriptors::D3Q19<descriptors::FORCE>>& fVel)
{
  for (auto pS : _pSystems) {
    pS->setVelToFluidVel(fVel);
  }
};

template<>
template<>
void SuperParticleSystem3D<double, MagneticParticle3D>::
setVelToFluidVel<descriptors::D3Q19<descriptors::FORCE>>(
  SuperLatticeInterpPhysVelocity3D<double, descriptors::D3Q19<descriptors::FORCE>>& fVel)
{
  for (auto pS : _pSystems) {
    pS->setVelToFluidVel(fVel);
  }
};
#endif

}  // namespace olb
