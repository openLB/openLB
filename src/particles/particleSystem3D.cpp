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

#include "particle3D.h"
#include "particle3D.hh"
#include "particleSystem3D.h"
#include "particleSystem3D.hh"
#include "dynamics/latticeDescriptors.h"
 
#include "functors/lattice/superLatticeLocalF3D.h"

namespace olb {


template class SimulateParticles<double,Particle3D>;
template class ParticleSystem3D<double,Particle3D>;

template class SimulateParticles<double, MagneticParticle3D>;
template class ParticleSystem3D<double, MagneticParticle3D>;

template<>
template<>
void ParticleSystem3D<double,Particle3D>::
setVelToFluidVel<descriptors::D3Q19<>>(
  SuperLatticeInterpPhysVelocity3D<double, descriptors::D3Q19<>>& fVel)
{
  for (auto& p : _particles) {
    if (p.getActive()) {
      fVel(&p.getVel()[0], &p.getPos()[0], p.getCuboid());
    }
  }
};

template<>
template<>
void ParticleSystem3D<double,MagneticParticle3D>::
setVelToFluidVel<descriptors::D3Q19<>>(
  SuperLatticeInterpPhysVelocity3D<double, descriptors::D3Q19<>>& fVel)
{
  for (auto& p : _particles) {
    if (p.getActive()) {
      fVel(&p.getVel()[0], &p.getPos()[0], p.getCuboid());
    }
  }
};

#ifndef OLB_PRECOMPILED
template<>
template<>
void ParticleSystem3D<double,Particle3D>::
setVelToFluidVel<descriptors::D3Q19<descriptors::FORCE>>(
  SuperLatticeInterpPhysVelocity3D<double, descriptors::D3Q19<descriptors::FORCE>>& fVel)
{
  for (auto& p : _particles) {
    if (p.getActive()) {
      fVel(&p.getVel()[0], &p.getPos()[0], p.getCuboid());
    }
  }
};

template<>
template<>
void ParticleSystem3D<double,MagneticParticle3D>::
setVelToFluidVel<descriptors::D3Q19<descriptors::FORCE>>(
  SuperLatticeInterpPhysVelocity3D<double, descriptors::D3Q19<descriptors::FORCE>>& fVel)
{
  for (auto& p : _particles) {
    if (p.getActive()) {
      fVel(&p.getVel()[0], &p.getPos()[0], p.getCuboid());
    }
  }
};
#endif

}  // namespace olb
