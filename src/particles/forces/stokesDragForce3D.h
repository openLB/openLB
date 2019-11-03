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

/* Formulation of Stokes drag force (6 pi mu r_p v_rel)
 * for actual time step with particle velocity of the last time step
 * with implicit Euler method. For more information see
 * paper of Thomas Henn: "Parallel dilute particulate flow simulations
 * in the human nasal cavity", Computers & Fluids, 2016.
 **/

#ifndef STOKESDRAGFORCE_3D_H
#define STOKESDRAGFORCE_3D_H

#include "functors/lattice/superLatticeLocalF3D.h"
#include "particles/particleSystem3D.h"
#include "force3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D;

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
class StokesDragForce3D : public Force3D<T, PARTICLETYPE> {
  
public:
  /// Constructor, FluidVelocity, physicalTimeStep, physicalDynamicViscosity
  StokesDragForce3D(SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& getVel, T dT, T mu);
  /// Constructor using values from converter
  StokesDragForce3D(SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& getVel, UnitConverter<T, DESCRIPTOR> const& converter);
  /// Constructor using values from converter
  StokesDragForce3D(SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& getVel, UnitConverter<T, DESCRIPTOR> const& converter, T scaleFactor);
  /// Destructor
  ~StokesDragForce3D() override {}
  void applyForce(typename std::deque<PARTICLETYPE<T> >::iterator p,
                  int pInt, ParticleSystem3D<T, PARTICLETYPE>& psSys) override;

  /// Compute Force for subgrid scale particles
  void computeForce(int pInt, ParticleSystem3D<T, PARTICLETYPE>* psSys,
                    T force[3]);
private:
  SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& _getVel;
  T _C1; // helperVariable
  T _mu; // dynamic viscosity
  T _dTinv; // inverse timestep
  T _scaleFactor;
};

}

#endif /* STOKESDRAGFORCE_3D_H */
