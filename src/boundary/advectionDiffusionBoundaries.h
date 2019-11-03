/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
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

#ifndef ADVECTION_DIFFUSION_BOUNDARIES_H
#define ADVECTION_DIFFUSION_BOUNDARIES_H

#include "dynamics/latticeDescriptors.h"
#include "dynamics/advectionDiffusionDynamics.h"
#include "dynamics/dynamics.h"

namespace olb {

//===================================================================================
//================= AdvectionDiffusionDynamcison Flat Boundaries =========
//===================================================================================
template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
class AdvectionDiffusionBoundariesDynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  AdvectionDiffusionBoundariesDynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step for flat boundary and given rho
  /* Working principle:
   * 1. Compute rho_current by summing up all known f_i
   * 2. Get difference (rho - rho_current) and initialise the unknown f_i
   */
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega_) override;
private:
  Dynamics boundaryDynamics;
};

//===================================================================================
//================= AdvectionDiffusionDynamcis On Edges =========
//===================================================================================
template<typename T, typename DESCRIPTOR, typename Dynamics, int plane, int normal1, int normal2>
class AdvectionDiffusionEdgesDynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  AdvectionDiffusionEdgesDynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega_) override;
private:
  Dynamics boundaryDynamics;
};


//===================================================================================
//================= AdvectionDiffusionDynamics on  Corners for 2D Boundaries =========
//===================================================================================
template<typename T, typename DESCRIPTOR, typename Dynamics, int xNormal, int yNormal>
class AdvectionDiffusionCornerDynamics2D : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  AdvectionDiffusionCornerDynamics2D(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega_) override;
private:
  Dynamics boundaryDynamics;
};

//===================================================================================
//================= AdvectionDiffusionDynamics on  Corners for 3D Boundaries =========
//===================================================================================
template<typename T, typename DESCRIPTOR, typename Dynamics, int xNormal, int yNormal, int zNormal>
class AdvectionDiffusionCornerDynamics3D : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  AdvectionDiffusionCornerDynamics3D(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega_) override;
private:
  Dynamics boundaryDynamics;
};



}  // namespace olb

#endif
