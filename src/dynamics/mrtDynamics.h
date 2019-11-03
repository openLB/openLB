/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

/** \file
 * This object is a MRT LB dynamics as described in D.Yu et al. in
 * Progress in Aerospace Sciences 39 (2003) 329-367
 */
#ifndef MRT_DYNAMICS_H
#define MRT_DYNAMICS_H

#include "dynamics/dynamics.h"

namespace olb {

template<typename T, typename DESCRIPTOR> class Cell;


/// Implementation of the entropic collision step
template<typename T, typename DESCRIPTOR>
class MRTdynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  MRTdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Compute all equilibrium moments
  void computeAllEquilibrium(T momentaEq[DESCRIPTOR::q], T rho,
                          const T u[DESCRIPTOR::d], const T uSqr);
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega_) override;
  /// Get local relaxation parameter of the dynamics
  T getLambda() const;
  /// Set local relaxation parameter of the dynamics
  void setLambda(T lambda_);
protected:
  T invM_S[DESCRIPTOR::q][DESCRIPTOR::q]; // relaxation times matrix.
  T omega; // the shear viscosity relaxation time
  T lambda;// the bulk viscosity relaxation time
};

/// Implementation of the entropic collision step
template<typename T, typename DESCRIPTOR>
class ForcedMRTdynamics : public MRTdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ForcedMRTdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Clone the object on its dynamic type.
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);

};

/// Implementation of the entropic collision step
template<typename T, typename DESCRIPTOR>
class MRTdynamics2 : public MRTdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  MRTdynamics2(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Clone the object on its dynamic type.
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);

protected:
  T invM_S_2[DESCRIPTOR::q][DESCRIPTOR::q]; // relaxation times matrix.
  T omega;
};



}

#endif
