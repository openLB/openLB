/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007, 2017 Orestis Malaspinas, Jonas Latt, Mathias J. Krause
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
 * A collection of entropic dynamics classes (e.g. EntropicEq,
 * ForcedEntropicEq, Entropic, ForcedEntropic) with which a Cell object
 * can be instantiated -- header file.
 *
 * Entropic Modell:
 * Ansumali, Santosh, Iliya V. Karlin, and Hans Christian Öttinger
 * Minimal entropic kinetic models for hydrodynamics
 * EPL (Europhysics Letters) 63.6 (2003): 798
 */

#ifndef ENTROPIC_LB_DYNAMICS_H
#define ENTROPIC_LB_DYNAMICS_H

#include "dynamics/dynamics.h"

namespace olb {

template<typename T, typename DESCRIPTOR> class Cell;


/// Implementation of the entropic collision step
template<typename T, typename DESCRIPTOR>
class EntropicEqDynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  EntropicEqDynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega_) override;
private:
  T omega;  ///< relaxation parameter
};

/// Implementation of the forced entropic collision step
template<typename T, typename DESCRIPTOR>
class ForcedEntropicEqDynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ForcedEntropicEqDynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const;
  /// Collision step
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);
  /// Get local relaxation parameter of the dynamics
  virtual T getOmega() const;
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
private:
  T omega;  ///< relaxation parameter
};

/// Implementation of the entropic collision step

template<typename T, typename DESCRIPTOR>
class EntropicDynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  EntropicDynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega_) override;
private:
  /// computes the entropy function H(f)=sum_i f_i*ln(f_i/t_i)
  T computeEntropy(const T f[]);
  /// computes the entropy growth H(f)-H(f-alpha*fNeq)
  T computeEntropyGrowth(const T f[], const T fNeq[], const T &alpha);
  /// computes the entropy growth derivative
  /// dH/dalpha=-sum_i fNeq_i*ln((f_i-alpha*fNeq_i)/t_i)
  T computeEntropyGrowthDerivative(const T f[], const T fNeq[], const T &alpha);
  /// Get the alpha parameter
  bool getAlpha(T &alpha, const T f[], const T fNeq[]);

  T omega;  ///< relaxation parameter
};

/// Implementation of the forced entropic collision step
template<typename T, typename DESCRIPTOR>
class ForcedEntropicDynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ForcedEntropicDynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const;
  /// Collision step
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);
  /// Get local relaxation parameter of the dynamics
  virtual T getOmega() const;
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
private:
  /// computes the entropy function H(f)=sum_i f_i*ln(f_i/t_i)
  T computeEntropy(const T f[]);
  /// computes the entropy growth H(f)-H(f-alpha*fNeq)
  T computeEntropyGrowth(const T f[], const T fNeq[], const T &alpha);
  /// computes the entropy growth derivative
  /// dH/dalpha=-sum_i fNeq_i*ln((f_i-alpha*fNeq_i)/t_i)
  T computeEntropyGrowthDerivative(const T f[], const T fNeq[], const T &alpha);
  /// Get the alpha parameter
  bool getAlpha(T &alpha, const T f[], const T fNeq[]);

  T omega;  ///< relaxation parameter
};
}

#endif
