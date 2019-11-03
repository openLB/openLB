/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Jonas Latt
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
 * BGK Dynamics for porous media -- header file.
 */
#ifndef POROUS_BGK_DYNAMICS_H
#define POROUS_BGK_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "core/cell.h"

namespace olb {

/// Implementation of the BGK collision step for a porosity model
template<typename T, typename DESCRIPTOR>
class PorousBGKdynamics : public BGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  PorousBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);

  /// Collision step
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);

  /// get relaxation parameter
  T    getOmega() const;
  /// set relaxation parameter
  void setOmega(T omega_);


private:
  T omega;      ///< relaxation parameter
};

/// Implementation of the BGK collision step for a porosity model enabling
/// drag computation
template<typename T, typename DESCRIPTOR>
class ExtendedPorousBGKdynamics : public BGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ExtendedPorousBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// extended Collision step, computes local drag in a given direction
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);
  /// get relaxation parameter
  T    getOmega() const;
  /// set relaxation parameter
  void setOmega(T omega_);


private:
  T omega;      ///< relaxation parameter
};

/// Implementation of the BGK collision step for subgridscale particles
template<typename T, typename DESCRIPTOR>
class SubgridParticleBGKdynamics : public BGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  SubgridParticleBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// extended Collision step, computes local drag in a given direction
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);
  /// get relaxation parameter
  T    getOmega() const;
  /// set relaxation parameter
  void setOmega(T omega_);


private:
  T omega;      ///< relaxation parameter
  T _fieldTmp[4];
};

/* Implementation of the BGK collision for moving porous media (HLBM approach).
 * As this scheme requires additionla data stored in an external field, 
 * it is meant to be used along with a PorousParticle descriptor.
 * \param omega Lattice relaxation frequency
 * \param momenta A standard object for the momenta computation
 */
template<typename T, typename DESCRIPTOR, bool isStatic=false>
class PorousParticleBGKdynamics : public BGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  PorousParticleBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// extended Collision step, computes local drag in a given direction
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);
  /// get relaxation parameter
  T    getOmega() const;
  /// set relaxation parameter
  void setOmega(T omega_);


private:
  T omega;      ///< relaxation parameter
  /// This structure is used to emulate a "static if" to switch between static and
  /// dynamic case. It can be replaced by a "constexpr if" when switching to C++17
  /// standard.
  template<bool isStaticStruct, bool dummy = true> struct effectiveVelocity {
      static void calculate(T* pExternal, T* pVelocity);
  };
};

/// Implementation of the HBGK collision step for a porosity model enabling
/// drag computation for many particles
/// including the Krause turbulence modell

template<typename T, typename DESCRIPTOR>
class KrauseHBGKdynamics : public BGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  KrauseHBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_, T dx_ = 1, T dt_ = 1);
  /// extended Collision step, computes local drag in a given direction
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);
  /// get relaxation parameter
  T    getOmega() const;
  /// set relaxation parameter
  void setOmega(T omega_);


private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_);
  /// Computes the local smagorinsky relaxation parameter
  void computeOmega(T omega0_, Cell<T,DESCRIPTOR>& cell, T preFactor_, T rho_,
                    T u[DESCRIPTOR::d],
                    T newOmega[DESCRIPTOR::q] );

private:

  T omega;      ///< relaxation parameter
  /// effective collision time based upon Smagorisnky approach
  T tau_eff;
  /// Smagorinsky constant
  T smagoConst;
  /// Precomputed constant which speeeds up the computation
  T preFactor;
  T dx;
  T dt;

  T _fieldTmp[4];

};

/// Implementation of the BGK collision step for a porosity model enabling
/// drag computation
template<typename T, typename DESCRIPTOR>
class ParticlePorousBGKdynamics : public BGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ParticlePorousBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// extended Collision step, computes local drag in a given direction
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);
  /// get relaxation parameter
  T    getOmega() const;
  /// set relaxation parameter
  void setOmega(T omega_);


private:
  T omega;      ///< relaxation parameter
};

/// Implementation of the BGK collision step for a small particles enabling
/// two way coupling
template<typename T, typename DESCRIPTOR>
class SmallParticleBGKdynamics : public BGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  SmallParticleBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// extended Collision step, computes local drag in a given direction
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);
  /// get relaxation parameter
  T    getOmega() const;
  /// set relaxation parameter
  void setOmega(T omega_);


private:
  T omega;      ///< relaxation parameter
};


enum Mode {M2, M3};

/// Implementation of the Partially Saturated Method (PSM),
/// see Krüger, Timm, et al. The Lattice Boltzmann Method. Springer, 2017. (p.447-451)
template<typename T, typename DESCRIPTOR>
class PSMBGKdynamics : public BGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  PSMBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, int mode_=0);
  ///  Compute fluid velocity on the cell.
  void computeU (
    Cell<T,DESCRIPTOR> const& cell,
    T u[DESCRIPTOR::d] ) const override;
  /// Compute fluid velocity and particle density on the cell.
  void computeRhoU (
    Cell<T,DESCRIPTOR> const& cell,
    T& rho, T u[DESCRIPTOR::d]) const override;
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// get relaxation parameter
  T    getOmega() const;
  /// set relaxation parameter
  void setOmega(T omega_);


private:
  T omega;      ///< relaxation parameter
  T paramA;      /// speed up parameter
  Mode mode;
};

/// Implementation of the Partially Saturated Method (PSM),
/// see Krüger, Timm, et al. The Lattice Boltzmann Method. Springer, 2017. (p.447-451)
template<typename T, typename DESCRIPTOR>
class ForcedPSMBGKdynamics : public ForcedBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ForcedPSMBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, int mode_=0);
  ///  Compute fluid velocity on the cell.
  void computeU (
    Cell<T,DESCRIPTOR> const& cell,
    T u[DESCRIPTOR::d] ) const override;
  /// Compute fluid velocity and particle density on the cell.
  void computeRhoU (
    Cell<T,DESCRIPTOR> const& cell,
    T& rho, T u[DESCRIPTOR::d]) const override;
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// get relaxation parameter
  T    getOmega() const;
  /// set relaxation parameter
  void setOmega(T omega_);


private:
  T omega;      ///< relaxation parameter
  T paramA;      /// speed up parameter
  Mode mode;
};

} // olb

#endif
