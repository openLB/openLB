/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2015 Mathias J. Krause, Jonas Latt, Patrick Nathen
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
 * BGK Dynamics with adjusted omega -- header file.
 */
#ifndef SMAGORINSKY_BGK_DYNAMICS_H
#define SMAGORINSKY_BGK_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "core/cell.h"

#include<complex> // For shear kalman Smagorinsky - Populations

namespace olb {

/// Interface for the Large-Eddy-Simulation dynamics classes
template<typename T, typename DESCRIPTOR>
struct LESDynamics {
  /// Destructor: virtual to enable inheritance
  virtual ~LESDynamics() { }
  /// Get local effective relaxation parameter of the dynamics
  virtual T getEffectiveOmega(Cell<T,DESCRIPTOR>& cell_) =0;

};

/// Implementation of Smagorinsky Dynamics
template<typename T, typename DESCRIPTOR>
class SmagorinskyDynamics : public LESDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  SmagorinskyDynamics(T smagoConst_);
  /// get the constant preFactor variable used to speed up calculations
  virtual T getPreFactor();

private:
  /// Smagorinsky constant
  T smagoConst;

protected:
  /// get the Smagorinsky constant
  virtual T getSmagoConst();
  /// Compute constant prefactor variable in order to speed up the computation
  virtual T computePreFactor();
  /// Precomputed constant which speeeds up the computation
  T preFactor;
};

/// Implementation of the Smagorinsky BGK collision step
template<typename T, typename DESCRIPTOR>
class SmagorinskyBGKdynamics : public SmagorinskyDynamics<T,DESCRIPTOR>, public BGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  SmagorinskyBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_,
                         T smagoConst_);
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics_) override;
  /// Get local smagorinsky relaxation parameter of the dynamics
  T getEffectiveOmega(Cell<T,DESCRIPTOR>& cell_) override;

protected:
  /// Computes the local smagorinsky relaxation parameter
  virtual T computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell);
};

/// Implementation of the ForcedBGK collision step
template<typename T, typename DESCRIPTOR>
class SmagorinskyForcedBGKdynamics : public SmagorinskyDynamics<T,DESCRIPTOR>, public ForcedBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  SmagorinskyForcedBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_);
  /// Collision step
  virtual void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics_) override;
  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getEffectiveOmega(Cell<T,DESCRIPTOR>& cell_) override;

protected:
  /// Computes the local smagorinsky relaxation parameter
  virtual T computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell_);
};

/// Implementation of a LES BGK with non local effective tau calculation through external field
template<typename T, typename DESCRIPTOR>
class ExternalTauEffLESBGKdynamics : public SmagorinskyBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ExternalTauEffLESBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_ = (T)0);
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics_) override;
};

/// Implementation of a LES ForcedBGK with non local effective tau calculation through external field
template<typename T, typename DESCRIPTOR>
class ExternalTauEffLESForcedBGKdynamics : public SmagorinskyForcedBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ExternalTauEffLESForcedBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_ = (T)0);
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics_) override;
};

/// Implementation of the consistent Strain Smagorinsky BGK collision step
///
/// Consistent subgrid scale modelling for lattice Boltzmann methods
/// Orestis Malaspinas and Pierre Sagaut
/// Journal of Fluid Mechanics / Volume / June 2012, pp 514-542
/// DOI: http://dx.doi.org/10.1017/jfm.2012.155

template<typename T, typename DESCRIPTOR>
class ConStrainSmagorinskyBGKdynamics : public SmagorinskyBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ConStrainSmagorinskyBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_,
                                  T smagoConst_=T(.1));
protected:
  /// Computes the local smagorinsky relaxation parameter
  T computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell_);
};

/// Implementation of the consistent Smagorinsky BGK collision step
///
/// Consistent subgrid scale modelling for lattice Boltzmann methods
/// Orestis Malaspinas and Pierre Sagaut
/// Journal of Fluid Mechanics / Volume / June 2012, pp 514-542
/// DOI: http://dx.doi.org/10.1017/jfm.2012.155

template<typename T, typename DESCRIPTOR>
class ConSmagorinskyBGKdynamics : public SmagorinskyBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ConSmagorinskyBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_);
protected:
  /// should be remove --> David
  T computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell_);

};

/// Implementation of a the dynamic Smarorinsky BGK collision step
template<typename T, typename DESCRIPTOR>
class DynSmagorinskyBGKdynamics : public SmagorinskyBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  DynSmagorinskyBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);

protected:
  /// Computes the local smagorinsky relaxation parameter
  T computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell);
};

/// Implementation of the ADM BGK collision step

/*template<typename T, typename DESCRIPTOR>
class ADMBGKdynamics : public BGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ADMBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Collision step
  virtual void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics_);
private:
  T omega;
};*/

/// Implementation of the ForcedADMBGK collision step
template<typename T, typename DESCRIPTOR>
class ForcedADMBGKdynamics : public BGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ForcedADMBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);

  /// Collision step
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);
private:
  T omega;
};

/// Implementation of a Shear Smarorinsky BGK collision step
/// Shown good results for wall-bounded flows
/// Leveque et al.: Shear-Improved Smagorinsky Model for Large-Eddy Simulation
/// of Wall-Bounded Turbulent Flows
/// DOI: http://dx.doi.org/10.1017/S0022112006003429

template<typename T, typename DESCRIPTOR>
class ShearSmagorinskyBGKdynamics : public SmagorinskyBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ShearSmagorinskyBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_);
  /// Collision step
  virtual void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics_);
  /// Get Effective Omega stored in a external field
  virtual T getEffectiveOmega(Cell<T,DESCRIPTOR>& cell);
protected:
  /// Computes the local smagorinsky relaxation parameter
  T computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell, int iT);
  /// The external field variables' positions

};

/// Implementation of the ForcedBGK collision step
template<typename T, typename DESCRIPTOR>
class ShearSmagorinskyForcedBGKdynamics : public SmagorinskyForcedBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ShearSmagorinskyForcedBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_);
  /// Collision step
  virtual void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics_);
  /// Get Effective Omega stored in a external field
  virtual T getEffectiveOmega(Cell<T,DESCRIPTOR>& cell);
protected:
  /// Computes the local smagorinsky relaxation parameter
  T computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell, int iT);
  // Define current time step
  /// Smagorinsky constant
};

/// Implementation of the ForcedBGK collision step
template<typename T, typename DESCRIPTOR>
class SmagorinskyLinearVelocityForcedBGKdynamics : public SmagorinskyForcedBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  SmagorinskyLinearVelocityForcedBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_,
      T smagoConst_);
  /// Collision step
  virtual void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics_);
};

/// Implementation of the BGK collision step
template<typename T, typename DESCRIPTOR>
class KrauseBGKdynamics : public SmagorinskyBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  KrauseBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_);
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics_) override;
  /// Get local smagorinsky relaxation parameter of the dynamics
  T getEffectiveOmega(Cell<T,DESCRIPTOR>& cell_) override;

private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor() override;
  /// Computes the local smagorinsky relaxation parameter
  void computeEffectiveOmega(T omega0, Cell<T,DESCRIPTOR>& cell, T preFactor_, T rho,
                             T u[DESCRIPTOR::d], T newOmega[DESCRIPTOR::q]);
  T preFactor;
};


/// Implementation of the BGK collision step
template<typename T, typename DESCRIPTOR>
class WALEBGKdynamics : public SmagorinskyBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  WALEBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_);

protected:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor() override;
  /// Computes the local smagorinsky relaxation parameter
  T computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell_) override;
};

/// Implementation of the BGK collision step
template<typename T, typename DESCRIPTOR>
class WALEForcedBGKdynamics : public SmagorinskyForcedBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  WALEForcedBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_);

protected:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor() override;
  /// Computes the local smagorinsky relaxation parameter
  T computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell_) override;
};

/// Implementation of the BGK collision step
template<typename T, typename DESCRIPTOR>
class FDKalmanShearSmagorinskyBGKdynamics : public SmagorinskyBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  FDKalmanShearSmagorinskyBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_, T u_char_lat, T f_char_lat);
  /// Get local effective relaxation parameter of the dynamics
  virtual T getEffectiveOmega(Cell<T,DESCRIPTOR>& cell_);

protected:
  /// Computes a constant prefactor in order to speed up the computation
  virtual T computePreFactor();
  /// Computes the local smagorinsky relaxation parameter
  virtual T computeOmega(Cell<T,DESCRIPTOR>& cell_);

  // The variance of increment of kalman filtered velocity
  T VarInVelKal;
  T UCharLat;
private:
  void computeNormStrainRate(Cell<T,DESCRIPTOR>& cell, T& NormStrainRate);
  void KalmanStep(Cell<T,DESCRIPTOR>& cell);
};



////////////////////////////////////////////////////////////////////////////////
/// Implementation of a Shear Smarorinsky BGK collision step with Kalman Filter
//
/// Leveque et al.: Shear-Improved Smagorinsky Model for Large-Eddy Simulation
/// of Wall-Bounded Turbulent Flows
///
/// Boudet et al. (2016) A Kalman filter adapted of the estimation of mean gradients
//   in the a large-eddy simulation of unsteady turbulent flows.

template<typename T, typename DESCRIPTOR>
class ShearKalmanSmagorinskyBGKdynamics : public SmagorinskyBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ShearKalmanSmagorinskyBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_,
                                     T smagoConst_, T u_char_lat, T f_char_lat);
  /// Get local effective relaxation parameter of the dynamics
  virtual T getEffectiveOmega(Cell<T,DESCRIPTOR>& cell_);

protected:
  /// Computes the local smagorinsky relaxation parameter
  T computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell_);
private:
  /// Updates the filtered velocity with a Kalman procedure
  void KalmanStep(Cell<T,DESCRIPTOR>& cell);
  /// Computes the kalman filtered velocity and strain rate using the filtered population stored in a externa field
  void computeKalmanUStress(Cell<T,DESCRIPTOR>& cell, T (&KalmanU)[DESCRIPTOR::d], T (&KalmanPi)[util::TensorVal<DESCRIPTOR >::n]);
  /// Computes The Kalman filtered velocity using the filtered populations stored in a external field
  void computeKalmanU(Cell<T,DESCRIPTOR>& cell, T (&KalmanU)[DESCRIPTOR::d]);
  /// Computes the Kalman filtered strain rate using the filtered populations stored in a external field
  void computeKalmanStress(Cell<T,DESCRIPTOR>& cell, T (&KalmanU)[DESCRIPTOR::d], T (&KalmanPi)[util::TensorVal<DESCRIPTOR >::n]);
  /// Computes instantaneous tau_sgs and update kalman tau_sgs
  void computeAndupdateTauSgs(Cell<T,DESCRIPTOR>& cell, T rho, T pi[util::TensorVal<DESCRIPTOR >::n],
                              T KalmanPiNeqN[util::TensorVal<DESCRIPTOR >::n], T KalmanPiNeqN1[util::TensorVal<DESCRIPTOR >::n],
                              T K, T &tau_sgs);
  /// Methods to compute the square Norm of second order moment non-quilibrium distribution function
  void computeNormSOM(T pi[util::TensorVal<DESCRIPTOR >::n], T &piNorm);
  void computeNormSOM(T pi1[util::TensorVal<DESCRIPTOR >::n], T pi2[util::TensorVal<DESCRIPTOR >::n], T rho, T &piNorm);
  void computeNormSOM(T pi[util::TensorVal<DESCRIPTOR >::n], T rho, T &piNorm);
  /// Compute the instantaneous tau_sgs
  void computeTauSgs(Cell<T,DESCRIPTOR>& cell, T rho, T KalmanPiNeqNormSqr, T KalmanInstPiNeqNormSqr, T PiNeqNormSqr, T K, T &tau_sgs);
  void computeRoots4thPoly(T A, T B, T C, T D, T E, std::complex<T> (&Roots)[4]);
  // Update the local kalman tau_sgs stored in a external field
  void updateTauSgsKalman(Cell<T,DESCRIPTOR>& cell, T NN, T Nn1, T n1n1, T N1N1, T K, T tau_sgs_n1);

  // The variance of increment of kalman filtered velocity
  T VarInVelKal;
  T UCharLat;
};



}

#endif
