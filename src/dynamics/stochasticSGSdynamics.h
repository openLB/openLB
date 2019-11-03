/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Patrick Nathen
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
 * MRT Dynamics with adjusted omega -- header file.
 */
#ifndef STOCHASTIC_SGS_DYNAMICS_H
#define STOCHASTIC_SGS_DYNAMICS_H

#include "mrtDynamics.h"
#include "core/cell.h"


namespace olb {

/// Implementation of the MRT collision step with stochastic relaxation based on
/// " A stochastic subgrid model with application to turbulent flow and scalar mixing"; Phys. of Fluids 19; 2007
template<typename T, typename DESCRIPTOR>
class StochasticSGSdynamics : public MRTdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  StochasticSGSdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T turbulenceInt_, T charU_, T smagoConst_, T dx_ = 1, T dt_ = 1 );


  // Collide
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);

  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);

  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getSmagorinskyOmega(Cell<T,DESCRIPTOR>& cell_, T X_lang_n_);

  /// Get local Random number of BoxMüllertransform -> returns randBM
  virtual T getRandBMTrans(Cell<T,DESCRIPTOR>& cell_,  T turbulenceInt_, T charU_);

  /// Get local Random number of BoxMüllertransform -> returns randBM
  // virtual void setRandomWalk(Cell<T,DESCRIPTOR>& cell_, T CharU, T drift_, T result_ );
  virtual T getRandomWalk(Cell<T,DESCRIPTOR>& cell_, T drift_, T result_);


private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_);

  /// Computes the local smagorinsky relaxation parameter
  T computeOmega(T omega0_, T preFactor_, T rho_, T pi_[util::TensorVal<DESCRIPTOR >::n] , T X_lang_n_);

  /// Computes the local time scale from SGS dissipation rate for BMtransform
  T computeTimeScale(T preFactor_, T rho_, T pi_[util::TensorVal<DESCRIPTOR >::n], T smagoConst_, T X_lang_n_);
  // virtual void setTimeScale(T preFactor_, T rho_, T pi_[util::TensorVal<DESCRIPTOR >::n], T smagoConst_ ,T X_lang_n_);

private:
  /// effective collision time based upon Smagorisnky approach
  T tau_eff;
  /// Initial turbulence intensity for random number generator
  T turbulenceInt;
  /// Smagorinsky Constant
  /// Precomputed constant which speeeds up the computation
  T smagoConst;
  T preFactor;

  T dx;
  T dt;

  T omega; // the shear viscosity relaxatin time
  T lambda;// the bulk viscosity relaxatin time

  //T result;
  T charU;
  //T drift;
  T X_lang_n;

  // Relaxation Time Matrix for
  T invM_S_SGS[DESCRIPTOR::q][DESCRIPTOR::q];

};

}

#endif
