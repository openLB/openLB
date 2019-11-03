/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Patrick Nathen, Mathias J. Krause
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
#ifndef SMAGORINSKY_MRT_DYNAMICS_H
#define SMAGORINSKY_MRT_DYNAMICS_H

#include "mrtDynamics.h"
#include "core/cell.h"


namespace olb {

/// Implementation of the MRT collision step
template<typename T, typename DESCRIPTOR>
class SmagorinskyMRTdynamics : public MRTdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  SmagorinskyMRTdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_, T dx_ = 1, T dt_ = 1);


  // Collide
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);

  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);

  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getSmagorinskyOmega(Cell<T,DESCRIPTOR>& cell_);

private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_);

  /// Computes the local smagorinsky relaxation parameter
  T computeOmega(T omega0_, T preFactor_, T rho_, T pi_[util::TensorVal<DESCRIPTOR >::n] );

protected:
  /// Smagorinsky constant
  T smagoConst;
  /// Precomputed constant which speeeds up the computation
  T preFactor;

  T omega; // the shear viscosity relaxatin time
  T lambda;// the bulk viscosity relaxatin time

  /// effective collision time based upon Smagorisnky approach
  T tau_eff;


  T dx;
  T dt;
  // Relaxation Time Matrix for
  T invM_S_SGS[DESCRIPTOR::q][DESCRIPTOR::q];

};


/// Implementation of the MRT collision step
template<typename T, typename DESCRIPTOR>
class SmagorinskyForcedMRTdynamics : public SmagorinskyMRTdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  SmagorinskyForcedMRTdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_, T dx_, T dt_ ) : SmagorinskyMRTdynamics<T,DESCRIPTOR>(omega_, momenta_, smagoConst_, dx_, dt_ ) {};

  // Collide
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);
};
}

#endif
