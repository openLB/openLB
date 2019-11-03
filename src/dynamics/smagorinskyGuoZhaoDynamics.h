/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016-2017 Davide Dapelo, Mathias J. Krause
 *  OpenLB e-mail contact: info@openlb.net
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
 * Specific dynamics classes for Guo and Zhao (2002) porous model
 * with a Smagorinsky LES turbulence model, with
 * which a Cell object can be instantiated -- header file.
 */
#ifndef LB_SMAGO_GUOZHAO_DYNAMICS_H
#define LB_SMAGO_GUOZHAO_DYNAMICS_H

#include "dynamics/guoZhaoLatticeDescriptors.h"
#include "core/util.h"
#include "core/postProcessing.h"
#include "core/latticeStatistics.h"

namespace olb {

/// Implementation of the BGK collision step with porous force according to
/// Guo and Zhao (2012), described as an external force
template<typename T, typename DESCRIPTOR>
class SmagorinskyGuoZhaoBGKdynamics : public GuoZhaoBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor.
  //Passing default value for smagoConst_ because 2D boundary conditions accept only
  //two-argument constructor for dynamics class.
  SmagorinskyGuoZhaoBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_ = 0.14,
                                T dx_ = 1, T dt_ = 1);
  /// Collision step
  virtual void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_);
  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getSmagorinskyOmega(Cell<T,DESCRIPTOR>& cell_);
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
protected:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_);
  /// Computes the local smagorinsky relaxation parameter
  T computeOmega(T omega0, T preFactor_, T rho,
                 T pi[util::TensorVal<DESCRIPTOR >::n] );

  /// effective collision time based upon Smagorisnky approach
  T tau_eff;
  /// Smagorinsky constant
  T smagoConst;
  /// Precomputed constant which speeeds up the computation
  T preFactor;
  T dx;
  T dt;
};

}

#endif
