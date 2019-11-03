/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Orestis Malaspinas, Jonas Latt
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

#ifndef INAMURO_NEWTON_RAPHSON_DYNAMICS_H
#define INAMURO_NEWTON_RAPHSON_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "io/ostreamManager.h"

namespace olb {

/**
* This class computes the inamuro BC with general dynamics. It uses the formula from the
 * paper by Inamuro et al. but since there is no explict solution
 * for a lattice different from the D2Q9 and for a speed of sound
 * c_s=q/sqrt(3), we have to use a Newton-Raphson algorithm to
 * implement these boundary conditions.
*/
template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
class InamuroNewtonRaphsonDynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  InamuroNewtonRaphsonDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega) override;
private:
  void computeApproxMomentum(T approxMomentum[DESCRIPTOR::d],
                             const Cell<T,DESCRIPTOR> &cell,
                             const T &rho, const T u[DESCRIPTOR::d], const T xi[DESCRIPTOR::d],
                             const std::vector<int> knownIndexes,const std::vector<int> missingIndexes);

  /// compute the error (L^2 norm of (u-uApprox))
  T computeError(const T &rho,const T u[DESCRIPTOR::d], const T approxMomentum[DESCRIPTOR::d]);

  void computeGradGradError(T gradGradError[DESCRIPTOR::d][DESCRIPTOR::d],
                            T gradError[DESCRIPTOR::d],
                            const T &rho, const T u[DESCRIPTOR::d],const T xi[DESCRIPTOR::d],
                            const T approxMomentum[DESCRIPTOR::d],
                            const std::vector<int> missingIndexes);

  /// compute the new xi with the newton raphson algorithm
  bool newtonRaphson(T xi[DESCRIPTOR::d],
                     const T gradError[DESCRIPTOR::d],
                     const T gradGradError[DESCRIPTOR::d][DESCRIPTOR::d]);

  bool invert(const T a[2][2],T b[2][2]);

  bool invert(const T a[3][3],T b[3][3]);

  Dynamics _boundaryDynamics;
  T _xi[DESCRIPTOR::d];
  mutable OstreamManager clout;
};

}

#endif
