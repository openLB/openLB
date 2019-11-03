/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Albert Mink
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

#ifndef RTLBM_BOUNDARY_DYNAMICS_H
#define RTLBM_BOUNDARY_DYNAMICS_H


#include "dynamics/dynamics.h"

namespace olb {

/** Defines incoming (axis parallel) directions on flat walls.
 *
 * \param directions   takes values 0, 1, 2; which corresponds to x, y, z plane normal
 * \param orientation  takes values -1, 1; which is the sign of plane normal
 *
 *  Note:
 *  Dynamics computes the arriving density on walls and redistributes it equally to
 *  all inward pointing directions.
 */
template<typename T, typename DESCRIPTOR, int direction, int orientation>
class RtlbmDiffuseBoundaryDynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  RtlbmDiffuseBoundaryDynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step for flat boundary
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
  /// place holder
  T getOmega() const override;
  /// place holder
  void setOmega(T omega_) override;
};

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
class RtlbmDiffuseEdgeBoundaryDynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  RtlbmDiffuseEdgeBoundaryDynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step for flat boundary
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
  /// place holder
  T getOmega() const override;
  /// place holder
  void setOmega(T omega_) override;
};

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
class RtlbmDiffuseCornerBoundaryDynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  RtlbmDiffuseCornerBoundaryDynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step for corner
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
  /// place holder
  T getOmega() const override;
  /// place holder
  void setOmega(T omega_) override;
};


/** Defines incoming directions on flat walls.
 *  Emposed dirichlet density is distributed equivalently to all incoming/unkown directions.
 *
 * \param directions   takes values 0, 1, 2; which corresponds to x, y, z plane normal
 * \param orientation  takes values -1, 1; which is the sign of plane normal
 *
 *  Note:
 *  Dynamics ensure that cell has constant 0th moment.
 *  To set the constant 0th moment defineRho() must be called.
 */
template<typename T, typename DESCRIPTOR, int direction, int orientation>
class RtlbmDiffuseConstBoundaryDynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  RtlbmDiffuseConstBoundaryDynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step for flat boundary
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
  /// place holder
  T getOmega() const override;
  /// place holder
  void setOmega(T omega_) override;
};



/** Defines incoming directions on edge boundaries.
*   Emposed dirichlet density is distributed equivalently to all incoming/unkown directions.
*
* \param plane    takes values 0, 1, 2; which denotes the alignment of edge either along x, y or z axis
* \param normal1  takes values -1 or 1; normal orientation of a surface side 'right'
* \param normal2  takes values -1 or 1; normal orientation of a surface side 'left'
*/
template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
class RtlbmDiffuseConstEdgeBoundaryDynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  RtlbmDiffuseConstEdgeBoundaryDynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step for edges
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
  /// place holder
  T getOmega() const override;
  /// place holder
  void setOmega(T omega_) override;
};


/** Defines incoming directions on corner boundaries.
*   Emposed dirichlet density is distributed equivalently to all incoming/unkown directions.
*
* \param xNormal  takes value -1 or 1
* \param yNormal  takes value -1 or 1
* \param zNormal  takes value -1 or 1
*/
template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
class RtlbmDiffuseConstCornerBoundaryDynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  RtlbmDiffuseConstCornerBoundaryDynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step for corner
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
  /// place holder
  T getOmega() const override;
  /// place holder
  void setOmega(T omega_) override;
};



template<typename T, typename DESCRIPTOR, int direction, int orientation>
class RtlbmDirectedBoundaryDynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  RtlbmDirectedBoundaryDynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step for directed boundary walls
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
  /// place holder
  T getOmega() const override;
  /// place holder
  void setOmega(T omega_) override;
};
}  // namespace olb

#endif
