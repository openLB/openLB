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
 * Local boundary cell 3D dynamics -- header file.
 */
#ifndef MOMENTA_ON_BOUNDARIES_3D_H
#define MOMENTA_ON_BOUNDARIES_3D_H

#include "momentaOnBoundaries.h"

namespace olb {

template<typename T, typename DESCRIPTOR,
         int plane, int normal1, int normal2>
class InnerEdgeVelBM3D : public DirichletBoundaryMomenta<T,DESCRIPTOR> {
public:
  enum { direction1 = (plane+1)%3, direction2 = (plane+2)%3 };
public:
  /// Default Constructor: initialization to zero
  InnerEdgeVelBM3D();
  /// Constructor with boundary initialization
  InnerEdgeVelBM3D(const T u_[DESCRIPTOR::d]);

  T computeRho(Cell<T,DESCRIPTOR> const& cell) const override;
  void computeU (
    Cell<T,DESCRIPTOR> const& cell,
    T u[DESCRIPTOR::d] ) const override;
  void computeJ (
    Cell<T,DESCRIPTOR> const& cell,
    T j[DESCRIPTOR::d] ) const override;
  void computeU(T u[DESCRIPTOR::d]) const;
  void defineRho(Cell<T,DESCRIPTOR>& cell, T rho) override ;
  void defineU(Cell<T,DESCRIPTOR>& cell,
                       const T u[DESCRIPTOR::d]) override ;
  void defineU(const T u[DESCRIPTOR::d]);
  void defineAllMomenta (
    Cell<T,DESCRIPTOR>& cell,
    T rho, const T u[DESCRIPTOR::d],
    const T pi[util::TensorVal<DESCRIPTOR >::n] ) override;
  /// Stress tensor
  void computeStress (
    Cell<T,DESCRIPTOR> const& cell,
    T rho, const T u[DESCRIPTOR::d],
    T pi[util::TensorVal<DESCRIPTOR >::n] ) const override;
private:
  T _u[DESCRIPTOR::d];   ///< value of the velocity on the boundary
};


template<typename T, typename DESCRIPTOR,
         int normalX, int normalY, int normalZ>
class InnerCornerVelBM3D : public DirichletBoundaryMomenta<T,DESCRIPTOR> {
public:
  /// Default Constructor: initialization to zero
  InnerCornerVelBM3D();
  /// Constructor with boundary initialization
  InnerCornerVelBM3D(const T u_[DESCRIPTOR::d]);

  T computeRho(Cell<T,DESCRIPTOR> const& cell) const override;
  void computeU (
    Cell<T,DESCRIPTOR> const& cell,
    T u[DESCRIPTOR::d] ) const override;
  void computeJ (
    Cell<T,DESCRIPTOR> const& cell,
    T j[DESCRIPTOR::d] ) const override;
  void computeU(T u[DESCRIPTOR::d]) const;
  void defineRho(Cell<T,DESCRIPTOR>& cell, T rho) override ;
  void defineU(Cell<T,DESCRIPTOR>& cell,
                       const T u[DESCRIPTOR::d]) override ;
  void defineU(const T u[DESCRIPTOR::d]);
  void defineAllMomenta (
    Cell<T,DESCRIPTOR>& cell,
    T rho, const T u[DESCRIPTOR::d],
    const T pi[util::TensorVal<DESCRIPTOR >::n] ) override;
  /// Stress tensor
  void computeStress (
    Cell<T,DESCRIPTOR> const& cell,
    T rho, const T u[DESCRIPTOR::d],
    T pi[util::TensorVal<DESCRIPTOR >::n] ) const override;
private:
  T _u[DESCRIPTOR::d];   ///< value of the velocity on the boundary
};

}

#endif
