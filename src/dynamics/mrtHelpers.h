/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007, 2013 Jonas Latt, Mathias J. Krause, Geng Liu
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
 * Helper functions for the implementation of LB dynamics. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef MRT_HELPERS_H
#define MRT_HELPERS_H

#include "lbHelpers.h"
#include "mrtLatticeDescriptors.h"

namespace olb {

 

/// All helper functions are inside this structure
template<typename T, typename DESCRIPTOR>
struct mrtHelpers {
  static_assert(
    std::is_same<typename DESCRIPTOR::category_tag, descriptors::tag::MRT>::value,
    "DESCRIPTOR is tagged as MRT");

  /// Computation of equilibrium distribution (in momenta space)
  static T equilibrium( int iPop, T rho, const T u[DESCRIPTOR::d],
                        const T uSqr )
  {
    T equ = T();
    for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
      equ += descriptors::m<T,DESCRIPTOR>(iPop,jPop) *
             (lbHelpers<T,DESCRIPTOR>::equilibrium(jPop,rho,u,uSqr) +
              descriptors::t<T,DESCRIPTOR>(jPop));
    }

    return equ;
  }

  /// Computation of all equilibrium distribution (in momenta space)
  static void computeEquilibrium( T momentaEq[DESCRIPTOR::q],
                                  T rho, const T u[DESCRIPTOR::d],
                                  const T uSqr )
  {
    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      momentaEq[iPop] = T();
      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
        momentaEq[iPop] += descriptors::m<T,DESCRIPTOR>(iPop,jPop) *
                           (lbHelpers<T,DESCRIPTOR>::equilibrium(jPop,rho,u,uSqr) +
                            descriptors::t<T,DESCRIPTOR>(jPop));
      }
    }
  }

  static void computeMomenta(T momenta[DESCRIPTOR::q], Cell<T,DESCRIPTOR> &cell)
  {
    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      momenta[iPop] = T();
      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
        momenta[iPop] += descriptors::m<T,DESCRIPTOR>(iPop,jPop) *
                         (cell[jPop] + descriptors::t<T,DESCRIPTOR>(jPop));
      }
    }
  }

  /// MRT collision step
  static T mrtCollision( Cell<T,DESCRIPTOR>& cell,
                         T rho, const T u[DESCRIPTOR::d],
                         T invM_S[DESCRIPTOR::q][DESCRIPTOR::q])
  {
    T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    T momenta[DESCRIPTOR::q];
    T momentaEq[DESCRIPTOR::q];

    computeMomenta(momenta,cell);
    computeEquilibrium(momentaEq,rho,u,uSqr);

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      T collisionTerm = T();
      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
        collisionTerm += invM_S[iPop][jPop] *
                         (momenta[jPop] - momentaEq[jPop]);
      }
      cell[iPop] -= collisionTerm;
    }

    return uSqr;
  }

  /// MRT SGS collision step
  static T mrtSGSCollision( Cell<T,DESCRIPTOR>& cell,
                            T rho, const T u[DESCRIPTOR::d],
                            T omega,
                            T invM_S_SGS[DESCRIPTOR::q][DESCRIPTOR::q])
  {
    T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    T momenta[DESCRIPTOR::q];
    T momentaEq[DESCRIPTOR::q];


    computeMomenta(momenta,cell);
    computeEquilibrium(momentaEq,rho,u,uSqr);

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      T collisionTerm = T();
      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
        collisionTerm += invM_S_SGS[iPop][jPop] *
                         (momenta[jPop] - momentaEq[jPop]);
      }
      cell[iPop] -= collisionTerm;
    }

    return uSqr;
  }


  /// Ladd-Verberg-I body force model for MRT
  /// A.Ladd, R. Verberg, DESCRIPTOR-Boltzmann simulations of particle-fluid suspensions, Journal of Statistical Physics 104(2001)
  static void addExternalForce( Cell<T,DESCRIPTOR>& cell,
                                T rho,
                                const T u[DESCRIPTOR::d],
                                T invM_S[DESCRIPTOR::q][DESCRIPTOR::q])
  {
    T* force = cell.template getFieldPointer<descriptors::FORCE>();
    T f_u = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      f_u += force[iD]*u[iD];
    }

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      T c_u = T();
      T c_f = T();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
        c_f += descriptors::c<DESCRIPTOR>(iPop,iD)*force[iD];
      }
      T f1 = descriptors::t<T,DESCRIPTOR>(iPop)*rho*c_f*descriptors::invCs2<T,DESCRIPTOR>();
      T f2 = descriptors::t<T,DESCRIPTOR>(iPop)*rho*(c_u*c_f*descriptors::invCs2<T,DESCRIPTOR>()-f_u)*descriptors::invCs2<T,DESCRIPTOR>();

      T invMsM = T();
      for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
        invMsM += invM_S[iPop][jPop]*descriptors::m<T,DESCRIPTOR>(jPop,iPop);
      }
      cell[iPop] += f1 + f2 - invMsM*f2/2.;
    }

    return;
  }

};  // struct mrtHelpers

}  // namespace olb

// The specialized code is directly included. That is because we never want
// it to be precompiled so that in both the precompiled and the
// "include-everything" version, the compiler can apply all the
// optimizations it wants.
#include "mrtHelpers2D.h"
#include "mrtHelpers3D.h"

#endif
