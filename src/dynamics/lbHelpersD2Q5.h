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
 * Template specializations for some computationally intensive LB
 * functions of the header file lbHelpers.h, for the D2Q5 grid.
 */

#ifndef LB_HELPERS_D2Q5_H
#define LB_HELPERS_D2Q5_H

#include "mrtLatticeDescriptors.h"

namespace olb {

// Efficient specialization for D2Q5 lattice
template<typename T, typename... FIELDS>
struct lbDynamicsHelpers<T, descriptors::D2Q5<FIELDS...> > {
  using SpecializedCellBase = CellBase<T,descriptors::D2Q5<FIELDS...>>;

  static T equilibrium( int iPop, T rho, const T u[2], T uSqr )
  {
    typedef descriptors::D2Q5<> L;
    T c_u = descriptors::c<L>(iPop,0)*u[0] + descriptors::c<L>(iPop,1)*u[1];
    return rho * descriptors::t<T,L>(iPop) * (
             1. + 3.*c_u + 4.5*c_u*c_u - 1.5*uSqr )
           - descriptors::t<T,L>(iPop);
  }

  /// equilibrium distribution
  static T equilibriumFirstOrder( int iPop, T rho, const T u[2] )
  {
    typedef descriptors::D2Q5<> L;
    T c_u = descriptors::c<L>(iPop,0) * u[0] + descriptors::c<L>(iPop,1) * u[1];

    return rho * descriptors::t<T,L>(iPop) * ( ( T )1 + c_u * descriptors::invCs2<T,L>() ) - descriptors::t<T,L>(iPop);
  }

  /// RLB advection diffusion collision step
  static T rlbCollision( SpecializedCellBase& cell, T rho, const T u[2], T omega )
  {
    typedef descriptors::D2Q5<> L;
    const T uSqr = u[0] * u[0] + u[1] * u[1];

    const T Cs2 = ( T )1 / descriptors::invCs2<T,L>();

    T rho_1 = rho - ( T )1;
    cell[0] = ( ( T )1 - ( T )2 * Cs2 ) * rho_1; //f[0]=(1-2c_s^2)(rho-1)

    const T omega_ = ( T )1 - omega;

    const T f1_3 = ( T )0.5 * omega_ * ( cell[1] - cell[3] );
    const T f2_4 = ( T )0.5 * omega_ * ( cell[2] - cell[4] );

    rho_1 *= ( T )0.5 * Cs2;

    const T ux_ = ( T )0.5 * omega * rho * u[0];

    cell[1] = rho_1 + f1_3 - ux_; //f[1]=1/2*(c_s^2(rho-1)+(1-omega)*(f[1]-f[3])-omega*rho*u[x])
    cell[3] = rho_1 - f1_3 + ux_; //f[3]=1/2*(c_s^2(rho-1)-(1-omega)*(f[1]-f[3])+omega*rho*u[x])

    const T uy_ = ( T )0.5 * omega * rho * u[1];

    cell[2] = rho_1 + f2_4 - uy_; //f[2]=1/2*(c_s^2(rho-1)+(1-omega)*(f[2]-f[4])-omega*rho*u[y])
    cell[4] = rho_1 - f2_4 + uy_; //f[4]=1/2*(c_s^2(rho-1)-(1-omega)*(f[2]-f[4])+omega*rho*u[y])

    return uSqr;
  }

  // BGK advection diffusion collision step
  static T bgkCollision(  SpecializedCellBase& cell, T rho, const T u[2], T omega )
  {
    typedef descriptors::D2Q5<> L;

    const T Cs2 = ( T )1 / descriptors::invCs2<T,L>();
    const T uSqr = u[0] * u[0] + u[1] * u[1];

    const T omega_ = ( T )1 - omega;
    const T omega_2 = ( T )0.5 * omega;
    T rho_ = ( rho - ( T )1 );

    cell[0] = omega_ * cell[0] + omega * ( ( T )1 - ( T )2 * Cs2 ) * rho_;

    const T jx = rho * u[0];
    const T jy = rho * u[1];

    rho_ *= Cs2;
    cell[1] = omega_ * cell[1] + omega_2 * ( rho_ - jx );
    cell[2] = omega_ * cell[2] + omega_2 * ( rho_ - jy );
    cell[3] = omega_ * cell[3] + omega_2 * ( rho_ + jx );
    cell[4] = omega_ * cell[4] + omega_2 * ( rho_ + jy );

    return uSqr;
  }

  static void computeMomentaEquilibrium( T momentaEq[5], T rho, const T u[2], T uSqr )
  {
    momentaEq[0] = rho;
    momentaEq[1] = rho*u[0];
    momentaEq[2] = rho*u[1];
    momentaEq[3] = -rho*(T)2/(T)3;
    momentaEq[4] = T();
  }

  static void computeMomenta( T momenta[5], SpecializedCellBase const& cell )
  {
    momenta[0] = cell[0] + cell[1] + cell[2] + cell[3] + cell[4];
    momenta[1] = -cell[1] + cell[3];
    momenta[2] = -cell[2] + cell[4];
    momenta[3] = -(T)4*cell[0] + cell[1] + cell[2] + cell[3] + cell[4] - (T)2 / (T)3;
    momenta[4] = cell[1] - cell[2] + cell[3] - cell[4];
  }

  static T mrtCollision( SpecializedCellBase& cell, T rho, const T u[2], T invM_S[5][5] )
  {
    T uSqr = util::normSqr<T,2>(u);
    T momenta[5];
    T momentaEq[5];

    computeMomenta(momenta, cell);
    computeMomentaEquilibrium(momentaEq, rho, u, uSqr);

    T mom1 = momenta[1] - momentaEq[1];
    T mom2 = momenta[2] - momentaEq[2];
    T mom3 = momenta[3] - momentaEq[3];
    T mom4 = momenta[4] - momentaEq[4];

    cell[0] -= invM_S[0][1]*mom1 + invM_S[0][2]*mom2 +
               invM_S[0][3]*mom3 + invM_S[0][4]*mom4;

    cell[1] -= invM_S[1][1]*mom1 + invM_S[1][2]*mom2 +
               invM_S[1][3]*mom3 + invM_S[1][4]*mom4;

    cell[2] -= invM_S[2][1]*mom1 + invM_S[2][2]*mom2 +
               invM_S[2][3]*mom3 + invM_S[2][4]*mom4;

    cell[3] -= invM_S[3][1]*mom1 + invM_S[3][2]*mom2 +
               invM_S[3][3]*mom3 + invM_S[3][4]*mom4;

    cell[4] -= invM_S[4][1]*mom1 + invM_S[4][2]*mom2 +
               invM_S[4][3]*mom3 + invM_S[4][4]*mom4;

    return uSqr;
  }

  static void computeFneq ( SpecializedCellBase const& cell, T fNeq[5], T rho, const T u[2] )
  {
    for (int iPop=0; iPop < 5; ++iPop) {
      // printf("WARNING: First order equilibrium function is used for the calculation of the non-equilibrium parts!\n");
      fNeq[iPop] = cell[iPop] - equilibriumFirstOrder(iPop, rho, u);
    }
  }

  static T computeRho(SpecializedCellBase const& cell)
  {
    T rho = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + (T)1;
    return rho;
  }

  static void computeRhoU(SpecializedCellBase const& cell, T& rho, T u[2])
  {
    rho = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + (T)1;
    T invRho= 1./rho;
    u[0]  = (cell[3] - cell[1])*invRho;
    u[1]  = (cell[4] - cell[2])*invRho;
  }

   static void computeRhoJ(SpecializedCellBase const& cell, T& rho, T j[2])
   {
     rho = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + (T)1;

     j[0]  = (cell[3] - cell[1]);
     j[1]  = (cell[4] - cell[2]);
   }

   static void computeJ(SpecializedCellBase const& cell, T j[2])
   {
     j[0]  = (cell[3] - cell[1]);
     j[1]  = (cell[4] - cell[2]);
   }

   static void computeStress(SpecializedCellBase const& cell, T rho, const T u[2], T pi[6])
   {
     // printf("ERROR: Stress tensor not defined in D2Q5!\n");
   }

   static void computeAllMomenta(SpecializedCellBase const& cell, T& rho, T u[2], T pi[6])
   {
     rho = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + (T)1;
     T invRho= 1./rho;
     u[0]  = (cell[3] - cell[1])*invRho;
     u[1]  = (cell[4] - cell[2])*invRho;
     assert(false);
   }

};


}  // namespace olb

#endif

