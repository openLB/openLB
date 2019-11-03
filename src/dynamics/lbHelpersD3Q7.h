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
 * functions of the header file lbHelpers.h, for the D3Q7 grid.
 */

#ifndef LB_HELPERS_D3Q7_H
#define LB_HELPERS_D3Q7_H

#include "rtlbmDescriptors.h"

namespace olb {

/** partial template specialization for D3Q7DescriptorBase.
 */
template<typename T, typename... FIELDS>
struct lbDynamicsHelpers<T, descriptors::D3Q7<FIELDS...> > {
  using SpecializedCellBase = CellBase<T,descriptors::D3Q7<FIELDS...>>;

  template <typename>
  using SpecializedDescriptor = descriptors::D3Q7<FIELDS...>;

  static T equilibrium( int iPop, T rho, const T u[3], T uSqr )
  {
    typedef descriptors::D3Q7<> L;
    T c_u = descriptors::c<L>(iPop,0)*u[0] + descriptors::c<L>(iPop,1)*u[1];
    return rho * descriptors::t<T,L>(iPop) * (
             1. + 3.*c_u + 4.5*c_u*c_u - 1.5*uSqr )
           - descriptors::t<T,L>(iPop);
  }

  static T equilibriumFirstOrder( int iPop, T rho, const T u[3] )
  {
    typedef descriptors::D3Q7<> L;
    T c_u = descriptors::c<L>(iPop,0) * u[0] + descriptors::c<L>(iPop,1) * u[1] + descriptors::c<L>(iPop,2) * u[2];
    return rho*descriptors::t<T,L>(iPop)*((T)1 + c_u*descriptors::invCs2<T,L>())-descriptors::t<T,L>(iPop);
  }

  /// RLB advection diffusion collision step
  static T rlbCollision( SpecializedCellBase& cell, T rho, const T u[3], T omega )
  {
    typedef descriptors::D3Q7<> L;
    const T uSqr = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];

    const T Cs2 = (T)1 / descriptors::invCs2<T,L>();

    T rho_1 = rho - (T)1;
    cell[0] = ( (T)1 - (T)3 * Cs2 ) * rho_1; //f[0]=(1-3c_s^2)(rho-1)

    const T omega_ = (T)1 - omega;

    const T f1_4 = omega_ * ( cell[1] - cell[4] );
    const T f2_5 = omega_ * ( cell[2] - cell[5] );
    const T f3_6 = omega_ * ( cell[3] - cell[6] );

    rho_1 *= Cs2;

    const T ux_ = omega * rho * u[0];

    cell[1] = (T)0.5 * ( rho_1 + f1_4 - ux_ ); //f[1]=1/2*(c_s^2(rho-1)+(1-omega)*(f[1]-f[4])-omega*rho*u[x])
    cell[4] = (T)0.5 * ( rho_1 - f1_4 + ux_ ); //f[4]=1/2*(c_s^2(rho-1)-(1-omega)*(f[1]-f[4])+omega*rho*u[x])

    const T uy_ = omega * rho * u[1];

    cell[2] = (T)0.5 * ( rho_1 + f2_5 - uy_ ); //f[2]=1/2*(c_s^2(rho-1)+(1-omega)*(f[2]-f[5])-omega*rho*u[y])
    cell[5] = (T)0.5 * ( rho_1 - f2_5 + uy_ ); //f[5]=1/2*(c_s^2(rho-1)-(1-omega)*(f[2]-f[5])+omega*rho*u[y])

    const T uz_ = omega * rho * u[2];

    cell[3] = (T)0.5 * ( rho_1 + f3_6 - uz_ ); //f[3]=1/2*(c_s^2(rho-1)+(1-omega)*(f[3]-f[6])-omega*rho*u[y])
    cell[6] = (T)0.5 * ( rho_1 - f3_6 + uz_ ); //f[6]=1/2*(c_s^2(rho-1)-(1-omega)*(f[3]-f[6])+omega*rho*u[y])

    return uSqr;
  }

  // BGK advection diffusion collision step
  static T bgkCollision( SpecializedCellBase& cell, T rho, const T u[3], T omega )
  {
    const T Cs2 = 0.25;
    const T uSqr = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];

    const T omega_ = (T)1 - omega;
    const T omega_2 = (T)0.5 * omega;
    T rho_ = ( rho - (T)1 );
    cell[0] = omega_ * cell[0] + omega * ( (T)1 - (T)3 * Cs2 ) * rho_;
    const T jx = rho * u[0];
    const T jy = rho * u[1];
    const T jz = rho * u[2];

    rho_ *= Cs2;
    cell[1] = omega_ * cell[1] + omega_2 * ( rho_ - jx );
    cell[2] = omega_ * cell[2] + omega_2 * ( rho_ - jy );
    cell[3] = omega_ * cell[3] + omega_2 * ( rho_ - jz );
    cell[4] = omega_ * cell[4] + omega_2 * ( rho_ + jx );
    cell[5] = omega_ * cell[5] + omega_2 * ( rho_ + jy );
    cell[6] = omega_ * cell[6] + omega_2 * ( rho_ + jz );

    return uSqr;
  }

  /// Paper: Mink et al. 2016 DOI: 10.1016/j.jocs.2016.03.014
  /// omega is expected to be one
  static T sinkCollision( SpecializedCellBase& cell, T intensity, T omega, T sink )
  {
    // fi = (1-w) fi + w fi^{eq} - absorption fi;
    // where fi^{eq} = rho ti and t0 = 1/4, t1=...t6=1/8
    // note: cell[i] = fi - ti

    const T omega_ = (T)1 - omega;
    T ti = 0.25;
    cell[0] =  omega_ * ( cell[0] + ti ) + omega * intensity * ti - sink * ( cell[0] + ti ) - ti;
    ti = 0.125;
    cell[1] = omega_ * ( cell[1] + ti ) + omega * intensity * ti - sink * ( cell[1] + ti ) - ti;
    cell[2] = omega_ * ( cell[2] + ti ) + omega * intensity * ti - sink * ( cell[2] + ti ) - ti;
    cell[3] = omega_ * ( cell[3] + ti ) + omega * intensity * ti - sink * ( cell[3] + ti ) - ti;
    cell[4] = omega_ * ( cell[4] + ti ) + omega * intensity * ti - sink * ( cell[4] + ti ) - ti;
    cell[5] = omega_ * ( cell[5] + ti ) + omega * intensity * ti - sink * ( cell[5] + ti ) - ti;
    cell[6] = omega_ * ( cell[6] + ti ) + omega * intensity * ti - sink * ( cell[6] + ti ) - ti;
    return 0.0;
  }

// MRT advection-diffusion functions
static void computeMomentaEquilibrium(T momentaEq[7], T rho, const T u[3], T uSqr) {
//   (Maria Lloret)
//    momentaEq[0] = rho;
//    momentaEq[1] = rho * u[0];
//    momentaEq[2] = rho * u[1];
//    momentaEq[3] = rho * u[2];
//    momentaEq[4] = -rho * (T) 3 / (T) 4;
//    momentaEq[5] = T();
//    momentaEq[6] = T();

//  Li, Yang et al 2016: The directions are modified for the OpenLB definition
    momentaEq[0] = rho;
    momentaEq[1] = rho * u[0];
    momentaEq[2] = rho * u[1];
    momentaEq[3] = rho * u[2];
    momentaEq[4] = rho * (T) 3 / (T) 4;    //(Maria Lloret) -rho*(T)3/(T)4;
    momentaEq[5] = T();
    momentaEq[6] = T();
  }

  static void computeMomenta(T momenta[7], SpecializedCellBase const& cell)
  {
    /* The implementation is based on the "Passive heat transfer
     * in a turbulent channel flow simulation using large eddy
     * simulation based on the lattice Boltzmann method framework"
     * by Hong Wu, Jiao Wang, Zhi Tao
     * (Maria Lloret)
     */
//    momenta[0] = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6];
//    momenta[1] = -cell[1] + cell[4];
//    momenta[2] = -cell[2] + cell[5];
//    momenta[3] = -cell[3] + cell[6];
//    momenta[4] = -(T) 6 * cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] - (T) 3 / (T) 4;
//    momenta[5] = cell[1] - cell[2] + cell[4] - cell[5];
//    momenta[6] = cell[1] + cell[2] - (T) 2 * cell[3] + cell[4] + cell[5] - (T) 2 * cell[6];

//  Li, Yang et al 2016: The directions are modified for the OpenLB definition
    momenta[0] = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6];
    momenta[1] = cell[1] - cell[3];
    momenta[2] = -cell[2] + cell[5];
    momenta[3] = cell[4] - cell[6];
    momenta[4] = (T) 6 * cell[0] - cell[1] - cell[2] - cell[3] - cell[4] - cell[5] - cell[6] + (T) 3 / (T) 4;
    momenta[5] = (T) 2 * cell[1] - cell[2] + (T) 2 * cell[3] - cell[4] - cell[5] - cell[6];
    momenta[6] = cell[2] - cell[4] + cell[5] - cell[6];
  }

  static T mrtCollision( SpecializedCellBase& cell, const T rho, const T u[3], const T invM_S[7][7])
  {
    T uSqr = util::normSqr<T, 3>(u);
    T momenta[7];
    T momentaEq[7];

    computeMomenta(momenta, cell);
    computeMomentaEquilibrium(momentaEq, rho, u, uSqr);

//    std::cout << "momenta = ";
//    for (int i=0; i < 7; ++i) {
//        std::cout << momenta[i] << ", ";
//    }
//    std::cout << std::endl;

//    std::cout << "momentaEq = ";
//    for (int i=0; i < 7; ++i) {
//        std::cout << momentaEq[i] << ", ";
//    }
//    std::cout << std::endl;

    T mom1 = momenta[1] - momentaEq[1];
    T mom2 = momenta[2] - momentaEq[2];
    T mom3 = momenta[3] - momentaEq[3];
    T mom4 = momenta[4] - momentaEq[4];
    T mom5 = momenta[5] - momentaEq[5];
    T mom6 = momenta[6] - momentaEq[6];

    cell[0] -= invM_S[0][1]*mom1 + invM_S[0][2]*mom2 +
               invM_S[0][3]*mom3 + invM_S[0][4]*mom4 +
               invM_S[0][5]*mom5 + invM_S[0][6]*mom6;

    cell[1] -= invM_S[1][1]*mom1 + invM_S[1][2]*mom2 +
               invM_S[1][3]*mom3 + invM_S[1][4]*mom4 +
               invM_S[1][5]*mom5 + invM_S[1][6]*mom6;

    cell[2] -= invM_S[2][1]*mom1 + invM_S[2][2]*mom2 +
               invM_S[2][3]*mom3 + invM_S[2][4]*mom4 +
               invM_S[2][5]*mom5 + invM_S[2][6]*mom6;

    cell[3] -= invM_S[3][1]*mom1 + invM_S[3][2]*mom2 +
               invM_S[3][3]*mom3 + invM_S[3][4]*mom4 +
               invM_S[3][5]*mom5 + invM_S[3][6]*mom6;

    cell[4] -= invM_S[4][1]*mom1 + invM_S[4][2]*mom2 +
               invM_S[4][3]*mom3 + invM_S[4][4]*mom4 +
               invM_S[4][5]*mom5 + invM_S[4][6]*mom6;

    cell[5] -= invM_S[5][1]*mom1 + invM_S[5][2]*mom2 +
               invM_S[5][3]*mom3 + invM_S[5][4]*mom4 +
               invM_S[5][5]*mom5 + invM_S[5][6]*mom6;

    cell[6] -= invM_S[6][1]*mom1 + invM_S[6][2]*mom2 +
               invM_S[6][3]*mom3 + invM_S[6][4]*mom4 +
               invM_S[6][5]*mom5 + invM_S[6][6]*mom6;
    return uSqr;
  }

    static void computeFneq (
      SpecializedCellBase const& cell, T fNeq[7], T rho, const T u[3] )
    {
      // std::cout << "WARNING: First order equilibrium function is used for the calculation of the non-equilibrium parts!" << std::endl;
      for (int iPop=0; iPop < 7; ++iPop) {
        fNeq[iPop] = cell[iPop] - equilibriumFirstOrder(iPop, rho, u);
      }
    }

    static T computeRho(SpecializedCellBase const& cell)
    {
      T rho = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + (T)1;
      return rho;
    }

    static void computeRhoU(SpecializedCellBase const& cell, T& rho, T u[3])
    {
       rho =  cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + (T)1;
       T invRho= 1./rho;

       u[0]  = ( cell[4] - cell[1] )*invRho;
       u[1]  = ( cell[5] - cell[2] )*invRho;
       u[2]  = ( cell[6] - cell[3] )*invRho;
    }

    static void computeRhoJ(SpecializedCellBase const& cell, T& rho, T j[3])
    {
      rho =  cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + (T)1;

      j[0]  = ( cell[4] - cell[1] );
      j[1]  = ( cell[5] - cell[2] );
      j[2]  = ( cell[6] - cell[3] );
    }

    static void computeJ(SpecializedCellBase const& cell, T j[3])
    {
      j[0]  = ( cell[4] - cell[1] );
      j[1]  = ( cell[5] - cell[2] );
      j[2]  = ( cell[6] - cell[3] );
    }

    static void computeStress(SpecializedCellBase const& cell, T rho, const T u[3], T pi[6])
    {
      // printf("ERROR: Stress tensor not defined in D3Q7!\n");
    }

    static void computeAllMomenta(SpecializedCellBase const& cell, T& rho, T u[3], T pi[6])
    {
      assert(false);
    }

};


}  // namespace olb

#endif
