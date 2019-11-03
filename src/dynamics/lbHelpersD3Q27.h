/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2015 Jonas Latt, Mathias J. Krause
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
 * functions of the header file lbHelpers.h, for some D3Q27 grids.
 */

#ifndef LB_HELPERS_D3Q27_H
#define LB_HELPERS_D3Q27_H

namespace olb {

// Efficient specialization for D3Q27 lattice
template<typename T, typename... FIELDS>
struct lbDynamicsHelpers<T, descriptors::D3Q27<FIELDS...> > {
  using SpecializedCellBase = CellBase<T,descriptors::D3Q27<FIELDS...>>;

  template <typename>
  using SpecializedDescriptor = descriptors::D3Q27<FIELDS...>;

  static T equilibrium( int iPop, T rho, const T u[3], const T uSqr )
  {
    typedef descriptors::D3Q27<> L;
    T c_u = descriptors::c<L>(iPop,0)*u[0] + descriptors::c<L>(iPop,1)*u[1] + descriptors::c<L>(iPop,2)*u[2];
    return rho * descriptors::t<T,L>(iPop) * ( 1. + 3.*c_u + 4.5*c_u*c_u - 1.5*uSqr ) - descriptors::t<T,L>(iPop);
  }

  static T equilibriumFirstOrder( int iPop, T rho, const T u[3] )
  {
    typedef descriptors::D3Q27<> L;
    T c_u = descriptors::c<L>(iPop,0) * u[0] + descriptors::c<L>(iPop,1) * u[1] + descriptors::c<L>(iPop,2) * u[2];

    return rho * descriptors::t<T,L>(iPop) * ( ( T )1 + c_u * descriptors::invCs2<T,L>() ) - descriptors::t<T,L>(iPop);
  }

  static T incEquilibrium(int iPop, const T j[3], const T jSqr, const T pressure)
  {
    typedef descriptors::D3Q27<> L;
    T c_j = descriptors::c<L>(iPop,0)*j[0] + descriptors::c<L>(iPop,1)*j[1] + descriptors::c<L>(iPop,2)*j[2];
    return descriptors::t<T,L>(iPop) * ( 3.*pressure + 3.*c_j + 4.5*c_j*c_j - 1.5*jSqr ) - descriptors::t<T,L>(iPop);
  }

  static void computeFneq (
    SpecializedCellBase const& cell, T fNeq[27], T rho, const T u[3] )
  {
    const T uSqr = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    for (int iPop=0; iPop < 27; ++iPop) {
      fNeq[iPop] = cell[iPop] - equilibrium(iPop, rho, u, uSqr);
    }
  }

  static T bgkCollision(SpecializedCellBase& cell, T const& rho, const T u[3], T const& omega)
  {
    typedef descriptors::D3Q27<> L;

    T one_m_omega = (T)1 - omega;
    T t0_omega = descriptors::t<T,L>(0)*omega; // weight for i=0
    T t1_omega = descriptors::t<T,L>(1)*omega; // weight for i=1,2,3,14,15,16
    T t4_omega = descriptors::t<T,L>(4)*omega; // weight for i=4,5,6,7,8,9,17,18,19,20,21,22
    T t10_omega = descriptors::t<T,L>(10)*omega; //weight for i=10,11,12,13,23,24,25,26

    T uSqr     = u[0]*u[0] + u[1]*u[1] + u[2]*u[2]; // compute of usqr

    T u3x      = (T)3*u[0]; // compute of 3*ux
    T u3y      = (T)3*u[1]; // compute of 3*uy
    T u3z      = (T)3*u[2]; // compute of 3*uz

    T u3xSqr_  = .5*u3x*u3x; // compute 9/2*(ux)²
    T u3ySqr_  = .5*u3y*u3y; // compute 9/2*(uy)²
    T u3zSqr_  = .5*u3z*u3z; // compute 9/2*(uz)²

    T u3xu3y_  = u3x*u3y; // compute 9(ux*uy)
    T u3xu3z_  = u3x*u3z; // compute 9(ux*uz)
    T u3yu3z_  = u3y*u3z; // compute 9(uy*uz)

    T C1 = (T)1 + (T)3*uSqr; // compute 1+3*((ux)²+(uy)²+(uz)²)
    T C2, C3;

    //**************************case i=0
    C3 = -u3xSqr_ - u3ySqr_ - u3zSqr_; // compute -9/2*((ux)²+(uy)²+(uz)²)

    // compute of feq0=rho*t0*(c1+c3)=rho*t0*(1-3/2*((ux)²+(uy)²+(uz)²))
    cell[0] *= one_m_omega;
    cell[0] += t0_omega*(rho*(C1 + C3) - (T)1);

    //**************************case i=1 and i=14
    C2 = -u3x; // compute -3*ux
    C3 = -u3ySqr_ - u3zSqr_; // compute -9/2*((uy)²+(uz)²)

    // compute of feq1=rho*t1*(c1+c2+c3)=rho*t0*(1-3*ux+9/2*(ux)²-3/2*((ux)²+(uy)²+(uz)²))
    cell[1]  *= one_m_omega;
    cell[1]  += t1_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq10=rho*t1*(c1-c2+c3)=rho*t0*(1+3*ux+9/2*(ux)²-3/2*((ux)²+(uy)²+(uz)²))
    cell[14] *= one_m_omega;
    cell[14] += t1_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=2 and i=15
    C2 = -u3y; // compute -3*uy
    C3 = -u3xSqr_ - u3zSqr_; // compute -9/2*((ux)²+(uz)²)

    // compute of feq2=rho*t1*(c1-c2+c3)=rho*t1*(1-3*uy+9/2*(uy)²-3/2*((ux)²+(uy)²+(uz)²))
    cell[2]  *= one_m_omega;
    cell[2]  += t1_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq2=rho*t1*(c1+c2+c3)=rho*t1*(1+3*uy+9/2*(uy)²-3/2*((ux)²+(uy)²+(uz)²))
    cell[15] *= one_m_omega;
    cell[15] += t1_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=3 and i=16
    C2 = -u3z; // compute -3*uz
    C3 = -u3xSqr_ - u3ySqr_; // compute -9/2*((ux)²+(uy)²)

    // compute of feq3=rho*t1*(c1+c2+c3)=rho*t1*(1-3*uz+9/2*(uz)²-3/2*((ux)²+(uy)²+(uz)²))
    cell[3]  *= one_m_omega;
    cell[3]  += t1_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq12=rho*t1*(c1-c2+c3)=rho*t1*(1+3*uz+9/2*(uz)²-3/2*((ux)²+(uy)²+(uz)²))
    cell[16] *= one_m_omega;
    cell[16] += t1_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=4 and i=17
    C2 = -u3x - u3y; // compute -3*(uz+uy)
    C3 = u3xu3y_ - u3zSqr_; // compute 9*(ux*uy)-9/2*(uz)²

    // compute of feq4=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(ux+uy)+9/2*((ux)²+(uy)²+2*ux*uy)-3/2*((ux)²+(uy)²+(uz)²))
    cell[4]  *= one_m_omega;
    cell[4]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq13=rho*t4*(c1-c2+c3)=rho*t4*(1+3*(ux+uy)+9/2*((ux)²+(uy)²+2*ux*uy)-3/2*((ux)²+(uy)²+(uz)²))
    cell[17] *= one_m_omega;
    cell[17] += t4_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=5 and i=18
    C2 = -u3x + u3y; // compute -3*(ux+uy)
    C3 = -u3xu3y_ - u3zSqr_; // compute -9*(ux*uy)-9/2*(uz)²

    // compute of feq5=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(ux-uy)+9/2*((ux)²+(uy)²-2*ux*uy)-3/2*((ux)²+(uy)²+(uz)²))
    cell[5]  *= one_m_omega;
    cell[5]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq14=rho*t4*(c1-c2+c3)=rho*t4*(1+3*(ux-uy)+9/2*((ux)²+(uy)²-2*ux*uy)-3/2*((ux)²+(uy)²+(uz)²))
    cell[18] *= one_m_omega;
    cell[18] +=t4_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=6 and i=19
    C2 = -u3x - u3z; // compute -3*(ux+uz)
    C3 = u3xu3z_ - u3ySqr_; // compute 9*(ux*uz)-9/2*(uy)²

    // compute of feq6=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(ux+uz)+9/2*((ux)²+(uy)²+2*ux*uz)-3/2*((ux)²+(uy)²+(uz)²))
    cell[6]  *= one_m_omega;
    cell[6]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq15=rho*t4*(c1-c2+c3)=rho*t15*(1+3*(ux+uz)+9/2*((ux)²+(uy)²+2*ux*uz)-3/2*((ux)²+(uy)²+(uz)²))
    cell[19] *= one_m_omega;
    cell[19] += t4_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=7 and i=20
    C2 = -u3x + u3z; // compute -3*(ux-uz)
    C3 = -u3xu3z_ - u3ySqr_; // compute -9*(ux*uz)-9/2*(uy)²

    // compute of feq7=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(ux-uz)+9/2*((ux)²+(uz)²-2*ux*uz)-3/2*((ux)²+(uy)²+(uz)²))
    cell[7]  *= one_m_omega;
    cell[7]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq16=rho*t4*(c1-c2+c3)=rho*t4*(1+3*(ux-uz)+9/2*((ux)²+(uz)²-2*ux*uz)-3/2*((ux)²+(uy)²+(uz)²))
    cell[20] *= one_m_omega;
    cell[20] += t4_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=8 and i=21
    C2 = -u3y - u3z; // compute -3*(uy+uz)
    C3 = u3yu3z_ - u3xSqr_; // compute 9*(uy*uz)-9/2*(ux)²

    // compute of feq8=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(uy+uz)+9/2*((uz)²+(uy)²+2*uy*uz)-3/2*((ux)²+(uy)²+(uz)²))
    cell[8]  *= one_m_omega;
    cell[8]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq17=rho*t4*(c1-c2+c3)=rho*t4*(1+3*(uy+uz)+9/2*((uz)²+(uy)²+2*uy*uz)-3/2*((ux)²+(uy)²+(uz)²))
    cell[21] *= one_m_omega;
    cell[21] += t4_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=9 and i=22
    C2 = -u3y + u3z; // compute -3*(uy-uz)
    C3 = -u3yu3z_ - u3xSqr_; // compute -9*(uy*uz)-9/2*(ux)²

    // compute of feq9=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(uy-uz)+9/2*((uy)²+(uz)²-2*uy*uz)-3/2*((ux)²+(uy)²+(uz)²))
    cell[9]  *= one_m_omega;
    cell[9]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq18=rho*t4*(c1-c2+c3)=rho*t4*(1+3*(uy-uz)+9/2*((uy)²+(uz)²-2*uy*uz)-3/2*((ux)²+(uy)²+(uz)²))
    cell[22] *= one_m_omega;
    cell[22] += t4_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=10 and i=23
    C2 = -u3x - u3y - u3z; //compute of -3*(ux+uy+uz)
    C3 = u3xu3y_ + u3xu3z_ + u3yu3z_; //compute of 9*(ux*uy+ux*uz+uy*uz)

    //compute of feq10=rho*t10*(c1+c2+c3)=rho*t10*(1-3*(ux+uy+uz)+9*(ux*uy+ux*uz+uy*uz)+3*((ux)²+(uy)²+(uz)²))
    cell[10] *= one_m_omega;
    cell[10] += t10_omega*(rho*(C1 + C2 + C3) - (T)1);

    //compute of feq23=rho*t10*(c1-c2+c3)=rho*t10*(1+3*(ux+uy+uz)+9*(ux*uy+ux*uz+uy*uz)+3*((ux)²+(uy)²+(uz)²))
    cell[23] *= one_m_omega;
    cell[23] += t10_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=11 and i=24
    C2 = -u3x - u3y + u3z; //compute of -3*(ux+uy-uz)
    C3 = u3xu3y_ - u3xu3z_ - u3yu3z_; //compute of 9*(ux*uy-ux*uz-uy*uz)

    //compute of feq11=rho*t10*(c1+c2+c3)=rho*t10*(1-3*(ux+uy-uz)+9*(ux*uy-ux*uz-uy*uz)+3*((ux)²+(uy)²+(uz)²))
    cell[11] *= one_m_omega;
    cell[11] += t10_omega*(rho*(C1 + C2 + C3) - (T)1);

    //compute feq24=rho*t10*(c1-c2+c3)=rho*t10*(1+3*(ux+uy-uz)+9*(ux*uy-ux*uz-uy*uz)+3*((ux)²+(uy)²+(uz)²))
    cell[24] *= one_m_omega;
    cell[24] += t10_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=12 and i=25
    C2 = -u3x + u3y - u3z; //compute of -3*(ux-uy+uz)
    C3 = -u3xu3y_ + u3xu3z_ - u3yu3z_; //compute of 9*(-ux*uy+ux*uz-uy*uz)

    //compute of feq11=rho*t10*(c1+c2+c3)=rho*t10*(1-3*(ux-uy+uz)+9*(-ux*uy+ux*uz-uy*uz)+3*((ux)²+(uy)²+(uz)²))
    cell[12] *= one_m_omega;
    cell[12] += t10_omega*(rho*(C1 + C2 + C3) - (T)1);

    //compute feq24=rho*t10*(c1-c2+c3)=rho*t10*(1+3*(ux-uy+uz)+9*(-ux*uy+ux*uz-uy*uz)+3*((ux)²+(uy)²+(uz)²))
    cell[25] *= one_m_omega;
    cell[25] += t10_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=13 and i=26
    C2 = -u3x + u3y + u3z; //compute of -3*(ux-uy-uz)
    C3 = -u3xu3y_ - u3xu3z_ + u3yu3z_; //compute of 9*(-ux*uy-ux*uz+uy*uz)

    //compute of feq11=rho*t10*(c1+c2+c3)=rho*t10*(1-3*(ux-uy-uz)+9*(-ux*uy-ux*uz+uy*uz)+3*((ux)²+(uy)²+(uz)²))
    cell[13] *= one_m_omega;
    cell[13] += t10_omega*(rho*(C1 + C2 + C3) - (T)1);

    //compute feq24=rho*t10*(c1-c2+c3)=rho*t10*(1+3*(ux-uy-uz)+9*(-ux*uy-ux*uz+uy*uz)+3*((ux)²+(uy)²+(uz)²))
    cell[26] *= one_m_omega;
    cell[26] += t10_omega*(rho*(C1 - C2 + C3) - (T)1);

    return uSqr;
  }

  static T incBgkCollision(SpecializedCellBase& cell, T pressure, const T j[3], T omega)
  {
    const T jSqr = util::normSqr<T,descriptors::D3Q27<>::d>(j);
    for (int iPop=0; iPop < descriptors::D3Q27<>::q; ++iPop) {
      cell[iPop] *= (T)1-omega;
      cell[iPop] += omega * lbDynamicsHelpers<T,descriptors::D3Q27<> >
                    ::incEquilibrium(iPop, j, jSqr, pressure);
    }
    return jSqr;
  }

  static T constRhoBgkCollision(SpecializedCellBase& cell, T rho, const T u[3], T ratioRho, T omega)
  {
    const T uSqr = util::normSqr<T,descriptors::D3Q27<>::d>(u);
    for (int iPop=0; iPop < descriptors::D3Q27<>::q; ++iPop) {
      T feq = lbDynamicsHelpers<T,descriptors::D3Q27<> >::
              equilibrium(iPop, rho, u, uSqr );
      cell[iPop] =
        ratioRho*(feq+descriptors::t<T,descriptors::D3Q27<>>(iPop))
        -descriptors::t<T,descriptors::D3Q27<>>(iPop) +
        ((T)1-omega)*(cell[iPop]-feq);
    }
    return uSqr;
  }

  static void partial_rho ( SpecializedCellBase const& cell,
                            T& surfX_M1, T& surfX_0, T& surfX_P1,
                            T& surfY_M1, T& surfY_P1, T& surfZ_M1, T& surfZ_P1 )
  {
    surfX_M1 = cell[1] + cell[4] + cell[5] + cell[6] + cell[7] +
               cell[10] + cell[11] + cell[12] + cell[13];
    surfX_0  = cell[0] + cell[2] + cell[3] + cell[8] + cell[9] +
               cell[15] + cell[16] + cell[21] + cell[22];
    surfX_P1 = cell[14] + cell[17] + cell[18] + cell[19] + cell[20] +
               cell[23] + cell[24] + cell[25] + cell[26];

    surfY_M1 = cell[2] + cell[4] + cell[8] + cell[9] + cell[10] +
               cell[11] + cell[18] + cell[25] + cell[26];
    surfY_P1 = cell[5] + cell[12] + cell[13] + cell[15] + cell[17] +
               cell[21] + cell[22] + cell[23] + cell[24];

    surfZ_M1 = cell[3] + cell[6] + cell[8] + cell[10] + cell[12] +
               cell[20] + cell[22] + cell[24] + cell[26];
    surfZ_P1 = cell[7] + cell[9] + cell[11] + cell[13] + cell[16] +
               cell[19] + cell[21] + cell[23] + cell[25];
  }

  static void computeRhoU(SpecializedCellBase const& cell, T& rho, T u[3])
  {
    T surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

    surfX_M1 = cell[1] + cell[4] + cell[5] + cell[6] + cell[7] +
               cell[10] + cell[11] + cell[12] + cell[13];
    surfX_0  = cell[0] + cell[2] + cell[3] + cell[8] + cell[9] +
               cell[15] + cell[16] + cell[21] + cell[22];
    surfX_P1 = cell[14] + cell[17] + cell[18] + cell[19] + cell[20] +
               cell[23] + cell[24] + cell[25] + cell[26];

    surfY_M1 = cell[2] + cell[4] + cell[8] + cell[9] + cell[10] +
               cell[11] + cell[18] + cell[25] + cell[26];
    surfY_P1 = cell[5] + cell[12] + cell[13] + cell[15] + cell[17] +
               cell[21] + cell[22] + cell[23] + cell[24];

    surfZ_M1 = cell[3] + cell[6] + cell[8] + cell[10] + cell[12] +
               cell[20] + cell[22] + cell[24] + cell[26];
    surfZ_P1 = cell[7] + cell[9] + cell[11] + cell[13] + cell[16] +
               cell[19] + cell[21] + cell[23] + cell[25];

    rho = surfX_M1 + surfX_0 + surfX_P1 + (T)1;
    T invRho= 1./rho;

    u[0]  = ( surfX_P1 - surfX_M1 )*invRho;
    u[1]  = ( surfY_P1 - surfY_M1 )*invRho;
    u[2]  = ( surfZ_P1 - surfZ_M1 )*invRho;
  }

  static void computeRhoJ(SpecializedCellBase const& cell, T& rho, T j[3])
  {
    T surfX_M1, surfX_0, surfX_P1,
    surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

    partial_rho(cell, surfX_M1, surfX_0, surfX_P1,
                surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

    rho = surfX_M1 + surfX_0 + surfX_P1 + (T)1;

    j[0]  = ( surfX_P1 - surfX_M1 );
    j[1]  = ( surfY_P1 - surfY_M1 );
    j[2]  = ( surfZ_P1 - surfZ_M1 );
  }

  static void computeJ(SpecializedCellBase const& cell, T j[3])
  {
    T surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

    surfX_M1 = cell[1] + cell[4] + cell[5] + cell[6] + cell[7] +
               cell[10] + cell[11] + cell[12] + cell[13];
    surfX_P1 = cell[14] + cell[17] + cell[18] + cell[19] + cell[20] +
               cell[23] + cell[24] + cell[25] + cell[26];

    surfY_M1 = cell[2] + cell[4] + cell[8] + cell[9] + cell[10] +
               cell[11] + cell[18] + cell[25] + cell[26];
    surfY_P1 = cell[5] + cell[12] + cell[13] + cell[15] + cell[17] +
               cell[21] + cell[22] + cell[23] + cell[24];

    surfZ_M1 = cell[3] + cell[6] + cell[8] + cell[10] + cell[12] +
               cell[20] + cell[22] + cell[24] + cell[26];
    surfZ_P1 = cell[7] + cell[9] + cell[11] + cell[13] + cell[16] +
               cell[19] + cell[21] + cell[23] + cell[25];

    j[0]  = ( surfX_P1 - surfX_M1 );
    j[1]  = ( surfY_P1 - surfY_M1 );
    j[2]  = ( surfZ_P1 - surfZ_M1 );
  }

  static void computeStress(SpecializedCellBase const& cell, T rho, const T u[3], T pi[6])
  {
    typedef descriptors::D3Q27<> L;
    // Workaround for Intel(r) compiler 9.1;
    // "using namespace util::tensorIndices3D" is not sufficient
    using util::tensorIndices3D::xx;
    using util::tensorIndices3D::yy;
    using util::tensorIndices3D::zz;
    using util::tensorIndices3D::xy;
    using util::tensorIndices3D::xz;
    using util::tensorIndices3D::yz;

    T surfX_M1, surfX_0, surfX_P1,
    surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

    partial_rho(cell, surfX_M1, surfX_0, surfX_P1,
                surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

    pi[xx] = surfX_P1+surfX_M1 - 1./descriptors::invCs2<T,L>()*(rho-(T)1) - rho*u[0]*u[0];
    pi[yy] = surfY_P1+surfY_M1 - 1./descriptors::invCs2<T,L>()*(rho-(T)1) - rho*u[1]*u[1];
    pi[zz] = surfZ_P1+surfZ_M1 - 1./descriptors::invCs2<T,L>()*(rho-(T)1) - rho*u[2]*u[2];

    pi[xy] = cell[4] + cell[10] + cell[11] - cell[5] - cell[12] - cell[13]
             + cell[17] + cell[23] + cell[24] - cell[18] - cell[25] - cell[26] - rho*u[0]*u[1];
    pi[xz] = cell[6] + cell[10] + cell[12] - cell[7] - cell[11] - cell[13]
             + cell[19] + cell[23] + cell[25] - cell[20] - cell[24] - cell[26] - rho*u[0]*u[2];
    pi[yz] = cell[8] + cell[10] + cell[26] - cell[9] - cell[11] - cell[25]
             + cell[21] + cell[13] + cell[23] - cell[22] - cell[12] - cell[24] - rho*u[1]*u[2];
  }

  static void computeAllMomenta(SpecializedCellBase const& cell, T& rho, T u[3], T pi[6])
  {
    typedef descriptors::D3Q27<> L;
    // Workaround for Intel(r) compiler 9.1;
    // "using namespace util::tensorIndices3D" is not sufficient
    using util::tensorIndices3D::xx;
    using util::tensorIndices3D::yy;
    using util::tensorIndices3D::zz;
    using util::tensorIndices3D::xy;
    using util::tensorIndices3D::xz;
    using util::tensorIndices3D::yz;

    T surfX_M1, surfX_0, surfX_P1,
    surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

    partial_rho(cell, surfX_M1, surfX_0, surfX_P1,
                surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

    rho = surfX_M1 + surfX_0 + surfX_P1 + (T)1;
    T invRho = 1. / rho;

    u[0]  = ( surfX_P1 - surfX_M1 ) * invRho;
    u[1]  = ( surfY_P1 - surfY_M1 ) * invRho;
    u[2]  = ( surfZ_P1 - surfZ_M1 ) * invRho;

    pi[xx] = surfX_P1+surfX_M1 - 1./descriptors::invCs2<T,L>()*(rho-(T)1) - rho*u[0]*u[0];
    pi[yy] = surfY_P1+surfY_M1 - 1./descriptors::invCs2<T,L>()*(rho-(T)1) - rho*u[1]*u[1];
    pi[zz] = surfZ_P1+surfZ_M1 - 1./descriptors::invCs2<T,L>()*(rho-(T)1) - rho*u[2]*u[2];

    pi[xy] = cell[4] + cell[10] + cell[11] - cell[5] - cell[12] - cell[13]
             + cell[17] + cell[23] + cell[24] - cell[18] - cell[25] - cell[26] - rho*u[0]*u[1];
    pi[xz] = cell[6] + cell[10] + cell[12] - cell[7] - cell[11] - cell[13]
             + cell[19] + cell[23] + cell[25] - cell[20] - cell[24] - cell[26] - rho*u[0]*u[2];
    pi[yz] = cell[8] + cell[10] + cell[26] - cell[9] - cell[11] - cell[25]
             + cell[21] + cell[13] + cell[23] - cell[22] - cell[12] - cell[24] - rho*u[1]*u[2];
  }

  static T computeRho(SpecializedCellBase const& cell)
  {
    T rho = cell[0] + cell[1] + cell[2] + cell[3] + cell[4]
            + cell[5] + cell[6] + cell[7] + cell[8]
            + cell[9] + cell[10] + cell[11] + cell[12]
            + cell[13] + cell[14] + cell[15] + cell[16]
            + cell[17] + cell[18] + cell[19] + cell[20]
            + cell[21] + cell[22] + cell[23] + cell[24]
            + cell[25] + cell[26] + (T)1;
    return rho;
  }

  static void modifyVelocity(SpecializedCellBase const& cell, const T newU[3])
  {
    T rho, oldU[3];
    computeRhoU(cell, rho, oldU);
    const T oldUSqr = util::normSqr<T,3>(oldU);
    const T newUSqr = util::normSqr<T,3>(newU);
    for (int iPop=0; iPop<27; ++iPop) {
      cell[iPop] = cell[iPop]
                   - equilibrium(iPop, rho, oldU, oldUSqr)
                   + equilibrium(iPop, rho, newU, newUSqr);
    }
  }

};  //struct lbDynamicsHelpers<D3Q27<>>

//Efficient specialization for D3Q27 lattice and for forced D3Q27 lattice
//  (operations applying to the whole lattice)

template<typename T>
struct lbLatticeHelpers<T, descriptors::D3Q27<>> {

  static void swapAndStreamCell (
    Cell<T,descriptors::D3Q27<>> ***grid,
    int iX, int iY, int iZ, int nX, int nY, int nZ, int iPop, T& fTmp )
  {
    fTmp                      = grid[iX][iY][iZ][iPop];
    grid[iX][iY][iZ][iPop]    = grid[iX][iY][iZ][iPop+13];
    grid[iX][iY][iZ][iPop+13] = grid[nX][nY][nZ][iPop];
    grid[nX][nY][nZ][iPop]    = fTmp;
  }

  static void swapAndStream3D(Cell<T,descriptors::D3Q27<>> ***grid,
                              int iX, int iY, int iZ)
  {
    T fTmp;
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY,   iZ,   1,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY-1, iZ,   2,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY,   iZ-1, 3,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ,   4,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ,   5,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY,   iZ-1, 6,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY,   iZ+1, 7,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY-1, iZ-1, 8,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY-1, iZ+1, 9,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ-1, 10, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ+1, 11, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ-1, 12, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ+1, 13, fTmp);
  }

};

template<typename T>
struct lbLatticeHelpers<T, descriptors::D3Q27<descriptors::FORCE>> {

  static void swapAndStreamCell (
    Cell<T,descriptors::D3Q27<descriptors::FORCE>> ***grid,
    int iX, int iY, int iZ, int nX, int nY, int nZ, int iPop, T& fTmp )
  {
    fTmp                      = grid[iX][iY][iZ][iPop];
    grid[iX][iY][iZ][iPop]    = grid[iX][iY][iZ][iPop+13];
    grid[iX][iY][iZ][iPop+13] = grid[nX][nY][nZ][iPop];
    grid[nX][nY][nZ][iPop]    = fTmp;
  }

  static void swapAndStream3D(Cell<T,descriptors::D3Q27<descriptors::FORCE>> ***grid,
                              int iX, int iY, int iZ)
  {
    T fTmp;
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY,   iZ,   1,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY-1, iZ,   2,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY,   iZ-1, 3,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ,   4,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ,   5,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY,   iZ-1, 6,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY,   iZ+1, 7,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY-1, iZ-1, 8,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY-1, iZ+1, 9,  fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ-1, 10, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ+1, 11, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ-1, 12, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ+1, 13, fTmp);
  }

};

}  // namespace olb

#endif
