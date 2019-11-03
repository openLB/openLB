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
 * functions of the header file lbHelpers.h, for some D3Q15 grids.
 */

#ifndef LB_HELPERS_D3Q15_H
#define LB_HELPERS_D3Q15_H

namespace olb {

// Efficient specialization for D3Q15 lattice
template<typename T, typename... FIELDS>
struct lbDynamicsHelpers<T, descriptors::D3Q15<FIELDS...> > {
  using SpecializedCellBase   = CellBase<T,descriptors::D3Q15<FIELDS...>>;
  using SpecializedDescriptor = descriptors::D3Q15<FIELDS...>;

  static T equilibrium( int iPop, T rho, const T u[3], const T uSqr )
  {
    typedef descriptors::D3Q15<> L;
    T c_u = descriptors::c<L>(iPop,0)*u[0] + descriptors::c<L>(iPop,1)*u[1] + descriptors::c<L>(iPop,2)*u[2];
    return rho * descriptors::t<T,L>(iPop) * ( 1. + 3.*c_u + 4.5*c_u*c_u - 1.5*uSqr ) - descriptors::t<T,L>(iPop);
  }

  static T incEquilibrium(int iPop, const T j[3], const T jSqr, const T pressure)
  {
    typedef descriptors::D3Q15<> L;
    T c_j = descriptors::c<L>(iPop,0)*j[0] + descriptors::c<L>(iPop,1)*j[1] + descriptors::c<L>(iPop,2)*j[2];
    return descriptors::t<T,L>(iPop) * ( 3.*pressure + 3.*c_j + 4.5*c_j*c_j - 1.5*jSqr ) - descriptors::t<T,L>(iPop);
  }

  static void computeFneq(SpecializedCellBase const& cell, T fNeq[15], T rho, const T u[3] )
  {
    const T uSqr = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    for (int iPop=0; iPop < 15; ++iPop) {
      fNeq[iPop] = cell[iPop] - equilibrium(iPop, rho, u, uSqr);
    }
  }

  static T bgkCollision(SpecializedCellBase& cell, T const& rho, const T u[3], T const& omega)
  {
    const T uSqr = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    for (int iPop=0; iPop < 15; ++iPop) {
      cell[iPop] *= (T)1-omega;
      cell[iPop] += omega *
                    lbDynamicsHelpers<T,SpecializedDescriptor >::equilibrium(iPop, rho, u, uSqr);
    }
    return uSqr;
  }

  static T incBgkCollision(SpecializedCellBase& cell, T pressure, const T j[3], T omega)
  {
    const T jSqr = util::normSqr<T,SpecializedDescriptor::d>(j);
    for (int iPop=0; iPop < SpecializedDescriptor::q; ++iPop) {
      cell[iPop] *= (T)1-omega;
      cell[iPop] += omega * lbDynamicsHelpers<T,SpecializedDescriptor >::incEquilibrium (
                      iPop, j, jSqr, pressure );
    }
    return jSqr;
  }

  static T constRhoBgkCollision(SpecializedCellBase& cell, T rho, const T u[3], T ratioRho, T omega)
  {
    const T uSqr = util::normSqr<T,SpecializedDescriptor::d>(u);
    for (int iPop=0; iPop < SpecializedDescriptor::q; ++iPop) {
      T feq = lbHelpers<T,descriptors::D3Q15<>>::
              equilibrium(iPop, rho, u, uSqr );
      cell[iPop] =
        ratioRho*(feq+descriptors::t<T,SpecializedDescriptor>(iPop))
        -descriptors::t<T,SpecializedDescriptor>(iPop) +
        ((T)1-omega)*(cell[iPop]-feq);
    }
    return uSqr;
  }

  static void partial_rho(SpecializedCellBase const& cell,
                          T& surfX_M1, T& surfX_0, T& surfX_P1,
                          T& surfY_M1, T& surfY_P1, T& surfZ_M1, T& surfZ_P1 )
  {
    surfX_M1 = cell[1] + cell[4] + cell[5] + cell[6] + cell[7];
    surfX_0  = cell[0] + cell[2] + cell[3] + cell[9] + cell[10];
    surfX_P1 = cell[8] + cell[11] + cell[12] + cell[13] + cell[14];

    surfY_M1 = cell[2] + cell[4] + cell[5] + cell[13] + cell[14];
    surfY_P1 = cell[6] + cell[7] + cell[9] + cell[11] + cell[12];

    surfZ_M1 = cell[3] + cell[4] + cell[6] + cell[12] + cell[14];
    surfZ_P1 = cell[5] + cell[7] + cell[10] + cell[11] + cell[13];
  }

  static T computeRho(SpecializedCellBase const& cell)
  {
    T rho = cell[0] + cell[1] + cell[2] + cell[3] + cell[4]
            + cell[5] + cell[6] + cell[7] + cell[8]
            + cell[9] + cell[10] + cell[11] + cell[12]
            + cell[13] + cell[14] + (T)1;
    return rho;
  }

  static void computeRhoU(SpecializedCellBase const& cell, T& rho, T u[3])
  {
    T surfX_M1, surfX_0, surfX_P1,
    surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

    partial_rho(cell, surfX_M1, surfX_0, surfX_P1,
                surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

    rho = surfX_M1 + surfX_0 + surfX_P1 + (T)1;

    u[0]  = ( surfX_P1 - surfX_M1 ) / rho;
    u[1]  = ( surfY_P1 - surfY_M1 ) / rho;
    u[2]  = ( surfZ_P1 - surfZ_M1 ) / rho;
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

    surfX_M1 = cell[1] + cell[4] + cell[5] + cell[6] + cell[7];
    surfX_P1 = cell[8] + cell[11] + cell[12] + cell[13] + cell[14];

    surfY_M1 = cell[2] + cell[4] + cell[5] + cell[13] + cell[14];
    surfY_P1 = cell[6] + cell[7] + cell[9] + cell[11] + cell[12];

    surfZ_M1 = cell[3] + cell[4] + cell[6] + cell[12] + cell[14];
    surfZ_P1 = cell[5] + cell[7] + cell[10] + cell[11] + cell[13];

    j[0]  = ( surfX_P1 - surfX_M1 );
    j[1]  = ( surfY_P1 - surfY_M1 );
    j[2]  = ( surfZ_P1 - surfZ_M1 );
  }

  static void computeStress(SpecializedCellBase const& cell, T rho, const T u[3], T pi[6])
  {
    typedef descriptors::D3Q15<> L;
    using namespace util::tensorIndices3D;

    T surfX_M1, surfX_0, surfX_P1,
    surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

    partial_rho(cell, surfX_M1, surfX_0, surfX_P1,
                surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

    pi[xx] = surfX_P1+surfX_M1 - 1./descriptors::invCs2<T,L>()*(rho-(T)1) - rho*u[0]*u[0];
    pi[yy] = surfY_P1+surfY_M1 - 1./descriptors::invCs2<T,L>()*(rho-(T)1) - rho*u[1]*u[1];
    pi[zz] = surfZ_P1+surfZ_M1 - 1./descriptors::invCs2<T,L>()*(rho-(T)1) - rho*u[2]*u[2];

    pi[xy] =   cell[4] + cell[5] - cell[6]  - cell[7]
               + cell[11] + cell[12] - cell[13] - cell[14] - rho*u[0]*u[1];
    pi[xz] =   cell[4] - cell[5] + cell[6]  - cell[7]
               + cell[11] - cell[12] + cell[13] - cell[14] - rho*u[0]*u[2];
    pi[yz] =   cell[4] - cell[5] - cell[6]  + cell[7]
               + cell[11] - cell[12] - cell[13] + cell[14] - rho*u[1]*u[2];

  }

  static void computeAllMomenta(SpecializedCellBase const& cell, T& rho, T u[3], T pi[6])
  {
    typedef descriptors::D3Q15<> L;
    using namespace util::tensorIndices3D;

    T surfX_M1, surfX_0, surfX_P1,
    surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

    partial_rho(cell, surfX_M1, surfX_0, surfX_P1,
                surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

    rho = surfX_M1 + surfX_0 + surfX_P1 + (T)1;

    T rhoU0  = ( surfX_P1 - surfX_M1 ) / rho;
    T rhoU1  = ( surfY_P1 - surfY_M1 ) / rho;
    T rhoU2  = ( surfZ_P1 - surfZ_M1 ) / rho;
    u[0] = rhoU0 / rho;
    u[1] = rhoU1 / rho;
    u[2] = rhoU2 / rho;

    pi[xx] = surfX_P1+surfX_M1 - 1./descriptors::invCs2<T,L>()*(rho-(T)1) - rhoU0*u[0];
    pi[yy] = surfY_P1+surfY_M1 - 1./descriptors::invCs2<T,L>()*(rho-(T)1) - rhoU1*u[1];
    pi[zz] = surfZ_P1+surfZ_M1 - 1./descriptors::invCs2<T,L>()*(rho-(T)1) - rhoU2*u[2];

    pi[xy] =   cell[4] + cell[5] - cell[6]  - cell[7]
               + cell[11] + cell[12] - cell[13] - cell[14] - rhoU0*u[1];
    pi[xz] =   cell[4] - cell[5] + cell[6]  - cell[7]
               + cell[11] - cell[12] + cell[13] - cell[14] - rhoU0*u[2];
    pi[yz] =   cell[4] - cell[5] - cell[6]  + cell[7]
               + cell[11] - cell[12] - cell[13] + cell[14] - rhoU1*u[2];
  }

  static void modifyVelocity(SpecializedCellBase& cell, const T newU[3])
  {
    T rho, oldU[3];
    computeRhoU(cell, rho, oldU);
    const T oldUSqr = util::normSqr<T,3>(oldU);
    const T newUSqr = util::normSqr<T,3>(newU);
    for (int iPop=0; iPop<15; ++iPop) {
      cell[iPop] = cell[iPop]
                   - equilibrium(iPop, rho, oldU, oldUSqr)
                   + equilibrium(iPop, rho, newU, newUSqr);
    }
  }

};  //struct lbDynamicsHelpers<D3Q15DescriptorBase>


// Efficient specialization for D3Q15 lattice and for forced D3Q15 lattice
//   (operations applying to the whole lattice)

template<typename T>
struct lbLatticeHelpers<T, descriptors::D3Q15<>> {
  static void swapAndStreamCell (
    Cell<T,descriptors::D3Q15<>> ***grid,
    int iX, int iY, int iZ, int nX, int nY, int nZ, int iPop, T& fTmp )
  {
    fTmp                     = grid[iX][iY][iZ][iPop];
    grid[iX][iY][iZ][iPop]   = grid[iX][iY][iZ][iPop+7];
    grid[iX][iY][iZ][iPop+7] = grid[nX][nY][nZ][iPop];
    grid[nX][nY][nZ][iPop]   = fTmp;
  }

  static void swapAndStream3D(Cell<T,descriptors::D3Q15<>> ***grid,
                              int iX, int iY, int iZ)
  {
    T fTmp;
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY,   iZ,   1, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY-1, iZ,   2, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY  , iZ-1, 3, fTmp);

    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ-1, 4, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ+1, 5, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ-1, 6, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ+1, 7, fTmp);
  }

};

template<typename T>
struct lbLatticeHelpers<T, descriptors::D3Q15<descriptors::FORCE>> {

  static void swapAndStreamCell (
    Cell<T,descriptors::D3Q15<>> ***grid,
    int iX, int iY, int iZ, int nX, int nY, int nZ, int iPop, T& fTmp )
  {
    fTmp                     = grid[iX][iY][iZ][iPop];
    grid[iX][iY][iZ][iPop]   = grid[iX][iY][iZ][iPop+7];
    grid[iX][iY][iZ][iPop+7] = grid[nX][nY][nZ][iPop];
    grid[nX][nY][nZ][iPop]   = fTmp;
  }

  static void swapAndStream3D(Cell<T,descriptors::D3Q15<>> ***grid,
                              int iX, int iY, int iZ)
  {
    T fTmp;
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY,   iZ,   1, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY-1, iZ,   2, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX,   iY  , iZ-1, 3, fTmp);

    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ-1, 4, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY-1, iZ+1, 5, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ-1, 6, fTmp);
    swapAndStreamCell(grid, iX, iY, iZ, iX-1, iY+1, iZ+1, 7, fTmp);
  }

};


}  // namespace olb

#endif
