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
 * functions of the header file lbHelpers.h, for the D2Q9 grid.
 */

#ifndef MRT_HELPERS_2D_H
#define MRT_HELPERS_2D_H

namespace olb {

 

// Efficient specialization for D2Q9 lattice
template<typename T>
struct mrtHelpers<T,descriptors::MRTD2Q9Descriptor> {
  using DESCRIPTOR = descriptors::MRTD2Q9Descriptor;

  /// Computation of all equilibrium distribution (in momenta space)
  static void computeEquilibrium( T momentaEq[9],
                                  T rho, const T u[2],
                                  const T uSqr )
  {
    // momentaEq[0] = rho;
    momentaEq[1] = rho*(-(T)2 + (T)3*uSqr);
    momentaEq[2] = rho*((T)1 - (T)3*uSqr);
    // momentaEq[3] = rho*u[0];
    momentaEq[4] = rho*-u[0];
    // momentaEq[5] = rho*u[1];
    momentaEq[6] = rho*-u[1];
    momentaEq[7] = rho*(u[0]*u[0] - u[1]*u[1]);
    momentaEq[8] = rho*u[0]*u[1];
  }

  /// Computation of all momenta (specialized for d2q9)
  static void computeMomenta(T momenta[9], Cell<T,DESCRIPTOR> &cell)
  {
    //         momenta[0] = cell[0] + cell[1] + cell[2] + cell[3] +
    //                 cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + (T)1;

    momenta[1] = -(T)4*cell[0] +(T)2*cell[1] - cell[2] + (T)2*cell[3] - cell[4] +
                 (T)2*cell[5] - cell[6] + (T)2*cell[7] - cell[8] - (T)2;

    momenta[2] = (T)4*cell[0] + cell[1] - (T)2*cell[2] + cell[3] - (T)2*cell[4] +
                 cell[5] - (T)2*cell[6] + cell[7] - (T)2*cell[8] + (T)1;

    //         momenta[3] = - cell[1] - cell[2] - cell[3] +
    //                 cell[5] + cell[6] + cell[7];

    momenta[4] = - cell[1] + (T)2*cell[2] - cell[3] +
                 cell[5] - (T)2*cell[6] + cell[7];

    //         momenta[5] = cell[1] - cell[3] - cell[4] -
    //                 cell[5] + cell[7] + cell[8];

    momenta[6] = cell[1] - cell[3] + (T)2*cell[4] -
                 cell[5] + cell[7] - (T)2*cell[8];

    momenta[7] = cell[2] - cell[4] + cell[6] - cell[8];

    momenta[8] = - cell[1] + cell[3] - cell[5] + cell[7];
  }

  /// MRT collision step
  static T mrtCollision( Cell<T,DESCRIPTOR>& cell,
                         T rho, const T u[2],
                         T invM_S[9][9] )
  {
    T uSqr = util::normSqr<T,2>(u);
    T momenta[9];
    T momentaEq[9];

    computeMomenta(momenta,cell);
    computeEquilibrium(momentaEq,rho,u,uSqr);

    T mom1 = momenta[1] - momentaEq[1];
    T mom2 = momenta[2] - momentaEq[2];
    T mom4 = momenta[4] - momentaEq[4];
    T mom6 = momenta[6] - momentaEq[6];
    T mom7 = momenta[7] - momentaEq[7];
    T mom8 = momenta[8] - momentaEq[8];

    cell[0] -= invM_S[0][1]*mom1 +
               invM_S[0][2]*mom2;

    cell[1] -= invM_S[1][1]*mom1 +
               invM_S[1][2]*mom2 +
               invM_S[1][4]*mom4 +
               invM_S[1][6]*mom6 +
               invM_S[1][8]*mom8;

    cell[2] -= invM_S[2][1]*mom1 +
               invM_S[2][2]*mom2 +
               invM_S[2][4]*mom4 +
               invM_S[2][7]*mom7;

    cell[3] -= invM_S[3][1]*mom1 +
               invM_S[3][2]*mom2 +
               invM_S[3][4]*mom4 +
               invM_S[3][6]*mom6 +
               invM_S[3][8]*mom8;

    cell[4] -= invM_S[4][1]*mom1 +
               invM_S[4][2]*mom2 +
               invM_S[4][6]*mom6 +
               invM_S[4][7]*mom7;

    cell[5] -= invM_S[5][1]*mom1 +
               invM_S[5][2]*mom2 +
               invM_S[5][4]*mom4 +
               invM_S[5][6]*mom6 +
               invM_S[5][8]*mom8;

    cell[6] -= invM_S[6][1]*mom1 +
               invM_S[6][2]*mom2 +
               invM_S[6][4]*mom4 +
               invM_S[6][7]*mom7;

    cell[7] -= invM_S[7][1]*mom1 +
               invM_S[7][2]*mom2 +
               invM_S[7][4]*mom4 +
               invM_S[7][6]*mom6 +
               invM_S[7][8]*mom8;

    cell[8] -= invM_S[8][1]*mom1 +
               invM_S[8][2]*mom2 +
               invM_S[8][6]*mom6 +
               invM_S[8][7]*mom7;

    return uSqr;
  }

  /// MRT collision step
  static T mrtSGSCollision( Cell<T,DESCRIPTOR>& cell,
                            T rho, const T u[2], T omega,
                            T invM_S_SGS[9][9] )
  {
    T uSqr = util::normSqr<T,2>(u);
    T momenta[9];
    T momentaEq[9];

    computeMomenta(momenta,cell);
    computeEquilibrium(momentaEq,rho,u,uSqr);

    T mom1 = momenta[1] - momentaEq[1];
    T mom2 = momenta[2] - momentaEq[2];
    T mom4 = momenta[4] - momentaEq[4];
    T mom6 = momenta[6] - momentaEq[6];
    T mom7 = momenta[7] - momentaEq[7];
    T mom8 = momenta[8] - momentaEq[8];

    cell[0] -= invM_S_SGS[0][1]*mom1 +
               invM_S_SGS[0][2]*mom2;

    cell[1] -= invM_S_SGS[1][1]*mom1/ +
               invM_S_SGS[1][2]*mom2 +
               invM_S_SGS[1][4]*mom4 +
               invM_S_SGS[1][6]*mom6 +
               invM_S_SGS[1][8]*mom8;

    cell[2] -= invM_S_SGS[2][1]*mom1 +
               invM_S_SGS[2][2]*mom2 +
               invM_S_SGS[2][4]*mom4 +
               invM_S_SGS[2][7]*mom7;

    cell[3] -= invM_S_SGS[3][1]*mom1 +
               invM_S_SGS[3][2]*mom2 +
               invM_S_SGS[3][4]*mom4 +
               invM_S_SGS[3][6]*mom6 +
               invM_S_SGS[3][8]*mom8;

    cell[4] -= invM_S_SGS[4][1]*mom1 +
               invM_S_SGS[4][2]*mom2 +
               invM_S_SGS[4][6]*mom6 +
               invM_S_SGS[4][7]*mom7;

    cell[5] -= invM_S_SGS[5][1]*mom1 +
               invM_S_SGS[5][2]*mom2 +
               invM_S_SGS[5][4]*mom4 +
               invM_S_SGS[5][6]*mom6 +
               invM_S_SGS[5][8]*mom8;

    cell[6] -= invM_S_SGS[6][1]*mom1 +
               invM_S_SGS[6][2]*mom2 +
               invM_S_SGS[6][4]*mom4 +
               invM_S_SGS[6][7]*mom7;

    cell[7] -= invM_S_SGS[7][1]*mom1 +
               invM_S_SGS[7][2]*mom2 +
               invM_S_SGS[7][4]*mom4 +
               invM_S_SGS[7][6]*mom6 +
               invM_S_SGS[7][8]*mom8;

    cell[8] -= invM_S_SGS[8][1]*mom1 +
               invM_S_SGS[8][2]*mom2 +
               invM_S_SGS[8][6]*mom6 +
               invM_S_SGS[8][7]*mom7;

    return uSqr;
  }

};


}  // namespace olb

#endif
