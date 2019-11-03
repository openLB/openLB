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

#ifndef FINITE_DIFFERENCE_3D_H
#define FINITE_DIFFERENCE_3D_H

#include "finiteDifference.h"

namespace olb {

namespace fd {

template<typename T, typename DESCRIPTOR,
         int direction, int orientation, int deriveDirection,
         bool orthogonal>
struct DirectedGradients3D {
  static void interpolateVector( T velDeriv[DESCRIPTOR::d],
                                 BlockLattice3D<T,DESCRIPTOR> const& blockLattice,
                                 int iX, int iY, int iZ );
  static void interpolateScalar( T& rhoDeriv,
                                 BlockLattice3D<T,DESCRIPTOR> const& blockLattice,
                                 int iX, int iY, int iZ );
};

// Implementation for orthogonal==true; i.e. the derivative is along the boundary normal.
template<typename T, typename DESCRIPTOR,
         int direction, int orientation, int deriveDirection>
struct DirectedGradients3D<T, DESCRIPTOR, direction, orientation,
         deriveDirection, true> {
  static void interpolateVector(T velDeriv[DESCRIPTOR::d],
                                BlockLattice3D<T,DESCRIPTOR> const& blockLattice,
                                int iX, int iY, int iZ)
  {
    using namespace fd;

    T u0[DESCRIPTOR::d], u1[DESCRIPTOR::d], u2[DESCRIPTOR::d];

    blockLattice.get(iX,iY,iZ).computeU(u0);
    blockLattice.get (
      iX+(direction==0 ? (-orientation):0),
      iY+(direction==1 ? (-orientation):0),
      iZ+(direction==2 ? (-orientation):0)  ).computeU(u1);
    blockLattice.get (
      iX+(direction==0 ? (-2*orientation):0),
      iY+(direction==1 ? (-2*orientation):0),
      iZ+(direction==2 ? (-2*orientation):0) ).computeU(u2);

    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      velDeriv[iD] = -orientation * boundaryGradient(u0[iD], u1[iD], u2[iD]);
    }
  }

  static void interpolateScalar(T& rhoDeriv,
                                BlockLattice3D<T,DESCRIPTOR> const& blockLattice,
                                int iX, int iY, int iZ)
  {
    using namespace fd;

    // note that the derivative runs along direction.
    T rho0 = blockLattice.get(iX,iY,iZ).computeRho();
    T rho1 = blockLattice.get (
               iX+(direction==0 ? (-orientation):0),
               iY+(direction==1 ? (-orientation):0),
               iZ+(direction==2 ? (-orientation):0)  ).computeRho();
    T rho2 = blockLattice.get (
               iX+(direction==0 ? (-2*orientation):0),
               iY+(direction==1 ? (-2*orientation):0),
               iZ+(direction==2 ? (-2*orientation):0) ).computeRho();

    rhoDeriv = -orientation * boundaryGradient(rho0, rho1, rho2);
  }
};

// Implementation for orthogonal==false; i.e. the derivative is aligned with the boundary.
template<typename T, typename DESCRIPTOR,
         int direction, int orientation, int deriveDirection>
struct DirectedGradients3D<T, DESCRIPTOR, direction, orientation,
         deriveDirection, false> {
  static void  interpolateVector(T velDeriv[DESCRIPTOR::d],
                                 BlockLattice3D<T,DESCRIPTOR> const& blockLattice,
                                 int iX, int iY, int iZ)
  {
    using namespace fd;

    T u_p1[DESCRIPTOR::d], u_m1[DESCRIPTOR::d];

    blockLattice.get (
      iX+(deriveDirection==0 ? 1:0),
      iY+(deriveDirection==1 ? 1:0),
      iZ+(deriveDirection==2 ? 1:0) ).computeU(u_p1);

    blockLattice.get (
      iX+(deriveDirection==0 ? (-1):0),
      iY+(deriveDirection==1 ? (-1):0),
      iZ+(deriveDirection==2 ? (-1):0) ).computeU(u_m1);

    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      velDeriv[iD] = centralGradient(u_p1[iD],u_m1[iD]);
    }
  }

  static void  interpolateScalar(T& rhoDeriv,
                                 BlockLattice3D<T,DESCRIPTOR> const& blockLattice,
                                 int iX, int iY, int iZ)
  {
    using namespace fd;

    T rho_p1 = blockLattice.get (
                 iX+(deriveDirection==0 ? 1:0),
                 iY+(deriveDirection==1 ? 1:0),
                 iZ+(deriveDirection==2 ? 1:0) ).computeRho();

    T rho_m1 = blockLattice.get (
                 iX+(deriveDirection==0 ? (-1):0),
                 iY+(deriveDirection==1 ? (-1):0),
                 iZ+(deriveDirection==2 ? (-1):0) ).computeRho();

    rhoDeriv = centralGradient(rho_p1, rho_m1);

  }
};

}  // namespace fd

}  // namespace olb


#endif
