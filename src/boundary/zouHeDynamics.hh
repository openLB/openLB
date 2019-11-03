/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006,2007 Orestis Malaspinas and Jonas Latt
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

#ifndef ZOU_HE_DYNAMICS_HH
#define ZOU_HE_DYNAMICS_HH

#include "zouHeDynamics.h"
#include "dynamics/latticeDescriptors.h"
#include "core/util.h"
#include "dynamics/lbHelpers.h"
#include <cmath>


namespace olb {

 

template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
ZouHeDynamics<T,DESCRIPTOR,Dynamics,direction,orientation>::ZouHeDynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_ )
  : BasicDynamics<T,DESCRIPTOR>(momenta_),
    boundaryDynamics(omega_, momenta_)
{ }

template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
T ZouHeDynamics<T,DESCRIPTOR, Dynamics, direction, orientation>::
computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return boundaryDynamics.computeEquilibrium(iPop, rho, u, uSqr);
}

template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
void ZouHeDynamics<T,DESCRIPTOR,Dynamics,direction,orientation>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  typedef lbHelpers<T,DESCRIPTOR> lbH;
  typedef DESCRIPTOR L;

  // Along all the commented parts of this code there will be an example based
  // on the situation where the wall's normal vector if (0,1) and the
  // numerotation of the velocites are done according to the D2Q9
  // lattice of the OpenLB library.

  // Find all the missing populations
  // (directions 3,4,5)
  std::vector<int> missingIndexes = util::subIndexOutgoing<L,direction,orientation>();

  // Will contain the missing poputations that are not normal to the wall.
  // (directions 3,5)
  std::vector<int> missingDiagonalIndexes = missingIndexes;
  for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
    int numOfNonNullComp = 0;
    for (int iDim = 0; iDim < L::d; ++iDim) {
      numOfNonNullComp += abs(descriptors::c<L>(missingIndexes[iPop],iDim));
    }

    if (numOfNonNullComp == 1) {
      missingDiagonalIndexes.erase(missingDiagonalIndexes.begin()+iPop);
      break;
    }
  }

  T rho, u[L::d];
  T falseRho, falseU[L::d];
  this->_momenta.computeRhoU(cell, rho, u);

  T uSqr = util::normSqr<T,L::d>(u);

  // The unknown non equilibrium populations are bounced back
  // (f[3] = feq[3] + fneq[7], f[4] = feq[4] + fneq[8],
  //  f[5] = feq[5] + fneq[1])
  for (unsigned iPop = 0; iPop < missingIndexes.size(); ++iPop) {
    cell[missingIndexes[iPop]] = cell[util::opposite<L>(missingIndexes[iPop])]
                                 - computeEquilibrium(util::opposite<L>(missingIndexes[iPop]), rho, u, uSqr)
                                 + computeEquilibrium(missingIndexes[iPop], rho, u, uSqr);
  }

  // We recompute rho and u in order to have the new momentum and density. Since
  // the momentum is not conserved from this scheme, we will corect it. By adding
  // a contribution to the missingDiagonalVelocities.
  lbH::computeRhoU(cell,falseRho,falseU);

  T diff[L::d];
  for (int iDim = 0; iDim < L::d; ++iDim) {
    diff[iDim] = (rho*u[iDim] - falseRho*falseU[iDim])/ (T)missingDiagonalIndexes.size();
  }

  for (unsigned iPop = 0; iPop < missingDiagonalIndexes.size(); ++iPop) {
    for (int iDim = 1; iDim < L::d; ++iDim) {
      cell[missingDiagonalIndexes[iPop]] +=
        descriptors::c<L>(missingDiagonalIndexes[iPop],(direction+iDim)%L::d) * diff[(direction+iDim)%L::d];
    }
  }

  boundaryDynamics.collide(cell, statistics);

  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
T ZouHeDynamics<T,DESCRIPTOR,Dynamics,direction,orientation>::getOmega() const
{
  return boundaryDynamics.getOmega();
}

template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
void ZouHeDynamics<T,DESCRIPTOR,Dynamics,direction,orientation>::setOmega(T omega_)
{
  boundaryDynamics.setOmega(omega_);
}


}  // namespace olb

#endif
