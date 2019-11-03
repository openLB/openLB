/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, Orestis Malaspinas and Jonas Latt
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

#ifndef INAMURO_ANALYTICAL_DYNAMICS_HH
#define INAMURO_ANALYTICAL_DYNAMICS_HH

#include "inamuroAnalyticalDynamics.h"
#include "dynamics/latticeDescriptors.h"
#include "core/util.h"
#include "dynamics/lbHelpers.h"
#include <cmath>

namespace olb {

template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
InamuroAnalyticalDynamics<T,DESCRIPTOR,Dynamics,direction,orientation>::InamuroAnalyticalDynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_ )
  : BasicDynamics<T,DESCRIPTOR>(momenta_),
    boundaryDynamics(omega_, momenta_)
{ }

template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
T InamuroAnalyticalDynamics<T,DESCRIPTOR, Dynamics, direction, orientation>::
computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return boundaryDynamics.computeEquilibrium(iPop, rho, u, uSqr);
}

template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
void InamuroAnalyticalDynamics<T,DESCRIPTOR,Dynamics,direction,orientation>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  typedef DESCRIPTOR L;

  // Along all the commented parts of this code there will be an example based
  // on the situation where the wall's normal vector if (0,1) and the
  // numerotation of the velocites are done according to the D2Q9
  // lattice of the OpenLB library.

  // Find all the missing populations
  // (directions 3,4,5)
  std::vector<int> missInd =
    util::subIndexOutgoing<L,direction,orientation>();

  // Will contain the missing poputations that are not normal to the wall.
  // (directions 3,5)
  std::vector<int> missDiagInd = missInd;

  for (unsigned iPop = 0; iPop < missInd.size(); ++iPop) {
    int numOfNonNullComp = 0;
    for (int iDim = 0; iDim < L::d; ++iDim) {
      numOfNonNullComp += abs(descriptors::c<L>(missInd[iPop],iDim));
    }

    if (numOfNonNullComp == 1) {
      missDiagInd.erase(missDiagInd.begin()+iPop);
      break;
    }
  }

  // Will contain the populations normal to the wall's normal vector.
  // (directions 2,6)
  std::vector<int> perpInd = util::subIndex<L,direction,0>();
  for (unsigned iPop = 0; iPop < perpInd.size(); ++iPop) {
    if (descriptors::c<L>(perpInd[iPop],0) == 0 && descriptors::c<L>(perpInd[iPop],1) == 0) {
      perpInd.erase(perpInd.begin() + iPop);
      break;
    }
  }

  T rho, u[L::d];
  this->_momenta.computeRhoU(cell, rho, u);

  T rhoCs = T();
  T uCs[L::d];
  for (int iDim = 0; iDim < L::d; ++iDim) {
    uCs[iDim] = T();
  }

  T fSum = T();
  for (unsigned iPop = 0; iPop < missInd.size(); ++iPop) {
    fSum += cell[util::opposite<L>(missInd[iPop])];
  }
  // do not forget the "+1" in the rhoCs equation in the numerator (it's
  // here because fEq = usualfEq - t[i]
  rhoCs = ((T)6 * (-orientation * rho * u[direction] + fSum) + (T)1) /
          ((T)3 * u[direction] * u[direction] - orientation * (T)3 * u[direction] + (T)1);

  T fDiffPerp = T();
  for (unsigned iPop = 0; iPop < perpInd.size(); ++iPop) {
    fDiffPerp += descriptors::c<L>(perpInd[iPop],(direction + 1)%2) * cell[perpInd[iPop]];
  }
  fDiffPerp *= orientation;

  T fDiffDiag = T();
  for (unsigned iPop = 0; iPop < missDiagInd.size(); ++iPop)
    fDiffDiag += descriptors::c<L>(util::opposite<L>(missDiagInd[iPop]),(direction + 1)%2)
                 * cell[util::opposite<L>(missDiagInd[iPop])];
  fDiffDiag *= orientation;

  uCs[(direction + 1)%L::d] = (
                                - orientation * (T)6 * rho * u[(direction+1)%L::d]
                                + orientation * rhoCs * u[(direction+1)%L::d]
                                - (T)3 * rhoCs * u[direction]*u[(direction+1)%L::d]
                                + (T)6*(fDiffPerp + fDiffDiag))
                              / (
                                rhoCs * (-orientation + (T)3 * u[direction]));

  for (int iDim = 0; iDim < L::d; ++iDim) {
    uCs[iDim] += u[iDim];
  }

  T uSqr = util::normSqr<T,L::d>(uCs);

  for (unsigned iPop = 0; iPop < missInd.size(); ++iPop) {
    cell[missInd[iPop]] = computeEquilibrium(missInd[iPop], rhoCs, uCs, uSqr);
  }

  boundaryDynamics.collide(cell, statistics);
}

template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
T InamuroAnalyticalDynamics<T,DESCRIPTOR,Dynamics,direction,orientation>::getOmega() const
{
  return boundaryDynamics.getOmega();
}

template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
void InamuroAnalyticalDynamics<T,DESCRIPTOR,Dynamics,direction,orientation>::setOmega(T omega_)
{
  boundaryDynamics.setOmega(omega_);
}


}  // namespace olb


#endif
