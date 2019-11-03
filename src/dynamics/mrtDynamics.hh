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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef MRT_DYNAMICS_HH
#define MRT_DYNAMICS_HH

#include <algorithm>
#include <limits>
#include "mrtHelpers.h"

namespace olb {

//==============================================================================//
/////////////////////////// Class MRTdynamics ///////////////////////////////
//==============================================================================//
/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param lambda_ will be used as an
 */

// Original implementation based on:
// D'Humieres et al., "Multiple-relaxation-time lattice Boltzmann models in three dimensions",
// Phil: Trans. R. soc. Lond. A (2002) 360, 437-451
// and
// Yu et al,, "LES of turbulent square jet flow using an MRT lattice Boltzmann model",
// Computers & Fluids 35 (2006), 957-965
template<typename T, typename DESCRIPTOR>
MRTdynamics<T,DESCRIPTOR>::MRTdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_ )
  : BasicDynamics<T,DESCRIPTOR>(momenta_), omega(omega_), lambda(omega_)
{
  T rt[DESCRIPTOR::q]; // relaxation times vector.
  for (int iPop  = 0; iPop < DESCRIPTOR::q; ++iPop) {
    rt[iPop] = descriptors::s<T,DESCRIPTOR>(iPop);
  }
  for (int iPop  = 0; iPop < descriptors::shearIndexes<DESCRIPTOR>(); ++iPop) {
    rt[descriptors::shearViscIndexes<DESCRIPTOR>(iPop)] = omega;
  }
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
      invM_S[iPop][jPop] = T();
      for (int kPop = 0; kPop < DESCRIPTOR::q; ++kPop) {
        if (kPop == jPop) {
          invM_S[iPop][jPop] += descriptors::invM<T,DESCRIPTOR>(iPop,kPop) *
                                rt[kPop];
        }
      }
    }
  }

}

template<typename T, typename DESCRIPTOR>
T MRTdynamics<T,DESCRIPTOR>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
}

template<typename T, typename DESCRIPTOR>
void MRTdynamics<T,DESCRIPTOR>::computeAllEquilibrium(T momentaEq[DESCRIPTOR::q],
                                                T rho, const T u[DESCRIPTOR::d],
                                                const T uSqr)
{
  mrtHelpers<T,DESCRIPTOR>::computeEquilibrium(momentaEq, rho, u, uSqr);
}

template<typename T, typename DESCRIPTOR>
void MRTdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  typedef DESCRIPTOR L;
  typedef mrtHelpers<T,DESCRIPTOR> mrtH;

  T rho, u[L::d];
  this->_momenta.computeRhoU(cell, rho, u);

  T uSqr = mrtH::mrtCollision(cell,rho,u,invM_S);

  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T MRTdynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void MRTdynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  omega = omega_;
}

template<typename T, typename DESCRIPTOR>
T MRTdynamics<T,DESCRIPTOR>::getLambda() const
{
  return lambda;
}

template<typename T, typename DESCRIPTOR>
void MRTdynamics<T,DESCRIPTOR>::setLambda(T lambda_)
{
  lambda = lambda_;
}




template<typename T, typename DESCRIPTOR>
MRTdynamics2<T,DESCRIPTOR>::MRTdynamics2 (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_ )
  : MRTdynamics<T,DESCRIPTOR>(omega_, momenta_)
{
  T rt[DESCRIPTOR::q]; // relaxation times vector.
  for (int iPop  = 0; iPop < DESCRIPTOR::q; ++iPop) {
    rt[iPop] = descriptors::s_2<T,DESCRIPTOR>(iPop);
  }
  for (int iPop  = 0; iPop < descriptors::shearIndexes<DESCRIPTOR>(); ++iPop) {
    rt[descriptors::shearViscIndexes<DESCRIPTOR>(iPop)] = omega;
  }
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
      invM_S_2[iPop][jPop] = T();
      for (int kPop = 0; kPop < DESCRIPTOR::q; ++kPop) {
        if (kPop == jPop) {
          invM_S_2[iPop][jPop] += descriptors::invM<T,DESCRIPTOR>(iPop,kPop) *
                                  rt[kPop];
        }
      }
    }
  }
}


// Stabalized MRT scheme with uniform relaxation times
template<typename T, typename DESCRIPTOR>
void MRTdynamics2<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  typedef mrtHelpers<T,DESCRIPTOR> mrtH;

  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);

  T uSqr = mrtH::mrtCollision(cell,rho,u,invM_S_2);

  statistics.incrementStats(rho, uSqr);
}



template<typename T, typename DESCRIPTOR>
ForcedMRTdynamics<T,DESCRIPTOR>::ForcedMRTdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_ )
  : MRTdynamics<T,DESCRIPTOR>(omega_, momenta_)
{
}

template<typename T, typename DESCRIPTOR>
void ForcedMRTdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);

  T uSqr = mrtHelpers<T,DESCRIPTOR>::mrtCollision(cell, rho, u, this->invM_S);
  mrtHelpers<T,DESCRIPTOR>::addExternalForce(cell, rho, u, this->invM_S);

  statistics.incrementStats(rho, uSqr);
}


} // end namespace

#endif

