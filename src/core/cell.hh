/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2007 Jonas Latt
 *                2015-2019 Mathias J. Krause, Adrian Kummerlaender
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
 * Definition of a LB cell -- generic implementation.
 */
#ifndef CELL_HH
#define CELL_HH

#include <algorithm>
#include "cell.h"
#include "util.h"
#include "dynamics/lbHelpers.h"

namespace olb {

////////////////////////// Class Cell /////////////////////////////

/** The possibility to default construct Cell objects facilitates
 * their use in various types of containers. However, they can not
 * be used directly after construction; the method defineDynamics()
 * must be used first.
 */
template<typename T, typename DESCRIPTOR>
Cell<T,DESCRIPTOR>::Cell()
  : dynamics(nullptr)
{
  initializeData();
}

/** This constructor initializes the dynamics, but not the values
 * of the distribution functions. Remember that the dynamics is not
 * owned by the Cell object, the user must ensure its proper
 * destruction and a sufficient life time.
 */
template<typename T, typename DESCRIPTOR>
Cell<T,DESCRIPTOR>::Cell(Dynamics<T,DESCRIPTOR>* dynamics_)
  : dynamics(dynamics_)
{
  initializeData();
}

template<typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::defineDynamics(Dynamics<T,DESCRIPTOR>* dynamics_)
{
  dynamics = dynamics_;
}

template<typename T, typename DESCRIPTOR>
Dynamics<T,DESCRIPTOR> const* Cell<T,DESCRIPTOR>::getDynamics() const
{
  OLB_PRECONDITION(dynamics);
  return dynamics;
}

template<typename T, typename DESCRIPTOR>
Dynamics<T,DESCRIPTOR>* Cell<T,DESCRIPTOR>::getDynamics()
{
  OLB_PRECONDITION(dynamics);
  return dynamics;
}

template<typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::computeFeq(T fEq[descriptors::q<DESCRIPTOR>()]) const
{
  T rho{};
  T u[2] {};
  computeRhoU(rho, u);
  const T uSqr = u[0]*u[0] + u[1]*u[1];
  for (int iPop=0; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
    fEq[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
  }
}

template <typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::computeFneq(T fNeq[descriptors::q<DESCRIPTOR>()]) const
{
  T rho{};
  T u[2] {};
  computeRhoU(rho, u);
  lbHelpers<T,DESCRIPTOR>::computeFneq(*this, fNeq, rho, u);
}

template<typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::revert()
{
  for (int iPop=1; iPop<=descriptors::q<DESCRIPTOR>()/2; ++iPop) {
    std::swap(this->data[iPop],this->data[iPop+descriptors::q<DESCRIPTOR>()/2]);
  }
}

template<typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::initializeData()
{
  for (int iPop=0; iPop<DESCRIPTOR::size(); ++iPop) {
    this->data[iPop] = T();
  }
}

template<typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::serialize(T* data) const
{
  for (int i=0; i < DESCRIPTOR::size(); ++i) {
    data[i] = this->data[i];
  }
}

template<typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::unSerialize(T const* data)
{
  for (int i=0; i < DESCRIPTOR::size(); ++i) {
    this->data[i] = data[i];
  }
}


template<typename T, typename DESCRIPTOR>
std::size_t Cell<T,DESCRIPTOR>::getSerializableSize() const
{
  return sizeof(T) * DESCRIPTOR::size();
}


template<typename T, typename DESCRIPTOR>
bool* Cell<T,DESCRIPTOR>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  this->registerVar(iBlock, sizeBlock, currentBlock, dataPtr, this->data[0], DESCRIPTOR::size());

  return dataPtr;
}

}  // namespace olb


#endif
