/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#ifndef SUPER_LATTICE_REFINEMENT_METRIC_F_2D_HH
#define SUPER_LATTICE_REFINEMENT_METRIC_F_2D_HH

#include "superLatticeRefinementMetricF2D.h"
#include "blockLatticeRefinementMetricF2D.h"


namespace olb {


template<typename T, typename DESCRIPTOR>
SuperLatticeKnudsen2D<T, DESCRIPTOR>::SuperLatticeKnudsen2D(
  SuperLattice2D<T, DESCRIPTOR>& lattice)
  : SuperLatticeF2D<T, DESCRIPTOR>(lattice, 1)
{
  this->getName() = "knudsen";

  for (int iC = 0; iC < lattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticeKnudsen2D<T, DESCRIPTOR>(lattice.getBlockLattice(iC))
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeKnudsen2D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  auto& load = this->_sLattice.getLoadBalancer();

  if (load.isLocal(input[0])) {
    const int loc = load.loc(input[0]);

    return this->getBlockF(loc)(output, &input[1]);
  }
  else {
    return false;
  }
}


template<typename T, typename DESCRIPTOR>
SuperLatticeRefinementMetricKnudsen2D<T, DESCRIPTOR>::SuperLatticeRefinementMetricKnudsen2D(
  SuperLattice2D<T, DESCRIPTOR>&      lattice,
  const UnitConverter<T, DESCRIPTOR>& converter)
  : SuperLatticeF2D<T, DESCRIPTOR>(lattice, 1)
{
  this->getName() = "refinementMetricKnudsen";

  for (int iC = 0; iC < lattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticeRefinementMetricKnudsen2D<T, DESCRIPTOR>(
        lattice.getBlockLattice(iC), converter)
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeRefinementMetricKnudsen2D<T, DESCRIPTOR>::operator()(T output[], int glob)
{
  auto& load = this->_sLattice.getLoadBalancer();

  if (load.isLocal(glob)) {
    const int loc = load.loc(glob);

    return this->getBlockF(loc)(output);
  }
  else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeRefinementMetricKnudsen2D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  auto& load = this->_sLattice.getLoadBalancer();

  if (load.isLocal(input[0])) {
    const int loc = load.loc(input[0]);

    return this->getBlockF(loc)(output, &input[1]);
  }
  else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLatticeRefinementMetricKnudsen2D<T, DESCRIPTOR>::print()
{
  const int nC = this->_sLattice.getCuboidGeometry().getNc();

  std::vector<T> factors(nC, T{});
  T output[1] = { };

  auto& load = this->_sLattice.getLoadBalancer();

  for (int iC = 0; iC < load.size(); ++iC) {
    this->getBlockF(iC)(output, iC);

    factors[load.glob(iC)] = output[0];
  }

  OstreamManager clout(std::cout, "refinement");

  for (int i = 0; i < nC; ++i) {
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().sendToMaster(&factors[i], 1, load.isLocal(i));
#endif

    clout << "factors[" << i << "]: " << factors[i] << std::endl;
  }
}


}

#endif
