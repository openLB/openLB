/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#ifndef REDUCTION_F_2D_HH
#define REDUCTION_F_2D_HH

#include "functors/lattice/reductionF2D.h"
#include "core/superLattice2D.h"
#include "dynamics/lbHelpers.h"

namespace olb {


template <typename T, typename DESCRIPTOR>
SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR>::SuperLatticeFfromAnalyticalF2D(
  FunctorPtr<AnalyticalF2D<T,T>>&& f,
  SuperLattice2D<T, DESCRIPTOR>&   sLattice)
  : SuperLatticeF2D<T, DESCRIPTOR>(sLattice, f->getTargetDim()),
    _f(std::move(f))
{
  this->getName() = "fromAnalyticalF(" + _f->getName() + ")";

  LoadBalancer<T>&     load   = sLattice.getLoadBalancer();
  CuboidGeometry2D<T>& cuboid = sLattice.getCuboidGeometry();

  for (int iC = 0; iC < load.size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticeFfromAnalyticalF2D<T,DESCRIPTOR>(
        *_f,
        sLattice.getExtendedBlockLattice(iC),
        cuboid.get(load.glob(iC)))
    );
  }
}

template <typename T, typename DESCRIPTOR>
bool SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR>::operator()(T output[], const int input[])
{
  T physR[2] = {};
  this->_sLattice.getCuboidGeometry().getPhysR(physR,input);
  return _f(output,physR);
}


template<typename T, typename DESCRIPTOR>
BlockLatticeFfromAnalyticalF2D<T, DESCRIPTOR>::BlockLatticeFfromAnalyticalF2D(
  AnalyticalF2D<T, T>&                    f,
  BlockLatticeStructure2D<T, DESCRIPTOR>& lattice,
  Cuboid2D<T>&                            cuboid)
  : BlockLatticeF2D<T, DESCRIPTOR>(lattice, f.getTargetDim()),
    _f(f),
    _cuboid(cuboid)
{
  this->getName() = "blockFfromAnalyticalF(" + _f.getName() + ")";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeFfromAnalyticalF2D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  T physR[2] = {};
  _cuboid.getPhysR(physR,input);
  return _f(output,physR);
}


} // end namespace olb

#endif
