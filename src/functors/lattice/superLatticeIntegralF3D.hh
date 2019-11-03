/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Fabian Klemens, Benjamin Förster, Adrian Kummerlaender
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

#ifndef SUPER_LATTICE_INTEGRAL_F_3D_HH
#define SUPER_LATTICE_INTEGRAL_F_3D_HH

#include <cmath>
#include <algorithm>

#include "superLatticeIntegralF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.hh"
#include "utilities/vectorHelpers.h"
#include "io/ostreamManager.h"
#include "utilities/functorPtr.hh"

using namespace olb::util;

namespace olb {


template<typename T, typename DESCRIPTOR>
SuperLatticePhysDrag3D<T, DESCRIPTOR>::SuperLatticePhysDrag3D(
  SuperLattice3D<T, DESCRIPTOR>&     sLattice,
  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3),
    _indicatorF(std::move(indicatorF)),
    _facesF(*_indicatorF, converter),
    _pBoundForceF(sLattice, *_indicatorF, converter),
    _sumF(_pBoundForceF, *_indicatorF),
    _factor(2./( converter.getPhysDensity()*converter.getCharPhysVelocity()*converter.getCharPhysVelocity() ))
{
  this->getName() = "physDrag";

  for (int iC = 0; iC < this->getSuperStructure().getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticePhysDrag3D<T,DESCRIPTOR>(
        sLattice.getBlockLattice(iC),
        indicatorF->getBlockIndicatorF(iC),
        converter)
    );
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysDrag3D<T, DESCRIPTOR>::SuperLatticePhysDrag3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice,
  SuperGeometry3D<T>& superGeometry, const int material,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysDrag3D(sLattice,
                           superGeometry.getMaterialIndicator(material),
                           converter)
{ }

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysDrag3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T faces[7] = { };
  T sum[4]   = { };
  _sumF(sum, input);
  _facesF(faces, input);

  output[0] = _factor * sum[0] / faces[0];
  output[1] = _factor * sum[1] / faces[1];
  output[2] = _factor * sum[2] / faces[2];

  return true;
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysCorrDrag3D<T, DESCRIPTOR>::SuperLatticePhysCorrDrag3D(
  SuperLattice3D<T, DESCRIPTOR>&     sLattice,
  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3),
    _indicatorF(std::move(indicatorF)),
    _facesF(*_indicatorF, converter),
    _pBoundForceF(sLattice, *_indicatorF, converter),
    _sumF(_pBoundForceF, *_indicatorF),
    _factor(2./( converter.getPhysDensity()*converter.getCharPhysVelocity()*converter.getCharPhysVelocity() ))
{
  this->getName() = "physCorrDrag";

  for (int iC = 0; iC < this->getSuperStructure().getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticePhysCorrDrag3D<T,DESCRIPTOR>(
        sLattice.getBlockLattice(iC),
        indicatorF->getBlockIndicatorF(iC),
        converter)
    );
  }
}

template<typename T, typename DESCRIPTOR>
SuperLatticePhysCorrDrag3D<T, DESCRIPTOR>::SuperLatticePhysCorrDrag3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice,
  SuperGeometry3D<T>& superGeometry, const int material,
  const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysCorrDrag3D(sLattice,
                               superGeometry.getMaterialIndicator(material),
                               converter)
{ }

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysCorrDrag3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T faces[7] = { };
  T sum[4]   = { };
  _sumF(sum, input);
  _facesF(faces, input);

  output[0] = _factor * sum[0] / faces[0];
  output[1] = _factor * sum[1] / faces[1];
  output[2] = _factor * sum[2] / faces[2];

  return true;
}


}  // end namespace olb

#endif
