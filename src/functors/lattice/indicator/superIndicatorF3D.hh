/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016-2018 Benjamin Förster, Adrian Kummerlaender
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

#ifndef SUPER_INDICATOR_F_3D_HH
#define SUPER_INDICATOR_F_3D_HH

#include <numeric>

#include "superIndicatorF3D.h"
#include "blockIndicatorF3D.h"
#include "core/util.h"

namespace olb {

template <typename T>
SuperIndicatorFfromIndicatorF3D<T>::SuperIndicatorFfromIndicatorF3D(
  FunctorPtr<IndicatorF3D<T>>&& indicatorF, SuperGeometry3D<T>& geometry)
  : SuperIndicatorF3D<T>(geometry),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "SuperIndicator_from_" + _indicatorF->getName();

  LoadBalancer<T>& load = this->getSuperStructure().getLoadBalancer();

  for (int iC = 0; iC < load.size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorFfromIndicatorF3D<T>(
        *_indicatorF, geometry.getBlockGeometry(iC))
    );
    this->_extendedBlockF.emplace_back(
      new BlockIndicatorFfromIndicatorF3D<T>(
        *_indicatorF, geometry.getExtendedBlockGeometry(iC))
    );
  }
}

template <typename T>
bool SuperIndicatorFfromIndicatorF3D<T>::operator() (bool output[], const int input[])
{
  T physR[3];
  this->_superStructure.getCuboidGeometry().getPhysR(physR, input);
  return _indicatorF(output, physR);
}


template <typename T, bool HLBM>
SuperIndicatorFfromSmoothIndicatorF3D<T, HLBM>::SuperIndicatorFfromSmoothIndicatorF3D(
  FunctorPtr<SmoothIndicatorF3D<T,T,HLBM>>&& indicatorF,
  SuperGeometry3D<T>&                     geometry)
  : SuperIndicatorF3D<T>(geometry),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "SuperIndicator_from_" + _indicatorF->getName();

  LoadBalancer<T>& load = this->getSuperStructure().getLoadBalancer();

  for (int iC = 0; iC < load.size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorFfromSmoothIndicatorF3D<T, HLBM>(
        *_indicatorF, geometry.getBlockGeometry(iC))
    );
    this->_extendedBlockF.emplace_back(
      new BlockIndicatorFfromSmoothIndicatorF3D<T, HLBM>(
        *_indicatorF, geometry.getExtendedBlockGeometry(iC))
    );
  }
}

template <typename T, bool HLBM>
bool SuperIndicatorFfromSmoothIndicatorF3D<T, HLBM>::operator() (bool output[], const int input[])
{
  T physR[3];
  T inside[1];
  this->_superStructure.getCuboidGeometry().getPhysR(physR, input);
  _indicatorF(inside, physR);
  return !util::nearZero(inside[0]);
}

template <typename T>
SuperIndicatorMaterial3D<T>::SuperIndicatorMaterial3D(
  SuperGeometry3D<T>& geometry, std::vector<int> materials)
  : SuperIndicatorF3D<T>(geometry)
{
  const std::string matString = std::accumulate(
                                  materials.begin()+1,
                                  materials.end(),
                                  std::to_string(materials[0]),
  [](const std::string& a, int b) {
    return a + '_' + std::to_string(b);
  });
  this->getName() = "SuperIndicator_on_Material_" + matString;

  for (int iC = 0; iC < this->_superGeometry.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorMaterial3D<T>(this->_superGeometry.getBlockGeometry(iC),
                                      materials)
    );
    this->_extendedBlockF.emplace_back(
      new BlockIndicatorMaterial3D<T>(this->_superGeometry.getExtendedBlockGeometry(iC),
                                      materials)
    );
  }
}

template <typename T>
SuperIndicatorMaterial3D<T>::SuperIndicatorMaterial3D(
  SuperGeometry3D<T>& geometry, std::list<int> materials)
  : SuperIndicatorMaterial3D(geometry,
                             std::vector<int>(materials.begin(), materials.end()))
{ }

template <typename T>
bool SuperIndicatorMaterial3D<T>::operator() (bool output[], const int input[])
{
  output[0] = false;

  LoadBalancer<T>& load = this->_superGeometry.getLoadBalancer();

  if (!this->_blockF.empty() && load.isLocal(input[0])) {
    // query material number of appropriate block indicator
    return this->getBlockF(load.loc(input[0]))(output,&input[1]);
  }
  else {
    return false;
  }
}

template <typename T>
SuperIndicatorIdentity3D<T>::SuperIndicatorIdentity3D(FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperIndicatorF3D<T>(indicatorF->getSuperGeometry()),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = _indicatorF->getName();

  for (int iC = 0; iC < _indicatorF->getBlockFSize(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorIdentity3D<T>(_indicatorF->getBlockIndicatorF(iC)));
    this->_extendedBlockF.emplace_back(
      new BlockIndicatorIdentity3D<T>(_indicatorF->getExtendedBlockIndicatorF(iC)));
  }
}

template <typename T>
bool SuperIndicatorIdentity3D<T>::operator()(bool output[], const int input[])
{
  return _indicatorF(output, input);
}

} // namespace olb

#endif
