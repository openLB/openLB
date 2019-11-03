/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerlaender
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

#ifndef BLOCK_INDICATOR_F_3D_HH
#define BLOCK_INDICATOR_F_3D_HH

#include <algorithm>

#include "blockIndicatorF3D.h"
#include "core/util.h"

namespace olb {

template <typename T>
BlockIndicatorFfromIndicatorF3D<T>::BlockIndicatorFfromIndicatorF3D(
  IndicatorF3D<T>& indicatorF, BlockGeometryStructure3D<T>& blockGeometry)
  : BlockIndicatorF3D<T>(blockGeometry),
    _indicatorF(indicatorF)
{ }

template <typename T>
bool BlockIndicatorFfromIndicatorF3D<T>::operator() (bool output[], const int input[])
{
  T physR[3];
  this->_blockGeometryStructure.getPhysR(physR,input);
  return _indicatorF(output,physR);
}

template <typename T>
Vector<int,3> BlockIndicatorFfromIndicatorF3D<T>::getMin()
{
  const Vector<T,3> min = _indicatorF.getMin();
  return Vector<int,3> {
    static_cast<int>(floor(min[0])),
    static_cast<int>(floor(min[1])),
    static_cast<int>(floor(min[2]))
  };
}

template <typename T>
Vector<int,3> BlockIndicatorFfromIndicatorF3D<T>::getMax()
{
  const Vector<T,3> max = _indicatorF.getMax();
  return Vector<int,3> {
    static_cast<int>(ceil(max[0])),
    static_cast<int>(ceil(max[1])),
    static_cast<int>(ceil(max[2]))
  };
}


template <typename T, bool HLBM>
BlockIndicatorFfromSmoothIndicatorF3D<T, HLBM>::BlockIndicatorFfromSmoothIndicatorF3D(
  SmoothIndicatorF3D<T,T,HLBM>& indicatorF, BlockGeometryStructure3D<T>& blockGeometry)
  : BlockIndicatorF3D<T>(blockGeometry),
    _indicatorF(indicatorF)
{ }

template <typename T, bool HLBM>
bool BlockIndicatorFfromSmoothIndicatorF3D<T,HLBM>::operator() (bool output[], const int input[])
{
  T physR[3];
  T inside[1];
  this->_blockGeometryStructure.getPhysR(physR,input);
  _indicatorF(inside, physR);
  return !util::nearZero(inside[0]);
}

template <typename T, bool HLBM>
Vector<int,3> BlockIndicatorFfromSmoothIndicatorF3D<T, HLBM>::getMin()
{
  const T min = -_indicatorF.getCircumRadius();
  return Vector<int,3> {
    static_cast<int>(floor(min)),
    static_cast<int>(floor(min)),
    static_cast<int>(floor(min))
  };
}

template <typename T, bool HLBM>
Vector<int,3> BlockIndicatorFfromSmoothIndicatorF3D<T, HLBM>::getMax()
{
  const T max = _indicatorF.getCircumRadius();
  return Vector<int,3> {
    static_cast<int>(ceil(max)),
    static_cast<int>(ceil(max)),
    static_cast<int>(ceil(max))
  };
}


template <typename T>
BlockIndicatorMaterial3D<T>::BlockIndicatorMaterial3D(
  BlockGeometryStructure3D<T>& blockGeometry, std::vector<int> materials)
  : BlockIndicatorF3D<T>(blockGeometry),
    _materials(materials)
{ }

template <typename T>
BlockIndicatorMaterial3D<T>::BlockIndicatorMaterial3D(
  BlockGeometryStructure3D<T>& blockGeometry, std::list<int> materials)
  : BlockIndicatorMaterial3D(blockGeometry,
                             std::vector<int>(materials.begin(), materials.end()))
{ }

template <typename T>
BlockIndicatorMaterial3D<T>::BlockIndicatorMaterial3D(
  BlockGeometryStructure3D<T>& blockGeometry, int material)
  : BlockIndicatorMaterial3D(blockGeometry, std::vector<int>(1,material))
{ }

template <typename T>
bool BlockIndicatorMaterial3D<T>::operator() (bool output[], const int input[])
{
  // read material number explicitly using the const version
  // of BlockGeometry3D<T>::get to avoid resetting geometry
  // statistics:
  const BlockGeometryStructure3D<T>& blockGeometry = this->_blockGeometryStructure;
  const int current = blockGeometry.getMaterial(input[0], input[1], input[2]);
  output[0] = std::any_of(_materials.cbegin(),
                          _materials.cend(),
                          [current](int material) { return current == material; });

  return true;
}

template <typename T>
bool BlockIndicatorMaterial3D<T>::isEmpty()
{
  auto& statistics = this->getBlockGeometryStructure().getStatistics();

  return std::none_of(_materials.cbegin(), _materials.cend(),
  [&statistics](int material) -> bool {
    return statistics.getNvoxel(material) > 0;
  });
}

template <typename T>
Vector<int,3> BlockIndicatorMaterial3D<T>::getMin()
{
  auto& blockGeometry = this->getBlockGeometryStructure();
  auto& statistics    = blockGeometry.getStatistics();

  Vector<int,3> globalMin{
    blockGeometry.getNx()-1,
    blockGeometry.getNy()-1,
    blockGeometry.getNz()-1,
  };

  for ( int material : _materials ) {
    if ( statistics.getNvoxel(material) > 0 ) {
      const Vector<int,3> localMin = statistics.getMinLatticeR(material);
      for ( int d = 0; d < 3; ++d ) {
        globalMin[d] = localMin[d] < globalMin[d] ? localMin[d] : globalMin[d];
      }
    }
  }

  return globalMin;
}

template <typename T>
Vector<int,3> BlockIndicatorMaterial3D<T>::getMax()
{
  auto& statistics = this->getBlockGeometryStructure().getStatistics();

  Vector<int,3> globalMax{ 0, 0, 0 };

  for ( int material : _materials ) {
    if ( statistics.getNvoxel(material) > 0 ) {
      const Vector<int,3> localMax = statistics.getMaxLatticeR(material);
      for ( int d = 0; d < 3; ++d ) {
        globalMax[d] = localMax[d] > globalMax[d] ? localMax[d] : globalMax[d];
      }
    }
  }

  return globalMax;
}


template <typename T>
BlockIndicatorIdentity3D<T>::BlockIndicatorIdentity3D(BlockIndicatorF3D<T>& indicatorF)
  : BlockIndicatorF3D<T>(indicatorF.getBlockGeometryStructure()),
    _indicatorF(indicatorF)
{ }

template <typename T>
bool BlockIndicatorIdentity3D<T>::operator() (bool output[], const int input[])
{
  return _indicatorF(output, input);
}

template <typename T>
Vector<int,3> BlockIndicatorIdentity3D<T>::getMin()
{
  return _indicatorF.getMin();
}

template <typename T>
Vector<int,3> BlockIndicatorIdentity3D<T>::getMax()
{
  return _indicatorF.getMax();
}

} // namespace olb

#endif
