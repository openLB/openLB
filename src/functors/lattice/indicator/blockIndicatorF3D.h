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

#ifndef BLOCK_INDICATOR_F_3D_H
#define BLOCK_INDICATOR_F_3D_H

#include "blockIndicatorBaseF3D.h"
#include "geometry/blockGeometryView3D.h"
#include "functors/analytical/indicator/smoothIndicatorBaseF3D.h"

namespace olb {

/// BlockIndicatorF3D from IndicatorF3D
template <typename T>
class BlockIndicatorFfromIndicatorF3D : public BlockIndicatorF3D<T> {
protected:
  IndicatorF3D<T>& _indicatorF;
public:
  /**
   * \param indicatorF    Indicator to be reduced to lattice space
   * \param blockGeometry Block geometry structure to be used for conversion
   *                      between lattice and physical coordinates.
   **/
  BlockIndicatorFfromIndicatorF3D(IndicatorF3D<T>&             indicatorF,
                                  BlockGeometryStructure3D<T>& blockGeometry);

  using BlockIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;

  /// Returns min lattice position of the indicated domain's bounding box
  Vector<int,3> getMin() override;
  /// Returns max lattice position of the indicated domain's bounding box
  Vector<int,3> getMax() override;
};


/// BlockIndicatorF3D from SmoothIndicatorF3D
template <typename T, bool HLBM>
class BlockIndicatorFfromSmoothIndicatorF3D : public BlockIndicatorF3D<T> {
protected:
  SmoothIndicatorF3D<T,T,HLBM>& _indicatorF;
public:
  /**
   * \param indicatorF    Smooth indicator to be reduced to lattice space
   * \param blockGeometry Block geometry structure to be used for conversion
   *                      between lattice and physical coordinates.
   **/
  BlockIndicatorFfromSmoothIndicatorF3D(SmoothIndicatorF3D<T,T,HLBM>&   indicatorF,
                                          BlockGeometryStructure3D<T>& blockGeometry);

  using BlockIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;

  /// Returns min lattice position of the indicated domain's bounding box
  Vector<int,3> getMin() override;
  /// Returns max lattice position of the indicated domain's bounding box
  Vector<int,3> getMax() override;
};


/// Block indicator functor from material numbers
template <typename T>
class BlockIndicatorMaterial3D : public BlockIndicatorF3D<T> {
protected:
  const std::vector<int> _materials;
public:
  /**
   * \param blockGeometry Block geometry structure to be queried
   * \param materials     Material number vector
   **/
  BlockIndicatorMaterial3D(BlockGeometryStructure3D<T>& blockGeometry,
                           std::vector<int>             materials);
  /**
   * \param blockGeometry Block geometry structure to be queried
   * \param materials     Material number list
   **/
  BlockIndicatorMaterial3D(BlockGeometryStructure3D<T>& blockGeometry,
                           std::list<int>             materials);
  /**
   * \param blockGeometry Block geometry structure to be queried
   * \param material      Material number
   **/
  BlockIndicatorMaterial3D(BlockGeometryStructure3D<T>& blockGeometry,
                           int                          material);

  using BlockIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;

  /// Returns true iff indicated domain subset is empty
  bool isEmpty() override;
  /// Returns min lattice position of the indicated domain's bounding box
  Vector<int,3> getMin() override;
  /// Returns max lattice position of the indicated domain's bounding box
  Vector<int,3> getMax() override;
};


/// Block indicator identity
template <typename T>
class BlockIndicatorIdentity3D : public BlockIndicatorF3D<T> {
protected:
  BlockIndicatorF3D<T>& _indicatorF;
public:
  /**
   * \param indicatorF Block indicator to be proxied
   **/
  BlockIndicatorIdentity3D(BlockIndicatorF3D<T>& indicatorF);

  using BlockIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;

  /// Returns min lattice position of the indicated domain's bounding box
  Vector<int,3> getMin() override;
  /// Returns max lattice position of the indicated domain's bounding box
  Vector<int,3> getMax() override;
};

} // namespace olb

#endif
