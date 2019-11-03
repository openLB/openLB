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

#ifndef BLOCK_INDICATOR_F_2D_H
#define BLOCK_INDICATOR_F_2D_H

#include "blockIndicatorBaseF2D.h"
#include "geometry/blockGeometryView2D.h"
#include "functors/analytical/indicator/smoothIndicatorBaseF2D.h"

namespace olb {

/// BlockIndicatorF2D from IndicatorF2D
template <typename T>
class BlockIndicatorFfromIndicatorF2D : public BlockIndicatorF2D<T> {
protected:
  IndicatorF2D<T>& _indicatorF;
public:
  /**
   * \param indicatorF    Indicator to be reduced to lattice space
   * \param blockGeometry Block geometry structure to be used for conversion
   *                      between lattice and physical coordinates.
   **/
  BlockIndicatorFfromIndicatorF2D(IndicatorF2D<T>&             indicatorF,
                                  BlockGeometryStructure2D<T>& blockGeometry);

  using BlockIndicatorF2D<T>::operator();
  bool operator() (bool output[], const int input[]) override;

  /// Returns min lattice position of the indicated domain's bounding box
  Vector<int,2> getMin() override;
  /// Returns max lattice position of the indicated domain's bounding box
  Vector<int,2> getMax() override;
};


/// BlockIndicatorF2D from SmoothIndicatorF2D
/**
 * Note on get(Min,Max): SmoothIndicatorF2D currently doesn't expose a bounding box
 *                       which is why these methods return the full block domain.
 **/
template <typename T, bool HLBM>
class BlockIndicatorFfromSmoothIndicatorF2D : public BlockIndicatorF2D<T> {
protected:
  SmoothIndicatorF2D<T,T,HLBM>& _indicatorF;
public:
  /**
   * \param indicatorF    Smooth indicator to be reduced to lattice space
   * \param blockGeometry Block geometry structure to be used for conversion
   *                      between lattice and physical coordinates.
   **/
  BlockIndicatorFfromSmoothIndicatorF2D(SmoothIndicatorF2D<T,T,HLBM>&   indicatorF,
                                          BlockGeometryStructure2D<T>& blockGeometry);

  using BlockIndicatorF2D<T>::operator();
  bool operator() (bool output[], const int input[]) override;

  /// Returns a min lattice position of the indicated domain's bounding box
  Vector<int,2> getMin() override;
  /// Returns a max lattice position of the indicated domain's bounding box
  Vector<int,2> getMax() override;
};


/// Block indicator functor from material numbers
template <typename T>
class BlockIndicatorMaterial2D : public BlockIndicatorF2D<T> {
protected:
  const std::vector<int> _materials;
public:
  /**
   * \param blockGeometry Block geometry structue to be queried
   * \param materials     Material number vector
   **/
  BlockIndicatorMaterial2D(BlockGeometryStructure2D<T>& blockGeometry,
                           std::vector<int>             materials);
  /**
   * \param blockGeometry Block geometry structure to be queried
   * \param materials     Material number list
   **/
  BlockIndicatorMaterial2D(BlockGeometryStructure2D<T>& blockGeometry,
                           std::list<int>             materials);
  /**
   * \param blockGeometry Block geometry structure to be queried
   * \param material      Material number
   **/
  BlockIndicatorMaterial2D(BlockGeometryStructure2D<T>& blockGeometry,
                           int                          material);

  using BlockIndicatorF2D<T>::operator();
  bool operator() (bool output[], const int input[]) override;

  /// Returns true iff indicated domain subset is empty
  bool isEmpty() override;
  /// Returns min lattice position of the indicated domain's bounding box
  Vector<int,2> getMin() override;
  /// Returns max lattice position of the indicated domain's bounding box
  Vector<int,2> getMax() override;
};


/// Block indicator identity
template <typename T>
class BlockIndicatorIdentity2D : public BlockIndicatorF2D<T> {
protected:
  BlockIndicatorF2D<T>& _indicatorF;
public:
  /**
   * \param indicatorF Block indicator to be proxied
   **/
  BlockIndicatorIdentity2D(BlockIndicatorF2D<T>& indicatorF);

  using BlockIndicatorF2D<T>::operator();
  bool operator() (bool output[], const int input[]) override;

  /// Returns min lattice position of the indicated domain's bounding box
  Vector<int,2> getMin() override;
  /// Returns max lattice position of the indicated domain's bounding box
  Vector<int,2> getMax() override;
};

} // namespace olb

#endif
