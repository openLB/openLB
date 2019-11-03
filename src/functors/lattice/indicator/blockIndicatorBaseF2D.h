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

#ifndef BLOCK_INDICATOR_BASE_F_2D_H
#define BLOCK_INDICATOR_BASE_F_2D_H

#include "functors/lattice/blockBaseF2D.h"
#include "core/blockData2D.h"
#include "core/blockStructure2D.h"

namespace olb {

/// Base block indicator functor (discrete)
template <typename T>
class BlockIndicatorF2D : public BlockF2D<bool> {
protected:
  BlockGeometryStructure2D<T>& _blockGeometryStructure;
  const BlockData2D<T,bool>*   _cachedData;
public:
  using BlockF2D<bool>::operator();

  BlockIndicatorF2D(BlockGeometryStructure2D<T>& geometry);
  /// Get underlying block geometry structure
  /**
   * \returns _blockGeometryStructure
   **/
  BlockGeometryStructure2D<T>& getBlockGeometryStructure();
  /// Block indicator specific function operator overload
  /**
   * The boolean return value of `operator()(T output[], S input[])` describes
   * the call's success and by convention must not describe the indicated domain.
   *
   * \return Domain indicator i.e. `true` iff the input lies within the described subset.
   **/
  bool operator() (const int input[]);
  /// Block indicator specific function operator overload
  /**
   * The boolean return value of `operator()(T output[], S input[])` describes
   * the call's success and by convention must not describe the indicated domain.
   *
   * \return Domain indicator i.e. `true` iff the input lies within the described subset.
   **/
  bool operator() (int iX, int iY);

  /// Set bool-mask cache to be used by indicator operator overloads
  void setCache(const BlockData2D<T,bool>& cache);

  /// Returns true only if the indicated domain subset is empty
  /**
   * May return false even if the indicated domain subset is in fact empty.
   * Primarily implemented to minimize block accesses if an empty domain can
   * be inferred by e.g. BlockGeometryStatistics2D data.
   *
   * i.e. only override this method if the domain can be checked for emptyness
   *      in an efficient fashion.
   **/
  virtual bool isEmpty();
  /// Returns min lattice position of the indicated subset's bounding box
  virtual Vector<int,2> getMin() = 0;
  /// Returns max lattice position of the indicated subset's bounding box
  virtual Vector<int,2> getMax() = 0;

};

} // namespace olb

#endif
