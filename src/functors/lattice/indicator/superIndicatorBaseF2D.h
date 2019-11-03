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

#ifndef SUPER_INDICATOR_BASE_F_2D_H
#define SUPER_INDICATOR_BASE_F_2D_H

#include "functors/genericF.h"
#include "functors/lattice/superBaseF2D.h"
#include "communication/superStructure2D.h"
#include "core/superData2D.h"
#include "blockIndicatorBaseF2D.h"

namespace olb {

template<typename T, typename W> class SuperF2D;
template<typename T> class SuperGeometry2D;
template<typename T> class SuperIndicatorIdentity2D;

template <typename T>
class SuperIndicatorF2D : public SuperF2D<T,bool> {
protected:
  SuperGeometry2D<T>& _superGeometry;

  std::unique_ptr<SuperData2D<T,bool>> _cachedData;
  std::vector<std::unique_ptr<BlockIndicatorF2D<T>>> _extendedBlockF;
public:
  using SuperF2D<T,bool>::operator();
  using identity_functor_type = SuperIndicatorIdentity2D<T>;

  SuperIndicatorF2D(SuperGeometry2D<T>& geometry);
  /**
   * Get block indicator
   *
   * \returns _blockF[iCloc] cast as BlockIndicatorF2D<T>&
   **/
  BlockIndicatorF2D<T>& getBlockIndicatorF(int iCloc);
  /**
   * Get extended block indicator
   *
   * \returns _extendedBlockF[iCloc] cast as BlockIndicatorF2D<T>&
   **/
  BlockIndicatorF2D<T>& getExtendedBlockIndicatorF(int iCloc);
  /**
   * Get underlying super geometry
   *
   * \returns _superGeometry
   **/
  SuperGeometry2D<T>& getSuperGeometry();
  /**
   * Indicator specific function operator overload
   *
   * The boolean return value of `operator()(T output[], S input[])` describes
   * the call's success and by convention must not describe the indicated domain.
   *
   * \return Domain indicator i.e. `true` iff the input lies within the described domain.
   **/
  bool operator() (const int input[]);
  /**
   * Indicator specific function operator overload
   *
   * The boolean return value of `operator()(T output[], S input[])` describes
   * the call's success and by convention must not describe the indicated domain.
   *
   * \return Domain indicator i.e. `true` iff the input lies within the described domain.
   **/
  bool operator() (int iC, int iX, int iY);

  /// Optional: initialize _cachedData for faster access
  void cache();
};


} // namespace olb

#endif
