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

#ifndef SUPER_CONST_2D_H
#define SUPER_CONST_2D_H

#include "superBaseF2D.h"
#include "core/vector.h"

namespace olb {


/// Functor returning a constant vector for all inputs
/**
 * Note that outside of SuperCalcF2D you should never be required to
 * instantiate this functor by itself. Functor arithmetic involving
 * constant values only makes sense in conjunction  with simulation
 * dependent functors.
 *
 * Instantiation of this functor is performed by SuperCalcF2D to compose
 * normal functors with constant vectors. e.g. scalar multiplication
 *
 * SuperCalcF2D handles compositions of blockified and non-blockified
 * functors which is why SuperConst2D is explicitly not blockified.
 **/
template <typename T, typename W = T>
class SuperConst2D final : public SuperF2D<T,W> {
protected:
  const std::vector<W> _c;
public:
  /// Constructor accepting std::vector
  /**
   * \param superStructure Only required to instantiate underlying SuperF2D
   * \param v              std::vector to be copied to output by operator.
   *                       Size determines target dimension.
   **/
  SuperConst2D(SuperStructure2D<T>& superStructure, std::vector<W> v);

  /// Constructor accepting single scalar
  SuperConst2D(SuperStructure2D<T>& superStructure, W scalar);
  /// Constructor template accepting vectors
  template <unsigned Size>
  SuperConst2D(SuperStructure2D<T>& superStructure, Vector<W,Size> v)
    : SuperConst2D(superStructure, v.toStdVector()) { };

  bool operator() (W output[], const int input[]) override;
};


}

#endif
