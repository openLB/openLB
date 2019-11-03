/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Mathias J. Krause,
 *                          Albert Mink, Adrian Kummerlaender
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

#ifndef BLOCK_CALC_F_3D_H
#define BLOCK_CALC_F_3D_H

#include "utilities/arithmetic.h"
#include "blockBaseF3D.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


/// Block level arithmetic operations for BlockF3D functors
/**
 * \tparam F Function object defining the arithmetic operation to be perfomed
 *         e.g. std::minus for substraction
 *
 * Functors BlockCalcF3D<T,F>::_f and BlockCalcF3D<T,F>::_f are stored as GenericF<T,int>
 * references to support arithmetic operations between block and super functors. For
 * details see BlockCalcF3D<T,F>::operator().
 **/
template <typename T, template<typename> class F>
class BlockCalcF3D : public BlockF3D<T> {
protected:
  GenericF<T,int>& _f;
  GenericF<T,int>& _g;

  /// Optional global cuboid ID for mixed super / block usage
  const int  _glob;
  const bool _fIsBlock;
  const bool _gIsBlock;
public:
  BlockCalcF3D(BlockF3D<T>& f, BlockF3D<T>& g);

  /**
   * \param f    Block functor
   * \param g    Generic functor to be restricted to block level
   * \param glob Global cuboid ID to be used to prefix calls to g
   **/
  BlockCalcF3D(BlockF3D<T>& f, GenericF<T,int>& g, int glob);

  /**
   * \param f    Generic functor to be restricted to block level
   * \param glob Global cuboid ID to be used to prefix calls to f
   * \param g    Block functor
   **/
  BlockCalcF3D(GenericF<T,int>& f, int glob, BlockF3D<T>& g);

  bool operator() (T output[], const int input[]) override;
};

/// Block level addition functor (T==bool: Union)
template <typename T>
using BlockCalcPlus3D = BlockCalcF3D<T,util::plus>;

/// Block level subtraction functor (T==bool: Without)
template <typename T>
using BlockCalcMinus3D = BlockCalcF3D<T,util::minus>;

/// Block level multiplication functor (T==bool: Intersection)
template <typename T>
using BlockCalcMultiplication3D = BlockCalcF3D<T,util::multiplies>;

/// Block level division functor
template <typename T>
using BlockCalcDivision3D = BlockCalcF3D<T,util::divides>;


} // end namespace olb

#endif
