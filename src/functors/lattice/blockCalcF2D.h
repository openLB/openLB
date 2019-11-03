/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013-2017 Albert Mink, Lukas Baron, Mathias J. Krause,
 *                          Adrian Kummerlaender
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

#ifndef BLOCK_CALC_F_2D_H
#define BLOCK_CALC_F_2D_H

#include "utilities/arithmetic.h"
#include "blockBaseF2D.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


/// Block level arithmetic operations for BlockF2D functors
/**
 * \tparam F Function object defining the arithmetic operation to be perfomed
 *         e.g. std::minus for substraction
 *
 * Functors BlockCalcF2D<T,F>::_f and BlockCalcF2D<T,F>::_f are stored as GenericF<T,int>
 * references to support arithmetic operations between block and super functors. For
 * details see BlockCalcF2D<T,F>::operator().
 **/
template <typename T, template<typename> class F>
class BlockCalcF2D : public BlockF2D<T> {
protected:
  GenericF<T,int>& _f;
  GenericF<T,int>& _g;

  /// Optional global cuboid ID for mixed super / block usage
  const int  _glob;
  const bool _fIsBlock;
  const bool _gIsBlock;
public:
  BlockCalcF2D(BlockF2D<T>& f, BlockF2D<T>& g);

  /**
   * \param f    Block functor
   * \param g    Generic functor to be restricted to block level
   * \param glob Global cuboid ID to be used to prefix calls to g
   **/
  BlockCalcF2D(BlockF2D<T>& f, GenericF<T,int>& g, int glob);

  /**
   * \param f    Generic functor to be restricted to block level
   * \param glob Global cuboid ID to be used to prefix calls to f
   * \param g    Block functor
   **/
  BlockCalcF2D(GenericF<T,int>& f, int glob, BlockF2D<T>& g);

  bool operator() (T output[], const int input[]) override;
};

/// Block level addition functor (T==bool: Union)
template <typename T>
using BlockCalcPlus2D = BlockCalcF2D<T,util::plus>;

/// Block level subtraction functor (T==bool: Without)
template <typename T>
using BlockCalcMinus2D = BlockCalcF2D<T,util::minus>;

/// Block level multiplication functor (T==bool: Intersection)
template <typename T>
using BlockCalcMultiplication2D = BlockCalcF2D<T,util::multiplies>;

/// Block level division functor
template <typename T>
using BlockCalcDivision2D = BlockCalcF2D<T,util::divides>;


} // end namespace olb

#endif
