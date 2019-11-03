/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2017 Albert Mink, Mathias J. Krause,
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

#ifndef SUPER_CALC_F_2D_H
#define SUPER_CALC_F_2D_H

#include "utilities/arithmetic.h"
#include "superBaseF2D.h"
#include "utilities/functorPtr.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


/// Arithmetic operations for SuperF2D functors
/**
 * \tparam F Function object defining the arithmetic operation to be perfomed
 *         e.g. std::minus for substraction
 *
 * Block level functors are instantiated for operations if at least one input
 * functor exposes block level functors. See BlockCalcF2D
 *
 * Global queries are not delegated to block level functors to prevent unnecessary
 * synchronization.
 *
 * Warning: Allocation error possible in functors that have multiple functor
 * evaluation like SuperSum2D
 **/
template <typename T, typename W, template<typename> class F>
class SuperCalcF2D : public SuperF2D<T,W> {
protected:
  FunctorPtr<SuperF2D<T,W>> _f;
  FunctorPtr<SuperF2D<T,W>> _g;
public:
  SuperCalcF2D(FunctorPtr<SuperF2D<T,W>>&& f,
               FunctorPtr<SuperF2D<T,W>>&& g);

  SuperCalcF2D(W scalar, FunctorPtr<SuperF2D<T,W>>&& g);
  SuperCalcF2D(FunctorPtr<SuperF2D<T,W>>&& f, W scalar);

  bool operator() (W output[], const int input[]) override;
};

/// Addition functor (W==bool: Union)
template <typename T, typename W>
using SuperCalcPlus2D = SuperCalcF2D<T,W,util::plus>;

/// Subtraction functor (W==bool: Without)
template <typename T, typename W>
using SuperCalcMinus2D = SuperCalcF2D<T,W,util::minus>;

/// Multiplication functor (W==bool: Intersection)
template <typename T, typename W>
using SuperCalcMultiplication2D = SuperCalcF2D<T,W,util::multiplies>;

/// Division functor
template <typename T, typename W>
using SuperCalcDivision2D = SuperCalcF2D<T,W,util::divides>;

/// Power functor
template <typename T, typename W = T>
using SuperCalcPower2D = SuperCalcF2D<T,W,util::power>;

/**
 * \name Arithmetic for functors managed by std::shared_ptr
 * \{
 **/

template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator+(std::shared_ptr<SuperF2D<T,W>> lhs, std::shared_ptr<SuperF2D<T,W>> rhs);
template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator+(std::shared_ptr<SuperF2D<T,W>> lhs, W rhs);
template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator+(W lhs, std::shared_ptr<SuperF2D<T,W>> rhs);

template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator-(std::shared_ptr<SuperF2D<T,W>> lhs, std::shared_ptr<SuperF2D<T,W>> rhs);
template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator-(std::shared_ptr<SuperF2D<T,W>> lhs, W rhs);
template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator-(W lhs, std::shared_ptr<SuperF2D<T,W>> rhs);

template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator*(std::shared_ptr<SuperF2D<T,W>> lhs, std::shared_ptr<SuperF2D<T,W>> rhs);
template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator*(std::shared_ptr<SuperF2D<T,W>> lhs, W rhs);
template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator*(W lhs, std::shared_ptr<SuperF2D<T,W>> rhs);

template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator/(std::shared_ptr<SuperF2D<T,W>> lhs, std::shared_ptr<SuperF2D<T,W>> rhs);
template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator/(std::shared_ptr<SuperF2D<T,W>> lhs, W rhs);
template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator/(W lhs, std::shared_ptr<SuperF2D<T,W>> rhs);

///\}


} // end namespace olb

#endif
