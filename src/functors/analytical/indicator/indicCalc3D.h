/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Albert Mink
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

#ifndef INDIC_CALC_3D_H
#define INDIC_CALC_3D_H

#include "utilities/arithmetic.h"
#include "indicatorBaseF3D.h"

namespace olb {


/*
 *  arithmetic helper classes for IndicatorF3D, smoothIndicator3D
 *  UNION         +
 *  WITHOUT       -
 *  INTERSECTION  *
*/

//////////////////////////////// indicCalc3D ////////////////////////////////
/// arithmetic helper class for Indicator 3d functors
template <typename S, template<typename U> class F>
class IndicCalc3D : public IndicatorF3D<S> {
protected:
  std::shared_ptr<IndicatorF3D<S>> _f;
  std::shared_ptr<IndicatorF3D<S>> _g;
public:
  IndicCalc3D( std::shared_ptr<IndicatorF3D<S>> f, std::shared_ptr<IndicatorF3D<S>> g );

  bool operator() (bool output[], const S input[3]) override;
};

/// Addition functor (W==bool: Union)
template <typename S>
using IndicPlus3D = IndicCalc3D<S,util::plus>;
template <typename S>
using IndicMinus3D = IndicCalc3D<S,util::minus>;
template <typename S>
using IndicMultiplication3D = IndicCalc3D<S,util::multiplies>;

/** Free function implements lhs+rhs, only for IndicaotrsF3D types through enable_if and is_base_of
 *
 * \tparam S usual type for source dimension of the functor
 * \tparam F1 lhs has to be derived from IndicatorF3D, otherwise function is disabled
 * \tparam F2 rhs
 */
template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorF3D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF3D<S>> operator+(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs);

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorF3D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF3D<S>> operator-(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs);

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorF3D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF3D<S>> operator*(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs);

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorIdentity3D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF3D<S>> operator+(F1<S> & lhs, std::shared_ptr<F2<S>> rhs);

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorIdentity3D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF3D<S>> operator-(F1<S> & lhs, std::shared_ptr<F2<S>> rhs);

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorIdentity3D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF3D<S>> operator*(F1<S> & lhs, std::shared_ptr<F2<S>> rhs);
} // end namespace olb

#endif
