/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Mathias J. Krause, Cyril Masquelier,
 *  Benjamin Förster, Albert Mink
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

#ifndef INDIC_CALC_2D_H
#define INDIC_CALC_2D_H

#include "indicatorBaseF2D.h"
#include "utilities/arithmetic.h"

namespace olb {


/*
 *  arithmetic helper classes for IndicatorF1D, IndicatorF2D, smoothIndicator2D
 *  UNION         +
 *  WITHOUT       -
 *  INTERSECTION  *
*/

//////////////////////////////// IndicCalc1D ////////////////////////////////
/// arithmetic helper class for Indicator 1d functors
template <typename S>
class IndicCalc1D : public IndicatorF1D<S> {
protected:
  IndicatorF1D<S>& _f;
  IndicatorF1D<S>& _g;
public:
  // set image/target dimensions of IndicCalc1D as well
  IndicCalc1D(IndicatorF1D<S>& f, IndicatorF1D<S>& g);
};

/// addition functor acts as union
template <typename S>
class IndicPlus1D : public IndicCalc1D<S> {
public:
  IndicPlus1D(IndicatorF1D<S>& f, IndicatorF1D<S>& g);
  bool operator() (bool output[], const S input[]) override;
};

/// subtraction functor acts as without
template <typename S>
class IndicMinus1D : public IndicCalc1D<S> {
public:
  IndicMinus1D(IndicatorF1D<S>& f, IndicatorF1D<S>& g);
  bool operator() (bool output[], const S input[]) override;
};

/// multiplication functor acts as intersection
template <typename S>
class IndicMultiplication1D : public IndicCalc1D<S> {
public:
  IndicMultiplication1D(IndicatorF1D<S>& f, IndicatorF1D<S>& g);
  bool operator() (bool output[], const S input[]) override;
};



//////////////////////////////// indicCalc2D ////////////////////////////////
/// arithmetic helper class for Indicator 2D functors
template <typename S, template<typename U> class F>
class IndicCalc2D : public IndicatorF2D<S> {
protected:
  std::shared_ptr<IndicatorF2D<S>> _f;
  std::shared_ptr<IndicatorF2D<S>> _g;
public:
  IndicCalc2D( std::shared_ptr<IndicatorF2D<S>> f, std::shared_ptr<IndicatorF2D<S>> g );

  bool operator() (bool output[], const S input[2]) override;
};

/// Addition functor (W==bool: Union)
template <typename S>
using IndicPlus2D = IndicCalc2D<S,util::plus>;
template <typename S>
using IndicMinus2D = IndicCalc2D<S,util::minus>;
template <typename S>
using IndicMultiplication2D = IndicCalc2D<S,util::multiplies>;

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorF2D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF2D<S>> operator+(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs);

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorF2D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF2D<S>> operator-(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs);

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorF2D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF2D<S>> operator*(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs);

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorIdentity2D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF2D<S>> operator+(F1<S> & lhs, std::shared_ptr<F2<S>> rhs);

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorIdentity2D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF2D<S>> operator-(F1<S> & lhs, std::shared_ptr<F2<S>> rhs);

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename=typename std::enable_if<std::is_base_of<IndicatorIdentity2D<S>, F1<S>>::value>::type>
std::shared_ptr<IndicatorF2D<S>> operator*(F1<S> & lhs, std::shared_ptr<F2<S>> rhs);


} // end namespace olb

#endif
