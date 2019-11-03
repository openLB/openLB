/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2018 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Adrian Kummerlaender
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

#ifndef ANALYTICAL_CALC_F_H
#define ANALYTICAL_CALC_F_H


#include "analyticalBaseF.h"
#include "utilities/functorPtr.h"
#include "utilities/arithmetic.h"


namespace olb {

/*
    arithmetic helper classes for AnalyticalF1D, AnalyticalF3D, AnalyticalF3D

    pointwise: difference, plus, multiplication, division

*/

//////////////////////////////// AnalyticCalcF1D ////////////////////////////////

/// arithmetic helper class for analytical 1D functors
template <typename T, typename S, template<typename> class F>
class AnalyticCalcF1D : public AnalyticalF1D<T,S> {
protected:
  FunctorPtr<AnalyticalF1D<T,S>> _f;
  FunctorPtr<AnalyticalF1D<T,S>> _g;
public:
  AnalyticCalcF1D(FunctorPtr<AnalyticalF1D<T,S>>&& f,
                  FunctorPtr<AnalyticalF1D<T,S>>&& g);

  AnalyticCalcF1D(T scalar, FunctorPtr<AnalyticalF1D<T,S>>&& g);
  AnalyticCalcF1D(FunctorPtr<AnalyticalF1D<T,S>>&& f, T scalar);

  bool operator() (T output[], const S input[]) override;
};

/// addition functor
template <typename T, typename S>
using AnalyticCalcPlus1D = AnalyticCalcF1D<T,S,util::plus>;

/// subtraction functor
template <typename T, typename S>
using AnalyticCalcMinus1D = AnalyticCalcF1D<T,S,util::minus>;

/// multiplication functor
template <typename T, typename S>
using AnalyticCalcMultiplication1D = AnalyticCalcF1D<T,S,util::multiplies>;

/// division functor
template <typename T, typename S>
using AnalyticCalcDivision1D = AnalyticCalcF1D<T,S,util::divides>;

//////////////////////////////// AnalyticCalcF2D ////////////////////////////////

/// arithmetic helper class for analytical 2D functors
template <typename T, typename S, template<typename> class F>
class AnalyticCalcF2D : public AnalyticalF2D<T,S> {
protected:
  FunctorPtr<AnalyticalF2D<T,S>> _f;
  FunctorPtr<AnalyticalF2D<T,S>> _g;
public:
  AnalyticCalcF2D(FunctorPtr<AnalyticalF2D<T,S>>&& f,
                  FunctorPtr<AnalyticalF2D<T,S>>&& g);

  AnalyticCalcF2D(T scalar, FunctorPtr<AnalyticalF2D<T,S>>&& g);
  AnalyticCalcF2D(FunctorPtr<AnalyticalF2D<T,S>>&& f, T scalar);

  bool operator() (T output[], const S input[]) override;
};

/// addition functor
template <typename T, typename S>
using AnalyticCalcPlus2D = AnalyticCalcF2D<T,S,util::plus>;

/// subtraction functor
template <typename T, typename S>
using AnalyticCalcMinus2D = AnalyticCalcF2D<T,S,util::minus>;

/// multiplication functor
template <typename T, typename S>
using AnalyticCalcMultiplication2D = AnalyticCalcF2D<T,S,util::multiplies>;

/// division functor
template <typename T, typename S>
using AnalyticCalcDivision2D = AnalyticCalcF2D<T,S,util::divides>;

//////////////////////////////// AnalyticCalcF3D ////////////////////////////////

/// arithmetic helper class for analytical 3D functors
template <typename T, typename S, template<typename> class F>
class AnalyticCalcF3D : public AnalyticalF3D<T,S> {
protected:
  FunctorPtr<AnalyticalF3D<T,S>> _f;
  FunctorPtr<AnalyticalF3D<T,S>> _g;
public:
  AnalyticCalcF3D(FunctorPtr<AnalyticalF3D<T,S>>&& f,
                  FunctorPtr<AnalyticalF3D<T,S>>&& g);

  AnalyticCalcF3D(T scalar, FunctorPtr<AnalyticalF3D<T,S>>&& g);
  AnalyticCalcF3D(FunctorPtr<AnalyticalF3D<T,S>>&& f, T scalar);

  bool operator() (T output[], const S input[]) override;
};

/// addition functor
template <typename T, typename S>
using AnalyticCalcPlus3D = AnalyticCalcF3D<T,S,util::plus>;

/// subtraction functor
template <typename T, typename S>
using AnalyticCalcMinus3D = AnalyticCalcF3D<T,S,util::minus>;

/// multiplication functor
template <typename T, typename S>
using AnalyticCalcMultiplication3D = AnalyticCalcF3D<T,S,util::multiplies>;

/// division functor
template <typename T, typename S>
using AnalyticCalcDivision3D = AnalyticCalcF3D<T,S,util::divides>;


/**
 * \name Arithmetic for functors managed by std::shared_ptr
 * \{
 **/

template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator+(std::shared_ptr<AnalyticalF1D<T,S>> lhs, std::shared_ptr<AnalyticalF1D<T,S>> rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator+(std::shared_ptr<AnalyticalF1D<T,S>> lhs, T rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator+(T lhs, std::shared_ptr<AnalyticalF1D<T,S>> rhs);

template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator-(std::shared_ptr<AnalyticalF1D<T,S>> lhs, std::shared_ptr<AnalyticalF1D<T,S>> rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator-(std::shared_ptr<AnalyticalF1D<T,S>> lhs, T rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator-(T lhs, std::shared_ptr<AnalyticalF1D<T,S>> rhs);

template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator*(std::shared_ptr<AnalyticalF1D<T,S>> lhs, std::shared_ptr<AnalyticalF1D<T,S>> rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator*(std::shared_ptr<AnalyticalF1D<T,S>> lhs, T rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator*(T lhs, std::shared_ptr<AnalyticalF1D<T,S>> rhs);

template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator/(std::shared_ptr<AnalyticalF1D<T,S>> lhs, std::shared_ptr<AnalyticalF1D<T,S>> rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator/(std::shared_ptr<AnalyticalF1D<T,S>> lhs, T rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator/(T lhs, std::shared_ptr<AnalyticalF1D<T,S>> rhs);


template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator+(std::shared_ptr<AnalyticalF2D<T,S>> lhs, std::shared_ptr<AnalyticalF2D<T,S>> rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator+(std::shared_ptr<AnalyticalF2D<T,S>> lhs, T rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator+(T lhs, std::shared_ptr<AnalyticalF2D<T,S>> rhs);

template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator-(std::shared_ptr<AnalyticalF2D<T,S>> lhs, std::shared_ptr<AnalyticalF2D<T,S>> rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator-(std::shared_ptr<AnalyticalF2D<T,S>> lhs, T rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator-(T lhs, std::shared_ptr<AnalyticalF2D<T,S>> rhs);

template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator*(std::shared_ptr<AnalyticalF2D<T,S>> lhs, std::shared_ptr<AnalyticalF2D<T,S>> rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator*(std::shared_ptr<AnalyticalF2D<T,S>> lhs, T rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator*(T lhs, std::shared_ptr<AnalyticalF2D<T,S>> rhs);

template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator/(std::shared_ptr<AnalyticalF2D<T,S>> lhs, std::shared_ptr<AnalyticalF2D<T,S>> rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator/(std::shared_ptr<AnalyticalF2D<T,S>> lhs, T rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator/(T lhs, std::shared_ptr<AnalyticalF2D<T,S>> rhs);


template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator+(std::shared_ptr<AnalyticalF3D<T,S>> lhs, std::shared_ptr<AnalyticalF3D<T,S>> rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator+(std::shared_ptr<AnalyticalF3D<T,S>> lhs, T rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator+(T lhs, std::shared_ptr<AnalyticalF3D<T,S>> rhs);

template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator-(std::shared_ptr<AnalyticalF3D<T,S>> lhs, std::shared_ptr<AnalyticalF3D<T,S>> rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator-(std::shared_ptr<AnalyticalF3D<T,S>> lhs, T rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator-(T lhs, std::shared_ptr<AnalyticalF3D<T,S>> rhs);

template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator*(std::shared_ptr<AnalyticalF3D<T,S>> lhs, std::shared_ptr<AnalyticalF3D<T,S>> rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator*(std::shared_ptr<AnalyticalF3D<T,S>> lhs, T rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator*(T lhs, std::shared_ptr<AnalyticalF3D<T,S>> rhs);

template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator/(std::shared_ptr<AnalyticalF3D<T,S>> lhs, std::shared_ptr<AnalyticalF3D<T,S>> rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator/(std::shared_ptr<AnalyticalF3D<T,S>> lhs, T rhs);
template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator/(T lhs, std::shared_ptr<AnalyticalF3D<T,S>> rhs);

///\}


} // end namespace olb

#endif
