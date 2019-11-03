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

#ifndef INDIC_CALC_3D_HH
#define INDIC_CALC_3D_HH

#include "indicCalc3D.h"

namespace olb {


template <typename S, template<typename U> class F>
IndicCalc3D<S,F>::IndicCalc3D(std::shared_ptr<IndicatorF3D<S>> f, std::shared_ptr<IndicatorF3D<S>> g)
  : _f(f), _g(g)
{
  for ( int i=0; i<3; i++) {
    this->_myMin[i] = std::min(_f->getMin()[i], _g->getMin()[i]);
    this->_myMax[i] = std::max(_f->getMax()[i], _g->getMax()[i]);
  }
}



template <typename S, template<typename> class F>
bool IndicCalc3D<S,F>::operator()( bool output[], const S input[3])
{
  // componentwise operation on equidimensional functors
  bool* outputF = output;
  _f->operator()(outputF, input);

  bool outputG[this->getTargetDim()];
  _g->operator()(outputG, input);

  for (int i = 0; i < this->getTargetDim(); i++) {
    output[i] = F<S>()(outputF[i], outputG[i]);
  }
  return output;
}




//// no association to a operator+ from a class is needed, thus we have these free functions
template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF3D<S>> operator+(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs)
{
  return std::make_shared<IndicPlus3D<S>>(lhs, rhs);
}

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF3D<S>> operator-(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs)
{
  return std::make_shared<IndicMinus3D<S>>(lhs, rhs);
}

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF3D<S>> operator*(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs)
{
  return std::make_shared<IndicMultiplication3D<S>>(lhs, rhs);
}

// template specialization for indicatorIdentity
template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF3D<S>> operator+(F1<S> & lhs, std::shared_ptr<F2<S>> rhs)
{
  return lhs._f + rhs;
}

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF3D<S>> operator-(F1<S> & lhs, std::shared_ptr<F2<S>> rhs)
{
  return lhs._f - rhs;
}

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF3D<S>> operator*(F1<S> & lhs, std::shared_ptr<F2<S>> rhs)
{
  return lhs._f * rhs;
}


} // end namespace olb

#endif
