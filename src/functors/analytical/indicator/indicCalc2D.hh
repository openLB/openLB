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

#ifndef INDIC_CALC_2D_HH
#define INDIC_CALC_2D_HH

#include "indicCalc2D.h"

namespace olb {


//////////////////////////////// IndicCalc1D ////////////////////////////////
template <typename S>
IndicCalc1D<S>::IndicCalc1D(IndicatorF1D<S>& f, IndicatorF1D<S>& g)
  : _f(f), _g(g)
{
  this->_myMin[0] = std::min(f.getMin()[0], g.getMin()[0]);
  this->_myMax[0] = std::max(f.getMax()[0], g.getMax()[0]);
  std::swap(f._ptrCalcC, this->_ptrCalcC);
}


template <typename S>
IndicPlus1D<S>::IndicPlus1D(IndicatorF1D<S>& f, IndicatorF1D<S>& g)
  : IndicCalc1D<S>(f, g)
{}

// returns 1 if( f==1 || g==1 ) UNION
template <typename S>
bool IndicPlus1D<S>::operator() (bool output[], const S input[])
{
  this->_f(output, input);
  bool tmp;
  this->_g(&tmp, input);
  output[0] |= tmp;
  return true;
}


template <typename S>
IndicMinus1D<S>::IndicMinus1D(IndicatorF1D<S>& f, IndicatorF1D<S>& g)
  : IndicCalc1D<S>(f, g)
{}

// returns 1 if( f==1 && g==0 ) WITHOUT
template <typename S>
bool IndicMinus1D<S>::operator()(bool output[], const S input[])
{
  this->_f(output, input);
  bool tmp;
  this->_g(&tmp, input);
  output[0] &= !tmp;
  return true;
}



template <typename S>
IndicMultiplication1D<S>::IndicMultiplication1D(IndicatorF1D<S>& f, IndicatorF1D<S>& g)
  : IndicCalc1D<S>(f, g)
{}

// returns 1 if( f==1 && g==1 ) INTERSECTION
template <typename S>
bool IndicMultiplication1D<S>::operator() (bool output[], const S input[])
{
  this->_f(output, input);
  bool tmp;
  this->_g(&tmp, input);
  output[0] &= tmp;
  return true;
}



template <typename S>
IndicatorF1D<S>& IndicatorF1D<S>::operator+(IndicatorF1D<S>& rhs)
{
  auto tmp = std::make_shared< IndicPlus1D<S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename S>
IndicatorF1D<S>& IndicatorF1D<S>::operator-(IndicatorF1D<S>& rhs)
{
  auto tmp = std::make_shared< IndicMinus1D<S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename S>
IndicatorF1D<S>& IndicatorF1D<S>::operator*(IndicatorF1D<S>& rhs)
{
  auto tmp = std::make_shared< IndicMultiplication1D<S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}




//////////////////////////////// IndicCalc2D ////////////////////////////////
template <typename S, template<typename U> class F>
IndicCalc2D<S,F>::IndicCalc2D(std::shared_ptr<IndicatorF2D<S>> f, std::shared_ptr<IndicatorF2D<S>> g)
  : _f(f), _g(g)
{
  for ( int i=0; i<2; i++) {
    this->_myMin[i] = std::min(_f->getMin()[i], _g->getMin()[i]);
    this->_myMax[i] = std::max(_f->getMax()[i], _g->getMax()[i]);
  }
}

template <typename S, template<typename U> class F>
bool IndicCalc2D<S,F>::operator()( bool output[], const S input[2])
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
std::shared_ptr<IndicatorF2D<S>> operator+(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs)
{
  return std::make_shared<IndicPlus2D<S>>(lhs, rhs);
}

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF2D<S>> operator-(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs)
{
  return std::make_shared<IndicMinus2D<S>>(lhs, rhs);
}

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF2D<S>> operator*(std::shared_ptr<F1<S>> lhs, std::shared_ptr<F2<S>> rhs)
{
  return std::make_shared<IndicMultiplication2D<S>>(lhs, rhs);
}

// template specialization for indicatorIdentity
template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF2D<S>> operator+(F1<S> & lhs, std::shared_ptr<F2<S>> rhs)
{
  return lhs._f + rhs;
}

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF2D<S>> operator-(F1<S> & lhs, std::shared_ptr<F2<S>> rhs)
{
  return lhs._f - rhs;
}

template<typename S, template <typename U> class F1, template <typename V> class F2,
         typename>
std::shared_ptr<IndicatorF2D<S>> operator*(F1<S> & lhs, std::shared_ptr<F2<S>> rhs)
{
  return lhs._f * rhs;
}


} // end namespace olb

#endif
