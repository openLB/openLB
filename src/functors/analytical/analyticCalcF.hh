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

#ifndef ANALYTICAL_CALC_F_HH
#define ANALYTICAL_CALC_F_HH


#include "analyticCalcF.h"
#include "analyticalF.h"
#include "core/olbDebug.h"

namespace olb {



//////////////////////////////// AnalyticCalcF1D ////////////////////////////////
template <typename T, typename S, template<typename> class F>
AnalyticCalcF1D<T,S,F>::AnalyticCalcF1D(FunctorPtr<AnalyticalF1D<T,S>>&& f,
                                        FunctorPtr<AnalyticalF1D<T,S>>&& g)
  : AnalyticalF1D<T,S>(f->getTargetDim()),
    _f(std::move(f)),
    _g(std::move(g))
{
  OLB_ASSERT(g->getTargetDim() == f->getTargetDim(),
             "the dimensions of both functors need to be equal");
  std::swap(f->_ptrCalcC, this->_ptrCalcC);
  this->getName() = "(" + _f->getName() + F<T>::symbol + _g->getName() + ")";
}

template <typename T, typename S, template<typename> class F>
AnalyticCalcF1D<T,S,F>::AnalyticCalcF1D(T scalar, FunctorPtr<AnalyticalF1D<T,S>>&& g)
  : AnalyticCalcF1D(
      std::unique_ptr<AnalyticalF1D<T,S>>(
        new AnalyticalConst1D<T,S>(std::vector<T>(g->getTargetDim(), scalar))),
      std::forward<decltype(g)>(g))
{ }

template <typename T, typename S, template<typename> class F>
AnalyticCalcF1D<T,S,F>::AnalyticCalcF1D(FunctorPtr<AnalyticalF1D<T,S>>&& f, T scalar)
  : AnalyticCalcF1D(
      std::forward<decltype(f)>(f),
      std::unique_ptr<AnalyticalF1D<T,S>>(
        new AnalyticalConst1D<T,S>(std::vector<T>(f->getTargetDim(), scalar))))
{ }

template <typename T, typename S, template<typename> class F>
bool AnalyticCalcF1D<T,S,F>::operator()(T output[], const S input[])
{
  T outputTmp[this->_g->getTargetDim()];
  this->_g(outputTmp, input);
  this->_f(output, input);
  for (int i = 0; i < this->_f->getTargetDim(); ++i) {
    output[i] = F<T>()(output[i], outputTmp[i]);
  }
  return true;
}


template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator+(std::shared_ptr<AnalyticalF1D<T,S>> lhs, std::shared_ptr<AnalyticalF1D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF1D<T,S>>(
           new AnalyticCalcPlus1D<T,S>(std::move(lhs), std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator+(std::shared_ptr<AnalyticalF1D<T,S>> lhs, T rhs)
{
  return std::shared_ptr<AnalyticalF1D<T,S>>(
           new AnalyticCalcPlus1D<T,S>(std::move(lhs), rhs));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator+(T lhs, std::shared_ptr<AnalyticalF1D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF1D<T,S>>(
           new AnalyticCalcPlus1D<T,S>(lhs, std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator-(std::shared_ptr<AnalyticalF1D<T,S>> lhs, std::shared_ptr<AnalyticalF1D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF1D<T,S>>(
           new AnalyticCalcMinus1D<T,S>(std::move(lhs), std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator-(std::shared_ptr<AnalyticalF1D<T,S>> lhs, T rhs)
{
  return std::shared_ptr<AnalyticalF1D<T,S>>(
           new AnalyticCalcMinus1D<T,S>(std::move(lhs), rhs));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator-(T lhs, std::shared_ptr<AnalyticalF1D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF1D<T,S>>(
           new AnalyticCalcMinus1D<T,S>(lhs, std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator*(std::shared_ptr<AnalyticalF1D<T,S>> lhs, std::shared_ptr<AnalyticalF1D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF1D<T,S>>(
           new AnalyticCalcMultiplication1D<T,S>(std::move(lhs), std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator*(std::shared_ptr<AnalyticalF1D<T,S>> lhs, T rhs)
{
  return std::shared_ptr<AnalyticalF1D<T,S>>(
           new AnalyticCalcMultiplication1D<T,S>(std::move(lhs), rhs));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator*(T lhs, std::shared_ptr<AnalyticalF1D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF1D<T,S>>(
           new AnalyticCalcMultiplication1D<T,S>(lhs, std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator/(std::shared_ptr<AnalyticalF1D<T,S>> lhs, std::shared_ptr<AnalyticalF1D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF1D<T,S>>(
           new AnalyticCalcDivision1D<T,S>(std::move(lhs), std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator/(std::shared_ptr<AnalyticalF1D<T,S>> lhs, T rhs)
{
  return std::shared_ptr<AnalyticalF1D<T,S>>(
           new AnalyticCalcDivision1D<T,S>(std::move(lhs), rhs));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF1D<T,S>> operator/(T lhs, std::shared_ptr<AnalyticalF1D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF1D<T,S>>(
           new AnalyticCalcDivision1D<T,S>(lhs, std::move(rhs)));
}

/////////////////////////////////operator()/// ////////////////////////////////
template <typename T, typename S>
AnalyticalF1D<T,S>& AnalyticalF1D<T,S>::operator+(AnalyticalF1D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticCalcPlus1D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF1D<T,S>& AnalyticalF1D<T,S>::operator-(AnalyticalF1D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticCalcMinus1D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF1D<T,S>& AnalyticalF1D<T,S>::operator*(AnalyticalF1D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticCalcMultiplication1D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF1D<T,S>& AnalyticalF1D<T,S>::operator/(AnalyticalF1D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticCalcDivision1D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}


//////////////////////////////// AnalyticCalcF2D ////////////////////////////////
template <typename T, typename S, template<typename> class F>
AnalyticCalcF2D<T,S,F>::AnalyticCalcF2D(FunctorPtr<AnalyticalF2D<T,S>>&& f,
                                        FunctorPtr<AnalyticalF2D<T,S>>&& g)
  : AnalyticalF2D<T,S>(f->getTargetDim()),
    _f(std::move(f)),
    _g(std::move(g))
{
  OLB_ASSERT(g->getTargetDim() == f->getTargetDim(),
             "the dimensions of both functors need to be equal");
  // pass through the shared_ptr from the first argument f to the arithmetic class itself.
  // used by secsessive calls: e.g. (functorA + functor B) followed by (functorA + functorC)
  // the result of the first operation is overwritten by the second.

  // equivalent operations
  //  std::swap(f._ptrCalcC, this->_ptrCalcC);
  //  this->_ptrCalcC = f._ptrCalcC;
  this->_ptrCalcC.swap(f->_ptrCalcC);

  this->getName() = "(" + _f->getName() + F<T>::symbol + _g->getName() + ")";
}

template <typename T, typename S, template<typename> class F>
AnalyticCalcF2D<T,S,F>::AnalyticCalcF2D(T scalar, FunctorPtr<AnalyticalF2D<T,S>>&& g)
  : AnalyticCalcF2D(
      std::unique_ptr<AnalyticalF2D<T,S>>(
        new AnalyticalConst2D<T,S>(std::vector<T>(g->getTargetDim(), scalar))),
      std::forward<decltype(g)>(g))
{ }

template <typename T, typename S, template<typename> class F>
AnalyticCalcF2D<T,S,F>::AnalyticCalcF2D(FunctorPtr<AnalyticalF2D<T,S>>&& f, T scalar)
  : AnalyticCalcF2D(
      std::forward<decltype(f)>(f),
      std::unique_ptr<AnalyticalF2D<T,S>>(
        new AnalyticalConst2D<T,S>(std::vector<T>(f->getTargetDim(), scalar))))
{ }

template <typename T, typename S, template<typename> class F>
bool AnalyticCalcF2D<T,S,F>::operator()(T output[], const S input[])
{
  T outputTmp[this->_g->getTargetDim()];
  this->_g(outputTmp, input);
  this->_f(output, input);
  for (int i = 0; i < this->_f->getTargetDim(); ++i) {
    output[i] = F<T>()(output[i], outputTmp[i]);
  }
  return true;
}


template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator+(std::shared_ptr<AnalyticalF2D<T,S>> lhs, std::shared_ptr<AnalyticalF2D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF2D<T,S>>(
           new AnalyticCalcPlus2D<T,S>(std::move(lhs), std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator+(std::shared_ptr<AnalyticalF2D<T,S>> lhs, T rhs)
{
  return std::shared_ptr<AnalyticalF2D<T,S>>(
           new AnalyticCalcPlus2D<T,S>(std::move(lhs), rhs));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator+(T lhs, std::shared_ptr<AnalyticalF2D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF2D<T,S>>(
           new AnalyticCalcPlus2D<T,S>(lhs, std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator-(std::shared_ptr<AnalyticalF2D<T,S>> lhs, std::shared_ptr<AnalyticalF2D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF2D<T,S>>(
           new AnalyticCalcMinus2D<T,S>(std::move(lhs), std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator-(std::shared_ptr<AnalyticalF2D<T,S>> lhs, T rhs)
{
  return std::shared_ptr<AnalyticalF2D<T,S>>(
           new AnalyticCalcMinus2D<T,S>(std::move(lhs), rhs));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator-(T lhs, std::shared_ptr<AnalyticalF2D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF2D<T,S>>(
           new AnalyticCalcMinus2D<T,S>(lhs, std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator*(std::shared_ptr<AnalyticalF2D<T,S>> lhs, std::shared_ptr<AnalyticalF2D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF2D<T,S>>(
           new AnalyticCalcMultiplication2D<T,S>(std::move(lhs), std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator*(std::shared_ptr<AnalyticalF2D<T,S>> lhs, T rhs)
{
  return std::shared_ptr<AnalyticalF2D<T,S>>(
           new AnalyticCalcMultiplication2D<T,S>(std::move(lhs), rhs));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator*(T lhs, std::shared_ptr<AnalyticalF2D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF2D<T,S>>(
           new AnalyticCalcMultiplication2D<T,S>(lhs, std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator/(std::shared_ptr<AnalyticalF2D<T,S>> lhs, std::shared_ptr<AnalyticalF2D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF2D<T,S>>(
           new AnalyticCalcDivision2D<T,S>(std::move(lhs), std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator/(std::shared_ptr<AnalyticalF2D<T,S>> lhs, T rhs)
{
  return std::shared_ptr<AnalyticalF2D<T,S>>(
           new AnalyticCalcDivision2D<T,S>(std::move(lhs), rhs));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> operator/(T lhs, std::shared_ptr<AnalyticalF2D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF2D<T,S>>(
           new AnalyticCalcDivision2D<T,S>(lhs, std::move(rhs)));
}

/////////////////////////////////operator()////////////////////////////////////
template <typename T, typename S>
AnalyticalF2D<T,S>& AnalyticalF2D<T,S>::operator+(AnalyticalF2D<T,S>& rhs)
{
  // version 1
  //    AnalyticalF2D<T,S>* tmp = new AnalyticCalcPlus2D<T,S>(*this,rhs);
  //    std::shared_ptr< GenericF<T,S> > ptr( tmp );
  //    this->_ptrCalcC = ptr;

  // version 2
  //  std::shared_ptr< AnalyticalF2D<T,S> > tmp = std::make_shared< AnalyticCalcPlus2D<T,S> >(*this,rhs);

  // version 2.5
  //  std::shared_ptr< AnalyticCalcPlus2D<T,S> > tmp( new AnalyticCalcPlus2D<T,S>(*this,rhs) );

  // version 3
  auto tmp = std::make_shared< AnalyticCalcPlus2D<T,S> >(*this,rhs);

  this->_ptrCalcC = tmp;

  //  std::cout << "operator+(): " << this->_ptrCalcC.get()->getName() << std::endl;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF2D<T,S>& AnalyticalF2D<T,S>::operator-(AnalyticalF2D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticCalcMinus2D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF2D<T,S>& AnalyticalF2D<T,S>::operator*(AnalyticalF2D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticCalcMultiplication2D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF2D<T,S>& AnalyticalF2D<T,S>::operator/(AnalyticalF2D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticCalcDivision2D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}


//////////////////////////////// AnalyticCalcF3D ////////////////////////////////
template <typename T, typename S, template<typename> class F>
AnalyticCalcF3D<T,S,F>::AnalyticCalcF3D(FunctorPtr<AnalyticalF3D<T,S>>&& f,
                                        FunctorPtr<AnalyticalF3D<T,S>>&& g):
  AnalyticalF3D<T,S>(f->getTargetDim()),
  _f(std::move(f)),
  _g(std::move(g))
{
  OLB_ASSERT(g->getTargetDim() == f->getTargetDim(),
             "the dimensions of both functors need to be equal");
  std::swap(f->_ptrCalcC, this->_ptrCalcC);
  this->getName() = "(" + _f->getName() + F<T>::symbol + _g->getName() + ")";
}

template <typename T, typename S, template<typename> class F>
AnalyticCalcF3D<T,S,F>::AnalyticCalcF3D(T scalar, FunctorPtr<AnalyticalF3D<T,S>>&& g)
  : AnalyticCalcF3D(
      std::unique_ptr<AnalyticalF3D<T,S>>(
        new AnalyticalConst3D<T,S>(std::vector<T>(g->getTargetDim(), scalar))),
      std::forward<decltype(g)>(g))
{ }

template <typename T, typename S, template<typename> class F>
AnalyticCalcF3D<T,S,F>::AnalyticCalcF3D(FunctorPtr<AnalyticalF3D<T,S>>&& f, T scalar)
  : AnalyticCalcF3D(
      std::forward<decltype(f)>(f),
      std::unique_ptr<AnalyticalF3D<T,S>>(
        new AnalyticalConst3D<T,S>(std::vector<T>(f->getTargetDim(), scalar))))
{ }

template <typename T, typename S, template<typename> class F>
bool AnalyticCalcF3D<T,S,F>::operator()(T output[], const S input[])
{
  T outputTmp[this->_g->getTargetDim()];
  this->_g(outputTmp, input);
  this->_f(output, input);
  for (int i = 0; i < this->_f->getTargetDim(); ++i) {
    output[i] = F<T>()(output[i], outputTmp[i]);
  }
  return true;
}


template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator+(std::shared_ptr<AnalyticalF3D<T,S>> lhs, std::shared_ptr<AnalyticalF3D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF3D<T,S>>(
           new AnalyticCalcPlus3D<T,S>(std::move(lhs), std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator+(std::shared_ptr<AnalyticalF3D<T,S>> lhs, T rhs)
{
  return std::shared_ptr<AnalyticalF3D<T,S>>(
           new AnalyticCalcPlus3D<T,S>(std::move(lhs), rhs));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator+(T lhs, std::shared_ptr<AnalyticalF3D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF3D<T,S>>(
           new AnalyticCalcPlus3D<T,S>(lhs, std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator-(std::shared_ptr<AnalyticalF3D<T,S>> lhs, std::shared_ptr<AnalyticalF3D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF3D<T,S>>(
           new AnalyticCalcMinus3D<T,S>(std::move(lhs), std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator-(std::shared_ptr<AnalyticalF3D<T,S>> lhs, T rhs)
{
  return std::shared_ptr<AnalyticalF3D<T,S>>(
           new AnalyticCalcMinus3D<T,S>(std::move(lhs), rhs));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator-(T lhs, std::shared_ptr<AnalyticalF3D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF3D<T,S>>(
           new AnalyticCalcMinus3D<T,S>(lhs, std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator*(std::shared_ptr<AnalyticalF3D<T,S>> lhs, std::shared_ptr<AnalyticalF3D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF3D<T,S>>(
           new AnalyticCalcMultiplication3D<T,S>(std::move(lhs), std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator*(std::shared_ptr<AnalyticalF3D<T,S>> lhs, T rhs)
{
  return std::shared_ptr<AnalyticalF3D<T,S>>(
           new AnalyticCalcMultiplication3D<T,S>(std::move(lhs), rhs));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator*(T lhs, std::shared_ptr<AnalyticalF3D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF3D<T,S>>(
           new AnalyticCalcMultiplication3D<T,S>(lhs, std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator/(std::shared_ptr<AnalyticalF3D<T,S>> lhs, std::shared_ptr<AnalyticalF3D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF3D<T,S>>(
           new AnalyticCalcDivision3D<T,S>(std::move(lhs), std::move(rhs)));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator/(std::shared_ptr<AnalyticalF3D<T,S>> lhs, T rhs)
{
  return std::shared_ptr<AnalyticalF3D<T,S>>(
           new AnalyticCalcDivision3D<T,S>(std::move(lhs), rhs));
}

template <typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> operator/(T lhs, std::shared_ptr<AnalyticalF3D<T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF3D<T,S>>(
           new AnalyticCalcDivision3D<T,S>(lhs, std::move(rhs)));
}

/////////////////////////////////operator()/// ////////////////////////////////
template <typename T, typename S>
AnalyticalF3D<T,S>& AnalyticalF3D<T,S>::operator+(AnalyticalF3D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticCalcPlus3D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF3D<T,S>& AnalyticalF3D<T,S>::operator-(AnalyticalF3D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticCalcMinus3D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF3D<T,S>& AnalyticalF3D<T,S>::operator*(AnalyticalF3D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticCalcMultiplication3D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename S>
AnalyticalF3D<T,S>& AnalyticalF3D<T,S>::operator/(AnalyticalF3D<T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticCalcDivision3D<T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

} // end namespace olb

#endif
