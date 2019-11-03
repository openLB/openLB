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

#ifndef SUPER_CALC_F_2D_HH
#define SUPER_CALC_F_2D_HH

#include "superCalcF2D.h"
#include "blockCalcF2D.h"
#include "superConst2D.h"

namespace olb {


template <typename T, typename W, template<typename> class F>
SuperCalcF2D<T,W,F>::SuperCalcF2D(FunctorPtr<SuperF2D<T,W>>&& f,
                                  FunctorPtr<SuperF2D<T,W>>&& g)
  : SuperF2D<T,W>(
      f->getSuperStructure(),
      f->getTargetDim() > g->getTargetDim() ? f->getTargetDim() : g->getTargetDim()),
    _f(std::move(f)),
    _g(std::move(g))
{
  OLB_ASSERT(
    _f->getTargetDim() == _g->getTargetDim() || _f->getTargetDim() == 1 || _g->getTargetDim() == 1,
    "Componentwise operation must be well defined.");

  this->getName() = "(" + _f->getName() + F<T>::symbol + _g->getName() + ")";

  std::swap(_f->_ptrCalcC, this->_ptrCalcC);

  LoadBalancer<T>& load = _f->getSuperStructure().getLoadBalancer();

  if ( _f->getBlockFSize() == load.size() ) {
    if ( _g->getBlockFSize() == load.size() ) {
      // both functors expose the correct count of block level functors
      for (int iC = 0; iC < load.size(); ++iC) {
        this->_blockF.emplace_back(
          new BlockCalcF2D<W,F>(_f->getBlockF(iC), _g->getBlockF(iC))
        );
      }
    }
    else {
      // operate on super functor `g` and block level functors provided by `f`
      for (int iC = 0; iC < load.size(); ++iC) {
        this->_blockF.emplace_back(
          new BlockCalcF2D<W,F>(_f->getBlockF(iC), *g, load.glob(iC))
        );
      }
    }
  }
  else if ( _g->getBlockFSize() == load.size() ) {
    // operate on block level functors provided by `f` and super functor `g`
    for (int iC = 0; iC < load.size(); ++iC) {
      this->_blockF.emplace_back(
        new BlockCalcF2D<W,F>(*f, load.glob(iC), _g->getBlockF(iC))
      );
    }
  }
}

template <typename T, typename W, template<typename> class F>
SuperCalcF2D<T,W,F>::SuperCalcF2D(W scalar, FunctorPtr<SuperF2D<T,W>>&& g)
  : SuperCalcF2D(
      std::unique_ptr<SuperF2D<T,W>>(new SuperConst2D<T,W>(g->getSuperStructure(), scalar)),
      std::forward<decltype(g)>(g))
{ }

template <typename T, typename W, template<typename> class F>
SuperCalcF2D<T,W,F>::SuperCalcF2D(FunctorPtr<SuperF2D<T,W>>&& f, W scalar)
  : SuperCalcF2D(
      std::forward<decltype(f)>(f),
      std::unique_ptr<SuperF2D<T,W>>(new SuperConst2D<T,W>(f->getSuperStructure(), scalar)))
{ }

template <typename T, typename W, template<typename> class F>
bool SuperCalcF2D<T,W,F>::operator()(W output[], const int input[])
{
  if ( _f->getTargetDim() == 1 || _g->getTargetDim() == 1 ) {
    // scalar operation
    W scalar;
    if ( _f->getTargetDim() == 1 ) {
      // apply the scalar f to possibly multidimensional g
      _f(&scalar, input);
      _g(output, input);

      for (int i = 0; i < this->getTargetDim(); i++) {
        output[i] = F<T>()(scalar, output[i]);
      }
    }
    else {
      // apply scalar g to possibly multidimensional f
      _f(output, input);
      _g(&scalar, input);

      for (int i = 0; i < this->getTargetDim(); i++) {
        output[i] = F<T>()(output[i], scalar);
      }
    }
  }
  else {
    // componentwise operation on equidimensional functors
    W* outputF = output;
    W outputG[this->getTargetDim()];

    _f(outputF, input);
    _g(outputG, input);

    for (int i = 0; i < this->getTargetDim(); i++) {
      output[i] = F<T>()(outputF[i], outputG[i]);
    }
  }
  return true;
}


template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator+(std::shared_ptr<SuperF2D<T,W>> lhs, std::shared_ptr<SuperF2D<T,W>> rhs)
{
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperCalcPlus2D<T,W>(std::move(lhs), std::move(rhs)));
}

template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator+(std::shared_ptr<SuperF2D<T,W>> lhs, W rhs)
{
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperCalcPlus2D<T,W>(std::move(lhs), rhs));
}

template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator+(W lhs, std::shared_ptr<SuperF2D<T,W>> rhs)
{
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperCalcPlus2D<T,W>(lhs, std::move(rhs)));
}

template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator-(std::shared_ptr<SuperF2D<T,W>> lhs, std::shared_ptr<SuperF2D<T,W>> rhs)
{
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperCalcMinus2D<T,W>(std::move(lhs), std::move(rhs)));
}

template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator-(std::shared_ptr<SuperF2D<T,W>> lhs, W rhs)
{
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperCalcMinus2D<T,W>(std::move(lhs), rhs));
}

template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator-(W lhs, std::shared_ptr<SuperF2D<T,W>> rhs)
{
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperCalcMinus2D<T,W>(lhs, std::move(rhs)));
}

template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator*(std::shared_ptr<SuperF2D<T,W>> lhs, std::shared_ptr<SuperF2D<T,W>> rhs)
{
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperCalcMultiplication2D<T,W>(std::move(lhs), std::move(rhs)));
}

template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator*(std::shared_ptr<SuperF2D<T,W>> lhs, W rhs)
{
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperCalcMultiplication2D<T,W>(std::move(lhs), rhs));
}

template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator*(W lhs, std::shared_ptr<SuperF2D<T,W>> rhs)
{
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperCalcMultiplication2D<T,W>(lhs, std::move(rhs)));
}

template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator/(std::shared_ptr<SuperF2D<T,W>> lhs, std::shared_ptr<SuperF2D<T,W>> rhs)
{
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperCalcDivision2D<T,W>(std::move(lhs), std::move(rhs)));
}

template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator/(std::shared_ptr<SuperF2D<T,W>> lhs, W rhs)
{
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperCalcDivision2D<T,W>(std::move(lhs), rhs));
}

template <typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> operator/(W lhs, std::shared_ptr<SuperF2D<T,W>> rhs)
{
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperCalcDivision2D<T,W>(lhs, std::move(rhs)));
}


template <typename T, typename W>
SuperF2D<T,W>& SuperF2D<T,W>::operator+(SuperF2D<T,W>& rhs)
{
  auto tmp = std::make_shared< SuperCalcPlus2D<T,W> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename W>
SuperF2D<T,W>& SuperF2D<T,W>::operator-(SuperF2D<T,W>& rhs)
{
  auto tmp = std::make_shared< SuperCalcMinus2D<T,W> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename W>
SuperF2D<T,W>& SuperF2D<T,W>::operator*(SuperF2D<T,W>& rhs)
{
  auto tmp = std::make_shared< SuperCalcMultiplication2D<T,W> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <typename T, typename W>
SuperF2D<T,W>& SuperF2D<T,W>::operator/(SuperF2D<T,W>& rhs)
{
  auto tmp = std::make_shared< SuperCalcDivision2D<T,W> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}


} // end namespace olb

#endif
