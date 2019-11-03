/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#ifndef FUNCTOR_DSL_2D_HH
#define FUNCTOR_DSL_2D_HH

#include "functorPtr.h"
#include "functors/lattice/superBaseF2D.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"
#include "functors/lattice/integral/superLpNorm2D.h"

namespace olb {

namespace functor_dsl {

template<typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> lift(SuperF2D<T,W>& f)
{
  return FunctorPtr<SuperF2D<T,W>>(f).toShared();
}

template<typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> lift(SuperF2D<T,W>* f)
{
  return FunctorPtr<SuperF2D<T,W>>(f).toShared();
}

template<typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> lift(AnalyticalF2D<T,S>& f)
{
  return FunctorPtr<AnalyticalF2D<T,S>>(f).toShared();
}

template<typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> lift(AnalyticalF2D<T,S>* f)
{
  return FunctorPtr<AnalyticalF2D<T,S>>(f).toShared();
}

template<typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> pow(std::shared_ptr<SuperF2D<T,W>> baseF,
                                   std::shared_ptr<SuperF2D<T,W>> exponentF)
{
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperCalcPower2D<T,W>(std::move(baseF),
                                     std::move(exponentF))
         );
}

template<typename T, typename W, typename E>
std::shared_ptr<SuperF2D<T,W>> pow(std::shared_ptr<SuperF2D<T,W>> baseF,
                                   E exponent)
{
  static_assert(std::is_arithmetic<E>::value,
                "Exponent must be an arithmetic value");
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperCalcPower2D<T,W>(std::move(baseF), exponent));
}

template<typename T, typename W, typename B>
std::shared_ptr<SuperF2D<T,W>> pow(B base,
                                   std::shared_ptr<SuperF2D<T,W>> exponentF)
{
  static_assert(std::is_arithmetic<B>::value,
                "Base must be an arithmetic value");
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperCalcPower2D<T,W>(base, std::move(exponentF)));
}

template<int P, typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> norm(std::shared_ptr<SuperF2D<T,W>>        f,
                                    std::shared_ptr<SuperIndicatorF2D<T>> indicatorF)
{
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperLpNorm2D<T,W,P>(std::move(f),
                                    std::move(indicatorF))
         );
}

template<typename T, typename W, typename DESCRIPTOR>
std::shared_ptr<SuperF2D<T,W>> restrict(std::shared_ptr<AnalyticalF2D<T,W>> f,
                                        SuperLattice2D<T, DESCRIPTOR>& sLattice)
{
  return std::shared_ptr<SuperF2D<T,W>>(
           new SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR>(std::move(f), sLattice));
}

}

}

#endif
