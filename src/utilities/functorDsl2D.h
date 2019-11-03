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

#ifndef FUNCTOR_DSL_2D_H
#define FUNCTOR_DSL_2D_H

#include "functors/lattice/superBaseF2D.h"
#include "functors/analytical/analyticalBaseF.h"

namespace olb {

namespace functor_dsl {

/// Lifts functor reference to std::shared_ptr functor arithmetic
template<typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> lift(SuperF2D<T,W>& f);

/// Lifts functor pointer to std::shared_ptr functor arithmetic
template<typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> lift(SuperF2D<T,W>* f);

/// Lifts functor reference to std::shared_ptr functor arithmetic
template<typename T, typename S>
std::shared_ptr<AnalyticalF2D<T,S>> lift(AnalyticalF2D<T,S>& f);

/// Lifts functor pointer to std::shared_ptr functor arithmetic
template<typename T, typename W>
std::shared_ptr<AnalyticalF2D<T,S>> lift(AnalyticalF2D<T,S>* f);

/// Returns baseF raised to the power of exponentF
template<typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> pow(std::shared_ptr<SuperF2D<T,W>> baseF,
                                   std::shared_ptr<SuperF2D<T,W>> exponentF);

/// Returns baseF raised to the power of exponent
template<typename T, typename W, typename E>
std::shared_ptr<SuperF2D<T,W>> pow(std::shared_ptr<SuperF2D<T,W>> baseF,
                                   E exponent);

/// Returns base raised to the power of exponentF
template<typename T, typename W, typename B>
std::shared_ptr<SuperF2D<T,W>> pow(B base,
                                   std::shared_ptr<SuperF2D<T,W>> exponentF);

/// Returns Lp norm for a functor f on the subset described by indicatorF
template<int P, typename T, typename W>
std::shared_ptr<SuperF2D<T,W>> norm(std::shared_ptr<SuperF2D<T,W>>        f,
                                    std::shared_ptr<SuperIndicatorF2D<T>> indicatorF);

/// Returns restriction of a analytical functor f to the lattice sLattice
template<typename T, typename W, typename DESCRIPTOR>
std::shared_ptr<SuperF2D<T,W>> restrict(std::shared_ptr<AnalyticalF2D<T,W>> f,
                                        SuperLattice2D<T, DESCRIPTOR>& sLattice);

}

}

#endif
