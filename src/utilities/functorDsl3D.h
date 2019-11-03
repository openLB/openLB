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

#ifndef FUNCTOR_DSL_3D_H
#define FUNCTOR_DSL_3D_H

#include "functors/lattice/superBaseF3D.h"
#include "functors/lattice/indicator/superIndicatorF3D.h"
#include "functors/analytical/analyticalBaseF.h"

namespace olb {

/// Helper functions for building functors via composition
/**
 * Certain types of functors such as error norms (see e.g. SuperRelativeErrorLpNorm3D)
 * lend themselves to being expressed as a composition of existing functors.
 *
 * std::shared_ptr based functor arithmetic offers a memory-safe way to perform
 * functor composition. However common operations such as lifting stack-allocated
 * functors into the std::shared_ptr context, constructing Lp norms or restricting
 * analytical functors on a lattice can become quite verbose.
 *
 * This is why the `functor_dsl` namespace contains a set of conveniently named
 * helper functions that aim to ease writing composed functors.
 *
 * e.g. a slightly modified version of SuperRelativeErrorLpNorm3D's constructor:
 *
 * \code{.cpp}
 * using namespace functor_dsl;
 *
 * // decltype(wantedF)    == std::shared_ptr<AnalyticalF3D<T,W>>
 * // decltype(f)          == std::shared_ptr<SuperF3D<T,W>>
 * // decltype(indicatorF) == std::shared_ptr<SuperIndicatorF3D<T>>
 *
 * auto wantedLatticeF = restrict(wantedF, sLattice);
 * auto relaticeErrorF = norm<P>(wantedLatticeF - f, indicatorF)
 *                     / norm<P>(wantedLatticeF, indicatorF);
 * \endcode
 *
 * Note that any of the functors managed by std::shared_ptr may be returned
 * and reused independently.
 *
 * i.e. it is not necessary to manually assure the lifetime of any functor
 * used by the composition.
 **/
namespace functor_dsl {

/// Lifts functor reference to std::shared_ptr functor arithmetic
template<typename T, typename W>
std::shared_ptr<SuperF3D<T,W>> lift(SuperF3D<T,W>& f);

/// Lifts functor pointer to std::shared_ptr functor arithmetic
template<typename T, typename W>
std::shared_ptr<SuperF3D<T,W>> lift(SuperF3D<T,W>* f);

/// Lifts functor reference to std::shared_ptr functor arithmetic
template<typename T, typename S>
std::shared_ptr<AnalyticalF3D<T,S>> lift(AnalyticalF3D<T,S>& f);

/// Lifts functor pointer to std::shared_ptr functor arithmetic
template<typename T, typename W>
std::shared_ptr<AnalyticalF3D<T,S>> lift(AnalyticalF3D<T,S>* f);

/// Returns baseF raised to the power of exponentF
template<typename T, typename W>
std::shared_ptr<SuperF3D<T,W>> pow(std::shared_ptr<SuperF3D<T,W>> baseF,
                                   std::shared_ptr<SuperF3D<T,W>> exponentF);

/// Returns baseF raised to the power of exponent
template<typename T, typename W, typename E>
std::shared_ptr<SuperF3D<T,W>> pow(std::shared_ptr<SuperF3D<T,W>> baseF,
                                   E exponent);

/// Returns base raised to the power of exponentF
template<typename T, typename W, typename B>
std::shared_ptr<SuperF3D<T,W>> pow(B base,
                                   std::shared_ptr<SuperF3D<T,W>> exponentF);

/// Returns Lp norm for a functor f on the subset described by indicatorF
template<int P, typename T, typename W>
std::shared_ptr<SuperF3D<T,W>> norm(std::shared_ptr<SuperF3D<T,W>>        f,
                                    std::shared_ptr<SuperIndicatorF3D<T>> indicatorF);

/// Returns restriction of a analytical functor f to the lattice sLattice
template<typename T, typename W, typename DESCRIPTOR>
std::shared_ptr<SuperF3D<T,W>> restrict(std::shared_ptr<AnalyticalF3D<T,W>> f,
                                        SuperLattice3D<T, DESCRIPTOR>& sLattice);

}

}

#endif
