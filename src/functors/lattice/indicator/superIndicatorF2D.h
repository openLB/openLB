/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerlaender
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

#ifndef SUPER_INDICATOR_F_2D_H
#define SUPER_INDICATOR_F_2D_H

#include "functors/analytical/indicator/indicatorBaseF2D.h"
#include "superIndicatorBaseF2D.h"
#include "geometry/superGeometry2D.h"
#include "functors/analytical/indicator/smoothIndicatorBaseF2D.h"

namespace olb {


/// SuperIndicatorF2D from IndicatorF2D
/**
 * Maintains block level BlockIndicatorFfromIndicatorF2D instances in SuperF2D::_blockF.
 **/
template <typename T>
class SuperIndicatorFfromIndicatorF2D : public SuperIndicatorF2D<T> {
protected:
  IndicatorF2D<T>& _indicatorF;
public:
  /**
   * \param indicatorF Indicator to be converted into a super indicator
   * \param geometry   Super geometry required for block indicator construction
   **/
  SuperIndicatorFfromIndicatorF2D(IndicatorF2D<T>&    indicatorF,
                                  SuperGeometry2D<T>& geometry);

  using SuperIndicatorF2D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
};


/// SuperIndicatorF2D from SmoothIndicatorF2D
/**
 * Returns true if SmoothIndicatorF2D output is not near zero.
 *
 * Exposes block level indicators via SuperIndicatorF2D::getBlockIndicatorF
 * respectively SuperIndicatorF2D::getExtendedBlockIndicatorF
 **/
template <typename T, bool HLBM>
class SuperIndicatorFfromSmoothIndicatorF2D : public SuperIndicatorF2D<T> {
protected:
  FunctorPtr<SmoothIndicatorF2D<T,T,HLBM>> _indicatorF;
public:
  /**
   * \param indicatorF Indicator to be converted into a super indicator
   * \param geometry   Super geometry required for block indicator construction
   **/
  SuperIndicatorFfromSmoothIndicatorF2D(FunctorPtr<SmoothIndicatorF2D<T,T,HLBM>>&& indicatorF,
                                          SuperGeometry2D<T>&                     geometry);

  using SuperIndicatorF2D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
};


/// Indicator functor from material numbers
/**
 * Maintains block level BlockIndicatorMaterial2D instances in SuperF2D::_blockF.
 **/
template <typename T>
class SuperIndicatorMaterial2D : public SuperIndicatorF2D<T> {
public:
  /**
   * \param geometry  Super geometry required for block indicator construction
   *                  and global material number queries
   * \param materials Vector of material numbers to be indicated
   **/
  SuperIndicatorMaterial2D(SuperGeometry2D<T>& geometry, std::vector<int> materials);

  using SuperIndicatorF2D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
};


/// Indicator identity functor
/**
 * Proxies a given indicator functor for simplified memory management.
 * i.e. mixing non-owning indicator references and owning indicator pointers
 *      in functor compositions.
 **/
template <typename T>
class SuperIndicatorIdentity2D : public SuperIndicatorF2D<T> {
protected:
  FunctorPtr<SuperIndicatorF2D<T>> _indicatorF;
public:
  SuperIndicatorIdentity2D(FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF);

  using SuperIndicatorF2D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
};


} // namespace olb

#endif
