/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016-2018 Benjamin Förster, Adrian Kummerlaender
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

#ifndef SUPER_INDICATOR_F_3D_H
#define SUPER_INDICATOR_F_3D_H

#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "superIndicatorBaseF3D.h"
#include "geometry/superGeometry3D.h"
#include "functors/analytical/indicator/smoothIndicatorBaseF3D.h"

namespace olb {


/// SuperIndicatorF3D from IndicatorF3D
/**
 * Exposes block level indicators via SuperIndicatorF3D::getBlockIndicatorF
 * respectively SuperIndicatorF3D::getExtendedBlockIndicatorF
 **/
template <typename T>
class SuperIndicatorFfromIndicatorF3D : public SuperIndicatorF3D<T> {
protected:
  FunctorPtr<IndicatorF3D<T>> _indicatorF;
public:
  /**
   * \param indicatorF Indicator to be converted into a super indicator
   * \param geometry   Super geometry required for block indicator construction
   **/
  SuperIndicatorFfromIndicatorF3D(FunctorPtr<IndicatorF3D<T>>&& indicatorF,
                                  SuperGeometry3D<T>&           geometry);

  using SuperIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
};


/// SuperIndicatorF3D from SmoothIndicatorF3D
/**
 * Returns true iff SmoothIndicatorF3D output is not near zero.
 *
 * Exposes block level indicators via SuperIndicatorF3D::getBlockIndicatorF
 * respectively SuperIndicatorF3D::getExtendedBlockIndicatorF
 **/
template <typename T, bool HLBM>
class SuperIndicatorFfromSmoothIndicatorF3D : public SuperIndicatorF3D<T> {
protected:
  FunctorPtr<SmoothIndicatorF3D<T,T,HLBM>> _indicatorF;
public:
  /**
   * \param indicatorF Indicator to be converted into a super indicator
   * \param geometry   Super geometry required for block indicator construction
   **/
  SuperIndicatorFfromSmoothIndicatorF3D(FunctorPtr<SmoothIndicatorF3D<T,T,HLBM>>&& indicatorF,
                                          SuperGeometry3D<T>&                     geometry);

  using SuperIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
};


/// Indicator functor from material numbers
/**
 * Exposes block level indicators via SuperIndicatorF3D::getBlockIndicatorF
 * respectively SuperIndicatorF3D::getExtendedBlockIndicatorF
 *
 * Material number comparsion is implemented in BlockIndicatorMaterial3D.
 **/
template <typename T>
class SuperIndicatorMaterial3D : public SuperIndicatorF3D<T> {
public:
  /**
   * \param geometry  Super geometry required for block indicator construction
   *                  and global material number queries
   * \param materials Vector of material numbers to be indicated
   **/
  SuperIndicatorMaterial3D(SuperGeometry3D<T>& geometry, std::vector<int> materials);
  /**
   * \param geometry  Super geometry required for block indicator construction
   *                  and global material number queries
   * \param materials List of material numbers to be indicated
   **/
  SuperIndicatorMaterial3D(SuperGeometry3D<T>& geometry, std::list<int> materials);

  using SuperIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
};


/// Indicator identity functor
/**
 * Proxies a given indicator functor for simplified memory management.
 * i.e. mixing non-owning indicator references and owning indicator pointers
 *      in functor compositions.
 **/
template <typename T>
class SuperIndicatorIdentity3D : public SuperIndicatorF3D<T> {
protected:
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
public:
  SuperIndicatorIdentity3D(FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);

  using SuperIndicatorF3D<T>::operator();
  bool operator() (bool output[], const int input[]) override;
};


} // namespace olb

#endif
