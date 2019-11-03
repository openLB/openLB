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

#ifndef SUPER_MAX_3D_H
#define SUPER_MAX_3D_H

#include "superBaseF3D.h"
#include "blockMax3D.h"
#include "indicator/superIndicatorBaseF3D.h"
#include "utilities/functorPtr.h"

namespace olb {


/// SuperMax3D returns the max in each component of f on a indicated subset
template <typename T, typename W = T>
class SuperMax3D final : public SuperF3D<T,W> {
private:
  FunctorPtr<SuperF3D<T,W>>        _f;
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
public:
  /// Constructor for determining the maximum of f on a indicated subset
  /**
   * \param f          functor of which the maximum is to be determined
   * \param indicatorF indicator describing the subset on which to evaluate f
   **/
  SuperMax3D(FunctorPtr<SuperF3D<T,W>>&&        f,
             FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
  /// Constructor for determining the maximum of f on a given material
  /**
   * \param f             functor of which the maximum is to be determined
   * \param superGeometry super geometry for constructing material indicator
   * \param material      number of the relevant material
   **/
  SuperMax3D(FunctorPtr<SuperF3D<T,W>>&& f,
             SuperGeometry3D<T>& superGeometry,
             const int material);

  bool operator() (W output[], const int input[]) override;
};


}

#endif
