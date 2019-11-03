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

#ifndef SMOOTH_INDICATOR_CALC_F_2D_H
#define SMOOTH_INDICATOR_CALC_F_2D_H

#include "smoothIndicatorBaseF2D.h"

namespace olb {


//////////////////////////////// IndicSmoothCalc2D ////////////////////////////////
/// arithmetic helper class for Indicator 2d functors
template <typename T, typename S>
class SmoothIndicCalc2D : public SmoothIndicatorF2D<T,S> {
protected:
  SmoothIndicatorF2D<T,S>& _f;
  SmoothIndicatorF2D<T,S>& _g;
public:
  SmoothIndicCalc2D(SmoothIndicatorF2D<T,S>& f, SmoothIndicatorF2D<T,S>& g);
};

/// addition functor acts as union
template <typename T, typename S>
class SmoothIndicPlus2D : public SmoothIndicCalc2D<T,S> {
public:
  SmoothIndicPlus2D(SmoothIndicatorF2D<T,S>& f, SmoothIndicatorF2D<T,S>& g);
  bool operator() (T output[], const S input[]) override;
};


} // end namespace olb

#endif
