/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Cyril Masquelier, Jan Marquardt, Mathias J. Krause,
 *  Davide Dapelo
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

#ifndef INDICATOR_F_2D_H
#define INDICATOR_F_2D_H

#include<vector>
#include "indicatorBaseF2D.h"
#include "io/xmlReader.h"

#include "core/blockData2D.h"
#include "core/unitConverter.h"
#include "indicatorBaseF3D.h"

/** \file
 * This file contains indicator functions. These return 1 if the given
 * coordinates are inside, and 0 if they are outside of the defined set.
 * Implemented are :
 - Cuboid
 - Circle

 * The smoothIndicator functors return values in [0,1]. In particular there is
 * an epsilon enclosure of the set, wherein the return values are smooth and do
 * not jump from 0 to 1.

 Boolean operators allow to create unions and intersections. They can be used
 for example for initialization of a SuperGeometry.
*/

namespace olb {

/// indicator function for a 2D-cuboid, parallel to the planes x=0, y=0;
/// theta rotates cuboid around its center, theta in radian measure
template <typename S>
class IndicatorF2DfromIndicatorF3D : public IndicatorF2D<S> {
private:
  IndicatorF3D<S>& _indicator3D;
public:
  IndicatorF2DfromIndicatorF3D(IndicatorF3D<S>& indicator3D);

  bool operator() (bool output[], const S input[]) override;
};

/// indicator function for a 2D-cuboid, parallel to the planes x=0, y=0;
/// theta rotates cuboid around its center, theta in radian measure
template <typename S>
class IndicatorCuboid2D : public IndicatorF2D<S> {
private:
  Vector<S,2> _center;
  S _xLength;
  S _yLength;
  S _theta;
public:
  /// constructs an cuboid with x axis dimension 0 to extend[0], ...
  IndicatorCuboid2D(Vector<S,2> extend, Vector<S,2> origin, S theta=0);
  /// constructs an cuboid with x axis dimension -xlength/2 to xlength/2
  IndicatorCuboid2D(S xlength, S ylength, Vector<S,2> center= {S(), S()}, S theta=0);
  /// returns true if input is inside, otherwise false
  bool operator() (bool output[], const S input[]) override;
};


/// indicator function for a 2D circle
template <typename S>
class IndicatorCircle2D : public IndicatorF2D<S> {
private:
  Vector<S,2> _center;
  S _radius2;
public:
  IndicatorCircle2D(Vector<S,2> center, S radius);
  bool operator() (bool output[], const S input[]) override;
  bool distance(S& distance, const Vector<S,2>& origin, const Vector<S,2>& direction,  int iC=-1) override;
  bool normal(Vector<S,2>& normal, const Vector<S,2>& origin, const Vector<S,2>& direction, int iC=-1) override;
};

/// Indicator function creating an layer around an input indicator (for positive \p layerSize) or
/// reducing the input indicator by a layer (for negative \p layerSize).
/// \param[in] indicatorF - some indicator (e.g. IndicatorCircle2D) 
/// \param[in] layerSize - size of the layer (can be negative) [physical units]
template <typename S>
class IndicatorLayer2D : public IndicatorF2D<S> {
private:
  IndicatorF2D<S>& _indicatorF;
  S _layerSize;
  bool _isPositive;
public:
  IndicatorLayer2D(IndicatorF2D<S>& indicatorF, S layerSize);
  bool operator() (bool output[], const S input[]) override;
};


/////////creatorFunctions//////////////////////
template <typename S>
IndicatorCuboid2D<S>* createIndicatorCuboid2D(XMLreader const& params, bool verbose=false);

template <typename S>
IndicatorCircle2D<S>* createIndicatorCircle2D(XMLreader const& params, bool verbose=false);


}

#endif
