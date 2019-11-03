/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Mathias J. Krause, Albert Mink
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

#ifndef SMOOTH_INDICATOR_F_3D_H
#define SMOOTH_INDICATOR_F_3D_H

#include <vector>

#include "smoothIndicatorBaseF3D.h"
#include "io/xmlReader.h"

#include "core/blockData3D.h"
#include "core/unitConverter.h"
#include "functors/analytical/indicator/indicatorBaseF2D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"

namespace olb {

/// implements a smooth particle cuboid in 3D with an _epsilon sector. 
template <typename T, typename S, bool HLBM=false>
class SmoothIndicatorCuboid3D final: public SmoothIndicatorF3D<T, S, HLBM> {
private:
  S _xLength;
  S _yLength;
  S _zLength;
public:
  SmoothIndicatorCuboid3D(Vector<S,3> center, S xLength, S yLength, S zLength, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.), S density=0, Vector<S,3> vel = Vector<S,3> (0.,0.,0.));
  bool operator()(T output[],const S input[]) override;
};

/// implements a smooth sphere in 3D with an _epsilon sector 
template <typename T, typename S, bool HLBM=false>
class SmoothIndicatorSphere3D final: public SmoothIndicatorF3D<T, S, HLBM> {
private:
  S _radius;
public:
  SmoothIndicatorSphere3D(Vector<S, 3> center, S radius, S epsilon, S density=0, Vector<S,3> vel = Vector<S,3> (0.,0.,0.));
  bool operator()(T output[], const S input[]) override;
};

/// implements a smooth particle cylinder in 3D with an _epsilon sector.
template <typename T, typename S, bool HLBM=false>
class SmoothIndicatorCylinder3D final: public SmoothIndicatorF3D<T, S, HLBM> {
private:
  S _radius;
  S _length;
  void initIndicatorCylinder3D(Vector<S,3> normal, Vector<S,3> theta, S density, Vector<S,3> vel);
public:
  SmoothIndicatorCylinder3D(Vector<S,3> pointA, Vector<S,3> pointB, S radius, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.), S density=0, Vector<S,3> vel = Vector<S,3> (0.,0.,0.));
  SmoothIndicatorCylinder3D(Vector<S,3> center, Vector<S,3> normal, S radius, S length, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.), S density=0, Vector<S,3> vel = Vector<S,3> (0.,0.,0.));
  bool operator()(T output[], const S input[]) override;
};

/// implements a smooth particle cone in 3D with an _epsilon sector
template <typename T, typename S, bool HLBM=false>
class SmoothIndicatorCone3D : public SmoothIndicatorF3D<T, S, HLBM> {
private:
  S _length;
  S _radiusA;
  S _radiusB;
  void initIndicatorCone3D(Vector<S,3> normal, Vector<S,3> theta, S density, Vector<S,3> vel);
public:
  SmoothIndicatorCone3D(Vector<S,3> pointA, Vector<S,3> pointB,
                        S radiusA, S radiusB, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.), S density=0, 
                        Vector<S,3> vel = Vector<S,3> (0.,0.,0.));
  SmoothIndicatorCone3D(Vector<S,3> center, Vector<S,3> normal, S lenght,
                        S radiusA, S radiusB, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.), 
                        S density=0, Vector<S,3> vel = Vector<S,3> (0.,0.,0.));
  bool operator() (T output[], const S input[]) override;
};

}

#endif

