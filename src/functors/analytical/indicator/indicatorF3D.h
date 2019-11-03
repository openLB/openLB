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

#ifndef INDICATOR_F_3D_H
#define INDICATOR_F_3D_H

#include<vector>
#include "indicatorBaseF3D.h"
#include "io/xmlReader.h"


/** \file
 * This file contains indicator functions. These return 1 if the given
 * coordinates are inside, and 0 if they are outside of the defined set.
 * Implemented are :
 - Sphere
 - Cylinder
 - Cone
 - Pipe (not yet)
 - Cube (not yet)
 - Cuboid
 - Circle

 * The smoothIndicator functors return values in [0,1]. In particular there is
 * an epsilon enclosure of the set, wherein the return values are smooth and do
 * not jump from 0 to 1.

 Boolean operators allow to create unions and intersections. They can be used
 for example for initialization of a SuperGeometry.
*/

namespace olb {


template <typename S>
class IndicatorTranslate3D : public IndicatorF3D<S> {
private:
  std::array<S,3> _translate;
  IndicatorF3D<S>& _indicator;
public:
  IndicatorTranslate3D(std::array<S,3> translate, IndicatorF3D<S>& indicator);
  bool operator() (bool output[], const S input[]) override;
};


/// indicator function for a 3D circle
template <typename S>
class IndicatorCircle3D : public IndicatorF3D<S> {
private:
  Vector<S,3> _center;
  Vector<S,3> _normal;
  S _radius2;
public:
  IndicatorCircle3D(Vector<S,3> center, Vector<S,3> normal, S radius);
  IndicatorCircle3D(S center0, S center1, S center2, S normal0, S normal1,
                    S normal2, S radius);
  bool operator() (bool output[], const S input[]) override;
  Vector<S,3> const& getCenter() const;
  Vector<S,3> const& getNormal() const;
  S getRadius() const;
  //virtual bool distance(S& distance, Vector<S,3> origin, Vector<S,3> direction, int iC=-1);
};



/// indicator function for a 3D-sphere
template <typename S>
class IndicatorSphere3D : public IndicatorF3D<S> {
private:
  Vector<S,3> _center;
  S _radius2;
public:
  IndicatorSphere3D(Vector<S,3> center, S radius);
  IndicatorSphere3D(const IndicatorSphere3D&);
  bool operator() (bool output[], const S input[]) override;
  bool distance(S& distance, const Vector<S,3>& origin,
                const Vector<S,3>& direction, int iC=-1) override;
  bool distance(S& distance, const Vector<S,3>& origin);
};

/// indicator function for a layer
template <typename S>
class IndicatorLayer3D : public IndicatorF3D<S> {
private:
  IndicatorF3D<S>& _indicatorF;
  S _layerSize;
public:
  IndicatorLayer3D(IndicatorF3D<S>& indicatorF, S layerSize);
  bool operator() (bool output[], const S input[]) override;
};

/// indicator function for the internal part of an input indicator
template <typename S>
class IndicatorInternal3D : public IndicatorF3D<S> {
private:
  IndicatorF3D<S>& _indicatorF;
  S _layerSize;
public:
  IndicatorInternal3D(IndicatorF3D<S>& indicatorF, S layerSize);
  bool operator() (bool output[], const S input[]) override;
};

/// indicator function for a 3d-cylinder
template <typename S>
class IndicatorCylinder3D : public IndicatorF3D<S> {
private:
  Vector<S,3> _center1;
  Vector<S,3> _center2;
  Vector<S,3> _I;
  Vector<S,3> _J;
  Vector<S,3> _K;
  S _length;
  S _radius2;
  void init();
public:
  IndicatorCylinder3D(Vector<S,3> center1, Vector<S,3> center2, S radius);
  // eps is length of cylinder
  IndicatorCylinder3D(Vector<S,3> center1, Vector<S,3> normal, S radius, S eps);
  IndicatorCylinder3D(IndicatorCircle3D<S> const& circleF, S eps);
  bool operator() (bool output[], const S input[]) override;
  Vector<S,3> const& getCenter1() const;
  Vector<S,3> const& getCenter2() const;
  S getRadius() const;
};

/// indicator function for a 3d frustum
template <typename S>
class IndicatorCone3D : public IndicatorF3D<S> {
private:
  Vector<S,3> _center1;
  Vector<S,3> _center2;
  Vector<S,3> _I;
  Vector<S,3> _J;
  Vector<S,3> _K;
  S _length;
  S _radius1;
  S _radius2; // The 2nd radius is optional: if not defined, _center2 is the vertex of the cone
public:
  IndicatorCone3D(Vector<S,3> center1, Vector<S,3> center2, S radius1, S radius2=0);
  bool operator() (bool output[], const S input[]) override;
};

/** indicator function for a 3d-cuboid, parallel to the planes x=0, y=0, z=0.
 * \param extend must have only positive elements
 * \param xLength must be positive
 */
template <typename S>
class IndicatorCuboid3D : public IndicatorF3D<S> {
private:
  const Vector<S,3> _center;
  const S _xLength;
  const S _yLength;
  const S _zLength;
public:
  /// constructs an cuboid with x axis dimension 0 to extend[0], ...
  IndicatorCuboid3D(Vector<S,3> extend, Vector<S,3> origin);
  /// constructs an cuboid with x axis dimension -xlength/2 to xlength/2
  IndicatorCuboid3D(S xlength, S ylength, S zlength, Vector<S,3> center);
  /// returns true if input is inside, otherwise false
  bool operator() (bool output[], const S input[]) override;
};


/** \brief indicator function for a 3d-cuboid, turned by an angle theta around an axis of rotation
 *
 * The cuboid is turned along the axis going through the point "centerRotation" and being orthogonal to either x=0 or y=0 or z=0 respectively
 * \param theta angle (rad) by which the cuboid is turned
 * \param plane the axis of rotation is orthogonal to the plane i=0, where 0 corresponds to the plane x=0, 1 to the plane y=0, and 2 to z=0
 * \param centerRotation vector of coordinates which defines, with "plane", the axis of rotation; centerRotation is a point on th axis of rotation
 *
 */
template <typename S>
class IndicatorCuboidRotate3D : public IndicatorCuboid3D<S> {
private:
    S _theta;
    int _plane;
    Vector<S,3> _centerRotation;
public:
    // constructs an cuboid turned by some angle theta around a given center of rotation
    IndicatorCuboidRotate3D(Vector<S,3> extend, Vector<S,3> origin, S theta, int plane, Vector<S,3> centerRotation);
    IndicatorCuboidRotate3D(S xlength, S ylength, S zlength, Vector<S,3> origin, S theta, int plane, Vector<S,3> centerRotation);
    bool operator() (bool output[], const S input[]);
};

/////////creatorFunctions//////////////////////
// creator function for geometric primitives
template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorCircle3D(XMLreader const& params, bool verbose=false);

template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorSphere3D(XMLreader const& params, bool verbose=false);

template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorCylinder3D(XMLreader const& params, bool verbose=false);

template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorCone3D(XMLreader const& params, bool verbose=false);

template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorCuboid3D(XMLreader const& params, bool verbose=false);

// arithmetic creator functions
template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorUnion3D(XMLreader const& params, bool verbose=false);

template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorWithout3D(XMLreader const& params, bool verbose=false);

template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorIntersection3D(XMLreader const&params, bool verbose=false);

// godfather
template <typename S>
std::shared_ptr<IndicatorF3D<S>> createIndicatorF3D(XMLreader const& params, bool verbose=false);



}

#endif

