/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Mathias J. Krause, Benjamin Förster
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

#ifndef INDICATOR_BASE_F_2D_HH
#define INDICATOR_BASE_F_2D_HH

#include <cmath>
#include "indicatorBaseF2D.h"
#include "math.h"

namespace olb {


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

template <typename S>
IndicatorF1D<S>::IndicatorF1D()
  : GenericF<bool,S>(1, 1)
{ }

template <typename S>
Vector<S,1>& IndicatorF1D<S>::getMin()
{
  return _myMin;
}

template <typename S>
Vector<S,1>& IndicatorF1D<S>::getMax()
{
  return _myMax;
}


template <typename S>
IndicatorF2D<S>::IndicatorF2D()
  : GenericF<bool,S>(1, 2)
{ }

template <typename S>
Vector<S,2>& IndicatorF2D<S>::getMin()
{
  return _myMin;
}

template <typename S>
Vector<S,2>& IndicatorF2D<S>::getMax()
{
  return _myMax;
}

template <typename S>
bool IndicatorF2D<S>::distance(S& distance, const Vector<S,2>& origin, const Vector<S,2>& direction, int iC)
{
  bool originValue;
  (*this)(&originValue, origin.data);
  Vector<S,2> currentPoint(origin);

  S precision = .0001;
  S pitch = 0.5;

  bool currentValue;
  (*this)(&currentValue, currentPoint.data);
  while (currentValue == originValue && isInsideBox(currentPoint)) {
    currentPoint += direction;
    (*this)(&currentValue, currentPoint.data);//changed until first point on the other side (inside/outside) is found
  }

  if (!isInsideBox(currentPoint) && !originValue) {
    return false;
  }

  while (pitch >= precision) {
    if (!isInsideBox(currentPoint) && originValue) {
      currentPoint -= pitch * direction;
      pitch /= 2.;
    } else {
      (*this)(&currentValue, currentPoint.data);
      if (currentValue == originValue) {
        currentPoint += pitch * direction;
        pitch /= 2.;
      } else {
        currentPoint-= pitch * direction;
        pitch /= 2.;
      }
    }
  }


  distance = (currentPoint - origin).norm();
  return true;
}

// given origin (inside) and direction first calculate distance to surface
// go -90� to direction using POS as origin and check if inside/outside
// if inside rotate outward, if outside rotate inward
// iterate until close enough to surface, then store point1
// repeat for 90� and store point2
// use point1 and point2 on surface to calculate normal
template <typename S>
bool IndicatorF2D<S>::normal(Vector<S,2>& normal, const Vector<S,2>& origin, const Vector<S,2>& direction, int iC)
{
  //OstreamManager clout(std::cout,"normal");
  //clout << "Calculating IndicatorF2D Normal " << endl;
  bool originValue;
  (*this)(&originValue, origin.data);
  Vector<S,2> currentPoint(origin);

  S precision = .0001;

  S dist;
  distance(dist, origin, direction, iC);


  Vector<S,2> POS(origin + dist*direction*(1/const_cast<Vector<S,2>&> (direction).norm())); //Point on Surface

  Vector<S,2> point1;
  Vector<S,2> point2;

  bool currentValue;

  for (int n: {
         -90,90
       }) { //std::range<int> n = {-90, 90};
    S rotate(n);
    S pitch(rotate/2.);
    while (std::abs(pitch) >= precision) {
      S theta(rotate*M_PI/180.);

      Vector<S,2> vec(std::cos(theta)*direction[0]+std::sin(theta)*direction[1],-std::sin(theta)*direction[0]+std::cos(theta)*direction[1]);
      currentPoint = POS + vec;
      (*this)(&currentValue, currentPoint.data);

      if (currentValue == originValue) {
        rotate -= pitch;
      }  else {
        rotate += pitch;
      }
      pitch /= 2.;
    }

    if (n == -90) {
      point1 = currentPoint;
    } else if (n == 90) {
      point2 = currentPoint;
    }
  }
  // Calculate Normal from point1 and point2
  normal = Vector<S,2>((point2[1] - point1[1]), (-1)*(point2[0] - point1[0]));


  //S dist;
  //Vector<S,2> dist;
  //distance(dist, origin, direction, iC);
  //normal = Vector<S,2>(dist);
  return true;
}

template <typename S>
bool IndicatorF2D<S>::isInsideBox(Vector<S,2> point)
{
  return point >= _myMin && point <= _myMax;
}

template <typename S>
bool IndicatorF2D<S>::operator() (const S input[])
{
  bool output;
  this->operator()(&output, input);
  return output;
}

template <typename S>
IndicatorIdentity2D<S>::IndicatorIdentity2D(std::shared_ptr<IndicatorF2D<S>> f)
  : _f(f)
{
  this->_myMin = _f->getMin();
  this->_myMax = _f->getMax();
}

template <typename S>
bool IndicatorIdentity2D<S>::operator() (bool output[1], const S input[2])
{
  return (_f)->operator()(output, input);
}

} // namespace olb

#endif
