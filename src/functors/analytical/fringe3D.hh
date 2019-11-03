/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause
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

#ifndef FRINGE_3D_HH
#define FRINGE_3D_HH

#include<vector>
#include<cmath>

#include "fringe3D.h"

namespace olb {

template <typename T, typename S>
Fringe3D<T,S>::Fringe3D(AnalyticalF3D<T,S>& wantedVelocity, T start, T end, int direction, T lambdaMax, T rise, T fall)
  : AnalyticalF3D<T,S>(12), _wantedVelocity(wantedVelocity), _start(start),
    _end(end), _direction(direction), _lambdaMax(lambdaMax), _rise(rise),
    _fall(fall)
{
  this->getName() = "fringe";
}

template <typename T, typename S>
bool Fringe3D<T,S>::operator()(T output[12], const S input[3])
{

  output[4] = T();
  output[5] = T();
  output[6] = T();
  output[8] = T();
  output[9] = T();
  output[10] = T();
  if ( _wantedVelocity(output, input) ) {
    T scaleLambda = lambda(input[_direction]);
    output[0] = -scaleLambda*output[0];
    output[1] = -scaleLambda*output[1];
    output[2] = -scaleLambda*output[2];
    output[3] = scaleLambda;
    output[7] = scaleLambda;
    output[11] = scaleLambda;
    return true;
  } else {
    output[0] = T();
    output[1] = T();
    output[2] = T();
    output[3] = T();
    output[7] = T();
    output[11] = T();
    return false;
  }
}

template <typename T, typename S>
T Fringe3D<T,S>::s(const T& x) const
{
  if ( x < std::numeric_limits<T>::epsilon() ) {
    return 0.;
  } else if (x > T(1)-std::numeric_limits<T>::epsilon() ) {
    return T(1);
  } else {
    return T(1)/(T(1) + exp( T(1)/(x - T(1) ) + T(1)/x ) );
  }
}

template <typename T, typename S>
T Fringe3D<T,S>::lambda(const T& x) const
{
  return _lambdaMax*(s((x - _start)/_rise*(_end - _start) ) - s((x - _end)/_fall*(_end - _start) + T(1.)) );
}

} // end namespace olb

#endif

