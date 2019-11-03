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

#ifndef FRINGE_2D_HH
#define FRINGE_2D_HH

#include<vector>
#include<cmath>

#include "fringe2D.h"

namespace olb {

template <typename T, typename S>
Fringe2D<T,S>::Fringe2D(AnalyticalF2D<T,S>& wantedVelocity, T start, T end, int direction, T lambdaMax, T rise, T fall)
  : AnalyticalF2D<T,S>(6), _wantedVelocity(wantedVelocity), _start(start),
    _end(end), _direction(direction), _lambdaMax(lambdaMax), _rise(rise),
    _fall(fall)
{
  this->getName() = "fringe";
}

template <typename T, typename S>
bool Fringe2D<T,S>::operator()(T output[6], const S input[2])
{
  output[3] = T();
  output[4] = T();
  if ( _wantedVelocity(output, input) ) {
    T scaleLambda = lambda(input[_direction]);
    output[0] = -scaleLambda*output[0];
    output[1] = -scaleLambda*output[1];
    output[2] = scaleLambda;
    output[5] = scaleLambda;
    return true;
  } else {
    output[0] = T();
    output[1] = T();
    output[2] = T();
    output[5] = T();
    return false;
  }
}

template <typename T, typename S>
T Fringe2D<T,S>::s(const T& x) const
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
T Fringe2D<T,S>::lambda(const T& x) const
{
  return _lambdaMax*(s((x - _start)/_rise*(_end - _start) ) - s((x - _end)/_fall*(_end - _start) + T(1.)) );
}

} // end namespace olb

#endif

