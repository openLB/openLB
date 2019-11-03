/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause, Albert Mink
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

#ifndef ANALYTICAL_BASE_F_HH
#define ANALYTICAL_BASE_F_HH

#include "analyticalBaseF.h"

namespace olb {


template <typename T, typename S>
AnalyticalF1D<T,S>::AnalyticalF1D(int n) : GenericF<T,S>(n,1) { }

template <typename T, typename S>
AnalyticalF2D<T,S>::AnalyticalF2D(int n) : GenericF<T,S>(n,2) { }

template <typename T, typename S>
AnalyticalF3D<T,S>::AnalyticalF3D(int n) : GenericF<T,S>(n,3) { }


// identity to "store results"
template <typename T, typename S>
AnalyticalIdentity1D<T,S>::AnalyticalIdentity1D(AnalyticalF1D<T,S>& f)
  : AnalyticalF1D<T,S>(f.getTargetDim()), _f(f)
{
  this->getName() = _f.getName();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename S>
bool AnalyticalIdentity1D<T,S>::operator()(T output[], const S input[])
{
  _f(output,input);
  return true;
}


// identity to "store results"
template <typename T, typename S>
AnalyticalIdentity2D<T,S>::AnalyticalIdentity2D(AnalyticalF2D<T,S>& f)
  : AnalyticalF2D<T,S>(f.getTargetDim()), _f(f)
{
  this->getName() = _f.getName();
  // pass through the shared_ptr from _f, e.g. an arithemticClass, to the identity
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename S>
bool AnalyticalIdentity2D<T,S>::operator()(T output[], const S input[])
{
  _f(output,input);
  return true;
}


// identity to "store results"
template <typename T, typename S>
AnalyticalIdentity3D<T,S>::AnalyticalIdentity3D(AnalyticalF3D<T,S>& f)
  : AnalyticalF3D<T,S>(f.getTargetDim()), _f(f)
{
  this->getName() = _f.getName();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename S>
bool AnalyticalIdentity3D<T,S>::operator()(T output[], const S input[])
{
  _f(output,input);
  return true;
}



template <typename T, typename S>
AnalyticalFfromIndicatorF3D<T, S>::AnalyticalFfromIndicatorF3D(IndicatorF3D<T>& indicatorF)
  : AnalyticalF3D<T,S>(1), _indicatorF(indicatorF)
{
  this->getName() = "IndicatorFfrom" + _indicatorF.getName();
}

template <typename T, typename S>
bool AnalyticalFfromIndicatorF3D<T, S>::operator() (T output[], const S input[])
{
  bool tmp = false;
  _indicatorF(&tmp, input);
  if ( tmp ) {
    output[0] = T(1);
  } else {
    output[0] = T(0);
  }
  return tmp;
}



} // end namespace olb

#endif
