/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2011-2013 Lukas Baron, Tim Dornieden, Mathias J. Krause
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

#include "analyticalF.h"
#include "analyticalF.hh"

namespace olb {


template class AnalyticalConst1D<double,int>;
template class AnalyticalConst1D<double,double>;

template class AnalyticalLinear1D<double,int>;
template class AnalyticalLinear1D<double,double>;

template class AnalyticalRandom1D<double,int>;
template class AnalyticalRandom1D<double,double>;

template class AnalyticalSquare1D<double,int>;
template class AnalyticalSquare1D<double,double>;

template class PolynomialStartScale<double,int>;
template class PolynomialStartScale<double,double>;

template class SinusStartScale<double,int>;
template class SinusStartScale<double,double>;

template class AnalyticalDiffFD1D<double>;


template class AnalyticalComposed2D<double,int>;
template class AnalyticalComposed2D<double,double>;

template class AnalyticalConst2D<double,int>;
template class AnalyticalConst2D<double,double>;

template class AnalyticalLinear2D<double,int>;
template class AnalyticalLinear2D<double,double>;

template class AnalyticalRandom2D<double,int>;
template class AnalyticalRandom2D<double,double>;


template class AnalyticalComposed3D<double,int>;
template class AnalyticalComposed3D<double,double>;

template class AnalyticalConst3D<double,int>;
template class AnalyticalConst3D<double,double>;

template class AnalyticalLinear3D<double,int>;
template class AnalyticalLinear3D<double,double>;

template class AnalyticalRandom3D<double,int>;
template class AnalyticalRandom3D<double,double>;

template class AnalyticalScaled3D<double,int>;
template class AnalyticalScaled3D<double,double>;

}

