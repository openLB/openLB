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

#include "analyticCalcF.h"
#include "analyticCalcF.hh"

namespace olb {


// arithmetic helper class for analytical 1d functors

template class AnalyticCalcF1D<double,int,util::plus>;
template class AnalyticCalcF1D<double,double,util::plus>;

template class AnalyticCalcF1D<double,int,util::minus>;
template class AnalyticCalcF1D<double,double,util::minus>;

template class AnalyticCalcF1D<double,int,util::multiplies>;
template class AnalyticCalcF1D<double,double,util::multiplies>;

template class AnalyticCalcF1D<double,int,util::divides>;
template class AnalyticCalcF1D<double,double,util::divides>;

// arithmetic helper class for analytical 2d functors

template class AnalyticCalcF2D<double,int,util::plus>;
template class AnalyticCalcF2D<double,double,util::plus>;

template class AnalyticCalcF2D<double,int,util::minus>;
template class AnalyticCalcF2D<double,double,util::minus>;

template class AnalyticCalcF2D<double,int,util::multiplies>;
template class AnalyticCalcF2D<double,double,util::multiplies>;

template class AnalyticCalcF2D<double,int,util::divides>;
template class AnalyticCalcF2D<double,double,util::divides>;

// arithmetic helper class for analytical 3d functors

template class AnalyticCalcF3D<double,int,util::plus>;
template class AnalyticCalcF3D<double,double,util::plus>;

template class AnalyticCalcF3D<double,int,util::minus>;
template class AnalyticCalcF3D<double,double,util::minus>;

template class AnalyticCalcF3D<double,int,util::multiplies>;
template class AnalyticCalcF3D<double,double,util::multiplies>;

template class AnalyticCalcF3D<double,int,util::divides>;
template class AnalyticCalcF3D<double,double,util::divides>;

}

