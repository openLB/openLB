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

#include "analyticalBaseF.h"
#include "analyticalBaseF.hh"
#include "analyticCalcF.hh"

namespace olb {

// 2nd level classes
template class AnalyticalF1D<double,int>;
template class AnalyticalF1D<double,double>;

template class AnalyticalIdentity1D<double,int>;
template class AnalyticalIdentity1D<double,double>;


template class AnalyticalF2D<double,int>;
template class AnalyticalF2D<double,double>;

template class AnalyticalIdentity2D<double,int>;
template class AnalyticalIdentity2D<double,double>;


template class AnalyticalF3D<double,int>;
template class AnalyticalF3D<double,double>;
template class AnalyticalF3D<bool,double>; // indicatorF

template class AnalyticalIdentity3D<double,int>;
template class AnalyticalIdentity3D<double,double>;

}

