/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#include "functorPtr.h"
#include "functorPtr.hh"

#include "functors/functors2D.h"
#include "functors/functors3D.h"

#include "dynamics/latticeDescriptors.h"

namespace olb {

template class FunctorPtr<AnalyticalF1D<double,double>>;
template class FunctorPtr<AnalyticalF1D<double,int>>;
template class FunctorPtr<AnalyticalF1D<bool,double>>;

template class FunctorPtr<AnalyticalF2D<double,double>>;
template class FunctorPtr<AnalyticalF2D<double,int>>;
template class FunctorPtr<AnalyticalF2D<bool,double>>;

template class FunctorPtr<AnalyticalF3D<double,double>>;
template class FunctorPtr<AnalyticalF3D<double,int>>;
template class FunctorPtr<AnalyticalF3D<bool,double>>;

template class FunctorPtr<SuperF2D<double>>;
template class FunctorPtr<SuperF2D<double,int>>;
template class FunctorPtr<SuperF2D<double,bool>>;
template class FunctorPtr<SuperLatticeF2D<double,descriptors::D2Q9<>>>;

template class FunctorPtr<SuperF3D<double>>;
template class FunctorPtr<SuperF3D<double,int>>;
template class FunctorPtr<SuperF3D<double,bool>>;
template class FunctorPtr<SuperLatticeF3D<double,descriptors::D3Q19<>>>;

template class FunctorPtr<SuperIndicatorF2D<double>>;
template class FunctorPtr<IndicatorF2D<double>>;
template class FunctorPtr<SmoothIndicatorF2D<double,double,false>>;

template class FunctorPtr<SuperIndicatorF3D<double>>;
template class FunctorPtr<IndicatorF3D<double>>;
template class FunctorPtr<SmoothIndicatorF3D<double,double,false>>;

}
