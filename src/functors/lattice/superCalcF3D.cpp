/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Tim Dornieden,
 *  Mathias J. Krause, Albert Mink, Adrian Kummerlaender
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

#include "superCalcF3D.h"
#include "superCalcF3D.hh"
#include "utilities/functorPtr.hh"

namespace olb {

template class SuperCalcF3D<double,double,util::plus>;
template class SuperCalcF3D<double,int,util::plus>;
template class SuperCalcF3D<double,bool,util::plus>;

template class SuperCalcF3D<double,double,util::minus>;
template class SuperCalcF3D<double,int,util::minus>;
template class SuperCalcF3D<double,bool,util::minus>;

template class SuperCalcF3D<double,double,util::multiplies>;
template class SuperCalcF3D<double,int,util::multiplies>;
template class SuperCalcF3D<double,bool,util::multiplies>;

template class SuperCalcF3D<double,double,util::divides>;
template class SuperCalcF3D<double,int,util::divides>;
template class SuperCalcF3D<double,bool,util::divides>;

template class SuperCalcF3D<double,double,util::power>;
template class SuperCalcF3D<double,int,util::power>;

} // end namespace olb
