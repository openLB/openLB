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

#include "io/stlReader.h"
#include "smoothIndicatorF3D.h"
#include "smoothIndicatorF3D.hh"

namespace olb {

// stl indicator functors
template class STLreader<double>;

// smoothIndicator functors
template class SmoothIndicatorCuboid3D<double,double,false>;
template class SmoothIndicatorSphere3D<double,double,false>;
template class SmoothIndicatorCylinder3D<double,double,false>;
//template class SmoothIndicatorCone3D<double,double,false>;
//template class SmoothIndicatorCustom3D<double,double,descriptors::D3Q19<>,false>;

}
