/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause, Benjamin Förster
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

/** \file
 * Reader Class for VTI Data -- template instantiation.
 */
#include "vtiReader.h"
#include "vtiReader.hh"

namespace olb {

template class BaseVTIreader<double>;
template class BaseVTIreader<int>;

template class BaseVTIreader3D<double,double>;
template class BaseVTIreader3D<double,int>;

template class BlockVTIreader3D<double,double>;
template class BlockVTIreader3D<double,int>;

template class SuperVTIreader3D<double,double>;
template class SuperVTIreader3D<double,int>;

}
