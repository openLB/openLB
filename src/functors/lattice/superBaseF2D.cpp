/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#include "superBaseF2D.h"
#include "superBaseF2D.hh"
#include "superCalcF2D.hh"
#include "dynamics/latticeDescriptors.h"

namespace olb {

template class SuperF2D<double,double>;
template class SuperF2D<double,int>;
template class SuperF2D<double,bool>;

template class SuperIdentity2D<double,double>;
template class SuperIdentity2D<double,int>;
template class SuperIdentity2D<double,bool>;

template class SuperLatticeF2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysF2D<double,descriptors::D2Q9<>>;

}

