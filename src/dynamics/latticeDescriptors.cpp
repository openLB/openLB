/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
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

#include "latticeDescriptors.h"
#include "mrtLatticeDescriptors.h"

namespace olb {

namespace descriptors {

template struct D2Q9<>;
template struct D2Q5<>;
template struct D2Q5<VELOCITY>;
template struct D3Q7<>;
template struct D3Q7<VELOCITY>;
template struct D3Q13<>;
template struct D3Q13<FORCE>;
template struct D3Q15<>;
template struct D3Q15<FORCE>;
template struct D3Q19<>;
template struct D3Q27<>;
template struct D3Q27<FORCE>;

template struct D2Q5<tag::MRT>;
template struct D2Q5<tag::MRT,VELOCITY>;
template struct D3Q7<tag::MRT>;
template struct D3Q7<tag::MRT,VELOCITY>;

}

}
