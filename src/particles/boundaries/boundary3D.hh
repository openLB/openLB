/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Thomas Henn, Mathias J. Krause
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


#ifndef BOUNDARY3D_HH
#define BOUNDARY3D_HH

#include "boundary3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
Boundary3D<T, PARTICLETYPE>::Boundary3D(Boundary3D<T, PARTICLETYPE>& f) : clout(f.clout)
{
}

template<typename T, template<typename U> class PARTICLETYPE>
Boundary3D<T, PARTICLETYPE>::Boundary3D(const Boundary3D<T, PARTICLETYPE>& f) : clout(f.clout)
{
}

template<typename T, template<typename U> class PARTICLETYPE>
Boundary3D<T, PARTICLETYPE>::Boundary3D(): clout(std::cout,"Boundary3D")
{
}

}
#endif /* BOUNDARY3D_HH */
