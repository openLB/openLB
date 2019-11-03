/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn
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


#ifndef FORCE_3D_HH
#define FORCE_3D_HH

#include "force3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
Force3D<T, PARTICLETYPE>::Force3D(Force3D<T, PARTICLETYPE>& f) :
  clout(f.clout)
{
}

template<typename T, template<typename U> class PARTICLETYPE>
Force3D<T, PARTICLETYPE>::Force3D(const Force3D<T, PARTICLETYPE>& f) :
  clout(f.clout)
{
}

template<typename T, template<typename U> class PARTICLETYPE>
Force3D<T, PARTICLETYPE>::Force3D():  clout(std::cout,"Force3D")
{
}

}
#endif /* FORCE3D_HH */
