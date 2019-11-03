/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Peter Weisbrod
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

#include "dynamics/interactionPotential.h"
#include "dynamics/interactionPotential.hh"

namespace olb {

// established -- original for both single- and multicomponent flow

template class ShanChen93<double,int>;
template class ShanChen93<double,double>;

// established -- only multicomponent flow

template class PsiEqualsRho<double,int>;
template class PsiEqualsRho<double,double>;

// established -- only singlecomponent flow

template class ShanChen94<double,int>;
template class ShanChen94<double,double>;

template class PengRobinson<double,int>;
template class PengRobinson<double,double>;

template class CarnahanStarling<double,int>;
template class CarnahanStarling<double,double>;

// under development -- for singlecomponent flow

template class Krause<double,int>;
template class Krause<double,double>;

template class WeisbrodKrause<double,int>;
template class WeisbrodKrause<double,double>;

template class Normal<double,int>;
template class Normal<double,double>;

}

