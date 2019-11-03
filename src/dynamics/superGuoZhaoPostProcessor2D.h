/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016-2017 Davide Dapelo, Mathias J. Krause
 *  OpenLB e-mail contact: info@openlb.net
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
 * Class to define the external fields involved
 * in the porous modelling (i.e., porosity, porous conductivity
 * with the addition of a generic body force)
 * -- header file
 */

#ifndef SUPER_GUO_ZAO_POST_PROCESSOR_2D_H
#define SUPER_GUO_ZAO_POST_PROCESSOR_2D_H

#include "dynamics/guoZhaoLbHelpers.h"

namespace olb {

template<typename T, typename DESCRIPTOR, class dynamicsManager>
class SuperGuoZhaoInstantiator2D {
public:
  SuperGuoZhaoInstantiator2D (SuperLattice2D<T, DESCRIPTOR>& sLattice_);
  void definePorousFields(AnalyticalF2D<T,T>& epsilon_, AnalyticalF2D<T,T>& K_);
  void defineEpsilon(SuperGeometry2D<T>& sGeometry, int material, AnalyticalF2D<T,T>& epsilon);
  void defineK(UnitConverter<T,DESCRIPTOR> const& converter, SuperGeometry2D<T>& sGeometry, int material, AnalyticalF2D<T,T>& K);
  void defineNu(UnitConverter<T,DESCRIPTOR> const& converter, SuperGeometry2D<T>& sGeometry, int material);
  void defineBodyForce(UnitConverter<T,DESCRIPTOR> const& converter, SuperGeometry2D<T>& sGeometry, int material, AnalyticalF2D<T,T>& bodyForce);

private:
  SuperLattice2D<T, DESCRIPTOR>& sLattice;
};

}

#endif
