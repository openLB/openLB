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
 * -- generic file
 */

#ifndef SUPER_GUO_ZAO_POST_PROCESSOR_2D_HH
#define SUPER_GUO_ZAO_POST_PROCESSOR_2D_HH

#include "dynamics/guoZhaoLbHelpers.h"

namespace olb {

template<typename T, typename DESCRIPTOR, class dynamicsManager>
SuperGuoZhaoInstantiator2D<T, DESCRIPTOR, dynamicsManager>::SuperGuoZhaoInstantiator2D (
  SuperLattice2D<T, DESCRIPTOR>& sLattice_) :
  sLattice(sLattice_)
{}

template<typename T, typename DESCRIPTOR, class dynamicsManager>
void SuperGuoZhaoInstantiator2D<T, DESCRIPTOR, dynamicsManager>::definePorousFields (
  AnalyticalF2D<T,T>& epsilon_, AnalyticalF2D<T,T>& K_)
{

}

template<typename T, typename DESCRIPTOR, class dynamicsManager>
void SuperGuoZhaoInstantiator2D<T, DESCRIPTOR, dynamicsManager>::defineEpsilon (
  SuperGeometry2D<T>& sGeometry, int material, AnalyticalF2D<T,T>& epsilon)
{

  sLattice.template defineField<descriptors::EPSILON>(sGeometry, material, epsilon);
}

template<typename T, typename DESCRIPTOR, class dynamicsManager>
void SuperGuoZhaoInstantiator2D<T, DESCRIPTOR, dynamicsManager>::defineK (
  UnitConverter<T,DESCRIPTOR> const& converter, SuperGeometry2D<T>& sGeometry, int material, AnalyticalF2D<T,T>& K)
{

  AnalyticalConst2D<T,T> normFactor(converter.getConversionFactorLength()*converter.getConversionFactorLength());
  AnalyticalIdentity2D<T,T> KLb(K / normFactor);
  sLattice.template defineField<descriptors::K>(sGeometry, material, KLb);
}

template<typename T, typename DESCRIPTOR, class dynamicsManager>
void SuperGuoZhaoInstantiator2D<T, DESCRIPTOR, dynamicsManager>::defineNu (
  UnitConverter<T,DESCRIPTOR> const& converter, SuperGeometry2D<T>& sGeometry, int material)
{

  AnalyticalConst2D<T,T> nu(converter.getLatticeViscosity());
  sLattice.template defineField<descriptors::NU>(sGeometry, material, nu);
}

template<typename T, typename DESCRIPTOR, class dynamicsManager>
void SuperGuoZhaoInstantiator2D<T, DESCRIPTOR, dynamicsManager>::defineBodyForce (
  UnitConverter<T,DESCRIPTOR> const& converter, SuperGeometry2D<T>& sGeometry, int material, AnalyticalF2D<T,T>& BodyForce)
{

  std::vector<T> normFactorValue ( 2,
                                   converter.getConversionFactorLength() / (converter.getConversionFactorTime()*converter.getConversionFactorTime()) );
  AnalyticalConst2D<T,T> normFactor(normFactorValue);
  AnalyticalIdentity2D<T,T> BodyForceLb(BodyForce / normFactor);
  sLattice.template defineField<descriptors::BODY_FORCE>(sGeometry, material, BodyForceLb);
}

}

#endif
