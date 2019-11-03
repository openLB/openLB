/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013-2018 Mathias Krause, Albert Mink, Adrian Kummerlaender
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

#ifndef SUPER_GEOMETRY_FACES_3D_H
#define SUPER_GEOMETRY_FACES_3D_H

#include "superBaseF3D.h"
#include "blockGeometryFaces3D.h"
#include "utilities/functorPtr.h"

namespace olb {


/// Accumulates the discrete surface of indicated cells facing unit directions
/// and returns the individual as well as the total surface in phys units.
/**
 * Unit directions: (1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1)
 **/
template <typename T>
class SuperGeometryFaces3D final : public SuperF3D<T> {
private:
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
public:
  /// Constructor accepting solid cell indicator and custom lattice length
  SuperGeometryFaces3D(FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
                       T _latticeL);
  /// Constructor accepting single solid cell material and custom lattice length
  SuperGeometryFaces3D(SuperGeometry3D<T>& superGeometry, const int material,
                       T _latticeL);

  /// Constructor accepting solid cell indicator and a unit converter to query lattice length
  template<typename DESCRIPTOR>
  SuperGeometryFaces3D(FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
                       const UnitConverter<T,DESCRIPTOR>& converter)
    : SuperGeometryFaces3D(std::forward<decltype(indicatorF)>(indicatorF),
                           converter.getConversionFactorLength()) { };
  /// Constructor accepting single solid cell material and a unit converter to query lattice length
  template<typename DESCRIPTOR>
  SuperGeometryFaces3D(SuperGeometry3D<T>& superGeometry, const int material,
                       const UnitConverter<T,DESCRIPTOR>& converter)
    : SuperGeometryFaces3D(superGeometry.getMaterialIndicator(material), converter) { };

  bool operator() (T output[], const int input[]) override;

};


}

#endif
