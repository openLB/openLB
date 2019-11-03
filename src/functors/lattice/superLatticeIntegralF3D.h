/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Benjamin Förster, Adrian Kummerlaender
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

#ifndef SUPER_LATTICE_INTEGRAL_F_3D_H
#define SUPER_LATTICE_INTEGRAL_F_3D_H

#include <vector>

#include "functors/genericF.h"
#include "blockLatticeIntegralF3D.h"
#include "superBaseF3D.h"
#include "indicator/superIndicatorBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "superLatticeLocalF3D.h"
#include "functors/analytical/interpolationF3D.h"
#include "functors/lattice/reductionF3D.h"
#include "integral/superIntegralF3D.h"
#include "core/superLattice3D.h"
#include "core/vector.h"
#include "io/ostreamManager.h"
#include "geometry/superGeometry3D.h"
#include "superGeometryFaces3D.h"
#include "utilities/functorPtr.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

////////////////////////////////////////////////////////////////////////////////
//////if globIC is not on the local processor, the returned vector is empty/////
////////////////////////////////////////////////////////////////////////////////


/// functor to get pointwise phys force acting on a indicated boundary on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysDrag3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  FunctorPtr<SuperIndicatorF3D<T>>              _indicatorF;
  SuperGeometryFaces3D<T>                       _facesF;
  SuperLatticePhysBoundaryForce3D<T,DESCRIPTOR> _pBoundForceF;
  SuperSum3D<T,T>                               _sumF;

  const T _factor;
public:
  SuperLatticePhysDrag3D(SuperLattice3D<T,DESCRIPTOR>&      sLattice,
                         FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
                         const UnitConverter<T,DESCRIPTOR>& converter);
  SuperLatticePhysDrag3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                         SuperGeometry3D<T>& superGeometry, const int material,
                         const UnitConverter<T,DESCRIPTOR>& converter);

  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys force acting on a indicated boundary on local lattice
/**
 *  see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
 */
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysCorrDrag3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  FunctorPtr<SuperIndicatorF3D<T>>                  _indicatorF;
  SuperGeometryFaces3D<T>                           _facesF;
  SuperLatticePhysCorrBoundaryForce3D<T,DESCRIPTOR> _pBoundForceF;
  SuperSum3D<T,T>                                   _sumF;

  const T _factor;
public:
  SuperLatticePhysCorrDrag3D(SuperLattice3D<T,DESCRIPTOR>&      sLattice,
                             FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
                             const UnitConverter<T,DESCRIPTOR>& converter);
  SuperLatticePhysCorrDrag3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                             SuperGeometry3D<T>& superGeometry, const int material,
                             const UnitConverter<T,DESCRIPTOR>& converter);

  bool operator() (T output[], const int input[]) override;
};


} // end namespace olb

#endif
