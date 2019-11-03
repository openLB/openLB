/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause, Adrian Kummerlaender
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

#ifndef SUPER_LATTICE_INTEGRAL_F_2D_H
#define SUPER_LATTICE_INTEGRAL_F_2D_H

#include<vector>
#include<cmath>

#include "functors/genericF.h"
#include "superBaseF2D.h"
#include "superLatticeLocalF2D.h"
#include "indicator/superIndicatorBaseF2D.h"
#include "functors/analytical/interpolationF2D.h"
#include "core/superLattice2D.h"
#include "core/vector.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

////////////////////////////////////////////////////////////////////////////////
//////if globIC is not on the local processor, the returned vector is empty/////
////////////////////////////////////////////////////////////////////////////////


/// functor that returns the max in each component of all points of a certain material
template <typename T>
class SuperMax2D final : public SuperF2D<T> {
private:
  SuperF2D<T>& _f;
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperMax2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry,
             const int material);
  bool operator() (T output[], const int input[]) override;
};


/// functor that returns the min in each component of all points of a certain material
template <typename T>
class SuperMin2D final : public SuperF2D<T> {
private:
  SuperF2D<T>& _f;
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperMin2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry,
             const int material);
  bool operator() (T output[], const int input[]) override;
};

/// sums over all cells of a certain material number
template <typename T>
class SuperSum2D final : public SuperF2D<T> {
private:
  SuperF2D<T>& _f;
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperSum2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry,
             const int material);
  bool operator() (T output[], const int input[]) override;
};


template <typename T>
class SuperIntegral2D final : public SuperF2D<T> {
private:
  SuperF2D<T>& _f;
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperIntegral2D(SuperF2D<T>& f, SuperGeometry2D<T>& superGeometry,
                  const int material);
  bool operator() (T output[], const int input[]) override;
};


/// functor counts to get the discrete surface for a material no. in direction (1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1) and total surface, then it converts it into phys units
template <typename T>
class SuperGeometryFaces2D final : public GenericF<T,int> {
private:
  SuperGeometry2D<T>&   _superGeometry;
  const int             _material;
  const T _latticeL;
public:
  template<typename DESCRIPTOR>
  SuperGeometryFaces2D(SuperGeometry2D<T>& superGeometry, const int material, const UnitConverter<T,DESCRIPTOR>& converter);
  SuperGeometryFaces2D(SuperGeometry2D<T>& superGeometry, const int material, T latticeL);
  bool operator() (T output[], const int input[]) override;
};


/// functor counts to get the discrete surface for a material no. in direction (1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1) and total surface, then it converts it into phys units
template <typename T, bool HLBM>
class SuperGeometryFacesIndicator2D final : public GenericF<T,int> {
private:
  SuperGeometry2D<T>&   _superGeometry;
  SmoothIndicatorF2D<T,T,HLBM>& _indicator;
  const int             _material;
  T _latticeL;
public:
  SuperGeometryFacesIndicator2D(SuperGeometry2D<T>& superGeometry, SmoothIndicatorF2D<T,T,HLBM>& indicator, const int material,
                                T deltaX);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise phys force acting on a boundary with a given material on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysDrag2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysDrag2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                         SuperGeometry2D<T>& superGeometry, const int material,
                         const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/**
 *  functor to get pointwise phys force acting on a boundary with a given material on local lattice
 *  see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
 */
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysCorrDrag2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysCorrDrag2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                             SuperGeometry2D<T>& superGeometry, const int material,
                             const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


} // end namespace olb

#endif
