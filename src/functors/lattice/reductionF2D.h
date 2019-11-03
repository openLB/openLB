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

#ifndef REDUCTION_F_2D_H
#define REDUCTION_F_2D_H

#include "functors/analytical/analyticalF.h"
#include "superBaseF2D.h"
#include "geometry/superGeometry2D.h"
#include "geometry/cuboidGeometry2D.h"

namespace olb {


template< typename T, typename DESCRIPTOR> class SuperLattice2D;

/// Functor used to convert analytical functions to lattice functions
/**
 *  Input functions are interpreted as SI->SI units, the resulting lattice
 *  function will map lattice->lattice units
 *
 *  Maintains block level BlockLatticeFfromAnalyticalF2D functors.
 */
template <typename T, typename DESCRIPTOR>
class SuperLatticeFfromAnalyticalF2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
protected:
  FunctorPtr<AnalyticalF2D<T,T>> _f;
public:
  /**
   * \param f        Analytical functor to be converted into a lattice functor
   * \param sLattice DESCRIPTOR reference required for conversion and block functor construction
   **/
  SuperLatticeFfromAnalyticalF2D(FunctorPtr<AnalyticalF2D<T,T>>&& f,
                                 SuperLattice2D<T,DESCRIPTOR>&    sLattice);
  bool operator() (T output[], const int input[]) override;
};


/// Block level functor for conversion of analytical to lattice functors.
/**
 * Instances are contained in SuperLatticeFfromAnalyticalF2D::_blockF.
 **/
template <typename T, typename DESCRIPTOR>
class BlockLatticeFfromAnalyticalF2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
protected:
  AnalyticalF2D<T,T>& _f;
  Cuboid2D<T>&        _cuboid;
public:
  /**
   * \param f       Analytical functor to be converted into a lattice functor
   * \param lattice Block lattice structure required for BlockLatticeF2D construction
   * \param cuboid  Cuboid reference required for input parameter conversion
   **/
  BlockLatticeFfromAnalyticalF2D(AnalyticalF2D<T,T>&                    f,
                                 BlockLatticeStructure2D<T,DESCRIPTOR>& lattice,
                                 Cuboid2D<T>&                           cuboid);
  bool operator() (T output[], const int input[]) override;
};


} // end namespace olb

#endif
