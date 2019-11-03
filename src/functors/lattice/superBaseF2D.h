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

#ifndef SUPER_BASE_F_2D_H
#define SUPER_BASE_F_2D_H

#include <memory>

#include "functors/genericF.h"
#include "blockBaseF2D.h"
#include "communication/superStructure2D.h"
#include "core/superData2D.h"
#include "core/superLattice2D.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

template<typename T, typename DESCRIPTOR> class SuperLattice2D;
template<typename T, typename BaseType> class SuperData2D;
template<typename T> class SuperStructure2D;
template<typename T, typename W> class SuperIdentity2D;
template<typename T> class BlockF2D;

/// represents all functors that operate on a SuperStructure2D<T> in general
template <typename T, typename W = T>
class SuperF2D : public GenericF<W, int> {
protected:
  SuperF2D(SuperStructure2D<T>& superStructure, int targetDim);

  SuperStructure2D<T>& _superStructure;
  /// Super functors may consist of several BlockF2D<W> derived functors
  /**
   * By convention: If block level functors are used at all they should
   * number exactly LoadBalancer<T>::size per process.
   **/
  std::vector<std::unique_ptr<BlockF2D<W>>> _blockF;
public:
  using identity_functor_type = SuperIdentity2D<T,W>;

  SuperF2D<T,W>& operator-(SuperF2D<T,W>& rhs);
  SuperF2D<T,W>& operator+(SuperF2D<T,W>& rhs);
  SuperF2D<T,W>& operator*(SuperF2D<T,W>& rhs);
  SuperF2D<T,W>& operator/(SuperF2D<T,W>& rhs);

  /// \return _superStructure
  SuperStructure2D<T>& getSuperStructure();
  /// \return Size of SuperF3D<T,W>::_blockF vector
  int getBlockFSize() const;
  /// \return _blockF[iCloc]
  BlockF2D<W>& getBlockF(int iCloc);
};


/// Functor from `SuperData2D`
template<typename T, typename BaseType>
class SuperDataF2D : public SuperF2D<T,BaseType> {
protected:
  /// `SuperData2D` object this functor was created from
  SuperData2D<T,BaseType>& _superData;
public:
  /// Constructor from `SuperData2D` - stores `_superData` reference
  SuperDataF2D(SuperData2D<T,BaseType>& superData);
  /// Operator for this functor - copies data from `_superData` object into output
  bool operator() (BaseType output[], const int input[]);
  /// Getter for `_superData`
  SuperData2D<T,BaseType>& getSuperData();
};

/// identity functor for memory management
template <typename T, typename W>
class SuperIdentity2D : public SuperF2D<T,W> {
protected:
  FunctorPtr<SuperF2D<T,W>> _f;
public:
  SuperIdentity2D(FunctorPtr<SuperF2D<T,W>>&& f);
  bool operator() (W output[], const int input[]) override;
};

/// represents all functors that operate on a SuperLattice in general, e.g. getVelocity(), getForce(), getPressure()
template <typename T, typename DESCRIPTOR>
class SuperLatticeF2D : public SuperF2D<T,T> {
protected:
  SuperLatticeF2D(SuperLattice2D<T,DESCRIPTOR>& superLattice, int targetDim);

  SuperLattice2D<T,DESCRIPTOR>& _sLattice;
public:
  SuperLattice2D<T,DESCRIPTOR>& getSuperLattice();
};

/// represents all functors that operate on a DESCRIPTOR with output in Phys, e.g. physVelocity(), physForce(), physPressure()
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysF2D : public SuperLatticeF2D<T,DESCRIPTOR> {
protected:
  SuperLatticePhysF2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                      const UnitConverter<T,DESCRIPTOR>& converter, int targetDim);
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  UnitConverter<T,DESCRIPTOR> const& getConverter() const;
};

/// represents all thermal functors that operate on a DESCRIPTOR with output in Phys, e.g. physTemperature(), physHeatFlux()
template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class SuperLatticeThermalPhysF2D : public SuperLatticeF2D<T,TDESCRIPTOR> {
protected:
  SuperLatticeThermalPhysF2D(SuperLattice2D<T,TDESCRIPTOR>& sLattice,
                             const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter, int targetDim);
  const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& _converter;
public:
  ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& getConverter() const;
};

} // end namespace olb

#endif
