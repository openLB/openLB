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

#ifndef SUPER_BASE_F_2D_HH
#define SUPER_BASE_F_2D_HH


#include "superBaseF2D.h"

namespace olb {

template<typename T, typename W>
SuperF2D<T,W>::SuperF2D(SuperStructure2D<T>& superStructure, int targetDim)
  : GenericF<W,int>(targetDim,3), _superStructure(superStructure) { }


template<typename T, typename W>
SuperStructure2D<T>& SuperF2D<T,W>::getSuperStructure()
{
  return _superStructure;
}

template <typename T, typename W>
int SuperF2D<T,W>::getBlockFSize() const
{
  OLB_ASSERT(_blockF.size() < INT32_MAX,
             "it is safe to cast std::size_t to int");
  return _blockF.size();
}

template <typename T, typename W>
BlockF2D<W>& SuperF2D<T,W>::getBlockF(int iCloc)
{
  OLB_ASSERT(size_t(iCloc) < _blockF.size() && iCloc >= 0,
             "block functor index within bounds");
  return *(_blockF[iCloc]);
}


template <typename T,typename BaseType>
SuperDataF2D<T,BaseType>::SuperDataF2D(SuperData2D<T,BaseType>& superData)
  : SuperF2D<T,BaseType>(superData, superData.getDataSize()),
    _superData(superData)
{
  for (int iC = 0; iC < _superData.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockDataViewF2D<T,BaseType>(_superData.get(iC),
                                       _superData.getOverlap())
    );
  }
}

template <typename T,typename BaseType>
bool SuperDataF2D<T,BaseType>::operator() (BaseType output[], const int input[])
{
  const auto& load = _superData.getLoadBalancer();
  if (load.rank(input[0]) == singleton::mpi().getRank()) {
    return this->getBlockF(load.loc(input[0]))(output, &input[1]);
  }
  else {
    return false;
  }
}

template <typename T,typename BaseType>
SuperData2D<T,BaseType>& SuperDataF2D<T,BaseType>::getSuperData()
{
  return _superData;
}


template <typename T, typename W>
SuperIdentity2D<T,W>::SuperIdentity2D(FunctorPtr<SuperF2D<T,W>>&& f)
  : SuperF2D<T,W>(f->getSuperStructure(), f->getTargetDim()),
    _f(std::move(f))
{
  this->getName() = "Id(" + _f->getName() + ")";

  for (int iC = 0; iC < _f->getBlockFSize(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIdentity2D<W>(_f->getBlockF(iC)));
  }
}

template <typename T, typename W>
bool SuperIdentity2D<T,W>::operator()(W output[], const int input[])
{
  return _f(output, input);
}


template <typename T, typename DESCRIPTOR>
SuperLatticeF2D<T,DESCRIPTOR>::SuperLatticeF2D(SuperLattice2D<T,DESCRIPTOR>& superLattice,
    int targetDim)
  : SuperF2D<T,T>(superLattice, targetDim), _sLattice(superLattice) { }

template <typename T, typename DESCRIPTOR>
SuperLattice2D<T,DESCRIPTOR>& SuperLatticeF2D<T,DESCRIPTOR>::getSuperLattice()
{
  return _sLattice;
}

template <typename T, typename DESCRIPTOR>
SuperLatticePhysF2D<T,DESCRIPTOR>::SuperLatticePhysF2D
(SuperLattice2D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter,
 int targetDim)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, targetDim), _converter(converter) { }

template <typename T, typename DESCRIPTOR>
UnitConverter<T,DESCRIPTOR> const& SuperLatticePhysF2D<T,DESCRIPTOR>::getConverter() const
{
  return this->_converter;
}

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
SuperLatticeThermalPhysF2D<T,DESCRIPTOR,TDESCRIPTOR>::SuperLatticeThermalPhysF2D
(SuperLattice2D<T,TDESCRIPTOR>& sLattice, const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter,
 int targetDim)
  : SuperLatticeF2D<T,TDESCRIPTOR>(sLattice, targetDim), _converter(converter) { }

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& SuperLatticeThermalPhysF2D<T,DESCRIPTOR,TDESCRIPTOR>::getConverter() const
{
  return this->_converter;
}

} // end namespace olb

#endif
