/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#ifndef SUPER_GEOMETRY_FACES_3D_HH
#define SUPER_GEOMETRY_FACES_3D_HH

#include "superGeometryFaces3D.h"
#include "geometry/superGeometry3D.h"

namespace olb {


template<typename T>
SuperGeometryFaces3D<T>::SuperGeometryFaces3D(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF, T latticeL)
  : SuperF3D<T>(indicatorF->getSuperStructure(), 7),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "superGeometryFaces";
  for (int iC = 0; iC < this->getSuperStructure().getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockGeometryFaces3D<T>(indicatorF->getBlockIndicatorF(iC), latticeL));
  }
}

template<typename T>
SuperGeometryFaces3D<T>::SuperGeometryFaces3D(
  SuperGeometry3D<T>& superGeometry, const int material, T latticeL)
  : SuperGeometryFaces3D(superGeometry.getMaterialIndicator(material), latticeL)
{ }

template<typename T>
bool SuperGeometryFaces3D<T>::operator()(T output[], const int input[])
{
  this->getSuperStructure().communicate();

  T blockOutput[7] = { };
  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
    output[iDim] = T();
  }

  for (int iC = 0; iC < this->getSuperStructure().getLoadBalancer().size(); ++iC) {
    this->getBlockF(iC)(blockOutput, input);
    for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
      output[iDim] += blockOutput[iDim];
    }
  }

#ifdef PARALLEL_MODE_MPI
  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
    singleton::mpi().reduceAndBcast(output[iDim], MPI_SUM);
  }
#endif
  return true;
}


}

#endif
