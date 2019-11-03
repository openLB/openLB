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

#ifndef SUPER_PLANE_INTEGRAL_F_2D_HH
#define SUPER_PLANE_INTEGRAL_F_2D_HH

#include "superPlaneIntegralF2D.h"
#include "utilities/vectorHelpers.h"
#include "utilities/functorPtr.hh"
#include "functors/lattice/indicator/indicator2D.hh"

namespace olb {


template<typename T>
bool SuperPlaneIntegralF2D<T>::isToBeIntegrated(const Vector<T,2>& physR, int iC)
{
  Vector<int,3> latticeR;
  //get nearest lattice point
  if ( _geometry.getCuboidGeometry().getFloorLatticeR(physR, latticeR) ) {
    const int& iX = latticeR[1];
    const int& iY = latticeR[2];

    // interpolation is possible iff all neighbours are within the indicated subset
    return _integrationIndicatorF->operator()(   iC, iX,   iY  )
           && _integrationIndicatorF->operator()(iC, iX,   iY+1)
           && _integrationIndicatorF->operator()(iC, iX+1, iY  )
           && _integrationIndicatorF->operator()(iC, iX+1, iY+1);
  }
  else {
    return false;
  }
}

template<typename T>
SuperPlaneIntegralF2D<T>::SuperPlaneIntegralF2D(
  FunctorPtr<SuperF2D<T>>&& f,
  SuperGeometry2D<T>&       geometry,
  const HyperplaneLattice2D<T>& hyperplaneLattice,
  FunctorPtr<SuperIndicatorF2D<T>>&& integrationIndicator,
  FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
  BlockDataReductionMode             mode)
  : SuperF2D<T>(f->getSuperStructure(), 2 + f->getTargetDim()),
    _geometry(geometry),
    _f(std::move(f)),
    _integrationIndicatorF(std::move(integrationIndicator)),
    _subplaneIndicatorF(std::move(subplaneIndicator)),
    _reductionF(*_f,
                hyperplaneLattice,
                BlockDataSyncMode::None,
                mode),
    _origin(hyperplaneLattice.getHyperplane().origin),
    _u(hyperplaneLattice.getVectorU()),
    _normal(hyperplaneLattice.getHyperplane().normal)
{
  this->getName() = "SuperPlaneIntegralF2D";

  _normal.normalize();
  _u.normalize();

  for ( const std::tuple<int,int>& pos : _reductionF.getRankLocalSubplane() ) {
    const int& i  = std::get<0>(pos);
    const int& iC = std::get<1>(pos);
    const Vector<T,2> physR = _reductionF.getPhysR(i);
    if (isToBeIntegrated(physR, iC)) {
      // check if interpolated hyperplane is to be restricted further
      // e.g. using IndicatorCircle2D
      if ( _subplaneIndicatorF ) {
        // determine physical coordinates relative to original hyperplane origin
        // [!] different from _reductionF._origin in the general case.
        const Vector<T,2> physRelativeToOrigin = physR - _origin;
        const T physOnHyperplane = physRelativeToOrigin * _u;

        if ( _subplaneIndicatorF->operator()(&physOnHyperplane) ) {
          _rankLocalSubplane.emplace_back(i);
        }
      }
      else {
        // plane is not restricted further
        _rankLocalSubplane.emplace_back(i);
      }
    }
  }
}

template<typename T>
SuperPlaneIntegralF2D<T>::SuperPlaneIntegralF2D(
  FunctorPtr<SuperF2D<T>>&& f,
  SuperGeometry2D<T>&       geometry,
  const Hyperplane2D<T>&    hyperplane,
  FunctorPtr<SuperIndicatorF2D<T>>&& integrationIndicator,
  FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
  BlockDataReductionMode             mode)
  : SuperPlaneIntegralF2D(
      std::forward<decltype(f)>(f),
      geometry,
      HyperplaneLattice2D<T>(geometry.getCuboidGeometry(), hyperplane),
      std::forward<decltype(integrationIndicator)>(integrationIndicator),
      std::forward<decltype(subplaneIndicator)>(subplaneIndicator),
      mode)
{ }

template<typename T>
SuperPlaneIntegralF2D<T>::SuperPlaneIntegralF2D(
  FunctorPtr<SuperF2D<T>>&& f,
  SuperGeometry2D<T>&       geometry,
  const Hyperplane2D<T>&    hyperplane,
  FunctorPtr<SuperIndicatorF2D<T>>&& integrationIndicator,
  BlockDataReductionMode             mode)
  : SuperPlaneIntegralF2D(
      std::forward<decltype(f)>(f),
      geometry,
      hyperplane,
      std::forward<decltype(integrationIndicator)>(integrationIndicator),
      nullptr,
      mode)
{ }

template<typename T>
SuperPlaneIntegralF2D<T>::SuperPlaneIntegralF2D(
  FunctorPtr<SuperF2D<T>>&& f,
  SuperGeometry2D<T>& geometry,
  const Vector<T,2>& origin, const Vector<T,2>& u,
  std::vector<int> materials,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF2D(
      std::forward<decltype(f)>(f),
      geometry,
      Hyperplane2D<T>().originAt(origin).parallelTo(u),
      geometry.getMaterialIndicator(std::forward<decltype(materials)>(materials)),
      mode)
{ }

template<typename T>
SuperPlaneIntegralF2D<T>::SuperPlaneIntegralF2D(
  FunctorPtr<SuperF2D<T>>&& f,
  SuperGeometry2D<T>& geometry,
  const Vector<T,2>& origin, const Vector<T,2>& u,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF2D(
      std::forward<decltype(f)>(f),
      geometry,
      origin, u,
      std::vector<int>(1,1),
      mode)
{ }


template<typename T>
bool SuperPlaneIntegralF2D<T>::operator()(T output[], const int input[])
{
  this->getSuperStructure().communicate();

  _reductionF.update();

  const int flowDim = _reductionF.getTargetDim();

  std::vector<T> flow(flowDim,0.);

  for ( int pos : _rankLocalSubplane ) {
    T outputTmp[flowDim];
    _reductionF(outputTmp, pos);

    for ( int j = 0; j < flowDim; j++ ) {
      flow[j] += outputTmp[j];
    }
  }

  int vox = _rankLocalSubplane.size();

#ifdef PARALLEL_MODE_MPI
  for ( int j = 0; j < flowDim; j++ ) {
    singleton::mpi().reduceAndBcast(flow[j], MPI_SUM);
  }
  singleton::mpi().reduceAndBcast(vox, MPI_SUM);
#endif

  const T h = _reductionF.getPhysSpacing();

  switch ( flowDim ) {
  case 1: {
    output[0] = flow[0] * h;
    break;
  }
  case 2: {
    output[0] = (h * Vector<T,2>(flow)) * _normal;
    break;
  }
  }

  // area
  output[1] = vox * h;
  // write flow to output[2..]
  std::copy_n(flow.cbegin(), flowDim, &output[2]);

  return true;
}


}

#endif
