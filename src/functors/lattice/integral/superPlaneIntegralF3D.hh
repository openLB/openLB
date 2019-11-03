/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerlaender
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

#ifndef SUPER_PLANE_INTEGRAL_F_3D_HH
#define SUPER_PLANE_INTEGRAL_F_3D_HH

#include "superPlaneIntegralF3D.h"
#include "utilities/vectorHelpers.h"
#include "functors/analytical/indicator/indicator2D.hh"

namespace olb {


template<typename T>
bool SuperPlaneIntegralF3D<T>::isToBeIntegrated(const Vector<T,3>& physR, int iC)
{
  Vector<int,4> latticeR;
  //get nearest lattice point
  if ( _geometry.getCuboidGeometry().getFloorLatticeR(physR, latticeR) ) {
    const int& iX = latticeR[1];
    const int& iY = latticeR[2];
    const int& iZ = latticeR[3];

    // interpolation is possible iff all neighbours are within the indicated subset
    return _integrationIndicatorF->operator()(   iC, iX,   iY,   iZ  )
           && _integrationIndicatorF->operator()(iC, iX,   iY,   iZ+1)
           && _integrationIndicatorF->operator()(iC, iX,   iY+1, iZ  )
           && _integrationIndicatorF->operator()(iC, iX,   iY+1, iZ+1)
           && _integrationIndicatorF->operator()(iC, iX+1, iY,   iZ  )
           && _integrationIndicatorF->operator()(iC, iX+1, iY,   iZ+1)
           && _integrationIndicatorF->operator()(iC, iX+1, iY+1, iZ  )
           && _integrationIndicatorF->operator()(iC, iX+1, iY+1, iZ+1);
  }
  else {
    return false;
  }
}

template<typename T>
SuperPlaneIntegralF3D<T>::SuperPlaneIntegralF3D(
  FunctorPtr<SuperF3D<T>>&& f,
  SuperGeometry3D<T>&       geometry,
  const HyperplaneLattice3D<T>& hyperplaneLattice,
  FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
  FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
  BlockDataReductionMode             mode)
  : SuperF3D<T>(f->getSuperStructure(), 2 + f->getTargetDim()),
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
    _v(hyperplaneLattice.getVectorV()),
    _normal(hyperplaneLattice.getHyperplane().normal)
{
  this->getName() = "SuperPlaneIntegralF3D";

  _normal.normalize();
  _u.normalize();
  _v.normalize();

  for ( const std::tuple<int,int,int>& pos : _reductionF.getRankLocalSubplane() ) {
    const int& i  = std::get<0>(pos);
    const int& j  = std::get<1>(pos);
    const int& iC = std::get<2>(pos);
    const Vector<T,3> physR = _reductionF.getPhysR(i, j);
    if (isToBeIntegrated(physR, iC)) {
      // check if interpolated hyperplane is to be restricted further
      // e.g. using IndicatorCircle2D
      if ( _subplaneIndicatorF ) {
        // determine physical coordinates relative to original hyperplane origin
        // [!] different from _reductionF._origin in the general case.
        const Vector<T,3> physRelativeToOrigin = physR - _origin;
        const T physOnHyperplane[2] {
          physRelativeToOrigin * _u,
          physRelativeToOrigin * _v
        };

        if ( _subplaneIndicatorF->operator()(physOnHyperplane) ) {
          _rankLocalSubplane.emplace_back(i, j);
        }
      }
      else {
        // plane is not restricted further
        _rankLocalSubplane.emplace_back(i, j);
      }
    }
  }
}

template<typename T>
SuperPlaneIntegralF3D<T>::SuperPlaneIntegralF3D(
  FunctorPtr<SuperF3D<T>>&& f,
  SuperGeometry3D<T>&       geometry,
  const Hyperplane3D<T>&    hyperplane,
  FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
  FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
  BlockDataReductionMode             mode)
  : SuperPlaneIntegralF3D(
      std::forward<decltype(f)>(f),
      geometry,
      HyperplaneLattice3D<T>(geometry.getCuboidGeometry(), hyperplane),
      std::forward<decltype(integrationIndicator)>(integrationIndicator),
      std::forward<decltype(subplaneIndicator)>(subplaneIndicator),
      mode)
{ }

template<typename T>
SuperPlaneIntegralF3D<T>::SuperPlaneIntegralF3D(
  FunctorPtr<SuperF3D<T>>&& f,
  SuperGeometry3D<T>&       geometry,
  const Hyperplane3D<T>&    hyperplane,
  FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
  BlockDataReductionMode             mode)
  : SuperPlaneIntegralF3D(
      std::forward<decltype(f)>(f),
      geometry,
      hyperplane,
      std::forward<decltype(integrationIndicator)>(integrationIndicator),
      nullptr,
      mode)
{ }

template<typename T>
SuperPlaneIntegralF3D<T>::SuperPlaneIntegralF3D(
  FunctorPtr<SuperF3D<T>>&& f,
  SuperGeometry3D<T>& geometry,
  const Vector<T,3>& origin, const Vector<T,3>& u, const Vector<T,3>& v,
  std::vector<int>       materials,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D(
      std::forward<decltype(f)>(f),
      geometry,
      Hyperplane3D<T>().originAt(origin).spannedBy(u, v),
      geometry.getMaterialIndicator(std::forward<decltype(materials)>(materials)),
      mode)
{ }

template<typename T>
SuperPlaneIntegralF3D<T>::SuperPlaneIntegralF3D(
  FunctorPtr<SuperF3D<T>>&& f,
  SuperGeometry3D<T>& geometry,
  const Vector<T,3>& origin, const Vector<T,3>& u, const Vector<T,3>& v,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D(
      std::forward<decltype(f)>(f),
      geometry,
      origin, u, v,
      std::vector<int>(1,1),
      mode)
{ }

template<typename T>
SuperPlaneIntegralF3D<T>::SuperPlaneIntegralF3D(
  FunctorPtr<SuperF3D<T>>&& f,
  SuperGeometry3D<T>& geometry,
  const Vector<T,3>& origin, const Vector<T,3>& normal,
  std::vector<int>       materials,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D(
      std::forward<decltype(f)>(f),
      geometry,
      Hyperplane3D<T>().originAt(origin).normalTo(normal),
      geometry.getMaterialIndicator(std::forward<decltype(materials)>(materials)),
      mode)
{ }

template<typename T>
SuperPlaneIntegralF3D<T>::SuperPlaneIntegralF3D(
  FunctorPtr<SuperF3D<T>>&& f,
  SuperGeometry3D<T>& geometry,
  const Vector<T,3>& origin, const Vector<T,3>& normal,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D(
      std::forward<decltype(f)>(f),
      geometry,
      origin, normal,
      std::vector<int>(1,1),
      mode)
{ }

template<typename T>
SuperPlaneIntegralF3D<T>::SuperPlaneIntegralF3D(
  FunctorPtr<SuperF3D<T>>&& f,
  SuperGeometry3D<T>& geometry,
  const Vector<T,3>&  normal,
  std::vector<int>    materials,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D(
      std::forward<decltype(f)>(f),
      geometry,
      Hyperplane3D<T>()
      .centeredIn(geometry.getCuboidGeometry().getMotherCuboid())
      .normalTo(normal),
      geometry.getMaterialIndicator(std::forward<decltype(materials)>(materials)),
      mode)
{ }

template<typename T>
SuperPlaneIntegralF3D<T>::SuperPlaneIntegralF3D(
  FunctorPtr<SuperF3D<T>>&& f,
  SuperGeometry3D<T>& geometry,
  const Vector<T,3>&  normal,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D(
      std::forward<decltype(f)>(f),
      geometry,
      normal,
      std::vector<int>(1,1),
      mode)
{ }

template<typename T>
SuperPlaneIntegralF3D<T>::SuperPlaneIntegralF3D(
  FunctorPtr<SuperF3D<T>>&& f,
  SuperGeometry3D<T>& geometry,
  const IndicatorCircle3D<T>& circle, std::vector<int> materials,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D(
      std::forward<decltype(f)>(f),
      geometry,
      Hyperplane3D<T>().originAt(circle.getCenter()).normalTo(circle.getNormal()),
      geometry.getMaterialIndicator(std::forward<std::vector<int>>(materials)),
      std::unique_ptr<IndicatorF2D<T>>(new IndicatorCircle2D<T>({0,0}, circle.getRadius())),
      mode)
{ }

template<typename T>
SuperPlaneIntegralF3D<T>::SuperPlaneIntegralF3D(
  FunctorPtr<SuperF3D<T>>&& f,
  SuperGeometry3D<T>& geometry,
  const IndicatorCircle3D<T>& circle,
  BlockDataReductionMode mode)
  : SuperPlaneIntegralF3D(
      std::forward<decltype(f)>(f),
      geometry,
      circle,
      std::vector<int>(1,1),
      mode)
{ }

template<typename T>
bool SuperPlaneIntegralF3D<T>::operator()(T output[], const int input[])
{
  this->getSuperStructure().communicate();

  _reductionF.update();

  const int flowDim = _reductionF.getTargetDim();

  std::vector<T> flow(flowDim,0.);

  for ( std::tuple<int,int>& pos : _rankLocalSubplane ) {
    T outputTmp[flowDim];
    const int inputTmp[2] { std::get<0>(pos), std::get<1>(pos) };

    _reductionF(outputTmp, inputTmp);

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
    output[0] = flow[0] * h * h;
    break;
  }
  case 3: {
    output[0] = (h*h * Vector<T,3>(flow)) * _normal;
    break;
  }
  }

  // area
  output[1] = vox * h * h;
  // write flow to output[2..]
  std::copy_n(flow.cbegin(), flowDim, &output[2]);

  return true;
}


}

#endif
