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

#ifndef SUPER_LP_NORM_3D_HH
#define SUPER_LP_NORM_3D_HH

#include "superLpNorm3D.h"
#include "blockLpNorm3D.h"
#include "functors/lattice/superBaseF3D.h"
#include "functors/lattice/indicator/superIndicatorF3D.h"
#include "geometry/superGeometry3D.h"
#include "latticeIntegralCommon.h"
#include "utilities/functorPtr.hh"

namespace olb {

template <typename T, typename W, int P>
SuperLpNorm3D<T,W,P>::SuperLpNorm3D(FunctorPtr<SuperF3D<T,W>>&&        f,
                                    FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperF3D<T,W>(f->getSuperStructure(),1),
    _f(std::move(f)),
    _indicatorF(std::move(indicatorF))
{
  OLB_ASSERT(_f->getSourceDim() == _indicatorF->getSourceDim(),
             "functor source dimension equals indicator source dimension");

  this->getName() = "L" + std::to_string(P) + "Norm(" + _f->getName() + ")";

  LoadBalancer<T>& load = _f->getSuperStructure().getLoadBalancer();

  if ( _f->getBlockFSize()          == load.size() &&
       _indicatorF->getBlockFSize() == load.size() ) {
    for (int iC = 0; iC < load.size(); ++iC) {
      this->_blockF.emplace_back(
        new BlockLpNorm3D<T,W,P>(_f->getBlockF(iC),
                                 _indicatorF->getBlockIndicatorF(iC))
      );
    }
  }
}

template <typename T, typename W, int P>
SuperLpNorm3D<T,W,P>::SuperLpNorm3D(FunctorPtr<SuperF3D<T,W>>&&        f,
                                    SuperGeometry3D<T>&                geometry,
                                    FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperLpNorm3D(std::forward<decltype(f)>(f),
                  std::forward<decltype(indicatorF)>(indicatorF))
{ }

template <typename T, typename W, int P>
SuperLpNorm3D<T,W,P>::SuperLpNorm3D(FunctorPtr<SuperF3D<T,W>>&& f,
                                    SuperGeometry3D<T>&         geometry,
                                    std::vector<int>            materials)
  : SuperLpNorm3D(std::forward<decltype(f)>(f),
                  geometry.getMaterialIndicator(std::move(materials)))
{ }

template <typename T, typename W, int P>
SuperLpNorm3D<T,W,P>::SuperLpNorm3D(FunctorPtr<SuperF3D<T,W>>&& f,
                                    SuperGeometry3D<T>&         geometry,
                                    int                         material)
  : SuperLpNorm3D(std::forward<decltype(f)>(f),
                  geometry.getMaterialIndicator(material))
{ }

template <typename T, typename W, int P>
bool SuperLpNorm3D<T,W,P>::operator() (W output[], const int input[])
{
  _f->getSuperStructure().communicate();
  CuboidGeometry3D<T>& cGeometry = _f->getSuperStructure().getCuboidGeometry();
  LoadBalancer<T>&     load      = _f->getSuperStructure().getLoadBalancer();

  output[0] = W(0);
  W outputTmp[_f->getTargetDim()];
  int inputTmp[_f->getSourceDim()];

  for (int iC = 0; iC < load.size(); ++iC) {
    Cuboid3D<T>& cuboid = cGeometry.get(load.glob(iC));

    const int nX = cuboid.getNx();
    const int nY = cuboid.getNy();
    const int nZ = cuboid.getNz();
    const T weight = pow(cuboid.getDeltaR(), 3);

    inputTmp[0] = load.glob(iC);

    for (inputTmp[1] = 0; inputTmp[1] < nX; ++inputTmp[1]) {
      for (inputTmp[2] = 0; inputTmp[2] < nY; ++inputTmp[2]) {
        for (inputTmp[3] = 0; inputTmp[3] < nZ; ++inputTmp[3]) {
          if (_indicatorF(inputTmp)) {
            _f(outputTmp, inputTmp);
            for (int iDim = 0; iDim < _f->getTargetDim(); ++iDim) {
              output[0] = LpNormImpl<T,W,P>()(output[0], outputTmp[iDim], weight);
            }
          }
        }
      }
    }
  }

#ifdef PARALLEL_MODE_MPI
  if (P == 0) {
    singleton::mpi().reduceAndBcast(output[0], MPI_MAX);
  }
  else {
    singleton::mpi().reduceAndBcast(output[0], MPI_SUM);
  }
#endif

  output[0] = LpNormImpl<T,W,P>().enclose(output[0]);

  return true;
}

}

#endif
