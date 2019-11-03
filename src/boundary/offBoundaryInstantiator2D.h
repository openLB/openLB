/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Jonas Kratzke, Mathias J. Krause
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

/** \file
 * A helper for initialising 2D boundaries -- header file.
 */

#ifndef OFF_BOUNDARY_INSTANTIATOR_2D_H
#define OFF_BOUNDARY_INSTANTIATOR_2D_H

#include "offBoundaryCondition2D.h"
#include "offBoundaryPostProcessors2D.h"
#include "geometry/blockGeometry2D.h"
#include "geometry/blockGeometryStatistics2D.h"
#include "dynamics/dynamics.h"
#include "core/cell.h"
#include "io/ostreamManager.h"
#include "core/util.h"
#include "io/stlReader.h"
#include "functors/lattice/indicator/blockIndicatorF2D.h"

namespace olb {

/**
* This class gets the needed processors from BoundaryManager and adds them
* to the Processor Vector of the DESCRIPTOR
*/

template<typename T, typename DESCRIPTOR, class BoundaryManager>
class OffBoundaryConditionInstantiator2D: public OffLatticeBoundaryCondition2D<T,
  DESCRIPTOR> {
public:
  OffBoundaryConditionInstantiator2D(BlockLatticeStructure2D<T, DESCRIPTOR>& block_, T epsFraction_ = 0.0001);
  ~OffBoundaryConditionInstantiator2D() override;

  void addOnePointZeroVelocityBoundary(int x, int y, int iPop, T dist) override;
  void addTwoPointZeroVelocityBoundary(int x, int y, int iPop, T dist) override;

  void addOnePointVelocityBoundary(int x, int y, int iPop, T dist) override;
  void addTwoPointVelocityBoundary(int x, int y, int iPop, T dist) override;

  void addOffDynamics(int x, int y, T location[DESCRIPTOR::d]) override;
  void addOffDynamics(int x, int y, T location[DESCRIPTOR::d], T distances[DESCRIPTOR::q]) override;
  void addOffDynamics(BlockIndicatorF2D<T>& indicator) override;

  void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, int iPop, T dist) override;
  void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, T distances[DESCRIPTOR::q]);
  void addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator) override;
  void addZeroVelocityBoundary(BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator) override;
  void addZeroVelocityBoundary(BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator) override;

  void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, int iPop, T dist) override;
  void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, T distances[DESCRIPTOR::q]);
  void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator) override;
  void addVelocityBoundary(BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator) override;
  void addVelocityBoundary(BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator) override;

  void addPressureBoundary(BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator) override;
  void addPressureBoundary(BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator) override;

  void setBoundaryIntersection(int iX, int iY, int iPop, T distance);
  bool getBoundaryIntersection(int iX, int iY, int iPop, T point[DESCRIPTOR::d]);

  void defineU(int iX, int iY, int iPop, const T u[DESCRIPTOR::d]) override;
  void defineU(BlockIndicatorF2D<T>& indicator, BlockIndicatorF2D<T>& bulkIndicator, AnalyticalF2D<T,T>& u) override;

  void defineRho(int iX, int iY, int iPop, const T rho) override;
  void defineRho(BlockIndicatorF2D<T>& indicator, BlockIndicatorF2D<T>& bulkIndicator, AnalyticalF2D<T,T>& rho) override;

  void outputOn() override;
  void outputOff() override;

  BlockLatticeStructure2D<T, DESCRIPTOR>& getBlock() override;
  BlockLatticeStructure2D<T, DESCRIPTOR> const& getBlock() const override;
private:
  BlockLatticeStructure2D<T, DESCRIPTOR>& block;
  //std::vector<Momenta<T, DESCRIPTOR>*> momentaVector;
  std::vector<Dynamics<T, DESCRIPTOR>*> dynamicsVector;
  bool _output;
  mutable OstreamManager clout;
};

///////// class OffBoundaryConditionInstantiator2D ////////////////////////

template<typename T, typename DESCRIPTOR, class BoundaryManager>
BlockLatticeStructure2D<T, DESCRIPTOR>& OffBoundaryConditionInstantiator2D<T, DESCRIPTOR,
                        BoundaryManager>::getBlock()
{
  return block;
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
BlockLatticeStructure2D<T, DESCRIPTOR> const& OffBoundaryConditionInstantiator2D<T, DESCRIPTOR,
                        BoundaryManager>::getBlock() const
{
  return block;
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::OffBoundaryConditionInstantiator2D(
  BlockLatticeStructure2D<T, DESCRIPTOR>& block_, T epsFraction_) :
  block(block_), _output(false), clout(std::cout,"BoundaryConditionInstantiator2D")
{
  this->_epsFraction = epsFraction_;
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::~OffBoundaryConditionInstantiator2D()
{
  for (unsigned iDynamics = 0; iDynamics < dynamicsVector.size(); ++iDynamics) {
    delete dynamicsVector[iDynamics];
  }
  /*
  for (unsigned iMomenta = 0; iMomenta < dynamicsVector.size(); ++iMomenta) {
    delete momentaVector[iMomenta];
  }*/
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addOnePointZeroVelocityBoundary(
  int x, int y, int iPop, T dist)
{
  PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor =
    BoundaryManager::getOnePointZeroVelocityBoundaryProcessor
    (x, y, iPop, dist);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addTwoPointZeroVelocityBoundary(
  int x, int y, int iPop, T dist)
{
  PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor =
    BoundaryManager::getTwoPointZeroVelocityBoundaryProcessor
    (x, y, iPop, dist);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addOnePointVelocityBoundary(
  int x, int y, int iPop, T dist)
{
  PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor =
    BoundaryManager::getOnePointVelocityBoundaryProcessor
    (x, y, iPop, dist);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addTwoPointVelocityBoundary(
  int x, int y, int iPop, T dist)
{
  PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor =
    BoundaryManager::getTwoPointVelocityBoundaryProcessor
    (x, y, iPop, dist);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addOffDynamics(
  int x, int y, T location[DESCRIPTOR::d])
{
  Dynamics<T,DESCRIPTOR>* dynamics
    = BoundaryManager::getOffDynamics(location);
  this->getBlock().defineDynamics(x,x,y,y, dynamics);
  dynamicsVector.push_back(dynamics);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addOffDynamics(
  int x, int y, T location[DESCRIPTOR::d], T distances[DESCRIPTOR::q])
{
  Dynamics<T,DESCRIPTOR>* dynamics
    = BoundaryManager::getOffDynamics(location, distances);
  this->getBlock().defineDynamics(x,x,y,y, dynamics);
  dynamicsVector.push_back(dynamics);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addOffDynamics(
  BlockIndicatorF2D<T>& indicator)
{
  if ( !indicator.isEmpty() ) {
    const Vector<int,2> min = indicator.getMin();
    const Vector<int,2> max = indicator.getMax();

    for (int iX = min[0]; iX <= max[0]; ++iX) {
      for (int iY = min[1]; iY <= max[1]; ++iY) {
        if (indicator(iX,iY)) {
          T location[DESCRIPTOR::d];
          indicator.getBlockGeometryStructure().getPhysR(location, iX,iY);
          addOffDynamics(iX, iY, location);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addZeroVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, int iPop, T dist)
{
  if (blockGeometryStructure.getMaterial(x-descriptors::c<DESCRIPTOR >(iPop,0), y-descriptors::c<DESCRIPTOR >(iPop,1)) != 1) {
    addOnePointZeroVelocityBoundary(x, y, iPop, dist);
  }
  else {
    addTwoPointZeroVelocityBoundary(x, y, iPop, dist);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::
addZeroVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, T distances[DESCRIPTOR::q])
{
  typedef DESCRIPTOR L;
  //T location[DESCRIPTOR::d];
  //location[0] = blockGeometryStructure.physCoordX(x);
  //location[1] = blockGeometryStructure.physCoordY(y);
  //location[2] = blockGeometryStructure.physCoordZ(z);
  //T distancesCopy[L::q];
  //T spacing = blockGeometryStructure.getDeltaR();
  //for (int iPop = 1; iPop < L::q ; ++iPop) {
  //  distancesCopy[iPop] = spacing*(1.-distances[iPop]);
  //  if (distances[iPop] == -1)
  //    distancesCopy[iPop] = -1;
  //}
  //addOffDynamics(x, y, z, location, distancesCopy);

  for (int iPop = 1; iPop < L::q ; ++iPop) {
    if ( !util::nearZero(distances[iPop]+1) ) {
      addZeroVelocityBoundary(blockGeometryStructure, x-descriptors::c<L>(iPop,0), y-descriptors::c<L>(iPop,1), iPop, distances[iPop]);
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addZeroVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY,
  BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator)
{
  T distances[DESCRIPTOR::q];
  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    distances[iPop] = -1;
  }

  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    const int iXn = iX + descriptors::c<DESCRIPTOR >(iPop,0);
    const int iYn = iY + descriptors::c<DESCRIPTOR >(iPop,1);
    if (bulkIndicator(iXn,iYn)) {
      T dist = -1;

      T physR[2];
      blockGeometryStructure.getPhysR(physR,iXn,iYn);
      T voxelSize=blockGeometryStructure.getDeltaR();

      Vector<T,2> physC(physR);

      Vector<T,2> direction(-voxelSize*descriptors::c<DESCRIPTOR >(iPop,0),-voxelSize*descriptors::c<DESCRIPTOR >(iPop,1)) ;
      T cPhysNorm = voxelSize*sqrt(descriptors::c<DESCRIPTOR >(iPop,0)*descriptors::c<DESCRIPTOR >(iPop,0)+descriptors::c<DESCRIPTOR >(iPop,1)*descriptors::c<DESCRIPTOR >(iPop,1));

      if (!geometryIndicator.distance(dist,physC,direction,blockGeometryStructure.getIcGlob() ) ) {
        T epsX = voxelSize*descriptors::c<DESCRIPTOR >(iPop,0)*this->_epsFraction;
        T epsY = voxelSize*descriptors::c<DESCRIPTOR >(iPop,1)*this->_epsFraction;

        Vector<T,2> physC2(physC);
        physC2[0] += epsX;
        physC2[1] += epsY;
        Vector<T,2> direction2(direction);
        direction2[0] -= 2.*epsX;
        direction2[1] -= 2.*epsY;

        if ( !geometryIndicator.distance(dist,physC2,direction2,blockGeometryStructure.getIcGlob())) {
          clout << "ERROR: no boundary found at (" << iXn << "," << iYn <<") ~ ("
                << physR[0] << "," << physR[1] << "), "
                << "in direction " << util::opposite<DESCRIPTOR >(iPop)
                << std::endl;
        }
        T distNew = (dist - sqrt(epsX*epsX+epsY*epsY) )/cPhysNorm;
        if (distNew < 0.5) {
          dist = 0;
        }
        else {
          dist = 0.5 * cPhysNorm;
          clout << "WARNING: distance at (" << iXn << "," << iYn << ") ~ ("
                << physR[0] << "," << physR[1] << "), "
                << "in direction " << util::opposite<DESCRIPTOR >(iPop) << ": "
                << distNew
                << " rounded to "
                << dist/cPhysNorm
                << std::endl;
        }
      }
      distances[util::opposite<DESCRIPTOR >(iPop)] = dist/cPhysNorm;
    } // bulk indicator if
  } // iPop loop
  addZeroVelocityBoundary(blockGeometryStructure, iX, iY, distances);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addZeroVelocityBoundary(
  BlockIndicatorF2D<T>& boundaryIndicator,
  BlockIndicatorF2D<T>& bulkIndicator,
  IndicatorF2D<T>&      geometryIndicator)
{
  if ( !boundaryIndicator.isEmpty() ) {
    const Vector<int,2> min = boundaryIndicator.getMin();
    const Vector<int,2> max = boundaryIndicator.getMax();

    for (int iX = min[0]; iX <= max[0]; ++iX) {
      for (int iY = min[1]; iY <= max[1]; ++iY) {
        if (boundaryIndicator(iX,iY)) {
          addZeroVelocityBoundary(bulkIndicator.getBlockGeometryStructure(), iX, iY,
                                  bulkIndicator, geometryIndicator);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addZeroVelocityBoundary(
  BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator)
{
  if ( boundaryIndicator.isEmpty() ) {
    return;
  }

  auto& blockGeometryStructure = boundaryIndicator.getBlockGeometryStructure();
  const Vector<int,2> min = boundaryIndicator.getMin();
  const Vector<int,2> max = boundaryIndicator.getMax();

  for (int iX = min[0]; iX <= max[0]; ++iX) {
    for (int iY = min[1]; iY <= max[1]; ++iY) {
      if (boundaryIndicator(iX, iY)) {
        bool streamDirections[DESCRIPTOR::q];
        // check if (ix,iY) is really a boundary node and compute the stream direction
        if (blockGeometryStructure.template findStreamDirections<T,DESCRIPTOR>(iX, iY, boundaryIndicator, bulkIndicator, streamDirections) ) {
          T distances[DESCRIPTOR::q];
          for (int iPop=0; iPop<DESCRIPTOR::q; iPop++) {
            if (streamDirections[iPop]) {
              distances[descriptors::opposite<DESCRIPTOR>(iPop)] = .5;
            }
            else {
              distances[descriptors::opposite<DESCRIPTOR>(iPop)] = -1;
            }
          }
          addZeroVelocityBoundary(blockGeometryStructure, iX, iY, distances);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, int iPop, T dist)
{
  if (blockGeometryStructure.getMaterial(x-descriptors::c<DESCRIPTOR >(iPop,0), y-descriptors::c<DESCRIPTOR >(iPop,1)) != 1) {
    addOnePointVelocityBoundary(x, y, iPop, dist);
  }
  else {
    addTwoPointVelocityBoundary(x, y, iPop, dist);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::
addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int x, int y, T distances[DESCRIPTOR::q])
{
  typedef DESCRIPTOR L;
  T location[DESCRIPTOR::d];
  blockGeometryStructure.getPhysR(location, x,y);

  T distancesCopy[L::q];
  T spacing = blockGeometryStructure.getDeltaR();
  for (int iPop = 1; iPop < L::q ; ++iPop) {
    distancesCopy[iPop] = spacing*(1.-distances[iPop]);
    if ( !util::nearZero(distances[iPop]+1) ) {
      distancesCopy[iPop] = -1;
    }
  }
  addOffDynamics(x, y, location, distancesCopy);

  for (int iPop = 1; iPop < L::q ; ++iPop) {
    if ( !util::nearZero(distances[iPop]+1) ) {
      addVelocityBoundary(blockGeometryStructure, x-descriptors::c<L>(iPop,0), y-descriptors::c<L>(iPop,1), iPop, distances[iPop]);
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addVelocityBoundary(
  BlockGeometryStructure2D<T>& blockGeometryStructure, int iX, int iY,
  BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator)
{
  T distances[DESCRIPTOR::q];
  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    distances[iPop] = -1;
  }

  for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
    const int iXn = iX + descriptors::c<DESCRIPTOR >(iPop,0);
    const int iYn = iY + descriptors::c<DESCRIPTOR >(iPop,1);
    if (bulkIndicator(iXn,iYn)) {
      T dist = -1;
      T physR[2];
      blockGeometryStructure.getPhysR(physR,iXn,iYn);
      T voxelSize=blockGeometryStructure.getDeltaR();
      Vector<T,2> physC(physR);

      Vector<T,2> direction(-voxelSize*descriptors::c<DESCRIPTOR >(iPop,0),-voxelSize*descriptors::c<DESCRIPTOR >(iPop,1));
      T cPhysNorm = voxelSize*sqrt(descriptors::c<DESCRIPTOR >(iPop,0)*descriptors::c<DESCRIPTOR >(iPop,0)+descriptors::c<DESCRIPTOR >(iPop,1)*descriptors::c<DESCRIPTOR >(iPop,1));

      if (!geometryIndicator.distance(dist,physC,direction,blockGeometryStructure.getIcGlob() ) ) {
        T epsX = voxelSize*descriptors::c<DESCRIPTOR >(iPop,0)*this->_epsFraction;
        T epsY = voxelSize*descriptors::c<DESCRIPTOR >(iPop,1)*this->_epsFraction;

        Vector<T,2> physC2(physC);
        physC2[0] += epsX;
        physC2[1] += epsY;
        Vector<T,2> direction2(direction);
        direction2[0] -= 2.*epsX;
        direction2[1] -= 2.*epsY;

        if ( !geometryIndicator.distance(dist,physC2,direction2,blockGeometryStructure.getIcGlob())) {
          clout << "ERROR: no boundary found at (" << iXn << "," << iYn <<") ~ ("
                << physR[0] << "," << physR[1] << "), "
                << "in direction " << util::opposite<DESCRIPTOR >(iPop)
                << std::endl;
        }
        T distNew = (dist - sqrt(epsX*epsX+epsY*epsY))/cPhysNorm;
        if (distNew < 0.5) {
          dist = 0;
        }
        else {
          dist = 0.5 * cPhysNorm;
          clout << "WARNING: distance at (" << iXn << "," << iYn <<") ~ ("
                << physR[0] << "," << physR[1] <<"), "
                << "in direction " << util::opposite<DESCRIPTOR >(iPop) << ": "
                << distNew
                << " rounded to "
                << dist/cPhysNorm
                << std::endl;
        }
      }
      distances[util::opposite<DESCRIPTOR >(iPop)] = dist/cPhysNorm;
    } // bulk indicator
  } // iPop loop
  addVelocityBoundary(blockGeometryStructure, iX, iY, distances);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addVelocityBoundary(
  BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator)
{
  if ( !boundaryIndicator.isEmpty() ) {
    const Vector<int,2> min = boundaryIndicator.getMin();
    const Vector<int,2> max = boundaryIndicator.getMax();

    for (int iX = min[0]; iX <= max[0]; ++iX) {
      for (int iY = min[1]; iY <= max[1]; ++iY) {
        if (boundaryIndicator(iX, iY)) {
          addVelocityBoundary(bulkIndicator.getBlockGeometryStructure(), iX, iY,
                              bulkIndicator, geometryIndicator);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addVelocityBoundary(
  BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator)
{
  if ( boundaryIndicator.isEmpty() ) {
    return;
  }

  auto& blockGeometryStructure = boundaryIndicator.getBlockGeometryStructure();
  const Vector<int,2> min = boundaryIndicator.getMin();
  const Vector<int,2> max = boundaryIndicator.getMax();

  for (int iX = min[0]; iX <= max[0]; ++iX) {
    for (int iY = min[1]; iY <= max[1]; ++iY) {
      if (boundaryIndicator(iX,iY)) {
        bool streamDirections[DESCRIPTOR::q];
        // check if (ix,iY) is really a boundary node and compute the stream direction
        if (blockGeometryStructure.template findStreamDirections<T,DESCRIPTOR>(iX, iY, boundaryIndicator, bulkIndicator, streamDirections) ) {
          T distances[DESCRIPTOR::q];
          for (int iPop=0; iPop<DESCRIPTOR::q; iPop++) {
            if (streamDirections[iPop]) {
              distances[descriptors::opposite<DESCRIPTOR>(iPop)] = .5;
            }
            else {
              distances[descriptors::opposite<DESCRIPTOR>(iPop)] = -1;
            }
          }
          addVelocityBoundary(blockGeometryStructure, iX, iY, distances);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addPressureBoundary(
  BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator, IndicatorF2D<T>& geometryIndicator)
{
  auto& blockGeometryStructure = boundaryIndicator.getBlockGeometryStructure();
  if ( boundaryIndicator.isEmpty() ) {
    const Vector<int,2> min = boundaryIndicator.getMin();
    const Vector<int,2> max = boundaryIndicator.getMax();

    for (int iX = min[0]; iX <= max[0]; ++iX) {
      for (int iY = min[1]; iY <= max[1]; ++iY) {
        bool streamDirections[DESCRIPTOR::q];
        // check if (ix,iY) is really a boundary node and compute the stream direction
        if (blockGeometryStructure.template findStreamDirections<T,DESCRIPTOR>(iX, iY, boundaryIndicator, bulkIndicator, streamDirections) ) {
          T location[DESCRIPTOR::d];
          blockGeometryStructure.getPhysR(location, iX, iX);
          block.defineDynamics(iX,iX,iY,iY,&instances::getBounceBackAnti<T, DESCRIPTOR>(1.));
        }
      }
    }
  }
  clout << "ERROR: OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addPressureBoundary NOT YET IMPLEMENTED" << std::endl;
  exit(-1);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addPressureBoundary(
  BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator)
{
  if ( boundaryIndicator.isEmpty() ) {
    return;
  }

  auto& blockGeometryStructure = boundaryIndicator.getBlockGeometryStructure();
  const Vector<int,2> min = boundaryIndicator.getMin();
  const Vector<int,2> max = boundaryIndicator.getMax();

  for (int iX = min[0]; iX <= max[0]; ++iX) {
    for (int iY = min[1]; iY <= max[1]; ++iY) {
      bool streamDirections[DESCRIPTOR::q];
      // check if (ix,iY) is really a boundary node and compute the stream direction
      if (blockGeometryStructure.template findStreamDirections<T,DESCRIPTOR>(iX, iY, boundaryIndicator, bulkIndicator, streamDirections) ) {
        T location[DESCRIPTOR::d];
        blockGeometryStructure.getPhysR(location, iX, iX);
        block.defineDynamics(iX,iX,iY,iY,&instances::getBounceBackAnti<T, DESCRIPTOR>(1.));
        for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
          if (streamDirections[iPop]) {
            PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor = new AntiBounceBackPostProcessorGenerator2D<T, DESCRIPTOR>(iX+descriptors::c<DESCRIPTOR>(iPop,0), iY+descriptors::c<DESCRIPTOR>(iPop,1), descriptors::opposite<DESCRIPTOR>(iPop));
            this->getBlock().addPostProcessor(*postProcessor);
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::
setBoundaryIntersection(int iX, int iY, int iPop, T distance)
{

  this->getBlock().getDynamics(iX, iY)->setBoundaryIntersection(iPop, distance);
  if (_output) {
    clout << "setBoundaryIntersection(" << iX << ", " << iY << " )" << std::endl;
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
bool OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::
getBoundaryIntersection(int iX, int iY, int iPop, T point[DESCRIPTOR::d])
{

  return this->getBlock().getDynamics(iX, iY)->getBoundaryIntersection(iPop, point);
}


template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::
defineU(int iX, int iY, int iPop, const T u[DESCRIPTOR::d])
{

  this->getBlock().getDynamics(iX, iY)->defineU(iPop, u);
  if (_output) {
    clout << "defineU(" << iX << ", " << iY << " )" << std::endl;
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::
defineU(BlockIndicatorF2D<T>& indicator,
        BlockIndicatorF2D<T>& bulkIndicator,
        AnalyticalF2D<T,T>& u)
{
  if ( indicator.isEmpty() ) {
    return;
  }

  const Vector<int,2> min = indicator.getMin();
  const Vector<int,2> max = indicator.getMax();

  for (int iX = min[0]; iX <= max[0]; ++iX) {
    for (int iY = min[1]; iY <= max[1]; ++iY) {
      if (indicator(iX, iY)) {
        for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
          // Get direction
          int iXn = iX + descriptors::c<DESCRIPTOR >(iPop,0);
          int iYn = iY + descriptors::c<DESCRIPTOR >(iPop,1);
          if (bulkIndicator(iXn, iYn)) {
            T intersection[] = { T(), T() }; // coord. of intersection
            int opp = util::opposite<DESCRIPTOR >(iPop);
            if (this->getBoundaryIntersection(iX, iY, opp, intersection) ) {
              T vel[]= {T(),T()};
              u(vel,intersection);
              this->defineU(iX, iY, opp, vel);
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::
defineRho(int iX, int iY, int iPop, const T rho)
{
  this->getBlock().getDynamics(iX, iY)->defineRho(iPop, rho);
  if (_output) {
    clout << "defineRho(" << iX << ", " << iY << " )" << std::endl;
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::
defineRho(BlockIndicatorF2D<T>& indicator,
          BlockIndicatorF2D<T>& bulkIndicator,
          AnalyticalF2D<T,T>&   rho)
{
  if ( indicator.isEmpty() ) {
    return;
  }

  const Vector<int,2> min = indicator.getMin();
  const Vector<int,2> max = indicator.getMax();

  for (int iX = min[0]; iX <= max[0]; ++iX) {
    for (int iY = min[1]; iY <= max[1]; ++iY) {
      if (indicator(iX, iY)) {
        for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
          // Get direction
          int iXn = iX + descriptors::c<DESCRIPTOR >(iPop,0);
          int iYn = iY + descriptors::c<DESCRIPTOR >(iPop,1);
          if (bulkIndicator(iXn, iYn)) {
            T intersection[] = { T(), T() }; // coord. of intersection
            int opp = util::opposite<DESCRIPTOR >(iPop);
            if (this->getBoundaryIntersection(iX, iY, opp, intersection) ) {
              T rhoLocal[]= {T(1)};
              rho(rhoLocal,intersection);
              this->defineRho(iX, iY, opp, rhoLocal[0]);
            }
          }
        }
      }
    }
  }
}


template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::outputOn()
{
  _output = true;
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void OffBoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::outputOff()
{
  _output = false;
}

}

#endif
