/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
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
 * A helper for initialising 3D boundaries -- header file.
 */

#ifndef BOUNDARY_INSTANTIATOR_3D_H
#define BOUNDARY_INSTANTIATOR_3D_H

#include "boundaryCondition3D.h"
#include "boundaryPostProcessors3D.h"
#include "wallFunctionBoundaryPostProcessors3D.h"
#include "io/ostreamManager.h"
#include "dynamics/wallFunctionLatticeDescriptors.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "dynamics/freeEnergyDynamics.h"

namespace olb {


template<typename T, typename DESCRIPTOR, class BoundaryManager>
class BoundaryConditionInstantiator3D : public OnLatticeBoundaryCondition3D<T,DESCRIPTOR> {
public:
  BoundaryConditionInstantiator3D( BlockLatticeStructure3D<T,DESCRIPTOR>& block_ );
  ~BoundaryConditionInstantiator3D() override;

  void addVelocityBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addVelocityBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addVelocityBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addVelocityBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addVelocityBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addVelocityBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;

  void addSlipBoundary(int x0, int x1, int y0, int y1, int z0, int z1, int discreteNormalX, int discreteNormalY, int discreteNormalZ);
  void addPartialSlipBoundary(T tuner, int x0, int x1, int y0, int y1, int z0, int z1, int discreteNormalX, int discreteNormalY, int discreteNormalZ);

  void addWallFunctionBoundary(int x0, int x1, int y0, int y1, int z0, int z1, BlockGeometryStructure3D<T>& blockGeometryStructure,
                               std::vector<int> discreteNormal, std::vector<int> missingIndices,
                               UnitConverter<T, DESCRIPTOR> const& converter, wallFunctionParam<T> const& wallFunctionParam,
                               IndicatorF3D<T>* geoIndicator);

  void addPressureBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addPressureBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addPressureBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addPressureBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addPressureBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addPressureBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;

  void addConvectionBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) override;

  void addExternalVelocityEdge0NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge0NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge0PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge0PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge1NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge1NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge1PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge1PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge2NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge2NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge2PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addExternalVelocityEdge2PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;

  void addInternalVelocityEdge0NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge0NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge0PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge0PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge1NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge1NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge1PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge1PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge2NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge2NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge2PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;
  void addInternalVelocityEdge2PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) override;

  void addExternalVelocityCornerNNN(int x, int y, int z, T omega) override;
  void addExternalVelocityCornerNNP(int x, int y, int z, T omega) override;
  void addExternalVelocityCornerNPN(int x, int y, int z, T omega) override;
  void addExternalVelocityCornerNPP(int x, int y, int z, T omega) override;
  void addExternalVelocityCornerPNN(int x, int y, int z, T omega) override;
  void addExternalVelocityCornerPNP(int x, int y, int z, T omega) override;
  void addExternalVelocityCornerPPN(int x, int y, int z, T omega) override;
  void addExternalVelocityCornerPPP(int x, int y, int z, T omega) override;

  void addInternalVelocityCornerNNN(int x, int y, int z, T omega) override;
  void addInternalVelocityCornerNNP(int x, int y, int z, T omega) override;
  void addInternalVelocityCornerNPN(int x, int y, int z, T omega) override;
  void addInternalVelocityCornerNPP(int x, int y, int z, T omega) override;
  void addInternalVelocityCornerPNN(int x, int y, int z, T omega) override;
  void addInternalVelocityCornerPNP(int x, int y, int z, T omega) override;
  void addInternalVelocityCornerPPN(int x, int y, int z, T omega) override;
  void addInternalVelocityCornerPPP(int x, int y, int z, T omega) override;

  void addVelocityBoundary(BlockIndicatorF3D<T>& indicator,
                           int x0, int x1, int y0, int y1, int z0, int z1,
                           T omega) override;
  void addSlipBoundary(BlockIndicatorF3D<T>& indicator,
                       int x0, int x1, int y0, int y1, int z0, int z1) override;
  void addPartialSlipBoundary(T tuner, BlockIndicatorF3D<T>& indicator,
                              int x0, int x1, int y0, int y1, int z0, int z1) override;
  void addPressureBoundary(BlockIndicatorF3D<T>& indicator,
                           int x0, int x1, int y0, int y1, int z0, int z1,
                           T omega) override;
  void addConvectionBoundary(BlockIndicatorF3D<T>& indicator,
                             int x0, int x1, int y0, int y1, int z0, int z1,
                             T omega, T* uAv=NULL) override;
  void addWallFunctionBoundary(BlockIndicatorF3D<T>& indicator,
                               int x0, int x1, int y0, int y1, int z0, int z1,
                               UnitConverter<T, DESCRIPTOR> const& converter,
                               wallFunctionParam<T> const& wallFunctionParam,
                               IndicatorF3D<T>* geoIndicator=NULL) override;
  void addFreeEnergyWallBoundary(BlockIndicatorF3D<T>& indicator,
                                 int x0, int x1, int y0, int y1, int z0, int z1,
                                 T addend, int latticeNumber) override;
  void addFreeEnergyInletBoundary(BlockIndicatorF3D<T>& indicator,
                                 int x0, int x1, int y0, int y1, int z0, int z1,
                                 T omega, std::string type, int latticeNumber) override;
  void addFreeEnergyOutletBoundary(BlockIndicatorF3D<T>& indicator,
                                 int x0, int x1, int y0, int y1, int z0, int z1,
                                 T omega, std::string type, int latticeNumber) override;

  void outputOn() override;
  void outputOff() override;

private:
  template<int direction, int orientation>
  void addVelocityBoundary(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  template<int direction, int orientation>
  void addPressureBoundary(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  template<int direction, int orientation>
  void addConvectionBoundary(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL);
  template<int plane, int normal1, int normal2>
  void addExternalVelocityEdge(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  template<int plane, int normal1, int normal2>
  void addInternalVelocityEdge(int x0, int x1, int y0, int y1, int z0, int z1, T omega);
  template<int normalX, int normalY, int normalZ>
  void addExternalVelocityCorner(int x, int y, int z, T omega);
  template<int normalX, int normalY, int normalZ>
  void addInternalVelocityCorner(int x, int y, int z, T omega);
private:
  BlockLatticeStructure3D<T,DESCRIPTOR>& _block;
  std::vector<Momenta<T,DESCRIPTOR>*>  momentaVector;
  std::vector<Dynamics<T,DESCRIPTOR>*> dynamicsVector;
  bool _output;
  mutable OstreamManager clout;
};


///////// class BoundaryConditionInstantiator3D ////////////////////////

template<typename T, typename DESCRIPTOR, class BoundaryManager>
BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::BoundaryConditionInstantiator3D (
  BlockLatticeStructure3D<T,DESCRIPTOR>& block)
  : _block(block), _output(false), clout(std::cout,"BoundaryConditionInstantiator3D")
{ }

template<typename T, typename DESCRIPTOR, class BoundaryManager>
BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::~BoundaryConditionInstantiator3D()
{
  for (auto &iDynamics : dynamicsVector) {
    delete iDynamics;
  }
  for (auto &iMomenta : momentaVector) {
    delete iMomenta;
  }
}

// Velocity BC

template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int direction, int orientation>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addVelocityBoundary(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  OLB_PRECONDITION( x0==x1 || y0==y1 || z0==z1 );

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        Momenta<T,DESCRIPTOR>* momenta
          = BoundaryManager::template getVelocityBoundaryMomenta<direction,orientation>();
        Dynamics<T,DESCRIPTOR>* dynamics
          = BoundaryManager::template getVelocityBoundaryDynamics<direction,orientation>(omega, *momenta);
        _block.defineDynamics(iX,iX,iY,iY,iZ,iZ, dynamics);
        momentaVector.push_back(momenta);
        dynamicsVector.push_back(dynamics);
        if (_output) {
          clout << "addVelocityBoundary<" << direction << ", " << orientation << ">(" << iX << ", " << iX << ", "<< iY << ", " << iY << ", " << iZ << ", " << iZ << ", "<< omega << " )" << std::endl;
        }
      }
    }
  }

  PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor
    = BoundaryManager::template getVelocityBoundaryProcessor<direction,orientation>(x0,x1, y0,y1, z0,z1);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}

// Slip BC

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addSlipBoundary(
  int x0, int x1, int y0, int y1, int z0, int z1, int discreteNormalX, int discreteNormalY, int discreteNormalZ)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);

  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if (_output) {
          clout << "addSlipBoundary<" << discreteNormalX << ","<< discreteNormalY << ","<< discreteNormalZ << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << z0 << ", " << z1 << " )" << std::endl;
        }
      }
    }
  }

  PostProcessorGenerator3D<T, DESCRIPTOR>* postProcessor = new SlipBoundaryProcessorGenerator3D<T, DESCRIPTOR>(x0, x1, y0, y1, z0, z1, discreteNormalX, discreteNormalY, discreteNormalZ);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}

// Partial slip BC

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addPartialSlipBoundary(
  T tuner, int x0, int x1, int y0, int y1, int z0, int z1, int discreteNormalX, int discreteNormalY, int discreteNormalZ)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);

  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if (_output) {
          clout << "addPartialSlipBoundary<" << discreteNormalX << ","<< discreteNormalY << ","<< discreteNormalZ << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << z0 << ", " << z1 << " )" << std::endl;
        }
      }
    }
  }

  PostProcessorGenerator3D<T, DESCRIPTOR>* postProcessor = new PartialSlipBoundaryProcessorGenerator3D<T, DESCRIPTOR>(tuner, x0, x1, y0, y1, z0, z1, discreteNormalX, discreteNormalY, discreteNormalZ);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}

// Wall Function BC

namespace {


template<typename T, typename DESCRIPTOR, class BoundaryManager, typename Enable = void>
struct WallFunctionBoundaryProcessorGenerator3DInstatiator {
  static PostProcessorGenerator3D<T, DESCRIPTOR>* getWallFunctionBoundaryProcessorGenerator(int x0, int x1, int y0, int y1, int z0, int z1, BlockGeometryStructure3D<T>&, std::vector<int>& discreteNormal, const std::vector<int>& missingIndices, UnitConverter<T, DESCRIPTOR> const&, wallFunctionParam<T> const& wallFunctionParam, IndicatorF3D<T>* geoIndicator)
  {
    std::abort();
  }
};

template<typename T, typename DESCRIPTOR, class BoundaryManager>
struct WallFunctionBoundaryProcessorGenerator3DInstatiator<
  T, DESCRIPTOR, BoundaryManager,
  typename std::enable_if<
    DESCRIPTOR::template provides<descriptors::AV_SHEAR>() &&
    DESCRIPTOR::template provides<descriptors::TAU_W>()    &&
    DESCRIPTOR::template provides<descriptors::TAU_EFF>()  &&
    DESCRIPTOR::template provides<descriptors::FORCE>()    &&
    DESCRIPTOR::template provides<descriptors::V12>()
  >::type
> {
  static PostProcessorGenerator3D<T, DESCRIPTOR>* getWallFunctionBoundaryProcessorGenerator(int x0, int x1, int y0, int y1, int z0, int z1, BlockGeometryStructure3D<T >& blockGeometryStructure, std::vector<int>& discreteNormal, const std::vector<int>& missingIndices, UnitConverter<T, DESCRIPTOR> const& converter, wallFunctionParam<T> const& wallFunctionParam, IndicatorF3D<T>* geoIndicator)
  {
    PostProcessorGenerator3D<T, DESCRIPTOR>* postProcessor = new WallFunctionBoundaryProcessorGenerator3D<T, DESCRIPTOR>(x0, x1, y0, y1, z0, z1, blockGeometryStructure, discreteNormal, missingIndices,
        converter, wallFunctionParam, geoIndicator);

    return postProcessor;
  }
};

}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addWallFunctionBoundary(
  int x0, int x1, int y0, int y1, int z0, int z1, BlockGeometryStructure3D<T>& blockGeometryStructure,
  std::vector<int> discreteNormal, std::vector<int> missingIndices,
  UnitConverter<T, DESCRIPTOR> const& converter, wallFunctionParam<T> const& wallFunctionParam, IndicatorF3D<T>* geoIndicator)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);

  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if (_output) {
          clout << "addWallFunctionBoundary<" << discreteNormal[0] << ","<< discreteNormal[1] << ","<< discreteNormal[2] << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << z0 << ", " << z1 << " )" << std::endl;
        }
      }
    }
  }

  PostProcessorGenerator3D<T, DESCRIPTOR>* postProcessor = WallFunctionBoundaryProcessorGenerator3DInstatiator<T, DESCRIPTOR, BoundaryManager>::getWallFunctionBoundaryProcessorGenerator(x0, x1, y0, y1, z0, z1, blockGeometryStructure, discreteNormal, missingIndices, converter, wallFunctionParam, geoIndicator);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}

// Pressure BC

template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int direction, int orientation>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addPressureBoundary(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  OLB_PRECONDITION( x0==x1 || y0==y1 || z0==z1 );

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        Momenta<T,DESCRIPTOR>* momenta
          = BoundaryManager::template getPressureBoundaryMomenta<direction,orientation>();
        Dynamics<T,DESCRIPTOR>* dynamics
          = BoundaryManager::template getPressureBoundaryDynamics<direction,orientation>(omega, *momenta);
        _block.defineDynamics(iX,iX,iY,iY,iZ,iZ, dynamics);
        momentaVector.push_back(momenta);
        dynamicsVector.push_back(dynamics);
        if (_output) {
          clout << "addPressureBoundary<" << direction << ", " << orientation << ">(" << iX << ", " << iX << ", "<< iY << ", " << iY << ", " << iZ << ", " << iZ << ", "<< omega << " )" << std::endl;
        }
      }
    }
  }

  PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor
    = BoundaryManager::template getPressureBoundaryProcessor<direction,orientation>(x0,x1, y0,y1, z0,z1);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}

// Convection BC

template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int direction, int orientation>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addConvectionBoundary(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv)
{
  OLB_PRECONDITION( x0==x1 || y0==y1 || z0==z1 );

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        if (_output) {
          clout << "addConvectionBoundary<" << direction << ", " << orientation << ">(" << iX << ", " << iX << ", "<< iY << ", " << iY << ", " << iZ << ", " << iZ << ", "<< omega << " )" << std::endl;
        }
      }
    }
  }

  PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor
    = BoundaryManager::template getConvectionBoundaryProcessor<direction,orientation>(x0,x1, y0,y1, z0,z1, uAv);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int plane, int normal1, int normal2>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityEdge(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  OLB_PRECONDITION(
    ( x0==x1 && y0==y1 ) ||
    ( x0==x1 && z0==z1 ) ||
    ( y0==y1 && z0==z1 ) );

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        Momenta<T,DESCRIPTOR>* momenta
          = BoundaryManager::template getExternalVelocityEdgeMomenta<plane,normal1,normal2>();
        Dynamics<T,DESCRIPTOR>* dynamics
          = BoundaryManager::template getExternalVelocityEdgeDynamics<plane,normal1,normal2>(omega, *momenta);
        _block.defineDynamics(iX,iX,iY,iY,iZ,iZ, dynamics);
        momentaVector.push_back(momenta);
        dynamicsVector.push_back(dynamics);
        if (_output) {
          clout << "addExternalVelocityEdge<" << plane << ", " << normal1 << ", " << normal2 << ">(" << iX << ", " << iX << ", "<< iY << ", " << iY << ", " << iZ << ", " << iZ << ", "<< omega << " )" << std::endl;
        }
      }
    }
  }

  PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor
    = BoundaryManager::template getExternalVelocityEdgeProcessor<plane,normal1,normal2>(x0,x1, y0,y1, z0,z1);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int plane, int normal1, int normal2>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityEdge(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  if (!(( x0==x1 && y0==y1 ) ||
        ( x0==x1 && z0==z1 ) ||
        ( y0==y1 && z0==z1 ) )) {
    clout << x0 <<" "<< x1 <<" "<< y0 <<" "<< y1 <<" "<< z0 <<" "<< z1 << std::endl;
  }

  OLB_PRECONDITION(
    ( x0==x1 && y0==y1 ) ||
    ( x0==x1 && z0==z1 ) ||
    ( y0==y1 && z0==z1 ) );

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        Momenta<T,DESCRIPTOR>* momenta
          = BoundaryManager::template getInternalVelocityEdgeMomenta<plane,normal1,normal2>();
        Dynamics<T,DESCRIPTOR>* dynamics
          = BoundaryManager::template getInternalVelocityEdgeDynamics<plane,normal1,normal2>(omega, *momenta);
        _block.defineDynamics(iX,iX,iY,iY,iZ,iZ, dynamics);
        momentaVector.push_back(momenta);
        dynamicsVector.push_back(dynamics);
        if (_output) {
          clout << "addInternalVelocityEdge<" << plane << ", " << normal1 << ", " << normal2 << ">(" << iX << ", " << iX << ", "<< iY << ", " << iY << ", " << iZ << ", " << iZ << ", "<< omega << " )" << std::endl;
        }
      }
    }
  }

  PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor
    = BoundaryManager::template getInternalVelocityEdgeProcessor<plane,normal1,normal2>(x0,x1, y0,y1, z0,z1);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}


template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int xNormal, int yNormal, int zNormal>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityCorner(int x, int y, int z, T omega)
{
  Momenta<T,DESCRIPTOR>* momenta
    = BoundaryManager::template getExternalVelocityCornerMomenta<xNormal,yNormal,zNormal>();
  Dynamics<T,DESCRIPTOR>* dynamics
    = BoundaryManager::template getExternalVelocityCornerDynamics<xNormal,yNormal,zNormal>(omega, *momenta);

  _block.defineDynamics(x,x,y,y,z,z, dynamics);

  momentaVector.push_back(momenta);
  dynamicsVector.push_back(dynamics);

  PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor
    = BoundaryManager::template getExternalVelocityCornerProcessor<xNormal,yNormal,zNormal>(x, y, z);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
  if (_output) {
    clout << "addExternalVelocityCorner<" << xNormal << ", " << yNormal << ", " << zNormal << ">(" << x << ", " << y << ", "<< z << omega << " )" << std::endl;
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int xNormal, int yNormal, int zNormal>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityCorner(int x, int y, int z, T omega)
{
  Momenta<T,DESCRIPTOR>* momenta
    = BoundaryManager::template getInternalVelocityCornerMomenta<xNormal,yNormal,zNormal>();
  Dynamics<T,DESCRIPTOR>* dynamics
    = BoundaryManager::template getInternalVelocityCornerDynamics<xNormal,yNormal,zNormal>(omega, *momenta);

  _block.defineDynamics(x,x,y,y,z,z, dynamics);

  momentaVector.push_back(momenta);
  dynamicsVector.push_back(dynamics);

  PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor
    = BoundaryManager::template getInternalVelocityCornerProcessor<xNormal,yNormal,zNormal>(x, y, z);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
  if (_output) {
    clout << "addInternalVelocityCorner<" << xNormal << ", " << yNormal << ", " << zNormal << ">(" << x << ", " << y << ", "<< z << omega << " )" << std::endl;
  }
}

// Velocity BC

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addVelocityBoundary(BlockIndicatorF3D<T>& indicator, int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  std::vector<int> discreteNormal(4,0);
  for (int iX = x0; iX <= x1; iX++) {
    for (int iY = y0; iY <= y1; iY++) {
      for (int iZ = z0; iZ <= z1; iZ++) {
        if (indicator(iX, iY, iZ)) {
          discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX, iY, iZ);
          if (discreteNormal[0] == 0) {
            if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {
              addVelocityBoundary<0,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
              addVelocityBoundary<0,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
              addVelocityBoundary<1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
              addVelocityBoundary<1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
              addVelocityBoundary<2,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
              addVelocityBoundary<2,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
          }

          else if (discreteNormal[0] == 1) {
            if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              addExternalVelocityCorner<1,1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              addExternalVelocityCorner<1,-1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              addExternalVelocityCorner<1,1,-1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              addExternalVelocityCorner<1,-1,-1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              addExternalVelocityCorner<-1,1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              addExternalVelocityCorner<-1,-1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              addExternalVelocityCorner<-1,1,-1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              addExternalVelocityCorner<-1,-1,-1>(iX,iY,iZ, omega);
            }
            ///                     addExternalVelocityCorner<discreteNormal[1],discreteNormal[2],discreteNormal[3]>(iX,iY,iZ, omega);
          }

          else if (discreteNormal[0] == 2) {
            if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              addInternalVelocityCorner<1,1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              addExternalVelocityCorner<1,-1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              addInternalVelocityCorner<1,1,-1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              addInternalVelocityCorner<1,-1,-1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              addInternalVelocityCorner<-1,1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              addInternalVelocityCorner<-1,-1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              addInternalVelocityCorner<-1,1,-1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              addInternalVelocityCorner<-1,-1,-1>(iX,iY,iZ, omega);
            }
            ///                     addInternalVelocityCorner<discreteNormal[1],discreteNormal[2],discreteNormal[3]>(iX,iY,iZ, omega);
          }

          else if (discreteNormal[0] == 3) {
            if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              addExternalVelocityEdge<0,1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              addExternalVelocityEdge<0,-1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              addExternalVelocityEdge<0,1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              addExternalVelocityEdge<0,-1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              addExternalVelocityEdge<1,1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              addExternalVelocityEdge<1,1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              addExternalVelocityEdge<1,-1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              addExternalVelocityEdge<1,-1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              addExternalVelocityEdge<2,1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              addExternalVelocityEdge<2,-1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              addExternalVelocityEdge<2,1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              addExternalVelocityEdge<2,-1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
          }

          else if (discreteNormal[0] == 4) {
            if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              addInternalVelocityEdge<0,1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              addInternalVelocityEdge<0,-1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              addInternalVelocityEdge<0,1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              addInternalVelocityEdge<0,-1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              addInternalVelocityEdge<1,1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              addInternalVelocityEdge<1,1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              addInternalVelocityEdge<1,-1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              addInternalVelocityEdge<1,-1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              addInternalVelocityEdge<2,1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              addInternalVelocityEdge<2,-1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              addInternalVelocityEdge<2,1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              addInternalVelocityEdge<2,-1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
          }
        }
      }
    }
  }
}

// Slip BC

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addSlipBoundary(
  BlockIndicatorF3D<T>& indicator,
  int x0, int x1, int y0, int y1, int z0, int z1)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  std::vector<int> discreteNormal(4, 0);
  for (int iX = x0; iX <= x1; iX++) {
    for (int iY = y0; iY <= y1; iY++) {
      for (int iZ = z0; iZ <= z1; iZ++) {
        if (indicator(iX, iY, iZ)) {
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);
          if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {
            addSlipBoundary(iX, iX, iY, iY, iZ, iZ,
                            discreteNormal[1], discreteNormal[2], discreteNormal[3]);
          }
          else {
            clout << "Warning: Could not addSlipBoundary (" << iX << ", " << iY << ", " << iZ << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] <<","<< discreteNormal[3] <<"), set to bounceBack" << std::endl;
            _block.defineDynamics(iX, iY, iZ, &instances::getBounceBack<T, DESCRIPTOR>());
          }
        }
      }
    }
  }
}

// Partial slip BC

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addPartialSlipBoundary(
  T tuner, BlockIndicatorF3D<T>& indicator,
  int x0, int x1, int y0, int y1, int z0, int z1)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  std::vector<int> discreteNormal(4, 0);
  for (int iX = x0; iX <= x1; iX++) {
    for (int iY = y0; iY <= y1; iY++) {
      for (int iZ = z0; iZ <= z1; iZ++) {
        if (indicator(iX, iY, iZ)) {
          if (tuner < 0. || tuner > 1.) {
            clout << "Warning: Could not addPartialSlipBoundary (" << iX << ", " << iY << ", " << iZ << "), tuner must be between 0.1 and instead is=" << tuner <<", set to bounceBack" << std::endl;
            _block.defineDynamics(iX, iY, iZ, &instances::getBounceBack<T, DESCRIPTOR>());
          } else {
            discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);
            if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {
              addPartialSlipBoundary(tuner, iX, iX, iY, iY, iZ, iZ,
                                     discreteNormal[1], discreteNormal[2], discreteNormal[3]);
            }
            else {
              clout << "Warning: Could not addPartialSlipBoundary (" << iX << ", " << iY << ", " << iZ << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] <<","<< discreteNormal[3] <<"), set to bounceBack" << std::endl;
              _block.defineDynamics(iX, iY, iZ, &instances::getBounceBack<T, DESCRIPTOR>());
            }
          }
        }
      }
    }
  }
}

// Pressure BC

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addPressureBoundary(BlockIndicatorF3D<T>& indicator,
                    int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  std::vector<int> discreteNormal(4,0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if (indicator(iX, iY, iZ)) {
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);

          if (discreteNormal[0] == 0) {
            if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {
              addPressureBoundary<0,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
              addPressureBoundary<0,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
              addPressureBoundary<1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
              addPressureBoundary<1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
              addPressureBoundary<2,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
              addPressureBoundary<2,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
          }
        }
      }
    }
  }
}

// Convection BC

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addConvectionBoundary(BlockIndicatorF3D<T>& indicator,
                      int x0, int x1, int y0, int y1, int z0, int z1,
                      T omega, T* uAv)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  std::vector<int> discreteNormal(4,0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if (indicator(iX, iY, iZ)) {
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);

          if (discreteNormal[0] == 0) {
            if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {
              addConvectionBoundary<0,-1>(iX,iX,iY,iY,iZ,iZ, omega, uAv);
            }
            else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
              addConvectionBoundary<0,1>(iX,iX,iY,iY,iZ,iZ, omega, uAv);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
              addConvectionBoundary<1,-1>(iX,iX,iY,iY,iZ,iZ, omega, uAv);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
              addConvectionBoundary<1,1>(iX,iX,iY,iY,iZ,iZ, omega, uAv);
            }

            else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
              addConvectionBoundary<2,-1>(iX,iX,iY,iY,iZ,iZ, omega, uAv);
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
              addConvectionBoundary<2,1>(iX,iX,iY,iY,iZ,iZ, omega, uAv);
            }
          }
        }
      }
    }
  }
}

// Wall Function BC

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addWallFunctionBoundary(
  BlockIndicatorF3D<T>& indicator,
  int x0, int x1, int y0, int y1, int z0, int z1,
  UnitConverter<T, DESCRIPTOR> const& converter,
  wallFunctionParam<T> const& wallFunctionParam,
  IndicatorF3D<T>* geoIndicator)
{
  const auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  std::vector<int> discreteNormal(4, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if (indicator(iX, iY, iZ)) {
          discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX, iY, iZ);
          std::vector<int> missingIndices;
          for (int x = -1 ; x < 2; ++x) {
            for (int y = -1 ; y < 2; ++y) {
              for (int z = -1 ; z < 2; ++z) {
                if (blockGeometryStructure.getMaterial(iX + x, iY + y, iZ + z) == 0) {
                  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
                    if (descriptors::c<DESCRIPTOR>(iPop,0) == x &&
                        descriptors::c<DESCRIPTOR>(iPop,1) == y &&
                        descriptors::c<DESCRIPTOR>(iPop,2) == z) {
                      missingIndices.push_back(descriptors::opposite<DESCRIPTOR>(iPop));
                    }
                  }
                }
              }
            }
          }
          if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {
            discreteNormal.erase(discreteNormal.begin());
            addWallFunctionBoundary(iX, iX, iY, iY, iZ, iZ,
                                    indicator.getBlockGeometryStructure(),
                                    discreteNormal, missingIndices,
                                    converter, wallFunctionParam, geoIndicator);
          }
          else {
            //clout << "Warning: Could not add WallFunction (" << iX << ", " << iY << ", " << iZ << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] <<","<< discreteNormal[3] <<"), set to bounceBack" << std::endl;
            _block.defineDynamics(iX, iY, iZ, &instances::getBounceBack<T, DESCRIPTOR>());
          }
        }
      }
    }
  }
}

// Free Energy BC

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addFreeEnergyWallBoundary(
  BlockIndicatorF3D<T>& indicator,
  int x0, int x1, int y0, int y1, int z0, int z1, T addend, int latticeNumber)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  std::vector<int> discreteNormal(4, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if(indicator(iX,iY,iZ)) {
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ, true);
          if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {

            Dynamics<T, DESCRIPTOR>* dynamics = NULL;
            if (latticeNumber == 1) {
              dynamics = &instances::getBounceBack<T, DESCRIPTOR>();
            } else {
              dynamics = new FreeEnergyWallDynamics<T, DESCRIPTOR>;
              dynamicsVector.push_back(dynamics);
            }
            _block.get(iX,iY,iZ).defineDynamics(dynamics);

            PostProcessorGenerator3D<T, DESCRIPTOR>* wettingPostProcessor =
                new FreeEnergyWallProcessorGenerator3D<T, DESCRIPTOR> ( iX, iX, iY, iY, iZ, iZ,
                    discreteNormal[1], discreteNormal[2], discreteNormal[3], addend );
            PostProcessorGenerator3D<T, DESCRIPTOR>* chemPotPostProcessor =
                new FreeEnergyChemPotBoundaryProcessorGenerator3D<T, DESCRIPTOR> ( iX, iX, iY, iY, iZ, iZ,
                    discreteNormal[1], discreteNormal[2], discreteNormal[3], latticeNumber );
            if (wettingPostProcessor) {
              _block.addPostProcessor(*wettingPostProcessor);
            }
            if (chemPotPostProcessor) {
              _block.addPostProcessor(*chemPotPostProcessor);
            }
          }
          if (_output) {
            clout << "addFreeEnergyWallBoundary<" << "," << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << z0 << ", " << z1 << ")" << std::endl;
          }
        }
      }
    }
  }

}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addFreeEnergyInletBoundary(
  BlockIndicatorF3D<T>& indicator,
  int x0, int x1, int y0, int y1, int z0, int z1,
  T omega, std::string type, int latticeNumber)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  std::vector<int> discreteNormal(4, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if(indicator(iX,iY,iZ)) {
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ, true);

          if (discreteNormal[0] == 0) {
            Momenta<T, DESCRIPTOR>* momenta = NULL;
            Dynamics<T, DESCRIPTOR>* dynamics = NULL;

            if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {
              if (latticeNumber == 1) {
                if (type == "density") {
                  addPressureBoundary<0,-1>(iX, iX, iY, iY, iZ, iZ, omega);
                } else {
                  addVelocityBoundary<0,-1>(iX, iX, iY, iY, iZ, iZ, omega);
                }
              } else {
                momenta = BoundaryManager::template 
                  getPressureBoundaryMomenta<0,-1>();
                dynamics = new FreeEnergyInletOutletDynamics<T,DESCRIPTOR,0,-1>(omega,*momenta);
              }
            }

            else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
              if (latticeNumber == 1) {
                if (type == "density") {
                  addPressureBoundary<0,1>(iX, iX, iY, iY, iZ, iZ, omega);
                } else {
                  addVelocityBoundary<0,1>(iX, iX, iY, iY, iZ, iZ, omega);
                }
              } else {
                momenta = BoundaryManager::template 
                  getPressureBoundaryMomenta<0,1>();
                dynamics = new FreeEnergyInletOutletDynamics<T,DESCRIPTOR,0,1>(omega,*momenta);
              }
            }

            else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
              if (latticeNumber == 1) {
                if (type == "density") {
                  addPressureBoundary<1,-1>(iX, iX, iY, iY, iZ, iZ, omega);
                } else {
                  addVelocityBoundary<1,-1>(iX, iX, iY, iY, iZ, iZ, omega);
                }
              } else {
                momenta = BoundaryManager::template 
                  getPressureBoundaryMomenta<1,-1>();
                dynamics = new FreeEnergyInletOutletDynamics<T,DESCRIPTOR,1,-1>(omega,*momenta);
              }
            }

            else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
              if (latticeNumber == 1) {
                if (type == "density") {
                  addPressureBoundary<1,1>(iX, iX, iY, iY, iZ, iZ, omega);
                } else {
                  addVelocityBoundary<1,1>(iX, iX, iY, iY, iZ, iZ, omega);
                }
              } else {
                momenta = BoundaryManager::template 
                  getPressureBoundaryMomenta<1,1>();
                dynamics = new FreeEnergyInletOutletDynamics<T,DESCRIPTOR,1,1>(omega,*momenta);
              }
              momenta = BoundaryManager::template 
                getPressureBoundaryMomenta<1,1>();
              if (latticeNumber == 1) {
                dynamics = BoundaryManager::template
                  getPressureBoundaryDynamics<1,1>(omega, *momenta);
              } else {
                dynamics = new FreeEnergyInletOutletDynamics<T,DESCRIPTOR,1,1>(omega,*momenta);
              }
            }

            else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
              if (latticeNumber == 1) {
                if (type == "density") {
                  addPressureBoundary<2,-1>(iX, iX, iY, iY, iZ, iZ, omega);
                } else {
                  addVelocityBoundary<2,-1>(iX, iX, iY, iY, iZ, iZ, omega);
                }
              } else {
                momenta = BoundaryManager::template 
                  getPressureBoundaryMomenta<2,-1>();
                dynamics = new FreeEnergyInletOutletDynamics<T,DESCRIPTOR,2,-1>(omega,*momenta);
              }
            }

            else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
              if (latticeNumber == 1) {
                if (type == "density") {
                  addPressureBoundary<2,1>(iX, iX, iY, iY, iZ, iZ, omega);
                } else {
                  addVelocityBoundary<2,1>(iX, iX, iY, iY, iZ, iZ, omega);
                }
              } else {
                momenta = BoundaryManager::template 
                  getPressureBoundaryMomenta<2,1>();
                dynamics = new FreeEnergyInletOutletDynamics<T,DESCRIPTOR,2,1>(omega,*momenta);
              }
            }

            if (latticeNumber != 1) {
              _block.defineDynamics(iX,iX,iY,iY,iZ,iZ, dynamics);
              momentaVector.push_back(momenta);
              dynamicsVector.push_back(dynamics);
            }
          }

          PostProcessorGenerator3D<T, DESCRIPTOR>* chemPotPostProcessor =
              new FreeEnergyChemPotBoundaryProcessorGenerator3D<T, DESCRIPTOR> ( iX, iX, iY, iY, iZ, iZ,
                  discreteNormal[1], discreteNormal[2], discreteNormal[3], latticeNumber );
          if (chemPotPostProcessor) {
            _block.addPostProcessor(*chemPotPostProcessor);
          }

          if (_output) {
            clout << "addFreeEnergyInletBoundary<" << "," << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << z0 << ", " << z1 << ")" << std::endl;
          }
        }
      }
    }
  }

}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T, DESCRIPTOR, BoundaryManager>::addFreeEnergyOutletBoundary(
  BlockIndicatorF3D<T>& indicator, int x0, int x1, int y0, int y1, int z0, int z1,
  T omega, std::string type, int latticeNumber)
{
  addFreeEnergyInletBoundary(indicator, x0, x1, y0, y1, z0, z1, omega, type, latticeNumber);

  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  std::vector<int> discreteNormal(4, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if(indicator(iX,iY,iZ)) {
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ, true);
          if (discreteNormal[0] == 0) {
            
            PostProcessorGenerator3D<T, DESCRIPTOR>* convectivePostProcessor =
              new FreeEnergyConvectiveProcessorGenerator3D<T, DESCRIPTOR> (
                iX, iX, iY, iY, iZ, iZ, discreteNormal[1], discreteNormal[2], discreteNormal[3] );
            if (convectivePostProcessor) {
              _block.addPostProcessor(*convectivePostProcessor);
            }

            if (_output) {
              clout << "addFreeEnergyOutletBoundary<" << "," << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << z0 << ", " << z1 << ")" << std::endl;
            }
          }
        }
      }
    }
  }

}


template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addVelocityBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addVelocityBoundary<0,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addVelocityBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addVelocityBoundary<0,1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addVelocityBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addVelocityBoundary<1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addVelocityBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addVelocityBoundary<1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addVelocityBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addVelocityBoundary<2,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addVelocityBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addVelocityBoundary<2, 1>(x0,x1,y0,y1,z0,z1, omega);
}

// Pressure BC

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addPressureBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addPressureBoundary<0,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addPressureBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addPressureBoundary<0,1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addPressureBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addPressureBoundary<1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addPressureBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addPressureBoundary<1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addPressureBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addPressureBoundary<2,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addPressureBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addPressureBoundary<2, 1>(x0,x1,y0,y1,z0,z1, omega);
}

// Convection BC

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addConvectionBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv)
{
  addConvectionBoundary<0,-1>(x0,x1,y0,y1,z0,z1, omega, uAv);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addConvectionBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv)
{
  addConvectionBoundary<0,1>(x0,x1,y0,y1,z0,z1, omega, uAv);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addConvectionBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv)
{
  addConvectionBoundary<1,-1>(x0,x1,y0,y1,z0,z1, omega, uAv);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addConvectionBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv)
{
  addConvectionBoundary<1, 1>(x0,x1,y0,y1,z0,z1, omega, uAv);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addConvectionBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv)
{
  addConvectionBoundary<2,-1>(x0,x1,y0,y1,z0,z1, omega, uAv);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addConvectionBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv)
{
  addConvectionBoundary<2, 1>(x0,x1,y0,y1,z0,z1, omega, uAv);
}


// Velocity BC

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityEdge0NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<0,-1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityEdge0NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<0,-1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityEdge0PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<0, 1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityEdge0PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<0, 1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityEdge1NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<1,-1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityEdge1NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<1,-1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityEdge1PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<1, 1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityEdge1PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<1, 1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityEdge2NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<2,-1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityEdge2NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<2,-1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityEdge2PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<2, 1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityEdge2PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addExternalVelocityEdge<2, 1, 1>(x0,x1,y0,y1,z0,z1, omega);
}



template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityEdge0NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<0,-1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityEdge0NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<0,-1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityEdge0PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<0, 1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityEdge0PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<0, 1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityEdge1NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<1,-1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityEdge1NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<1,-1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityEdge1PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<1, 1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityEdge1PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<1, 1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityEdge2NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<2,-1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityEdge2NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<2,-1, 1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityEdge2PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<2, 1,-1>(x0,x1,y0,y1,z0,z1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityEdge2PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  addInternalVelocityEdge<2, 1, 1>(x0,x1,y0,y1,z0,z1, omega);
}


template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityCornerNNN(int x, int y, int z, T omega)
{
  addExternalVelocityCorner<-1,-1,-1>(x,y,z, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityCornerNNP(int x, int y, int z, T omega)
{
  addExternalVelocityCorner<-1,-1, 1>(x,y,z, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityCornerNPN(int x, int y, int z, T omega)
{
  addExternalVelocityCorner<-1, 1,-1>(x,y,z, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityCornerNPP(int x, int y, int z, T omega)
{
  addExternalVelocityCorner<-1, 1, 1>(x,y,z, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityCornerPNN(int x, int y, int z, T omega)
{
  addExternalVelocityCorner< 1,-1,-1>(x,y,z, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityCornerPNP(int x, int y, int z, T omega)
{
  addExternalVelocityCorner< 1,-1, 1>(x,y,z, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityCornerPPN(int x, int y, int z, T omega)
{
  addExternalVelocityCorner< 1, 1,-1>(x,y,z, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExternalVelocityCornerPPP(int x, int y, int z, T omega)
{
  addExternalVelocityCorner< 1, 1, 1>(x,y,z, omega);
}


template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityCornerNNN(int x, int y, int z, T omega)
{
  addInternalVelocityCorner<-1,-1,-1>(x,y,z, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityCornerNNP(int x, int y, int z, T omega)
{
  addInternalVelocityCorner<-1,-1, 1>(x,y,z, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityCornerNPN(int x, int y, int z, T omega)
{
  addInternalVelocityCorner<-1, 1,-1>(x,y,z, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityCornerNPP(int x, int y, int z, T omega)
{
  addInternalVelocityCorner<-1, 1, 1>(x,y,z, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityCornerPNN(int x, int y, int z, T omega)
{
  addInternalVelocityCorner< 1,-1,-1>(x,y,z, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityCornerPNP(int x, int y, int z, T omega)
{
  addInternalVelocityCorner< 1,-1, 1>(x,y,z, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityCornerPPN(int x, int y, int z, T omega)
{
  addInternalVelocityCorner< 1, 1,-1>(x,y,z, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addInternalVelocityCornerPPP(int x, int y, int z, T omega)
{
  addInternalVelocityCorner< 1, 1, 1>(x,y,z, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::outputOn()
{
  _output = true;
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::outputOff()
{
  _output = false;
}


}  // namespace olb

#endif
