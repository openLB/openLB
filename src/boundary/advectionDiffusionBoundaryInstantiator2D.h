/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
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

#ifndef ADVECTION_DIFFUSION_BOUNDARY_INSTANTIATOR_2D_H
#define ADVECTION_DIFFUSION_BOUNDARY_INSTANTIATOR_2D_H

#include "advectionDiffusionBoundaryCondition2D.h"
#include "advectionDiffusionBoundaryCondition2D.hh"
#include "functors/lattice/indicator/blockIndicatorF2D.h"

namespace olb {

template<typename T, typename DESCRIPTOR, class BoundaryManager>
class AdvectionDiffusionBoundaryConditionInstantiator2D : public OnLatticeAdvectionDiffusionBoundaryCondition2D<T,DESCRIPTOR> {
public:
  AdvectionDiffusionBoundaryConditionInstantiator2D( BlockLatticeStructure2D<T,DESCRIPTOR>& block_ );
  ~AdvectionDiffusionBoundaryConditionInstantiator2D() override;

  void addTemperatureBoundary0N(int x0, int x1, int y0, int y1, T omega) override;
  void addTemperatureBoundary0P(int x0, int x1, int y0, int y1, T omega) override;
  void addTemperatureBoundary1N(int x0, int x1, int y0, int y1, T omega) override;
  void addTemperatureBoundary1P(int x0, int x1, int y0, int y1, T omega) override;

  void addTemperatureCornerNN(int x, int y, T omega) override;
  void addTemperatureCornerNP(int x, int y, T omega) override;
  void addTemperatureCornerPN(int x, int y, T omega) override;
  void addTemperatureCornerPP(int x, int y, T omega) override;

  virtual BlockLatticeStructure2D<T,DESCRIPTOR>& getBlock();
  virtual BlockLatticeStructure2D<T,DESCRIPTOR> const& getBlock() const;

  void addTemperatureBoundary(BlockIndicatorF2D<T>& indicator, int x0, int x1, int y0, int y1, T omega) override;

  void addRegularizedTemperatureBoundary(BlockIndicatorF2D<T>& indicator, int x0, int x1, int y0, int y1, T omega) override;
  void addRegularizedTemperatureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, T omega) override;
  void addRegularizedTemperatureBoundary(BlockIndicatorF2D<T>& indicator, T omega) override;
  void addRegularizedTemperatureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega) override;

  void addRegularizedHeatFluxBoundary(BlockIndicatorF2D<T>& indicator, int x0, int x1, int y0, int y1, T omega, T *heatFlux) override;
  void addRegularizedHeatFluxBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0, int x1, int y0, int y1, T omega, T *heatFlux) override;
  void addRegularizedHeatFluxBoundary(BlockIndicatorF2D<T>& indicator, T omega, T *heatFlux) override;
  void addRegularizedHeatFluxBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega, T *heatFlux) override;

private:
  template<int direction, int orientation>
  void addTemperatureBoundary(int x0, int x1, int y0, int y1, T omega);
  template<int normalX, int normalY>
  void addTemperatureCorner(int x, int y, T omega);

  template<int direction, int orientation>
  void addRegularizedTemperatureBoundary(int x0, int x1, int y0, int y1, T omega);
  template<int normalX, int normalY>
  void addRegularizedTemperatureCorner(int x, int y, T omega);

  template<int direction, int orientation>
  void addRegularizedHeatFluxBoundary(int x0, int x1, int y0, int y1, T omega, T *heatFlux);
  template<int normalX, int normalY>
  void addRegularizedHeatFluxCorner(int x, int y, T omega);

  BlockLatticeStructure2D<T,DESCRIPTOR>& block;
  std::vector<Momenta<T,DESCRIPTOR>*>  momentaVector;
  std::vector<Dynamics<T,DESCRIPTOR>*> dynamicsVector;
};

///////// class AdvectionDiffusionBoundaryConditionInstantiator2D ////////////////////////

template<typename T, typename DESCRIPTOR, class BoundaryManager>
AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::AdvectionDiffusionBoundaryConditionInstantiator2D (
  BlockLatticeStructure2D<T,DESCRIPTOR>& block_)
  : block(block_)
{ }

template<typename T, typename DESCRIPTOR, class BoundaryManager>
AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::~AdvectionDiffusionBoundaryConditionInstantiator2D()
{
  for (unsigned iDynamics=0; iDynamics<dynamicsVector.size(); ++iDynamics) {
    delete dynamicsVector[iDynamics];
  }
  for (unsigned iMomenta=0; iMomenta<dynamicsVector.size(); ++iMomenta) {
    delete momentaVector[iMomenta];
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int direction, int orientation>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addTemperatureBoundary(int x0, int x1, int y0, int y1, T omega)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      Momenta<T,DESCRIPTOR>* momenta
        = BoundaryManager::template getTemperatureBoundaryMomenta<direction,orientation>();
      Dynamics<T,DESCRIPTOR>* dynamics
        = BoundaryManager::template getTemperatureBoundaryDynamics<direction,orientation>(omega, *momenta);
      this->getBlock().defineDynamics(iX,iX,iY,iY, dynamics);
      momentaVector.push_back(momenta);
      dynamicsVector.push_back(dynamics);
    }
  }

  PostProcessorGenerator2D<T,DESCRIPTOR>* postProcessor
    = BoundaryManager::template getTemperatureBoundaryProcessor<direction,orientation>(x0,x1, y0,y1);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int xNormal, int yNormal>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addTemperatureCorner(int x, int y, T omega)
{
  Momenta<T,DESCRIPTOR>* momenta
    = BoundaryManager::template getTemperatureCornerMomenta<xNormal,yNormal>();
  Dynamics<T,DESCRIPTOR>* dynamics
    = BoundaryManager::template getTemperatureCornerDynamics<xNormal,yNormal>(omega, *momenta);

  this->getBlock().defineDynamics(x,x,y,y, dynamics);

  PostProcessorGenerator2D<T,DESCRIPTOR>* postProcessor
    = BoundaryManager::template getTemperatureCornerProcessor<xNormal,yNormal>(x, y);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}




template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int direction, int orientation>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addRegularizedTemperatureBoundary(int x0, int x1, int y0, int y1, T omega)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      Momenta<T,DESCRIPTOR>* momenta
        = BoundaryManager::template getRegularizedTemperatureBoundaryMomenta<direction,orientation>();
      Dynamics<T,DESCRIPTOR>* dynamics
        = BoundaryManager::template getRegularizedTemperatureBoundaryDynamics<direction,orientation>(omega, *momenta);
      this->getBlock().defineDynamics(iX,iX,iY,iY, dynamics);
      momentaVector.push_back(momenta);
      dynamicsVector.push_back(dynamics);
    }
  }

  PostProcessorGenerator2D<T,DESCRIPTOR>* postProcessor
    = BoundaryManager::template getRegularizedTemperatureBoundaryProcessor<direction,orientation>(x0,x1, y0,y1);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int xNormal, int yNormal>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addRegularizedTemperatureCorner(int x, int y, T omega)
{
  Momenta<T,DESCRIPTOR>* momenta
    = BoundaryManager::template getRegularizedTemperatureCornerMomenta<xNormal,yNormal>();
  Dynamics<T,DESCRIPTOR>* dynamics
    = BoundaryManager::template getRegularizedTemperatureCornerDynamics<xNormal,yNormal>(omega, *momenta);

  this->getBlock().defineDynamics(x,x,y,y, dynamics);

  PostProcessorGenerator2D<T,DESCRIPTOR>* postProcessor
    = BoundaryManager::template getRegularizedTemperatureCornerProcessor<xNormal,yNormal>(x, y);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}



template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int direction, int orientation>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addRegularizedHeatFluxBoundary(int x0, int x1, int y0, int y1, T omega, T *heatFlux)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      Momenta<T,DESCRIPTOR>* momenta
        = BoundaryManager::template getRegularizedHeatFluxBoundaryMomenta<direction,orientation>(heatFlux);
      Dynamics<T,DESCRIPTOR>* dynamics
        = BoundaryManager::template getRegularizedHeatFluxBoundaryDynamics<direction,orientation>(omega, *momenta);
      this->getBlock().defineDynamics(iX,iX,iY,iY, dynamics);
      momentaVector.push_back(momenta);
      dynamicsVector.push_back(dynamics);
    }
  }

  PostProcessorGenerator2D<T,DESCRIPTOR>* postProcessor
    = BoundaryManager::template getRegularizedHeatFluxBoundaryProcessor<direction,orientation>(x0,x1, y0,y1);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int xNormal, int yNormal>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addRegularizedHeatFluxCorner(int x, int y, T omega)
{
  Momenta<T,DESCRIPTOR>* momenta
    = BoundaryManager::template getRegularizedHeatFluxCornerMomenta<xNormal,yNormal>();
  Dynamics<T,DESCRIPTOR>* dynamics
    = BoundaryManager::template getRegularizedHeatFluxCornerDynamics<xNormal,yNormal>(omega, *momenta);

  this->getBlock().defineDynamics(x,x,y,y, dynamics);

  PostProcessorGenerator2D<T,DESCRIPTOR>* postProcessor
    = BoundaryManager::template getRegularizedHeatFluxCornerProcessor<xNormal,yNormal>(x, y);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}




template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addTemperatureBoundary(BlockIndicatorF2D<T>& indicator,
                       int x0, int x1, int y0, int y1, T omega)
{
  std::vector<int> discreteNormal(3,0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      if (indicator(iX, iY)) {
        discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX,iY);
        if (discreteNormal[0] == 0) {
          if (discreteNormal[1] == 1) {
            addTemperatureBoundary<0,1>(iX,iX,iY,iY, omega);
          }
          else if (discreteNormal[1] == -1) {
            addTemperatureBoundary<0,-1>(iX,iX,iY,iY, omega);
          }
          else if (discreteNormal[2] == 1) {
            addTemperatureBoundary<1,1>(iX,iX,iY,iY, omega);
          }
          else if (discreteNormal[2] == -1) {
            addTemperatureBoundary<1,-1>(iX,iX,iY,iY, omega);
          }
        }
        else if (discreteNormal[0] == 1) {
          if (discreteNormal[1] == 1 && discreteNormal[2] == 1) {
            addTemperatureCorner<1,1>(iX,iY, omega);
          }
          else if (discreteNormal[1] == 1 && discreteNormal[2] == -1) {
            addTemperatureCorner<1,-1>(iX,iY, omega);
          }
          else if (discreteNormal[1] == -1 && discreteNormal[2] == 1) {
            addTemperatureCorner<-1,1>(iX,iY, omega);
          }
          else if (discreteNormal[1] == -1 && discreteNormal[2] == -1) {
            addTemperatureCorner<-1,-1>(iX,iY, omega);
          }
        }
      }
    }
  }
}




template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addRegularizedTemperatureBoundary(BlockIndicatorF2D<T>& indicator,
                       int x0, int x1, int y0, int y1, T omega)
{
  std::vector<int> discreteNormal(3,0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      if (indicator(iX, iY)) {
        discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX,iY);
        // cout << discreteNormal[0] << " " << discreteNormal[1] << " " << discreteNormal[2] << endl;
        if (discreteNormal[0] == 0) {
          if (discreteNormal[1] == 1) {
            addRegularizedTemperatureBoundary<0,1>(iX,iX,iY,iY, omega);
          }
          else if (discreteNormal[1] == -1) {
            addRegularizedTemperatureBoundary<0,-1>(iX,iX,iY,iY, omega);
          }
          else if (discreteNormal[2] == 1) {
            addRegularizedTemperatureBoundary<1,1>(iX,iX,iY,iY, omega);
          }
          else if (discreteNormal[2] == -1) {
            addRegularizedTemperatureBoundary<1,-1>(iX,iX,iY,iY, omega);
          }
        }
        else if (discreteNormal[0] == 1) {
          if (discreteNormal[1] == 1 && discreteNormal[2] == 1) {
            addRegularizedTemperatureCorner<1,1>(iX,iY, omega);
          }
          else if (discreteNormal[1] == 1 && discreteNormal[2] == -1) {
            addRegularizedTemperatureCorner<1,-1>(iX,iY, omega);
          }
          else if (discreteNormal[1] == -1 && discreteNormal[2] == 1) {
            addRegularizedTemperatureCorner<-1,1>(iX,iY, omega);
          }
          else if (discreteNormal[1] == -1 && discreteNormal[2] == -1) {
            addRegularizedTemperatureCorner<-1,-1>(iX,iY, omega);
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addRegularizedTemperatureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
                       int x0, int x1, int y0, int y1, T omega)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addRegularizedTemperatureBoundary(indicator, x0, x1, y0, y1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addRegularizedTemperatureBoundary(BlockIndicatorF2D<T>& indicator, T omega)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  addRegularizedTemperatureBoundary(indicator,
                         0, blockGeometryStructure.getNx()-1,
                         0, blockGeometryStructure.getNy()-1,
                         omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addRegularizedTemperatureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addRegularizedTemperatureBoundary(indicator,
                         0, blockGeometryStructure.getNx()-1,
                         0, blockGeometryStructure.getNy()-1,
                         omega);
}



template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addRegularizedHeatFluxBoundary(BlockIndicatorF2D<T>& indicator,
                       int x0, int x1, int y0, int y1, T omega, T *heatFlux)
{
  std::vector<int> discreteNormal(3,0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      if (indicator(iX, iY)) {
        discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX,iY);
        // cout << discreteNormal[0] << " " << discreteNormal[1] << " " << discreteNormal[2] << endl;
        if (discreteNormal[0] == 0) {
          if (discreteNormal[1] == 1) {
            addRegularizedHeatFluxBoundary<0,1>(iX,iX,iY,iY, omega, heatFlux);
          }
          else if (discreteNormal[1] == -1) {
            addRegularizedHeatFluxBoundary<0,-1>(iX,iX,iY,iY, omega, heatFlux);
          }
          else if (discreteNormal[2] == 1) {
            addRegularizedHeatFluxBoundary<1,1>(iX,iX,iY,iY, omega, heatFlux);
          }
          else if (discreteNormal[2] == -1) {
            addRegularizedHeatFluxBoundary<1,-1>(iX,iX,iY,iY, omega, heatFlux);
          }
        }
        else if (discreteNormal[0] == 1) {
          if (discreteNormal[1] == 1 && discreteNormal[2] == 1) {
            addRegularizedHeatFluxCorner<1,1>(iX,iY, omega);
          }
          else if (discreteNormal[1] == 1 && discreteNormal[2] == -1) {
            addRegularizedHeatFluxCorner<1,-1>(iX,iY, omega);
          }
          else if (discreteNormal[1] == -1 && discreteNormal[2] == 1) {
            addRegularizedHeatFluxCorner<-1,1>(iX,iY, omega);
          }
          else if (discreteNormal[1] == -1 && discreteNormal[2] == -1) {
            addRegularizedHeatFluxCorner<-1,-1>(iX,iY, omega);
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addRegularizedHeatFluxBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
                       int x0, int x1, int y0, int y1, T omega, T *heatFlux)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addRegularizedHeatFluxBoundary(indicator, x0, x1, y0, y1, omega, heatFlux);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addRegularizedHeatFluxBoundary(BlockIndicatorF2D<T>& indicator, T omega, T *heatFlux)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  addRegularizedHeatFluxBoundary(indicator,
                         0, blockGeometryStructure.getNx()-1,
                         0, blockGeometryStructure.getNy()-1,
                         omega, heatFlux);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addRegularizedHeatFluxBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega, T *heatFlux)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometryStructure, material);
  addRegularizedHeatFluxBoundary(indicator,
                         0, blockGeometryStructure.getNx()-1,
                         0, blockGeometryStructure.getNy()-1,
                         omega, heatFlux);
}


template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addTemperatureBoundary0N(int x0, int x1, int y0, int y1, T omega)
{
  addTemperatureBoundary<0,-1>(x0,x1,y0,y1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addTemperatureBoundary0P(int x0, int x1, int y0, int y1, T omega)
{
  addTemperatureBoundary<0,1>(x0,x1,y0,y1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addTemperatureBoundary1N(int x0, int x1, int y0, int y1, T omega)
{
  addTemperatureBoundary<1,-1>(x0,x1,y0,y1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addTemperatureBoundary1P(int x0, int x1, int y0, int y1, T omega)
{
  addTemperatureBoundary<1,1>(x0,x1,y0,y1, omega);
}


template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addTemperatureCornerNN(int x, int y, T omega)
{
  addTemperatureCorner<-1,-1>(x,y, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addTemperatureCornerNP(int x, int y, T omega)
{
  addTemperatureCorner<-1,1>(x,y, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addTemperatureCornerPN(int x, int y, T omega)
{
  addTemperatureCorner<1,-1>(x,y, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::
addTemperatureCornerPP(int x, int y, T omega)
{
  addTemperatureCorner<1,1>(x,y, omega);
}


template<typename T, typename DESCRIPTOR, class BoundaryManager>
BlockLatticeStructure2D<T,DESCRIPTOR>& AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::getBlock()
{
  return block;
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
BlockLatticeStructure2D<T,DESCRIPTOR> const& AdvectionDiffusionBoundaryConditionInstantiator2D<T,DESCRIPTOR,BoundaryManager>::getBlock() const
{
  return block;
}

}


#endif
