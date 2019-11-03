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


#ifndef ADVECTION_DIFFUSION_BOUNDARY_INSTANTIATOR_3D_HH
#define ADVECTION_DIFFUSION_BOUNDARY_INSTANTIATOR_3D_HH

#include "advectionDiffusionBoundaryInstantiator3D.h"
#include "advectionDiffusionBoundaryCondition3D.h"
#include "advectionDiffusionBoundaryCondition3D.hh"
#include "advectionDiffusionBoundaryPostProcessor3D.hh"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "dynamics/rtlbmDescriptors.h"

namespace olb {

///////// class AdvectionDiffusionBoundaryConditionInstantiator3D ////////////////////////

template<typename T, typename DESCRIPTOR, class BoundaryManager>
AdvectionDiffusionBoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
AdvectionDiffusionBoundaryConditionInstantiator3D(BlockLatticeStructure3D<T,DESCRIPTOR>& block)
  : _block(block)
{
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
AdvectionDiffusionBoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
~AdvectionDiffusionBoundaryConditionInstantiator3D()
{
  for (auto &iDynamics : dynamicsVector) {
    delete iDynamics;
  }
  for (auto &iMomenta : momentaVector) {
    delete iMomenta;
  }
}

// ---- flat -------------
template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int direction, int orientation>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addTemperatureBoundary(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        Momenta<T,DESCRIPTOR>* momenta
          = BoundaryManager::template getTemperatureBoundaryMomenta<direction,orientation>();
        Dynamics<T,DESCRIPTOR>* dynamics
          = BoundaryManager::template getTemperatureBoundaryDynamics<direction,orientation>(omega, *momenta);
        _block.defineDynamics(iX,iX,iY,iY,iZ,iZ, dynamics);
        momentaVector.push_back(momenta);
        dynamicsVector.push_back(dynamics);
      }
    }
  }
  PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor
    = BoundaryManager::template getTemperatureBoundaryProcessor<direction,orientation>(x0,x1, y0,y1, z0,z1);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}

// ---- edges -------------
template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int plane, int normal1, int normal2>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addTemperatureBoundaryEdge(int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  OLB_PRECONDITION(
    ( x0==x1 && y0==y1 ) ||
    ( x0==x1 && z0==z1 ) ||
    ( y0==y1 && z0==z1 ) );

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        Momenta<T,DESCRIPTOR>* momenta
          = BoundaryManager::template getTemperatureBoundaryEdgeMomenta<plane,normal1,normal2>();
        Dynamics<T,DESCRIPTOR>* dynamics
          = BoundaryManager::template getTemperatureBoundaryEdgeDynamics<plane,normal1,normal2>(omega, *momenta);
        _block.defineDynamics(iX,iX,iY,iY,iZ,iZ, dynamics);
        momentaVector.push_back(momenta);
        dynamicsVector.push_back(dynamics);
      }
    }
  }

  PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor
    = BoundaryManager::template getTemperatureBoundaryEdgeProcessor<plane,normal1,normal2>(x0,x1, y0,y1, z0,z1);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}


// ---- corner -------------
template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int xNormal, int yNormal, int zNormal>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addTemperatureBoundaryCorner(int x, int y, int z, T omega)
{
  Momenta<T,DESCRIPTOR>* momenta
    = BoundaryManager::template getTemperatureBoundaryCornerMomenta<xNormal,yNormal,zNormal>();
  Dynamics<T,DESCRIPTOR>* dynamics
    = BoundaryManager::template getTemperatureBoundaryCornerDynamics<xNormal,yNormal,zNormal>(omega, *momenta);
  _block.defineDynamics(x,x,y,y,z,z, dynamics);
  PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor
    = BoundaryManager::template getTemperatureBoundaryCornerProcessor<xNormal,yNormal,zNormal>(x, y, z);
  if (postProcessor) {
    _block.addPostProcessor(*postProcessor);
  }
}


template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addTemperatureBoundary(BlockIndicatorF3D<T>& indicator,
                       int x0, int x1, int y0, int y1, int z0, int z1, T omega)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  std::vector<int> discreteNormal(4,0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if (indicator(iX, iY, iZ)) {
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);
          if (discreteNormal[0] == 0) { // flat
            if (discreteNormal[1] != 0 && discreteNormal[1] == -1) {
              addTemperatureBoundary<0,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] != 0 && discreteNormal[1] == 1) {
              addTemperatureBoundary<0,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == -1) {
              addTemperatureBoundary<1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[2] != 0 && discreteNormal[2] == 1) {
              addTemperatureBoundary<1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == -1) {
              addTemperatureBoundary<2,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[3] != 0 && discreteNormal[3] == 1) {
              addTemperatureBoundary<2,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
          }
          else if (discreteNormal[0] == 1) {    // corner
            if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              addTemperatureBoundaryCorner<1,1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              addTemperatureBoundaryCorner<1,-1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              addTemperatureBoundaryCorner<1,1,-1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              addTemperatureBoundaryCorner<1,-1,-1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              addTemperatureBoundaryCorner<-1,1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              addTemperatureBoundaryCorner<-1,-1,1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              addTemperatureBoundaryCorner<-1,1,-1>(iX,iY,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              addTemperatureBoundaryCorner<-1,-1,-1>(iX,iY,iZ, omega);
            }
          }
          else if (discreteNormal[0] == 3) {    // edge
            if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == 1) {
              addTemperatureBoundaryEdge<0,1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == 1) {
              addTemperatureBoundaryEdge<0,-1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == 1 && discreteNormal[3] == -1) {
              addTemperatureBoundaryEdge<0,1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 0 && discreteNormal[2] == -1 && discreteNormal[3] == -1) {
              addTemperatureBoundaryEdge<0,-1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              addTemperatureBoundaryEdge<1,1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == 1) {
              addTemperatureBoundaryEdge<1,1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              addTemperatureBoundaryEdge<1,-1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 0 && discreteNormal[3] == -1) {
              addTemperatureBoundaryEdge<1,-1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              addTemperatureBoundaryEdge<2,1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == 1 && discreteNormal[3] == 0) {
              addTemperatureBoundaryEdge<2,-1,1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == 1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              addTemperatureBoundaryEdge<2,1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
            else if (discreteNormal[1] == -1 && discreteNormal[2] == -1 && discreteNormal[3] == 0) {
              addTemperatureBoundaryEdge<2,-1,-1>(iX,iX,iY,iY,iZ,iZ, omega);
            }
          }

        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addConvectionBoundary(BlockIndicatorF3D<T>& indicator,
                      int x0, int x1, int y0, int y1, int z0, int z1)
{
  std::vector<int> discreteNormal(4, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if (indicator(iX, iY, iZ)) {
          discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX, iY, iZ);
          if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {
            PostProcessorGenerator3D<T, DESCRIPTOR>* postProcessor = new ConvectionBoundaryProcessorGenerator3D<T, DESCRIPTOR>(iX, iX, iY, iY, iZ, iZ, -discreteNormal[1], -discreteNormal[2], -discreteNormal[3]);
            if (postProcessor) {
              _block.addPostProcessor(*postProcessor);
            }
          }
          else {
//            cout << "Warning: Could not addConvectionBoundary (" << iX << ", " << iY << ", " << iZ << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] << ", " << discreteNormal[3] << "), set to bounceBack" << std::endl;
          }
        }

      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addZeroDistributionBoundary(BlockIndicatorF3D<T>& indicator,
                            int x0, int x1, int y0, int y1, int z0, int z1)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  std::vector<int> discreteNormal(4, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if (indicator(iX, iY, iZ)) {
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);
          if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {
            PostProcessorGenerator3D<T, DESCRIPTOR>* postProcessor = new ZeroDistributionBoundaryProcessorGenerator3D<T, DESCRIPTOR>(iX, iX, iY, iY, iZ, iZ, -discreteNormal[1], -discreteNormal[2], -discreteNormal[3]);
            if (postProcessor) {
              _block.addPostProcessor(*postProcessor);
            }
          }
          else {
//            cout << "Warning: Could not addZeroDistributionBoundary (" << iX << ", " << iY << ", " << iZ << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] << "," << discreteNormal[3] << ")" << std::endl;
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,DESCRIPTOR,BoundaryManager>::
addExtFieldBoundary(BlockIndicatorF3D<T>& indicator,
                    int offset, int x0, int x1, int y0, int y1, int z0, int z1)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  std::vector<int> discreteNormal(4, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        if (indicator(iX, iY, iZ)) {
          discreteNormal = blockGeometryStructure.getStatistics().getType(iX, iY, iZ);
          if (discreteNormal[1]!=0 || discreteNormal[2]!=0 || discreteNormal[3]!=0) {
            PostProcessorGenerator3D<T, DESCRIPTOR>* postProcessor = new ExtFieldBoundaryProcessorGenerator3D<T, DESCRIPTOR>(iX, iX, iY, iY, iZ, iZ, -discreteNormal[1], -discreteNormal[2], -discreteNormal[3], offset);
            if (postProcessor) {
              _block.addPostProcessor(*postProcessor);
            }
          }
          else {
//            cout << "Warning: Could not addZeroDistributionBoundary (" << iX << ", " << iY << ", " << iZ << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] << "," << discreteNormal[3] << ")" << std::endl;
          }
        }
      }
    }
  }
}


}

#endif
