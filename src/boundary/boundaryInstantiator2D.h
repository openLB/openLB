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
 * A helper for initialising 2D boundaries -- header file.
 */

#ifndef BOUNDARY_INSTANTIATOR_2D_H
#define BOUNDARY_INSTANTIATOR_2D_H

#include "boundaryCondition2D.h"
#include "boundaryPostProcessors2D.h"
#include "geometry/blockGeometry2D.h"
#include "geometry/blockGeometryStatistics2D.h"
#include "io/ostreamManager.h"
#include "functors/lattice/indicator/blockIndicatorF2D.h"
#include "dynamics/freeEnergyDynamics.h"


namespace olb {

template<typename T, typename DESCRIPTOR, class BoundaryManager>
class BoundaryConditionInstantiator2D: public OnLatticeBoundaryCondition2D<T,
  DESCRIPTOR> {
public:
  BoundaryConditionInstantiator2D(BlockLatticeStructure2D<T, DESCRIPTOR>& block_);
  ~BoundaryConditionInstantiator2D() override;

  void addVelocityBoundary0N(int x0, int x1, int y0, int y1, T omega) override;
  void addVelocityBoundary0P(int x0, int x1, int y0, int y1, T omega) override;
  void addVelocityBoundary1N(int x0, int x1, int y0, int y1, T omega) override;
  void addVelocityBoundary1P(int x0, int x1, int y0, int y1, T omega) override;

  void addSlipBoundary(int x0, int x1, int y0, int y1, int discreteNormalX, int discreteNormalY);

  void addPartialSlipBoundary(T tuner, int x0, int x1, int y0, int y1, int discreteNormalX, int discreteNormalY);

  void addConvectionBoundary0N(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary0P(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary1N(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) override;
  void addConvectionBoundary1P(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) override;

  void addExternalVelocityCornerNN(int x, int y, T omega) override;
  void addExternalVelocityCornerNP(int x, int y, T omega) override;
  void addExternalVelocityCornerPN(int x, int y, T omega) override;
  void addExternalVelocityCornerPP(int x, int y, T omega) override;

  void addInternalVelocityCornerNN(int x, int y, T omega) override;
  void addInternalVelocityCornerNP(int x, int y, T omega) override;
  void addInternalVelocityCornerPN(int x, int y, T omega) override;
  void addInternalVelocityCornerPP(int x, int y, T omega) override;

  void addVelocityBoundary(BlockIndicatorF2D<T>& indicator,
                           int x0, int x1, int y0, int y1,
                           T omega) override;
  void addSlipBoundary(BlockIndicatorF2D<T>& indicator,
                       int x0, int x1, int y0, int y1) override;
  void addPartialSlipBoundary(T tuner, BlockIndicatorF2D<T>& indicator,
                       int x0, int x1, int y0, int y1) override;
  void addPressureBoundary(BlockIndicatorF2D<T>& indicator,
                           int x0, int x1, int y0, int y1,
                           T omega) override;
  void addConvectionBoundary(BlockIndicatorF2D<T>& indicator,
                             int x0, int x1, int y0, int y1,
                             T omega, T* uAv=NULL) override;
  void addFreeEnergyWallBoundary(BlockIndicatorF2D<T>& indicator,
                           int x0, int x1, int y0, int y1,
                           T addend, int latticeNumber) override;
  void addFreeEnergyInletBoundary(BlockIndicatorF2D<T>& indicator,
                           int x0, int x1, int y0, int y1, T omega,
                           std::string type, int latticeNumber) override;
  void addFreeEnergyOutletBoundary(BlockIndicatorF2D<T>& indicator,
                           int x0, int x1, int y0, int y1, T omega,
                           std::string type, int latticeNumber) override;

  void outputOn() override;
  void outputOff() override;

  BlockLatticeStructure2D<T, DESCRIPTOR>& getBlock() override;
  BlockLatticeStructure2D<T, DESCRIPTOR> const& getBlock() const override;
private:
  template<int direction, int orientation>
  void addVelocityBoundary(int x0, int x1, int y0, int y1, T omega);
  template<int direction, int orientation>
  void addPressureBoundary(int x0, int x1, int y0, int y1, T omega);
  template<int direction, int orientation>
  void addConvectionBoundary(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL);
  template<int normalX, int normalY>
  void addExternalVelocityCorner(int x, int y, T omega);
  template<int normalX, int normalY>
  void addInternalVelocityCorner(int x, int y, T omega);
private:
  BlockLatticeStructure2D<T, DESCRIPTOR>& block;
  std::vector<Momenta<T, DESCRIPTOR>*> momentaVector;
  std::vector<Dynamics<T, DESCRIPTOR>*> dynamicsVector;
  bool _output;
  mutable OstreamManager clout;
};

///////// class BoundaryConditionInstantiator2D ////////////////////////

template<typename T, typename DESCRIPTOR, class BoundaryManager>
BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::BoundaryConditionInstantiator2D(
  BlockLatticeStructure2D<T, DESCRIPTOR>& block_) :
  block(block_), _output(false), clout(std::cout,"BoundaryConditionInstantiator2D")
{
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::~BoundaryConditionInstantiator2D()
{
  for (auto &iDynamics : dynamicsVector) {
    delete iDynamics;
  }
  for (auto &iMomenta : momentaVector) {
    delete iMomenta;
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int direction, int orientation>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addVelocityBoundary(
  int x0, int x1, int y0, int y1, T omega)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);

  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      Momenta<T, DESCRIPTOR>* momenta =
        BoundaryManager::template getVelocityBoundaryMomenta<
          direction, orientation>();
      Dynamics<T, DESCRIPTOR>* dynamics =
        BoundaryManager::template getVelocityBoundaryDynamics<
          direction, orientation>(omega, *momenta);
      this->getBlock().defineDynamics(iX, iX, iY, iY, dynamics);
      momentaVector.push_back(momenta);
      dynamicsVector.push_back(dynamics);
      if (_output) {
        clout << "addVelocityBoundary<" << direction << ","<< orientation << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << omega << " )" << std::endl;
      }
    }
  }

  PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor =
    BoundaryManager::template getVelocityBoundaryProcessor<direction,
        orientation>(x0, x1, y0, y1);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addSlipBoundary(
  int x0, int x1, int y0, int y1, int discreteNormalX, int discreteNormalY)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);

  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      if (_output) {
        clout << "addSlipBoundary<" << discreteNormalX << ","<< discreteNormalY << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << " )" << std::endl;
      }
    }
  }

  PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor = new SlipBoundaryProcessorGenerator2D<T, DESCRIPTOR>(x0, x1, y0, y1, discreteNormalX, discreteNormalY);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addPartialSlipBoundary(
  T tuner, int x0, int x1, int y0, int y1, int discreteNormalX, int discreteNormalY)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);

  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      if (_output) {
        clout << "addPartialSlipBoundary<" << discreteNormalX << ","<< discreteNormalY << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << " )" << std::endl;
      }
    }
  }

  PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor = new PartialSlipBoundaryProcessorGenerator2D<T, DESCRIPTOR>(tuner, x0, x1, y0, y1, discreteNormalX, discreteNormalY);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int direction, int orientation>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addPressureBoundary(
  int x0, int x1, int y0, int y1, T omega)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);

  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      Momenta<T, DESCRIPTOR>* momenta =
        BoundaryManager::template getPressureBoundaryMomenta<
          direction, orientation>();
      Dynamics<T, DESCRIPTOR>* dynamics =
        BoundaryManager::template getPressureBoundaryDynamics<
          direction, orientation>(omega, *momenta);
      this->getBlock().defineDynamics(iX, iX, iY, iY, dynamics);
      momentaVector.push_back(momenta);
      dynamicsVector.push_back(dynamics);
      if (_output) {
        clout << "addPressureBoundary<" << direction << ","<< orientation << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << omega << " )" << std::endl;
      }
    }
  }

  PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor =
    BoundaryManager::template getPressureBoundaryProcessor<direction,
        orientation>(x0, x1, y0, y1);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int direction, int orientation>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addConvectionBoundary(
  int x0, int x1, int y0, int y1, T omega, T* uAv)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);

  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      if (_output) {
        clout << "addConvectionBoundary<" << direction << ","<< orientation << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << omega << " )" << std::endl;
      }
    }
  }

  PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor =
    BoundaryManager::template getConvectionBoundaryProcessor<direction,orientation>(x0, x1, y0, y1, uAv);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addFreeEnergyWallBoundary(
  BlockIndicatorF2D<T>& indicator, int x0, int x1, int y0, int y1, T addend, int latticeNumber)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  std::vector<int> discreteNormal(3, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      if(indicator(iX,iY)) {
        discreteNormal = blockGeometryStructure.getStatistics().getType(iX,iY);
        if (discreteNormal[1]!=0 || discreteNormal[2]!=0) {

          Dynamics<T, DESCRIPTOR>* dynamics = NULL;
          if (latticeNumber == 1) {
            dynamics = &instances::getBounceBack<T, DESCRIPTOR>();
          } else {
            dynamics = new FreeEnergyWallDynamics<T, DESCRIPTOR>;
            dynamicsVector.push_back(dynamics);
          }
          this->getBlock().get(iX,iY).defineDynamics(dynamics);
          dynamicsVector.push_back(dynamics);

          PostProcessorGenerator2D<T, DESCRIPTOR>* wettingPostProcessor =
            new FreeEnergyWallProcessorGenerator2D<T, DESCRIPTOR> (
                iX, iX, iY, iY, discreteNormal[1], discreteNormal[2], addend );
          PostProcessorGenerator2D<T, DESCRIPTOR>* chemPotPostProcessor =
            new FreeEnergyChemPotBoundaryProcessorGenerator2D<T, DESCRIPTOR> (
              iX, iX, iY, iY, discreteNormal[1], discreteNormal[2], latticeNumber );
          if (wettingPostProcessor) {
            this->getBlock().addPostProcessor(*wettingPostProcessor);
          }
          if (chemPotPostProcessor) {
            this->getBlock().addPostProcessor(*chemPotPostProcessor);
          }
        }
        if (_output) {
          clout << "addFreeEnergyWallBoundary<" << "," << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << " )" << std::endl;
        }
      }
    }
  }

}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addFreeEnergyInletBoundary(
  BlockIndicatorF2D<T>& indicator, int x0, int x1, int y0, int y1,
  T omega, std::string type, int latticeNumber)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  std::vector<int> discreteNormal(3, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      if(indicator(iX,iY)) {
        discreteNormal = blockGeometryStructure.getStatistics().getType(iX,iY);
        if (discreteNormal[0] == 0) {
          Momenta<T, DESCRIPTOR>* momenta = NULL;
          Dynamics<T, DESCRIPTOR>* dynamics = NULL;

          if (discreteNormal[1] == -1) {
            if (latticeNumber == 1) {
              if (type == "density") {
                addPressureBoundary<0,-1>(iX, iX, iY, iY, omega);
              } else {
                addVelocityBoundary<0,-1>(iX, iX, iY, iY, omega);
              }
            } else {
              momenta = BoundaryManager::template 
                getPressureBoundaryMomenta<0,-1>();
              dynamics = new FreeEnergyInletOutletDynamics<T,DESCRIPTOR,0,-1>(omega,*momenta);
            }
          }

          else if (discreteNormal[1] == 1) {
            if (latticeNumber == 1) {
              if (type == "density") {
                addPressureBoundary<0,1>(iX, iX, iY, iY, omega);
              } else {
                addVelocityBoundary<0,1>(iX, iX, iY, iY, omega);
              }
            } else {
              momenta = BoundaryManager::template 
                getPressureBoundaryMomenta<0,1>();
              dynamics = new FreeEnergyInletOutletDynamics<T,DESCRIPTOR,0,1>(omega,*momenta);
            }
          }

          else if (discreteNormal[2] == -1) {
            if (latticeNumber == 1) {
              if (type == "density") {
                addPressureBoundary<1,-1>(iX, iX, iY, iY, omega);
              } else {
                addVelocityBoundary<1,-1>(iX, iX, iY, iY, omega);
              }
            } else {
              momenta = BoundaryManager::template 
                getPressureBoundaryMomenta<1,-1>();
              dynamics = new FreeEnergyInletOutletDynamics<T,DESCRIPTOR,1,-1>(omega,*momenta);
            }
          }

          else if (discreteNormal[2] == 1) {
            if (latticeNumber == 1) {
              if (type == "density") {
                addPressureBoundary<1,1>(iX, iX, iY, iY, omega);
              } else {
                addVelocityBoundary<1,1>(iX, iX, iY, iY, omega);
              }
            } else {
              momenta = BoundaryManager::template 
                getPressureBoundaryMomenta<1,1>();
              dynamics = new FreeEnergyInletOutletDynamics<T,DESCRIPTOR,1,1>(omega,*momenta);
            }
          }

          if (latticeNumber != 1) {
            this->getBlock().get(iX,iY).defineDynamics(dynamics);
            momentaVector.push_back(momenta);
            dynamicsVector.push_back(dynamics);
          }
          PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor =
            new FreeEnergyChemPotBoundaryProcessorGenerator2D<T, DESCRIPTOR> (
              iX, iX, iY, iY, discreteNormal[1], discreteNormal[2], latticeNumber );
          if (postProcessor) {
            this->getBlock().addPostProcessor(*postProcessor);
          }

          if (_output) {
            clout << "addFreeEnergyInletBoundary<" << "," << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << " )" << std::endl;
          }

        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addFreeEnergyOutletBoundary(
  BlockIndicatorF2D<T>& indicator, int x0, int x1, int y0, int y1,
  T omega, std::string type, int latticeNumber)
{
  addFreeEnergyInletBoundary(indicator, x0, x1, y0, y1, omega, type, latticeNumber);

  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  std::vector<int> discreteNormal(3, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      if(indicator(iX,iY)) {
        discreteNormal = blockGeometryStructure.getStatistics().getType(iX,iY);
        if (discreteNormal[0] == 0) {

          PostProcessorGenerator2D<T, DESCRIPTOR>* convectivePostProcessor =
            new FreeEnergyConvectiveProcessorGenerator2D<T, DESCRIPTOR> (
              iX, iX, iY, iY, discreteNormal[1], discreteNormal[2] );
          if (convectivePostProcessor) {
            this->getBlock().addPostProcessor(*convectivePostProcessor);
          }

          if (_output) {
            clout << "addFreeEnergyOutletBoundary<" << "," << ">("  << x0 << ", "<< x1 << ", " << y0 << ", " << y1 << ", " << " )" << std::endl;
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int xNormal, int yNormal>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addExternalVelocityCorner(
  int x, int y, T omega)
{
  Momenta<T, DESCRIPTOR>* momenta =
    BoundaryManager::template getExternalVelocityCornerMomenta<xNormal,
        yNormal>();
  Dynamics<T, DESCRIPTOR>* dynamics =
    BoundaryManager::template getExternalVelocityCornerDynamics<
      xNormal, yNormal>(omega, *momenta);

  this->getBlock().defineDynamics(x, x, y, y, dynamics);
  momentaVector.push_back(momenta);
  dynamicsVector.push_back(dynamics);

  PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor =
    BoundaryManager::template getExternalVelocityCornerProcessor<
      xNormal, yNormal>(x, y);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
  if (_output) {
    clout << "addExternalVelocityCorner<" << xNormal << ","<< yNormal << ">("  << x << ", "<< y << omega << " )" << std::endl;
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
template<int xNormal, int yNormal>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addInternalVelocityCorner(
  int x, int y, T omega)
{
  Momenta<T, DESCRIPTOR>* momenta =
    BoundaryManager::template getInternalVelocityCornerMomenta<xNormal,
        yNormal>();
  Dynamics<T, DESCRIPTOR>* dynamics =
    BoundaryManager::template getInternalVelocityCornerDynamics<
      xNormal, yNormal>(omega, *momenta);

  this->getBlock().defineDynamics(x, x, y, y, dynamics);
  momentaVector.push_back(momenta);
  dynamicsVector.push_back(dynamics);

  PostProcessorGenerator2D<T, DESCRIPTOR>* postProcessor =
    BoundaryManager::template getInternalVelocityCornerProcessor<
      xNormal, yNormal>(x, y);
  if (postProcessor) {
    this->getBlock().addPostProcessor(*postProcessor);
  }
  if (_output) {
    clout << "addInternalVelocityCorner<" << xNormal << ","<< yNormal << ">("  << x << ", "<< y << omega << " )" << std::endl;
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addVelocityBoundary(
  BlockIndicatorF2D<T>& indicator, int x0, int x1, int y0, int y1, T omega)
{
  std::vector<int> discreteNormal(3, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      if (indicator(iX, iY)) {
        discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX, iY);
        if (discreteNormal[0] == 0) {
          if (discreteNormal[1] == 1) {
            addVelocityBoundary<0, 1> (iX, iX, iY, iY, omega);
          }
          else if (discreteNormal[1] == -1) {
            addVelocityBoundary<0, -1> (iX, iX, iY, iY, omega);
          }
          else if (discreteNormal[2] == 1) {
            addVelocityBoundary<1, 1> (iX, iX, iY, iY, omega);
          }
          else if (discreteNormal[2] == -1) {
            addVelocityBoundary<1, -1> (iX, iX, iY, iY, omega);
          }
          else {
            clout << "Could not addVelocityBoundary (" << iX
                  << ", " << iY << ")" << std::endl;
          }
        }
        else if (discreteNormal[0] == 1) {
          if (discreteNormal[1] == 1) {
            if (discreteNormal[2] == 1) {
              addExternalVelocityCorner<1, 1> (iX, iY, omega);
            }
            else if (discreteNormal[2] == -1) {
              addExternalVelocityCorner<1, -1> (iX, iY, omega);
            }
            else {
              clout << "Could not addVelocityBoundary ("
                    << iX << ", " << iY << ")" << std::endl;
            }
          }
          else if (discreteNormal[1] == -1) {
            if (discreteNormal[2] == 1) {
              addExternalVelocityCorner<-1, 1> (iX, iY, omega);
            }
            else if (discreteNormal[2] == -1) {
              addExternalVelocityCorner<-1, -1> (iX, iY, omega);
            }
            else {
              clout << "Could not addVelocityBoundary ("
                    << iX << ", " << iY << ")" << std::endl;
            }
          }
        }
        else if (discreteNormal[0] == 2) {
          if (discreteNormal[1] == 1) {
            if (discreteNormal[2] == 1) {
              addInternalVelocityCorner<1, 1> (iX, iY, omega);
            }
            else if (discreteNormal[2] == -1) {
              addInternalVelocityCorner<1, -1> (iX, iY, omega);
            }
            else {
              clout << "Could not addVelocityBoundary ("
                    << iX << ", " << iY << ")" << std::endl;
            }
          }
          else if (discreteNormal[1] == -1) {
            if (discreteNormal[2] == 1) {
              addInternalVelocityCorner<-1, 1> (iX, iY, omega);
            }
            else if (discreteNormal[2] == -1) {
              addInternalVelocityCorner<-1, -1> (iX, iY, omega);
            }
            else {
              clout << "Could not addVelocityBoundary ("
                    << iX << ", " << iY << ")" << std::endl;
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addSlipBoundary(
  BlockIndicatorF2D<T>& indicator, int x0, int x1, int y0, int y1)
{
  std::vector<int> discreteNormal(3, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      if (indicator(iX, iY)) {
        discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX, iY);
        if (discreteNormal[1]!=0 || discreteNormal[2]!=0) {
          addSlipBoundary(iX, iX, iY, iY, discreteNormal[1], discreteNormal[2]);
        }
        else {
          clout << "Warning: Could not addSlipBoundary (" << iX << ", " << iY << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] <<"), set to bounceBack" << std::endl;
          this->getBlock().defineDynamics(iX, iY, &instances::getBounceBack<T, DESCRIPTOR>() );
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addPartialSlipBoundary(
  T tuner, BlockIndicatorF2D<T>& indicator, int x0, int x1, int y0, int y1)
{
  std::vector<int> discreteNormal(3, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      if (indicator(iX, iY)) {
        if (tuner < 0. || tuner > 1.) {
          clout << "Warning: Could not addPartialSlipBoundary (" << iX << ", " << iY << "), tuner must be between 0.1 and instead is=" << tuner <<", set to bounceBack" << std::endl;
         this->getBlock().defineDynamics(iX, iY, &instances::getBounceBack<T, DESCRIPTOR>() );
        } else {
          discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX, iY);
          if (discreteNormal[1]!=0 || discreteNormal[2]!=0) {
            addPartialSlipBoundary(tuner, iX, iX, iY, iY, discreteNormal[1], discreteNormal[2]);
          }
          else {
            clout << "Warning: Could not addPartialSlipBoundary (" << iX << ", " << iY << "), discreteNormal=(" << discreteNormal[0] <<","<< discreteNormal[1] <<","<< discreteNormal[2] <<"), set to bounceBack" << std::endl;
            this->getBlock().defineDynamics(iX, iY, &instances::getBounceBack<T, DESCRIPTOR>() );
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addPressureBoundary(
  BlockIndicatorF2D<T>& indicator, int x0, int x1, int y0, int y1, T omega)
{
  std::vector<int> discreteNormal(3, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      if (indicator(iX, iY)) {
        discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX, iY);
        if (discreteNormal[0] == 0) {
          if (discreteNormal[1] == -1) {
            addPressureBoundary<0, -1> (iX, iX, iY, iY, omega);
          }
          else if (discreteNormal[1] == 1) {
            addPressureBoundary<0, 1> (iX, iX, iY, iY, omega);
          }
          else if (discreteNormal[2] == -1) {
            addPressureBoundary<1, -1> (iX, iX, iY, iY, omega);
          }
          else if (discreteNormal[2] == 1) {
            addPressureBoundary<1, 1> (iX, iX, iY, iY, omega);
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addConvectionBoundary(
  BlockIndicatorF2D<T>& indicator, int x0, int x1, int y0, int y1, T omega, T* uAv)
{
  std::vector<int> discreteNormal(3, 0);
  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      if (indicator(iX, iY)) {
        discreteNormal = indicator.getBlockGeometryStructure().getStatistics().getType(iX, iY);
        if (discreteNormal[0] == 0) {
          if (discreteNormal[1] == -1) {
            addConvectionBoundary0N(iX, iX, iY, iY, omega, uAv);
          }
          else if (discreteNormal[1] == 1) {
            addConvectionBoundary0P(iX, iX, iY, iY, omega, uAv);
          }
          else if (discreteNormal[2] == -1) {
            addConvectionBoundary1N(iX, iX, iY, iY, omega, uAv);
          }
          else if (discreteNormal[2] == 1) {
            addConvectionBoundary1P(iX, iX, iY, iY, omega, uAv);
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addVelocityBoundary0N(
  int x0, int x1, int y0, int y1, T omega)
{
  addVelocityBoundary<0, -1> (x0, x1, y0, y1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addVelocityBoundary0P(
  int x0, int x1, int y0, int y1, T omega)
{
  addVelocityBoundary<0, 1> (x0, x1, y0, y1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addVelocityBoundary1N(
  int x0, int x1, int y0, int y1, T omega)
{
  addVelocityBoundary<1, -1> (x0, x1, y0, y1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addVelocityBoundary1P(
  int x0, int x1, int y0, int y1, T omega)
{
  addVelocityBoundary<1, 1> (x0, x1, y0, y1, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addConvectionBoundary0N(
  int x0, int x1, int y0, int y1, T omega, T* uAv)
{
  addConvectionBoundary<0, -1> (x0, x1, y0, y1, omega, uAv);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addConvectionBoundary0P(
  int x0, int x1, int y0, int y1, T omega, T* uAv)
{
  addConvectionBoundary<0, 1> (x0, x1, y0, y1, omega, uAv);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addConvectionBoundary1N(
  int x0, int x1, int y0, int y1, T omega, T* uAv)
{
  addConvectionBoundary<1, -1> (x0, x1, y0, y1, omega, uAv);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addConvectionBoundary1P(
  int x0, int x1, int y0, int y1, T omega, T* uAv)
{
  addConvectionBoundary<1, 1> (x0, x1, y0, y1, omega, uAv);
}


template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addExternalVelocityCornerNN(
  int x, int y, T omega)
{
  addExternalVelocityCorner<-1, -1> (x, y, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addExternalVelocityCornerNP(
  int x, int y, T omega)
{
  addExternalVelocityCorner<-1, 1> (x, y, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addExternalVelocityCornerPN(
  int x, int y, T omega)
{
  addExternalVelocityCorner<1, -1> (x, y, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addExternalVelocityCornerPP(
  int x, int y, T omega)
{
  addExternalVelocityCorner<1, 1> (x, y, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addInternalVelocityCornerNN(
  int x, int y, T omega)
{
  addInternalVelocityCorner<-1, -1> (x, y, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addInternalVelocityCornerNP(
  int x, int y, T omega)
{
  addInternalVelocityCorner<-1, 1> (x, y, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addInternalVelocityCornerPN(
  int x, int y, T omega)
{
  addInternalVelocityCorner<1, -1> (x, y, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::addInternalVelocityCornerPP(
  int x, int y, T omega)
{
  addInternalVelocityCorner<1, 1> (x, y, omega);
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
BlockLatticeStructure2D<T, DESCRIPTOR>& BoundaryConditionInstantiator2D<T, DESCRIPTOR,
                        BoundaryManager>::getBlock()
{
  return block;
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
BlockLatticeStructure2D<T, DESCRIPTOR> const& BoundaryConditionInstantiator2D<T, DESCRIPTOR,
                        BoundaryManager>::getBlock() const
{
  return block;
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::outputOn()
{
  _output = true;
}

template<typename T, typename DESCRIPTOR, class BoundaryManager>
void BoundaryConditionInstantiator2D<T, DESCRIPTOR, BoundaryManager>::outputOff()
{
  _output = false;
}

}

#endif
