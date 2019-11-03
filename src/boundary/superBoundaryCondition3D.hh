/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Mathias J. Krause
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
 * A helper for initialising 3D boundaries -- generic implementation.
 */

#ifndef SUPER_BOUNDARY_CONDITION_3D_HH
#define SUPER_BOUNDARY_CONDITION_3D_HH

#include <vector>
#include "boundaryCondition3D.h"
#include "geometry/superGeometry3D.h"
#include "extendedFiniteDifferenceBoundary3D.h"
#include "superBoundaryCondition3D.h"
#include "core/superLattice3D.h"
#include "functors/lattice/indicator/superIndicatorBaseF3D.h"
#include "advectionDiffusionBoundaryCondition3D.h"

namespace olb {

///////// class superBoundaryCondition3D ///////////////////////////////

template<typename T, typename DESCRIPTOR>
sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::sOnLatticeBoundaryCondition3D(
  SuperLattice3D<T, DESCRIPTOR>& sLattice) :
  clout(std::cout,"sOnLatticeBoundaryCondition3D"),
  _sLattice(sLattice),
  _output(false)
{
}

template<typename T, typename DESCRIPTOR>
sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::sOnLatticeBoundaryCondition3D(
  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> const& rhs) :
  clout(std::cout,"sOnLatticeBoundaryCondition3D"),
  _sLattice(rhs._sLattice),
  _output(false)
{
  _blockBCs = rhs._blockBCs;
  _ADblockBCs = rhs._ADblockBCs;
  _overlap = rhs._overlap;
}

template<typename T, typename DESCRIPTOR>
sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sOnLatticeBoundaryCondition3D<T,
                              DESCRIPTOR>::operator=(sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> rhs)
{
  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> tmp(rhs);
  return tmp;
}

template<typename T, typename DESCRIPTOR>
sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::~sOnLatticeBoundaryCondition3D()
{
  for (auto &iC : _blockBCs) {
    delete iC;
  }
  for (auto &iC : _ADblockBCs) {
    delete iC;
  }
}


template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addVelocityBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T omega)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iC = 0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
    _blockBCs[iC]->addVelocityBoundary(indicator->getExtendedBlockIndicatorF(iC),
                                       omega, includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addVelocityBoundary(
  SuperGeometry3D<T>& superGeometry, int material, T omega)
{
  addVelocityBoundary(superGeometry.getMaterialIndicator(material),
                      omega);
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addSlipBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->addSlipBoundary(
      indicator->getExtendedBlockIndicatorF(iCloc), includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addSlipBoundary(
  SuperGeometry3D<T>& superGeometry, int material)
{
  addSlipBoundary(superGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addPartialSlipBoundary(
  T tuner, FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->addPartialSlipBoundary(
      tuner, indicator->getExtendedBlockIndicatorF(iCloc), includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addPartialSlipBoundary(
  T tuner, SuperGeometry3D<T>& superGeometry, int material)
{
  addPartialSlipBoundary(tuner, superGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addWallFunctionBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
  UnitConverter<T, DESCRIPTOR> const& converter,
  wallFunctionParam<T> const& wallFunctionParam,
  IndicatorF3D<T>* geoIndicator)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->addWallFunctionBoundary(
      indicator->getExtendedBlockIndicatorF(iCloc),
      converter, wallFunctionParam, geoIndicator, includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addWallFunctionBoundary(
  SuperGeometry3D<T>& superGeometry, int material,
  UnitConverter<T, DESCRIPTOR> const& converter,
  wallFunctionParam<T> const& wallFunctionParam,
  IndicatorF3D<T>* geoIndicator)
{
  addWallFunctionBoundary(superGeometry.getMaterialIndicator(material),
                          converter,
                          wallFunctionParam,
                          geoIndicator);
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addPressureBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T omega)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->addPressureBoundary(indicator->getExtendedBlockIndicatorF(iCloc),
                                          omega, includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addPressureBoundary(
  SuperGeometry3D<T>& superGeometry, int material, T omega)
{
  addPressureBoundary(superGeometry.getMaterialIndicator(material), omega);
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addConvectionBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T omega, T* uAv)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->addConvectionBoundary(indicator->getExtendedBlockIndicatorF(iCloc),
                                            omega, uAv, includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addConvectionBoundary(
  SuperGeometry3D<T>& superGeometry, int material, T omega, T* uAv)
{
  addConvectionBoundary(superGeometry.getMaterialIndicator(material), omega, uAv);
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addConvectionBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _ADblockBCs[iCloc]->addConvectionBoundary(indicator->getExtendedBlockIndicatorF(iCloc), includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addConvectionBoundary(
  SuperGeometry3D<T>& superGeometry, int material)
{
  addConvectionBoundary(superGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addTemperatureBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T omega)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); iCloc++) {
    _ADblockBCs[iCloc]->addTemperatureBoundary(
      indicator->getExtendedBlockIndicatorF(iCloc), omega, includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addTemperatureBoundary(
  SuperGeometry3D<T>& superGeometry, int material, T omega)
{
  addTemperatureBoundary(superGeometry.getMaterialIndicator(material),
                         omega);
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addExtFieldBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, int offset)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _ADblockBCs[iCloc]->addExtFieldBoundary(indicator->getExtendedBlockIndicatorF(iCloc),
                                            offset, includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addExtFieldBoundary(
  SuperGeometry3D<T>& superGeometry, int material, int offset)
{
  addExtFieldBoundary(superGeometry.getMaterialIndicator(material),
                      offset);
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addZeroDistributionBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _ADblockBCs[iCloc]->addZeroDistributionBoundary(
      indicator->getExtendedBlockIndicatorF(iCloc), includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addZeroDistributionBoundary(
  SuperGeometry3D<T>& superGeometry, int material)
{
  addZeroDistributionBoundary(superGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyWallBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T alpha, T kappa1, T kappa2, T h1, T h2, int latticeNumber)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  T addend = 0;
  if(latticeNumber==1)
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (h2/kappa2) );
  else if(latticeNumber==2)
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (-h2/kappa2) );
  else if(latticeNumber==3)
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (h2/kappa2) );
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->addFreeEnergyWallBoundary(
      indicator->getExtendedBlockIndicatorF(iCloc), addend, latticeNumber, includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyWallBoundary(
  SuperGeometry3D<T>& superGeometry, int material, T alpha, T kappa1, T kappa2, T h1, T h2, int latticeNumber)
{
  addFreeEnergyWallBoundary(superGeometry.getMaterialIndicator(material),
		  alpha, kappa1, kappa2, h1, h2, latticeNumber);
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyWallBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T alpha,
  T kappa1, T kappa2, T kappa3, T h1, T h2, T h3, int latticeNumber)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  T addend = 0;
  if(latticeNumber==1)
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (h2/kappa2) + (h3/kappa3) );
  else if(latticeNumber==2)
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (-h2/kappa2) );
  else if(latticeNumber==3)
    addend = 1./(alpha*alpha) * ( (h3/kappa3) );
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->addFreeEnergyWallBoundary(
      indicator->getExtendedBlockIndicatorF(iCloc), addend, latticeNumber, includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyWallBoundary(
  SuperGeometry3D<T>& superGeometry, int material, T alpha,
  T kappa1, T kappa2, T kappa3, T h1, T h2, T h3, int latticeNumber)
{
  addFreeEnergyWallBoundary(superGeometry.getMaterialIndicator(material),
		  alpha, kappa1, kappa2, kappa3, h1, h2, h3, latticeNumber);
}


template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyInletBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T omega, std::string type, int latticeNumber)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->addFreeEnergyInletBoundary(
      indicator->getExtendedBlockIndicatorF(iCloc), omega, type, latticeNumber, includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyInletBoundary(
  SuperGeometry3D<T>& superGeometry, int material, T omega, std::string type, int latticeNumber)
{
  addFreeEnergyInletBoundary(superGeometry.getMaterialIndicator(material), omega, type, latticeNumber);
}


template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyOutletBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T omega, std::string type, int latticeNumber)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->addFreeEnergyOutletBoundary(
      indicator->getExtendedBlockIndicatorF(iCloc), omega, type, latticeNumber, includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addFreeEnergyOutletBoundary(
  SuperGeometry3D<T>& superGeometry, int material, T omega, std::string type, int latticeNumber)
{
  addFreeEnergyOutletBoundary(superGeometry.getMaterialIndicator(material), omega, type, latticeNumber);
}


template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addPoints2CommBC(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator)
{
  if (_overlap == 0) {
    return;
  }

  SuperGeometry3D<T>& superGeometry = indicator->getSuperGeometry();
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    const int nX = superGeometry.getBlockGeometry(iCloc).getNx();
    const int nY = superGeometry.getBlockGeometry(iCloc).getNy();
    const int nZ = superGeometry.getBlockGeometry(iCloc).getNz();

    for (int iX = -_overlap; iX < nX+_overlap; ++iX) {
      for (int iY = -_overlap; iY < nY+_overlap; ++iY) {
        for (int iZ = -_overlap; iZ < nZ+_overlap; ++iZ) {
          if (iX < 0 || iX > nX - 1 ||
              iY < 0 || iY > nY - 1 ||
              iZ < 0 || iZ > nZ - 1 ) { // if within overlap
            if (superGeometry.getBlockGeometry(iCloc).getMaterial(iX,iY,iZ) != 0) {
              bool found = false;
              for (int iXo = -_overlap; iXo <= _overlap && !found; ++iXo) {
                for (int iYo = -_overlap; iYo <= _overlap && !found; ++iYo) {
                  for (int iZo = -_overlap; iZo <= _overlap && !found; ++iZo) {
                    const int nextX = iXo + iX;
                    const int nextY = iYo + iY;
                    const int nextZ = iZo + iZ;
                    if (indicator->getBlockIndicatorF(iCloc)(nextX, nextY, nextZ)
                        && nextX >= -_overlap && nextX < nX+_overlap
                        && nextY >= -_overlap && nextY < nY+_overlap
                        && nextZ >= -_overlap && nextZ < nZ+_overlap) {
                      _sLattice.get_commBC().add_cell(iCloc, iX, iY, iZ);
                      found = true;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::addPoints2CommBC(SuperGeometry3D<T>& superGeometry, int material)
{
  addPoints2CommBC(superGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
SuperLattice3D<T, DESCRIPTOR>& sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::getSuperLattice()
{
  return _sLattice;
}

template<typename T, typename DESCRIPTOR>
std::vector<OnLatticeBoundaryCondition3D<T, DESCRIPTOR>*>& sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::getBlockBCs()
{
  return _blockBCs;
}

template<typename T, typename DESCRIPTOR>
std::vector<OnLatticeAdvectionDiffusionBoundaryCondition3D<T, DESCRIPTOR>*>& sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::getADblockBCs()
{
  return _ADblockBCs;
}

template<typename T, typename DESCRIPTOR>
int sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::getOverlap()
{
  return _overlap;
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::setOverlap(int overlap)
{
  _overlap = overlap;
}

//////////////// Output functions //////////////////////////////////
template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::outputOn()
{
  _output = true;
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->outputOn();
  }
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>::outputOff()
{
  _output = false;
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->outputOff();
  }
}


////////////////// Factory functions //////////////////////////////////

template<typename T, typename DESCRIPTOR>
void createLocalBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& sBC)
{
  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(0);
  for (int iC = 0; iC < nC; iC++) {
    OnLatticeBoundaryCondition3D<T, DESCRIPTOR>* blockBC =
      createLocalBoundaryCondition3D(
        sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getBlockBCs().push_back(blockBC);
  }
}

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void createInterpBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& sBC)
{
  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(1);
  for (int iC = 0; iC < nC; iC++) {
    OnLatticeBoundaryCondition3D<T, DESCRIPTOR>* blockBC =
      createInterpBoundaryCondition3D<T,DESCRIPTOR,MixinDynamics>(
        sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getBlockBCs().push_back(blockBC);
  }
}

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void createExtFdBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& sBC)
{
  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(1);
  for (int iC = 0; iC < nC; iC++) {
    OnLatticeBoundaryCondition3D<T, DESCRIPTOR>* blockBC =
      createExtendedFdBoundaryCondition3D<T,DESCRIPTOR,MixinDynamics>(
        sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getBlockBCs().push_back(blockBC);
  }
}



} // namespace olb

#endif
