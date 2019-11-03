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
 * A helper for initialising 2D boundaries -- generic implementation.
 */

#ifndef SUPER_BOUNDARY_CONDITION_2D_HH
#define SUPER_BOUNDARY_CONDITION_2D_HH

#include <vector>
#include "boundaryCondition2D.h"
#include "extendedFiniteDifferenceBoundary2D.h"
#include "superBoundaryCondition2D.h"
#include "core/superLattice2D.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"

namespace olb {

///////// class superBoundaryCondition2D ///////////////////////////////

template<typename T, typename DESCRIPTOR>
sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::sOnLatticeBoundaryCondition2D(
  SuperLattice2D<T, DESCRIPTOR>& sLattice) :
  clout(std::cout,"sOnLatticeBoundaryCondition2D"),
  _sLattice(sLattice),
  _output(false)
{
}

template<typename T, typename DESCRIPTOR>
sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::sOnLatticeBoundaryCondition2D(
  sOnLatticeBoundaryCondition2D<T, DESCRIPTOR> const& rhs) :
  clout(std::cout,"sOnLatticeBoundaryCondition2D"),
  _sLattice(rhs._sLattice),
  _output(false)
{

  _blockBCs = rhs._blockBCs;
  _overlap = rhs._overlap;
}

template<typename T, typename DESCRIPTOR>
sOnLatticeBoundaryCondition2D<T, DESCRIPTOR> sOnLatticeBoundaryCondition2D<T,
                              DESCRIPTOR>::operator=(sOnLatticeBoundaryCondition2D<T, DESCRIPTOR> rhs)
{

  sOnLatticeBoundaryCondition2D<T, DESCRIPTOR> tmp(rhs);
  return tmp;
}

template<typename T, typename DESCRIPTOR>
sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::~sOnLatticeBoundaryCondition2D()
{
  //for (unsigned iC = 0; iC < _blockBCs.size(); iC++) {
  //  delete _blockBCs[iC];
  //}
}


template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addVelocityBoundary(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T omega)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->addVelocityBoundary(
      indicator->getExtendedBlockIndicatorF(iCloc), omega, includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addVelocityBoundary(
  SuperGeometry2D<T>& superGeometry, int material, T omega)
{
  addVelocityBoundary(superGeometry.getMaterialIndicator(material), omega);
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addSlipBoundary(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->addSlipBoundary(indicator->getExtendedBlockIndicatorF(iCloc), includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addSlipBoundary(
  SuperGeometry2D<T>& superGeometry, int material)
{
  addSlipBoundary(superGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addPartialSlipBoundary(
  T tuner, FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
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
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addPartialSlipBoundary(
  T tuner, SuperGeometry2D<T>& superGeometry, int material)
{
  addPartialSlipBoundary(tuner, superGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addTemperatureBoundary(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T omega)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _ADblockBCs[iCloc]->addTemperatureBoundary(
      indicator->getExtendedBlockIndicatorF(iCloc), omega, includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addTemperatureBoundary(
  SuperGeometry2D<T>& superGeometry, int material, T omega)
{
  addTemperatureBoundary(superGeometry.getMaterialIndicator(material), omega);
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addRegularizedTemperatureBoundary(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T omega)
{
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _ADblockBCs[iCloc]->addRegularizedTemperatureBoundary(
      indicator->getExtendedBlockIndicatorF(iCloc), omega);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addRegularizedTemperatureBoundary(
  SuperGeometry2D<T>& superGeometry, int material, T omega)
{
  addRegularizedTemperatureBoundary(superGeometry.getMaterialIndicator(material), omega);
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addRegularizedHeatFluxBoundary(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T omega, T *heatFlux)
{
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _ADblockBCs[iCloc]->addRegularizedHeatFluxBoundary(
      indicator->getExtendedBlockIndicatorF(iCloc), omega, heatFlux);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addRegularizedHeatFluxBoundary(
  SuperGeometry2D<T>& superGeometry, int material, T omega, T *heatFlux)
{
  addRegularizedHeatFluxBoundary(superGeometry.getMaterialIndicator(material), omega, heatFlux);
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addPressureBoundary(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T omega)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->addPressureBoundary(
      indicator->getExtendedBlockIndicatorF(iCloc), omega, includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addPressureBoundary(
  SuperGeometry2D<T>& superGeometry, int material, T omega)
{
  addPressureBoundary(superGeometry.getMaterialIndicator(material), omega);
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addConvectionBoundary(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T omega, T* uAv)
{
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->addConvectionBoundary(
      indicator->getExtendedBlockIndicatorF(iCloc), omega, uAv, includeOuterCells);
  }
  addPoints2CommBC(std::forward<decltype(indicator)>(indicator));
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addConvectionBoundary(
  SuperGeometry2D<T>& superGeometry, int material, T omega, T* uAv)
{
  addConvectionBoundary(superGeometry.getMaterialIndicator(material),
                        omega, uAv);
}


template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyWallBoundary(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T alpha, T kappa1, T kappa2, T h1, T h2, int latticeNumber)
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
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyWallBoundary(
  SuperGeometry2D<T>& superGeometry, int material, T alpha, T kappa1, T kappa2, T h1, T h2, int latticeNumber)
{
  addFreeEnergyWallBoundary(superGeometry.getMaterialIndicator(material),
                            alpha, kappa1, kappa2, h1, h2, latticeNumber);
}


template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyWallBoundary(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T alpha,
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
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyWallBoundary(
  SuperGeometry2D<T>& superGeometry, int material, T alpha,
  T kappa1, T kappa2, T kappa3, T h1, T h2, T h3, int latticeNumber)
{
  addFreeEnergyWallBoundary(superGeometry.getMaterialIndicator(material),
                            alpha, kappa1, kappa2, kappa3, h1, h2, h3, latticeNumber);
}


template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyInletBoundary(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T omega, std::string type, int latticeNumber)
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
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyInletBoundary(
  SuperGeometry2D<T>& superGeometry, int material, T omega, std::string type, int latticeNumber)
{
  addFreeEnergyInletBoundary(superGeometry.getMaterialIndicator(material), omega, type, latticeNumber);
}


template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyOutletBoundary(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T omega, std::string type, int latticeNumber)
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
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addFreeEnergyOutletBoundary(
  SuperGeometry2D<T>& superGeometry, int material, T omega, std::string type, int latticeNumber)
{
  addFreeEnergyOutletBoundary(superGeometry.getMaterialIndicator(material), omega, type, latticeNumber);
}


template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addPoints2CommBC(
  FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
{
  if (_overlap == 0) {
    return;
  }

  SuperGeometry2D<T>& superGeometry = indicator->getSuperGeometry();
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    const int nX = superGeometry.getBlockGeometry(iCloc).getNx();
    const int nY = superGeometry.getBlockGeometry(iCloc).getNy();

    for (int iX = -_overlap; iX < nX+_overlap; ++iX) {
      for (int iY = -_overlap; iY < nY+_overlap; ++iY) {
        if (iX < 0 || iX > nX - 1 ||
            iY < 0 || iY > nY - 1 ) { // if within overlap
          if (superGeometry.getBlockGeometry(iCloc).getMaterial(iX,iY) != 0) {
            bool found = false;
            for (int iXo = -_overlap; iXo <= _overlap && !found; ++iXo) {
              for (int iYo = -_overlap; iYo <= _overlap && !found; ++iYo) {
                const int nextX = iXo + iX;
                const int nextY = iYo + iY;
                if (indicator->getBlockIndicatorF(iCloc)(nextX, nextY)) {
                  _sLattice.get_commBC().add_cell(iCloc, iX, iY);
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

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::addPoints2CommBC(
  SuperGeometry2D<T>& superGeometry, int material)
{
  addPoints2CommBC(superGeometry.getMaterialIndicator(material));
}

////////////////// Factory functions //////////////////////////////////

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void createLocalBoundaryCondition2D(
  sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>& sBC)
{

  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(0);
  for (int iC = 0; iC < nC; iC++) {
    OnLatticeBoundaryCondition2D<T, DESCRIPTOR>* blockBC =
      createLocalBoundaryCondition2D(
        sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getBlockBCs().push_back(blockBC);
  }
}

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void createInterpBoundaryCondition2D(
  sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>& sBC)
{

  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(1);
  for (int iC = 0; iC < nC; iC++) {
    OnLatticeBoundaryCondition2D<T, DESCRIPTOR>* blockBC =
      createInterpBoundaryCondition2D<T,DESCRIPTOR,MixinDynamics>(
        sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getBlockBCs().push_back(blockBC);
  }
}

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void createExtFdBoundaryCondition2D(
  sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>& sBC)
{

  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(1);
  for (int iC = 0; iC < nC; iC++) {
    OnLatticeBoundaryCondition2D<T, DESCRIPTOR>* blockBC =
      createExtendedFdBoundaryCondition2D<T,DESCRIPTOR,MixinDynamics>(
        sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getBlockBCs().push_back(blockBC);
  }
}

//////////////// Output functions //////////////////////////////////
template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::outputOn()
{
  _output = true;
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->outputOn();
  }
}

template<typename T, typename DESCRIPTOR>
void sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>::outputOff()
{
  _output = false;
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->outputOff();
  }
}

} // namespace olb

#endif
