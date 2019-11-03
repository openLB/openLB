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
 * A helper for initialising 3D boundaries -- generic implementation.
 */


#ifndef SUPER_OFF_BOUNDARY_CONDITION_3D_HH
#define SUPER_OFF_BOUNDARY_CONDITION_3D_HH

#include <vector>
#include <list>
#include "offBoundaryCondition3D.h"
#include "superOffBoundaryCondition3D.h"
#include "core/superLattice3D.h"
#include "core/util.h"
#include "functors/analytical/analyticalF.h"
#include "functors/lattice/indicator/superIndicatorF3D.h"

namespace olb {

///////// class superOffBoundaryCondition3D ///////////////////////////////

template<typename T, typename DESCRIPTOR>
sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::sOffLatticeBoundaryCondition3D(
  SuperLattice3D<T,DESCRIPTOR>& sLattice, T epsFraction )
  : clout(std::cout,"sOffLatticeBoundaryCondition3D"),
    _sLattice(sLattice),
    _epsFraction(epsFraction),
    _output(false)
{}

template<typename T, typename DESCRIPTOR>
sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::sOffLatticeBoundaryCondition3D(
  sOffLatticeBoundaryCondition3D<T,DESCRIPTOR> const& rhs)
  : clout(std::cout,"sOffLatticeBoundaryCondition3D"),
    _sLattice(rhs._sLattice),
    _epsFraction(rhs._epsFraction),
    _output(false)
{
  _blockBCs = rhs._blockBCs;
  _overlap = rhs._overlap;
}

template<typename T, typename DESCRIPTOR>
sOffLatticeBoundaryCondition3D<T,DESCRIPTOR> sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::operator=(
  sOffLatticeBoundaryCondition3D<T,DESCRIPTOR> rhs)
{
  sOffLatticeBoundaryCondition3D<T,DESCRIPTOR> tmp(rhs);
  return tmp;
}

template<typename T, typename DESCRIPTOR>
sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::~sOffLatticeBoundaryCondition3D()
{
  for (unsigned iC=0; iC<_blockBCs.size(); iC++) {
    delete _blockBCs[iC];
  }
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::addVelocityBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator,
  FunctorPtr<SuperIndicatorF3D<T>>&& bulkIndicator,
  IndicatorF3D<T>&                   geometryIndicator)
{
  if (_output) {
    clout << "epsFraction=" << _epsFraction << std::endl;
    clout.setMultiOutput(true);
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    if (_output) {
      clout << "Cuboid globiC " << _sLattice.getLoadBalancer().glob(iCloc)
            << " starts to read distances for Velocity Boundary..." << std::endl;
    }
    _blockBCs[iCloc]->addVelocityBoundary(
      boundaryIndicator->getExtendedBlockIndicatorF(iCloc),
      bulkIndicator->getExtendedBlockIndicatorF(iCloc),
      geometryIndicator);
    if (_output) {
      clout << "Cuboid globiC " << _sLattice.getLoadBalancer().glob(iCloc)
            << " finished reading distances for Velocity Boundary." << std::endl;
    }
  }
  if (_output) {
    clout.setMultiOutput(false);
  }
  addPoints2CommBC(std::forward<decltype(boundaryIndicator)>(boundaryIndicator));
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::addVelocityBoundary(
  SuperGeometry3D<T>& superGeometry, int material,
  IndicatorF3D<T>& geometryIndicator, std::vector<int> bulkMaterials)
{
  addVelocityBoundary(superGeometry.getMaterialIndicator(material),
                      superGeometry.getMaterialIndicator(std::move(bulkMaterials)),
                      geometryIndicator);
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::addZeroVelocityBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator,
  FunctorPtr<SuperIndicatorF3D<T>>&& bulkIndicator,
  IndicatorF3D<T>&                   geometryIndicator)
{
  if (_output) {
    clout << "epsFraction=" << _epsFraction << std::endl;
    clout.setMultiOutput(true);
  }
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    if (_output) {
      clout << "Cuboid globiC " << _sLattice.getLoadBalancer().glob(iCloc)
            << " starts to read distances for ZeroVelocity Boundary..." << std::endl;
    }
    _blockBCs[iCloc]->addZeroVelocityBoundary(
      boundaryIndicator->getExtendedBlockIndicatorF(iCloc),
      bulkIndicator->getExtendedBlockIndicatorF(iCloc),
      geometryIndicator);
    if (_output) {
      clout << "Cuboid globiC " << _sLattice.getLoadBalancer().glob(iCloc)
            << " finished reading distances for ZeroVelocity Boundary." << std::endl;
    }
  }
  if (_output) {
    clout.setMultiOutput(false);
  }
  addPoints2CommBC(std::forward<decltype(boundaryIndicator)>(boundaryIndicator));
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::addZeroVelocityBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator,
  IndicatorF3D<T>& geometryIndicator, std::vector<int> bulkMaterials)
{
  SuperGeometry3D<T>& superGeometry = boundaryIndicator->getSuperGeometry();
  addZeroVelocityBoundary(
    std::forward<decltype(boundaryIndicator)>(boundaryIndicator),
    superGeometry.getMaterialIndicator(std::move(bulkMaterials)),
    geometryIndicator);
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::addZeroVelocityBoundary(
  SuperGeometry3D<T>& superGeometry, int material,
  IndicatorF3D<T>& geometryIndicator, std::vector<int> bulkMaterials)
{
  addZeroVelocityBoundary(superGeometry.getMaterialIndicator(material),
                          geometryIndicator,
                          bulkMaterials);
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::defineU(
  FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator,
  FunctorPtr<SuperIndicatorF3D<T>>&& bulkIndicator,
  AnalyticalF3D<T,T>& u)
{
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->defineU(boundaryIndicator->getExtendedBlockIndicatorF(iCloc),
                              bulkIndicator->getExtendedBlockIndicatorF(iCloc),
                              u);
  }
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::defineU(
  FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator,
  AnalyticalF3D<T,T>& u, std::vector<int> bulkMaterials)
{
  defineU(std::forward<decltype(boundaryIndicator)>(boundaryIndicator),
          boundaryIndicator->getSuperGeometry().getMaterialIndicator(std::move(bulkMaterials)),
          u);
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::defineU(
  SuperGeometry3D<T>& superGeometry, int material,
  AnalyticalF3D<T,T>& u, std::vector<int> bulkMaterials)
{
  defineU(superGeometry.getMaterialIndicator(material), u, bulkMaterials);
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::addPoints2CommBC(
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

    for (int iX=-_overlap; iX<nX+_overlap; ++iX) {
      for (int iY=-_overlap; iY<nY+_overlap; ++iY) {
        for (int iZ=-_overlap; iZ<nZ+_overlap; ++iZ) {
          if (iX < 0 || iX > nX - 1 ||
              iY < 0 || iY > nY - 1 ||
              iZ < 0 || iZ > nZ - 1 ) { // is inside boundary
            int found = false;
            if (superGeometry.getBlockGeometry(iCloc).getMaterial(iX,iY,iZ) != 0) {
              for (int iXo=-_overlap; iXo<=_overlap && !found; ++iXo) {
                for (int iYo=-_overlap; iYo<=_overlap && !found; ++iYo) {
                  for (int iZo=-_overlap; iZo<=_overlap && !found; ++iZo) {
                    const int nextX = iXo + iX;
                    const int nextY = iYo + iY;
                    const int nextZ = iZo + iZ;
                    if (indicator->getBlockIndicatorF(iCloc)(nextX, nextY, nextZ)) {
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
void sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::addPoints2CommBC(
  SuperGeometry3D<T>& superGeometry, int material)
{
  addPoints2CommBC(superGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
SuperLattice3D<T,DESCRIPTOR>& sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::getSuperLattice()
{
  return _sLattice;
}

template<typename T, typename DESCRIPTOR>
std::vector<OffLatticeBoundaryCondition3D<T,DESCRIPTOR>* >& sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::getBlockBCs()
{
  return _blockBCs;
}

template<typename T, typename DESCRIPTOR>
int sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::getOverlap()
{
  return _overlap;
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::setOverlap(int overlap)
{
  _overlap = overlap;
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::outputOn()
{
  _output = true;
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->outputOn();
  }
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>::outputOff()
{
  _output = false;
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->outputOff();
  }
}

////////////////// Factory functions //////////////////////////////////

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void createBouzidiBoundaryCondition3D(sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& sBC)
{

  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(1);
  for (int iC=0; iC<nC; iC++) {
    OffLatticeBoundaryCondition3D<T,DESCRIPTOR>* blockBC
      = createBouzidiBoundaryCondition3D<T,DESCRIPTOR,MixinDynamics>(sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getBlockBCs().push_back(blockBC);
  }
}

}  // namespace olb

#endif
