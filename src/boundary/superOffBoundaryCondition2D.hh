/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2016 Jonas Kratzke, Mathias J. Krause
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


#ifndef SUPER_OFF_BOUNDARY_CONDITION_2D_HH
#define SUPER_OFF_BOUNDARY_CONDITION_2D_HH

#include <vector>
#include <list>

#include "offBoundaryCondition2D.h"
#include "superOffBoundaryCondition2D.h"
#include "core/superLattice2D.h"
#include "core/util.h"
#include "functors/analytical/analyticalF.h"
#include "functors/lattice/indicator/superIndicatorBaseF2D.h"

namespace olb {

///////// class superOffBoundaryCondition2D ///////////////////////////////

template<typename T, typename DESCRIPTOR>
sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
sOffLatticeBoundaryCondition2D (SuperLattice2D<T,DESCRIPTOR>& sLattice, T epsFraction )
  : clout(std::cout,"sOffLatticeBoundaryCondition2D"),
    _sLattice(sLattice),
    _epsFraction(epsFraction),
    _output(false)
{}

template<typename T, typename DESCRIPTOR>
sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
sOffLatticeBoundaryCondition2D(sOffLatticeBoundaryCondition2D<T,DESCRIPTOR> const& rhs)
  : clout(std::cout,"sOffLatticeBoundaryCondition2D"),
    _sLattice(rhs._sLattice),
    _epsFraction(rhs._epsFraction),
    _output(false)
{
  _blockBCs = rhs._blockBCs;
  _overlap = rhs._overlap;
}

template<typename T, typename DESCRIPTOR>
sOffLatticeBoundaryCondition2D<T,DESCRIPTOR> sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::operator=(
  sOffLatticeBoundaryCondition2D<T,DESCRIPTOR> rhs)
{

  sOffLatticeBoundaryCondition2D<T,DESCRIPTOR> tmp(rhs);
  return tmp;
}

template<typename T, typename DESCRIPTOR>
sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
~sOffLatticeBoundaryCondition2D()
{

  for (unsigned iC=0; iC<_blockBCs.size(); iC++) {
    delete _blockBCs[iC];
  }
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
addVelocityBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& boundaryIndicator,
                    FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
                    IndicatorF2D<T>&                   geometryIndicator)
{
  clout << "epsFraction=" << _epsFraction << std::endl;
  clout.setMultiOutput(true);
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    clout << "Cuboid globiC " << _sLattice.getLoadBalancer().glob(iCloc)
          << " starts to read distances for Velocity Boundary..." << std::endl;
    _blockBCs[iCloc]->addVelocityBoundary(boundaryIndicator->getExtendedBlockIndicatorF(iCloc),
                                          bulkIndicator->getExtendedBlockIndicatorF(iCloc),
                                          geometryIndicator);
    clout << "Cuboid globiC " << _sLattice.getLoadBalancer().glob(iCloc)
          << " finished reading distances for Velocity Boundary." << std::endl;
  }
  clout.setMultiOutput(false);
  addPoints2CommBC(std::forward<decltype(boundaryIndicator)>(boundaryIndicator));
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
addVelocityBoundary(SuperGeometry2D<T>& superGeometry, int material,
                    IndicatorF2D<T>& geometryIndicator,
                    std::vector<int> bulkMaterials)
{
  addVelocityBoundary(
    superGeometry.getMaterialIndicator(material),
    superGeometry.getMaterialIndicator(std::move(bulkMaterials)),
    geometryIndicator);
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
addVelocityBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& boundaryIndicator,
                    FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator)
{
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->addVelocityBoundary(boundaryIndicator->getExtendedBlockIndicatorF(iCloc),
                                          bulkIndicator->getExtendedBlockIndicatorF(iCloc));
  }
  addPoints2CommBC(std::forward<decltype(boundaryIndicator)>(boundaryIndicator));
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
addVelocityBoundary(SuperGeometry2D<T>& superGeometry, int material,
                    std::vector<int> bulkMaterials)
{
  addVelocityBoundary(superGeometry.getMaterialIndicator(material),
                      superGeometry.getMaterialIndicator(std::move(bulkMaterials)));
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
addZeroVelocityBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& boundaryIndicator,
                        FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
                        IndicatorF2D<T>&                   geometryIndicator)
{
  clout << "epsFraction=" << _epsFraction << std::endl;
  clout.setMultiOutput(true);
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    clout << "Cuboid globiC " << _sLattice.getLoadBalancer().glob(iCloc)
          << " starts to read distances for ZeroVelocity Boundary..." << std::endl;
    _blockBCs[iCloc]->addZeroVelocityBoundary(
      boundaryIndicator->getExtendedBlockIndicatorF(iCloc),
      bulkIndicator->getExtendedBlockIndicatorF(iCloc),
      geometryIndicator);
    clout << "Cuboid globiC " << _sLattice.getLoadBalancer().glob(iCloc)
          << " finished reading distances for ZeroVelocity Boundary." << std::endl;
  }
  clout.setMultiOutput(false);
  addPoints2CommBC(std::forward<decltype(boundaryIndicator)>(boundaryIndicator));
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
addZeroVelocityBoundary(SuperGeometry2D<T>& superGeometry, int material,
                        IndicatorF2D<T>& geometryIndicator,
                        std::vector<int> bulkMaterials)
{
  addZeroVelocityBoundary(superGeometry.getMaterialIndicator(material),
                          superGeometry.getMaterialIndicator(std::move(bulkMaterials)),
                          geometryIndicator);
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
addZeroVelocityBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& boundaryIndicator,
                        FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator)
{
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->addZeroVelocityBoundary(
      boundaryIndicator->getExtendedBlockIndicatorF(iCloc),
      bulkIndicator->getExtendedBlockIndicatorF(iCloc));
  }
  addPoints2CommBC(std::forward<decltype(boundaryIndicator)>(boundaryIndicator));
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
addZeroVelocityBoundary(SuperGeometry2D<T>& superGeometry, int material,
                        std::vector<int> bulkMaterials)
{
  addZeroVelocityBoundary(superGeometry.getMaterialIndicator(material),
                          superGeometry.getMaterialIndicator(std::move(bulkMaterials)));
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
addPressureBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& boundaryIndicator,
                    FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
                    IndicatorF2D<T>&                   geometryIndicator)
{
  clout << "epsFraction=" << _epsFraction << std::endl;
  clout.setMultiOutput(true);
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    clout << "Cuboid globiC " << _sLattice.getLoadBalancer().glob(iCloc)
          << " starts to read distances for Pressure Boundary..." << std::endl;
    _blockBCs[iCloc]->addPressureBoundary(boundaryIndicator->getExtendedBlockIndicatorF(iCloc),
                                          bulkIndicator->getExtendedBlockIndicatorF(iCloc),
                                          geometryIndicator);
    clout << "Cuboid globiC " << _sLattice.getLoadBalancer().glob(iCloc)
          << " finished reading distances for Pressure Boundary." << std::endl;
  }
  clout.setMultiOutput(false);
  addPoints2CommBC(std::forward<decltype(boundaryIndicator)>(boundaryIndicator));
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
addPressureBoundary(SuperGeometry2D<T>& superGeometry, int material,
                    IndicatorF2D<T>& geometryIndicator,
                    std::vector<int> bulkMaterials)
{
  addPressureBoundary(superGeometry.getMaterialIndicator(material),
                      superGeometry.getMaterialIndicator(std::move(bulkMaterials)),
                      geometryIndicator);
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
addPressureBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& boundaryIndicator,
                    FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator)
{
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->addPressureBoundary(boundaryIndicator->getExtendedBlockIndicatorF(iCloc),
                                          bulkIndicator->getExtendedBlockIndicatorF(iCloc));
  }
  addPoints2CommBC(std::forward<decltype(boundaryIndicator)>(boundaryIndicator));
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
addPressureBoundary(SuperGeometry2D<T>& superGeometry, int material,
                    std::vector<int> bulkMaterials)
{
  addPressureBoundary(superGeometry.getMaterialIndicator(material),
                      superGeometry.getMaterialIndicator(std::move(bulkMaterials)));
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
defineU(FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
        FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
        AnalyticalF2D<T,T>& u)
{
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->defineU(indicator->getExtendedBlockIndicatorF(iCloc),
                              bulkIndicator->getExtendedBlockIndicatorF(iCloc),
                              u);
  }
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
defineU(SuperGeometry2D<T>& superGeometry, int material,
        AnalyticalF2D<T,T>& u, std::vector<int> bulkMaterials)
{
  defineU(superGeometry.getMaterialIndicator(material),
          superGeometry.getMaterialIndicator(std::move(bulkMaterials)),
          u);
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
defineRho(FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
          FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
          AnalyticalF2D<T,T>&                rho)
{
  for (int iCloc = 0; iCloc < _sLattice.getLoadBalancer().size(); ++iCloc) {
    _blockBCs[iCloc]->defineRho(indicator->getExtendedBlockIndicatorF(iCloc),
                                bulkIndicator->getExtendedBlockIndicatorF(iCloc),
                                rho);
  }
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
defineRho(SuperGeometry2D<T>& superGeometry, int material,
          AnalyticalF2D<T,T>& rho, std::vector<int> bulkMaterials)
{
  defineRho(superGeometry.getMaterialIndicator(material),
            superGeometry.getMaterialIndicator(std::move(bulkMaterials)),
            rho);
}


template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
addPoints2CommBC(FunctorPtr<SuperIndicatorF2D<T>>&& indicator)
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
            iY < 0 || iY > nY - 1 ) { // is inside boundary
          bool found = false;
          if (superGeometry.getBlockGeometry(iCloc).getMaterial(iX,iY) != 0) {
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
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::
addPoints2CommBC(SuperGeometry2D<T>& superGeometry, int material)
{
  addPoints2CommBC(superGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::outputOn()
{
  _output = true;
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->outputOn();
  }
}

template<typename T, typename DESCRIPTOR>
void sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>::outputOff()
{
  _output = false;
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->outputOff();
  }
}

////////////////// Factory functions //////////////////////////////////

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void createBouzidiBoundaryCondition2D(sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>& sBC)
{

  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(1);
  for (int iC=0; iC<nC; iC++) {
    OffLatticeBoundaryCondition2D<T,DESCRIPTOR>* blockBC
      = createBouzidiBoundaryCondition2D<T,DESCRIPTOR,MixinDynamics>(sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getBlockBCs().push_back(blockBC);
  }
}

template<typename T, typename DESCRIPTOR>
void createBounceBackBoundaryCondition2D(sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>& sBC)
{

  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(1);
  for (int iC=0; iC<nC; iC++) {
    OffLatticeBoundaryCondition2D<T,DESCRIPTOR>* blockBC
      = createBouzidiBoundaryCondition2D<T,DESCRIPTOR>(sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getBlockBCs().push_back(blockBC);
  }
}

}  // namespace olb

#endif
