/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2014 Mathias J. Krause, Peter Weisbrod
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
 * Representation of the 3D geometry -- generic implementation.
 */

#ifndef SUPER_GEOMETRY_3D_HH
#define SUPER_GEOMETRY_3D_HH

#include <cmath>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>

#include "geometry/cuboid3D.h"
#include "geometry/cuboidGeometry3D.h"
#include "geometry/superGeometry3D.h"
#include "communication/superStructure3D.h"
#include "communication/loadBalancer.h"
#include "functors/analytical/indicator/indicatorF3D.h"
#include "functors/lattice/indicator/superIndicatorF3D.h"
#include "functors/analytical/indicator/indicCalc3D.h"
#include "io/ostreamManager.h"

namespace olb {

template<typename T>
SuperGeometry3D<T>::SuperGeometry3D(CuboidGeometry3D<T>& cuboidGeometry, LoadBalancer<T>& loadBalancer, int overlap)
  : SuperStructure3D<T>(cuboidGeometry, loadBalancer, overlap),
    _statistics(this),
    clout(std::cout,"SuperGeometry3D")
{
  // init communicator
  this->_communicator.init_nh();
  this->_communicator.add_cells(this->_overlap);
  this->_communicator.init();
  this->_communicationNeeded = true;

  // constructing the block and extended block geometries from the cuboid geometry
  _blockGeometries.clear();

  for (int iCloc=0; iCloc<this->getLoadBalancer().size(); iCloc++) {
    int iCglob = this->getLoadBalancer().glob(iCloc);
    Cuboid3D<T> extendedCuboid(cuboidGeometry.get(iCglob),overlap);
    _extendedBlockGeometries.emplace_back(extendedCuboid, iCglob);
  }

  for (int iCloc=0; iCloc<this->getLoadBalancer().size(); iCloc++) {
    int iCglob = this->getLoadBalancer().glob(iCloc);
    int nX = cuboidGeometry.get(iCglob).getNx();
    int nY = cuboidGeometry.get(iCglob).getNy();
    int nZ = cuboidGeometry.get(iCglob).getNz();
    _blockGeometries.emplace_back(_extendedBlockGeometries[iCloc],
                                  overlap, overlap+nX-1,
                                  overlap, overlap+nY-1,
                                  overlap, overlap+nZ-1);
  }
  _statistics.getStatisticsStatus() = true;
}

template<typename T>
SuperGeometry3D<T>::SuperGeometry3D(SuperGeometry3D const& rhs)
  : SuperStructure3D<T>(rhs._cuboidGeometry, rhs._loadBalancer, rhs._overlap),
    _statistics(this),
    clout(std::cout,"SuperGeometry3D")
{
  // init communicator
  this->_communicator.init_nh();
  this->_communicator.add_cells(this->_overlap);
  this->_communicator.init();
  this->_communicationNeeded = true;
  // copy block and extended block geometries
  _blockGeometries = rhs._blockGeometries;
  _extendedBlockGeometries = rhs._extendedBlockGeometries;
  _statistics.getStatisticsStatus() = true;
}

template<typename T>
SuperGeometry3D<T>& SuperGeometry3D<T>::operator=(SuperGeometry3D const& rhs)
{

  // init communicator
  this->_communicator.init_nh();
  this->_communicator.add_cells(this->_overlap);
  this->_communicator.init();
  this->_communicationNeeded = true;
  // copy mother data
  this->_cuboidGeometry = rhs._cuboidGeometry;
  this->_loadBalancer = rhs._loadBalancer;
  this->_overlap = rhs._overlap;
  // copy block and extended block geometrie
  _blockGeometries = rhs._blockGeometries;
  _extendedBlockGeometries = rhs._extendedBlockGeometries;
  _statistics = SuperGeometryStatistics3D<T>(this);
  return *this;
}

template<typename T>
bool* SuperGeometry3D<T>::operator() (int iCloc, int iX, int iY, int iZ, int iData)
{
  return (bool*)&getExtendedBlockGeometry(iCloc).get(iX+this->_overlap, iY+this->_overlap, iZ+this->_overlap);
}

template<typename T>
int SuperGeometry3D<T>::getDataSize() const
{
  return 1;
}

template<typename T>
int SuperGeometry3D<T>::getDataTypeSize() const
{
  return sizeof(int);
}


template<typename T>
int& SuperGeometry3D<T>::set(int iCglob, int iXloc, int iYloc, int iZloc)
{

  if ( this->getLoadBalancer().rank(iCglob) == singleton::mpi().getRank() ) {
    std::cout << "warning: read only access to data which is not available in";
    this->_communicationNeeded = true;
    _statistics.getStatisticsStatus() = true;
    return _extendedBlockGeometries[this->getLoadBalancer().loc(iCglob)].get(iXloc+this->_overlap, iYloc+this->_overlap, iZloc+this->_overlap);
  }
  else {
    std::cout << "error: write access to data which is not available in the any block geometry";
    exit(-1);
    //return 0;
  }
}

template<typename T>
int const& SuperGeometry3D<T>::get(int iCglob, int iXloc, int iYloc, int iZloc) const
{
  if ( this->getLoadBalancer().rank(iCglob) == singleton::mpi().getRank() ) {
    return _extendedBlockGeometries[this->getLoadBalancer().loc(iCglob)].get(iXloc+this->_overlap, iYloc+this->_overlap, iZloc+this->_overlap);
  }
  else {
    std::cout << "error: read only access to data which is not available in the any block geometry, returning 0 as default" << std::endl;
    exit(-1);
    //return 0;
  }
}

template<typename T>
int SuperGeometry3D<T>::getAndCommunicate(int iCglob, int iXloc, int iYloc, int iZloc) const
{
  int material = 0;
  if ( this->getLoadBalancer().rank(iCglob) == singleton::mpi().getRank() ) {
    material = _extendedBlockGeometries[this->getLoadBalancer().loc(iCglob)].get(iXloc + this->_overlap,
               iYloc + this->_overlap,
               iZloc + this->_overlap);
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().bCast(&material, 1, this->_loadBalancer.rank(iCglob));
#endif
  return material;
}

template<typename T>
int& SuperGeometry3D<T>::set(std::vector<int> latticeR)
{
  return set(latticeR[0], latticeR[1], latticeR[2], latticeR[3]);
}

template<typename T>
int const& SuperGeometry3D<T>::get(std::vector<int> latticeR) const
{
  return get(latticeR[0], latticeR[1], latticeR[2], latticeR[3]);
}

template<typename T>
int const& SuperGeometry3D<T>::get(const int latticeR[4]) const
{
  return get(latticeR[0], latticeR[1], latticeR[2], latticeR[3]);
}

template<typename T>
int SuperGeometry3D<T>::getAndCommunicate(std::vector<int> latticeR) const
{
  return getAndCommunicate(latticeR[0], latticeR[1], latticeR[2], latticeR[3]);
}

template<typename T>
std::vector<T> SuperGeometry3D<T>::getPhysR(int iCglob, int iX, int iY, int iZ) const
{
  T physRv[3];
  this->_cuboidGeometry.getPhysR(physRv, iCglob, iX, iY, iZ);
  std::vector<T> physR(physRv,physRv + 3);
  return physR;
}

template<typename T>
std::vector<T> SuperGeometry3D<T>::getPhysR(std::vector<int> latticeR) const
{
  T physRv[3];
  this->_cuboidGeometry.getPhysR(physRv, latticeR[0], latticeR[1], latticeR[2], latticeR[3]);
  std::vector<T> physR(physRv,physRv + 3);
  return physR;
}

template<typename T>
void  SuperGeometry3D<T>::getPhysR(T physR[3], const int& iCglob,  const int& iX, const int& iY, const int& iZ) const
{
  this->_cuboidGeometry.getPhysR(physR, iCglob, iX, iY, iZ);
}

template<typename T>
void  SuperGeometry3D<T>::getPhysR(T physR[3], const int latticeR[4]) const
{
  this->_cuboidGeometry.getPhysR(physR, latticeR);
}


template<typename T>
BlockGeometry3D<T>& SuperGeometry3D<T>::getExtendedBlockGeometry(int locIC)
{
  _statistics.getStatisticsStatus() = true;
  return _extendedBlockGeometries[locIC];
}

template<typename T>
BlockGeometry3D<T> const& SuperGeometry3D<T>::getExtendedBlockGeometry(int locIC) const
{
  return _extendedBlockGeometries[locIC];
}

template<typename T>
BlockGeometryView3D<T>& SuperGeometry3D<T>::getBlockGeometry(int locIC)
{
  _statistics.getStatisticsStatus() = true;
  return _blockGeometries[locIC];
}

template<typename T>
BlockGeometryView3D<T> const& SuperGeometry3D<T>::getBlockGeometry(int locIC) const
{
  return _blockGeometries[locIC];
}


template<typename T>
SuperGeometryStatistics3D<T>& SuperGeometry3D<T>::getStatistics()
{
  if (this->_communicationNeeded) {
    this->communicate();
    getStatisticsStatus()=true;
  }
  return _statistics;
}

template<typename T>
bool& SuperGeometry3D<T>::getStatisticsStatus()
{
  return _statistics.getStatisticsStatus();
}

template<typename T>
bool const& SuperGeometry3D<T>::getStatisticsStatus() const
{
  return _statistics.getStatisticsStatus();
}

template<typename T>
void SuperGeometry3D<T>::updateStatistics(bool verbose)
{
  if (this->_communicationNeeded) {
    this->communicate(verbose);
    getStatisticsStatus()=true;
  }
  _statistics.update(verbose);
  for (unsigned iC=0; iC<_extendedBlockGeometries.size(); iC++) {
    _extendedBlockGeometries[iC].getStatistics().update(verbose);
  }
}

template<typename T>
int SuperGeometry3D<T>::clean(bool verbose)
{
  this->communicate();
  int counter=0;
  for (unsigned iC=0; iC<_extendedBlockGeometries.size(); iC++) {
    counter+=_extendedBlockGeometries[iC].clean(false);
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(counter, MPI_SUM);
#endif

  if (verbose) {
    clout << "cleaned "<< counter << " outer boundary voxel(s)" << std::endl;
  }
  if (counter>0) {
    _statistics.getStatisticsStatus() = true;
    this->_communicationNeeded = true;
  }
  return counter;
}

template<typename T>
int SuperGeometry3D<T>::clean(int material, bool verbose)
{
  this->communicate();
  int counter=0;
  for (unsigned iC=0; iC<_extendedBlockGeometries.size(); iC++) {
    counter+=_extendedBlockGeometries[iC].clean(material, false);
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(counter, MPI_SUM);
#endif

  if (verbose) {
    clout << "cleaned "<< counter << " outer boundary voxel(s)" << std::endl;
  }
  if (counter>0) {
    _statistics.getStatisticsStatus() = true;
    this->_communicationNeeded = true;
  }
  return counter;
}

template<typename T>
int SuperGeometry3D<T>::outerClean(bool verbose)
{
  this->communicate();
  int counter=0;
  for (unsigned iC=0; iC<_extendedBlockGeometries.size(); iC++) {
    counter+=_extendedBlockGeometries[iC].outerClean(false);
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(counter, MPI_SUM);
#endif

  if (verbose) {
    clout << "cleaned "<< counter << " outer fluid voxel(s)" << std::endl;
  }
  if (counter>0) {
    _statistics.getStatisticsStatus() = true;
    this->_communicationNeeded = true;
  }
  return counter;
}

template<typename T>
int SuperGeometry3D<T>::innerClean(bool verbose)
{
  this->communicate();
  int counter=0;
  for (unsigned iC=0; iC<_extendedBlockGeometries.size(); iC++) {
    counter+=_extendedBlockGeometries[iC].innerClean(false);
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().barrier();
  singleton::mpi().reduceAndBcast(counter, MPI_SUM);
#endif

  if (verbose) {
    clout << "cleaned "<< counter << " inner boundary voxel(s)" << std::endl;
  }
  if (counter>0) {
    _statistics.getStatisticsStatus() = true;
    this->_communicationNeeded = true;
  }
  return counter;
}

template<typename T>
int SuperGeometry3D<T>::innerClean(int bcType, bool verbose)
{
  this->communicate();
  int counter=0;
  for (unsigned iC=0; iC<_extendedBlockGeometries.size(); iC++) {
    counter+=_extendedBlockGeometries[iC].innerClean(bcType,false);
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(counter, MPI_SUM);
#endif

  if (verbose) {
    clout << "cleaned "<< counter << " inner boundary voxel(s) of Type " << bcType << std::endl;
  }
  if (counter>0) {
    _statistics.getStatisticsStatus() = true;
    this->_communicationNeeded = true;
  }
  return counter;
}

template<typename T>
bool SuperGeometry3D<T>::checkForErrors(bool verbose)
{
  this->communicate();
  bool error = false;
  for (unsigned iC=0; iC<_blockGeometries.size(); iC++) {
    if (_blockGeometries[iC].checkForErrors(false)) {
      error = true;
    }
  }
  if (verbose) {
    if (error) {
      this->clout << "error!" << std::endl;
    }
    else {
      this->clout << "the model is correct!" << std::endl;
    }
  }
  return error;
}


template<typename T>
void SuperGeometry3D<T>::rename(int fromM, int toM)
{

  this->communicate();
  for (unsigned iC=0; iC<_extendedBlockGeometries.size(); iC++) {
    _blockGeometries[iC].rename(fromM,toM);
  }
  _statistics.getStatisticsStatus() = true;
  this->_communicationNeeded = true;
}

template<typename T>
void SuperGeometry3D<T>::rename(int fromM, int toM, FunctorPtr<IndicatorF3D<T>>&& condition)
{

  this->communicate();

  for (unsigned iC=0; iC<_extendedBlockGeometries.size(); iC++) {
    _extendedBlockGeometries[iC].rename(fromM,toM,*condition);
  }
  _statistics.getStatisticsStatus() = true;
}

template<typename T>
void SuperGeometry3D<T>::rename(int fromM, int toM, unsigned offsetX, unsigned offsetY, unsigned offsetZ)
{

  if ( offsetX<=unsigned(this->_overlap)
       && offsetY<=unsigned(this->_overlap)
       && offsetZ<=unsigned(this->_overlap) ) {
    this->communicate();
    for (unsigned iC=0; iC<_blockGeometries.size(); iC++) {
      _blockGeometries[iC].rename(fromM,toM,offsetX, offsetY, offsetZ);
    }
    _statistics.getStatisticsStatus() = true;
    this->_communicationNeeded = true;
  }
  else {
    clout << "error rename only implemented for offset<=overlap" << std::endl;
  }
}

template<typename T>
void SuperGeometry3D<T>::rename(int fromM, int toM, int testM, std::vector<int> testDirection)
{

  if ( testDirection[0]*testDirection[0]<=(this->_overlap)*(this->_overlap)
       && testDirection[1]*testDirection[1]<=(this->_overlap)*(this->_overlap)
       && testDirection[2]*testDirection[2]<=(this->_overlap)*(this->_overlap) ) {
    this->communicate();
    for (unsigned iC=0; iC<_blockGeometries.size(); iC++) {
      _blockGeometries[iC].rename(fromM,toM,testM,testDirection);
    }
    _statistics.getStatisticsStatus() = true;
    this->_communicationNeeded = true;
  }
  else {
    clout << "error rename only implemented for |testDirection[i]|<=overlap" << std::endl;
  }
}

template<typename T>
void SuperGeometry3D<T>::rename(int fromBcMat, int toBcMat, int fluidMat,
                                FunctorPtr<IndicatorF3D<T>>&& condition)
{
  if (this->_overlap>1) {
    this->communicate();
    rename(fromBcMat, toBcMat, *condition);
    std::vector<int> testDirection = this->getStatistics().computeDiscreteNormal(toBcMat);

    //std::cout << testDirection[0]<<testDirection[1]<<testDirection[2]<<std::endl;
    this->communicate();
    for (unsigned iC=0; iC<_blockGeometries.size(); iC++) {
      _blockGeometries[iC].rename(fromBcMat,toBcMat,fluidMat,*condition,testDirection);
    }
    _statistics.getStatisticsStatus() = true;
    this->_communicationNeeded = true;
  }
  else {
    clout << "error rename only implemented for overlap>=2" << std::endl;
  }
}

template<typename T>
void SuperGeometry3D<T>::copyMaterialLayer(IndicatorF3D<T>& condition, int discreteNormal[3], int numberOfLayers)
{

  this->communicate();

  for (unsigned iC=0; iC<_extendedBlockGeometries.size(); iC++) {
    _extendedBlockGeometries[iC].copyMaterialLayer(condition, discreteNormal, numberOfLayers);
  }
  _statistics.getStatisticsStatus() = true;

}


template<typename T>
void SuperGeometry3D<T>::print()
{
  this->_cuboidGeometry.print();
  getStatistics().print();
}

template<typename T>
std::unique_ptr<SuperIndicatorF3D<T>> SuperGeometry3D<T>::getMaterialIndicator(
                                     std::vector<int>&& materials)
{
  return this->getIndicator<SuperIndicatorMaterial3D>(
           std::forward<std::vector<int>>(materials));
}

template<typename T>
std::unique_ptr<SuperIndicatorF3D<T>> SuperGeometry3D<T>::getMaterialIndicator(
                                     int material)
{
  return this->getMaterialIndicator(std::vector<int> { material });
}


} // namespace olb

#endif
