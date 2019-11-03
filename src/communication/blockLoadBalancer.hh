/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Mathias Krause
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

#ifndef BLOCK_LOAD_BALANCER_HH
#define BLOCK_LOAD_BALANCER_HH

#include <vector>
#include <map>
#include "communication/mpiManager.h"
#include "communication/blockLoadBalancer.h"
#include "geometry/cuboidGeometry2D.h"
#include "geometry/cuboidGeometry3D.h"
#include "core/olbDebug.h"

namespace olb {

template<typename T>
BlockLoadBalancer<T>::BlockLoadBalancer(int rank, int size, int globChunkSize, int offset)
{
  init_chunkD(rank, size, globChunkSize, offset);
}

template<typename T>
BlockLoadBalancer<T>::BlockLoadBalancer(CuboidGeometry3D<T>& cGeometry)
{
  init_chunkD(singleton::mpi().getRank(), singleton::mpi().getSize(), cGeometry.getNc(), 0);
}

template<typename T>
BlockLoadBalancer<T>::BlockLoadBalancer(CuboidGeometry2D<T>& cGeometry)
{
  init_chunkD(singleton::mpi().getRank(), singleton::mpi().getSize(), cGeometry.getNc(), 0);
}

template<typename T>
void BlockLoadBalancer<T>::init_chunkD(int rank, int size, int globChunkSize, int offset)
{

  OLB_PRECONDITION(rank>=0 && size>=1 && offset>=0)
  OLB_PRECONDITION(size<=globChunkSize &&  rank<size);

  // nice way to calculate # of chunks per processor
  this->_locChunkSize = (globChunkSize+size-rank-1)/size;
  this->_size = _locChunkSize;
  if (rank+1 <= globChunkSize-(globChunkSize/size)*size) {
    this->_firstGlobNum = globChunkSize/size * rank + rank + offset;
    this->_lastGlobNum  = this->_firstGlobNum + this->_locChunkSize - 1;
  } else {
    this->_firstGlobNum = globChunkSize/size * rank + globChunkSize - (globChunkSize/size)*size + offset;
    this->_lastGlobNum  = this->_firstGlobNum + this->_locChunkSize - 1;
  }
  for (int i=0; i<this->_locChunkSize; i++) {
    this->_loc[this->_firstGlobNum + i] = i;
    this->_glob.push_back(this->_firstGlobNum + i);
  }
  int temp = offset;
  for (int iRank=0; iRank<size; iRank++) {
    int iLocChunkSize = (globChunkSize+size-iRank-1)/size;
    for (int i=0; i<iLocChunkSize; i++) {
      this->_rank[temp] = iRank;
      temp++;
    }
  }
}

template<typename T>
int BlockLoadBalancer<T>::locChunkSize() const
{
  return _locChunkSize;
}

template<typename T>
int BlockLoadBalancer<T>::firstGlobNum() const
{
  return _firstGlobNum;
}

template<typename T>
int BlockLoadBalancer<T>::lastGlobNum() const
{
  return _lastGlobNum;
}

}  // namespace olb
#endif
