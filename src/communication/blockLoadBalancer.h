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


#ifndef BLOCK_LOAD_BALANCER_H
#define BLOCK_LOAD_BALANCER_H


#include "communication/loadBalancer.h"



namespace olb {

template<typename T> class CuboidGeometry3D;
template<typename T> class CuboidGeometry2D;


template<typename T>
class BlockLoadBalancer : public LoadBalancer<T> {
private:
  int _locChunkSize;
  int _firstGlobNum;
  int _lastGlobNum;
public:
  BlockLoadBalancer() {}
  BlockLoadBalancer(int rank, int size, int globChunkSize, int offset);
  BlockLoadBalancer(CuboidGeometry2D<T>& cGeometry);
  BlockLoadBalancer(CuboidGeometry3D<T>& cGeometry);
  void init_chunkD(int rank, int size, int globChunkSize, int offset);
  int locChunkSize() const;
  int firstGlobNum() const;
  int lastGlobNum() const;
};
}  // namespace olb

#endif
