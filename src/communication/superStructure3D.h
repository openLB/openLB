/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Peter Weisbrod, Albert Mink, Mathias J. Krause
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
 * Mother class of SuperGeometry3D and SuperLattice3D -- header file.
 */

#ifndef SUPER_STRUCTURE_3D_H
#define SUPER_STRUCTURE_3D_H

#include "geometry/cuboidGeometry3D.h"
#include "communication/loadBalancer.h"
#include "communication/communicator3D.h"

/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T> class LoadBalancer;
template<typename T> class CuboidGeometry3D;
template<typename T> class Communicator3D;

template<typename T>
class SuperStructure3D {
protected:
  /// The grid structure is stored here
  CuboidGeometry3D<T>& _cuboidGeometry;
  /// Distribution of the cuboids of the cuboid structure
  LoadBalancer<T>& _loadBalancer;
  /// Size of ghost cell layer (must be greater than 1 and
  /// greater_overlapBC, default =1)
  int _overlap;
  /// This communicator handles the communication of the overlap
  Communicator3D<T> _communicator;
  /// Specifies if there has been some data updated which requires
  /// communication
  bool _communicationNeeded;
  /// class specific output stream
  mutable OstreamManager clout;
public:
  /// Virtual Destructor for inheritance
  virtual ~SuperStructure3D() {};
  /// Construction of a super structure
  SuperStructure3D(CuboidGeometry3D<T>& cuboidGeometry,
                   LoadBalancer<T>& loadBalancer, int overlap = 2);
  /// Default Constructor for empty SuperStructure
  SuperStructure3D(int overlap = 1);

  /// Write access to the memory of the data of the super structure where (iX, iY, iZ) is the point providing the data iData in the block iCloc
  virtual bool* operator() (int iCloc, int iX, int iY, int iZ, int iData) =0;
  /// Read only access to the dim of the data of the super structure
  virtual int getDataSize() const =0;
  /// Read only access to the data type dim of the data of the super structure
  virtual int getDataTypeSize() const=0;

  /// Read and write access to cuboid geometry
  CuboidGeometry3D<T>& getCuboidGeometry();
  /// Read only access to cuboid geometry
  CuboidGeometry3D<T> const& getCuboidGeometry() const;
  /// Read and write access to the overlap
  int getOverlap();
  /// Read only access to the overlap
  int getOverlap() const;
  /// Read and write access to the load balancer
  LoadBalancer<T>& getLoadBalancer();
  /// Read only access to the load balancer
  LoadBalancer<T> const& getLoadBalancer() const;
  /// Communicates the data in the overlap
  void communicate(bool verbose=false);

  void initCommunication();
};

} // namespace olb

#endif
