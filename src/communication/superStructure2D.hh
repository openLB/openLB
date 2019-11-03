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
 * Mother class of SuperGeometry2D and SuperLattice2D -- generic implementation.
 */

#ifndef SUPER_STRUCTURE_2D_HH
#define SUPER_STRUCTURE_2D_HH

#include "geometry/cuboidGeometry2D.h"
#include "communication/loadBalancer.h"
#include "communication/superStructure2D.h"
#include "io/ostreamManager.h"

namespace olb {

template<typename T>
SuperStructure2D<T>::SuperStructure2D(CuboidGeometry2D<T>& cuboidGeometry,
                                      LoadBalancer<T>& loadBalancer, int overlap)
  : _cuboidGeometry(cuboidGeometry), _loadBalancer(loadBalancer), _overlap(overlap),
    _communicator(*this), _communicationNeeded(true), clout(std::cout,"SuperGeometry2D")
{
}

template<typename T>
SuperStructure2D<T>::SuperStructure2D(int overlap)
  : SuperStructure2D( *(new CuboidGeometry2D<T> ()), *(new LoadBalancer<T> ()), overlap)
{
}

template<typename T>
CuboidGeometry2D<T>& SuperStructure2D<T>::getCuboidGeometry()
{
  return _cuboidGeometry;
}

template<typename T>
CuboidGeometry2D<T> const& SuperStructure2D<T>::getCuboidGeometry() const
{
  return _cuboidGeometry;
}

template<typename T>
int SuperStructure2D<T>::getOverlap()
{
  return _overlap;
}

template<typename T>
int SuperStructure2D<T>::getOverlap() const
{
  return _overlap;
}

template<typename T>
LoadBalancer<T>& SuperStructure2D<T>::getLoadBalancer()
{
  return _loadBalancer;
}

template<typename T>
LoadBalancer<T> const& SuperStructure2D<T>::getLoadBalancer() const
{
  return _loadBalancer;
}

template<typename T>
void SuperStructure2D<T>::communicate(bool verbose)
{
  if (_communicationNeeded) {
    if (verbose) {
      std::cout << "Communicate ..." << std::endl;
    }
    _communicator.send();
    _communicator.receive();
    _communicator.write();
    _communicationNeeded = false;
    if (verbose) {
      std::cout << "Communicate ...ok" << std::endl;
    }
  }
}

template<typename T>
void SuperStructure2D<T>::initCommunication()
{
  _communicator.init_nh();
  _communicator.init();
}

} // namespace olb

#endif
