/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias J. Krause
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
 * A communincator provides a cuboids with cells of other
 * cuboids -- header file.
 */


#ifndef COMMUNICATOR_3D_H
#define COMMUNICATOR_3D_H

#include <vector>
#include "communication/cuboidNeighbourhood3D.h"
#include "communication/superStructure3D.h"


/// All OpenLB code is contained in this namespace.
namespace olb {

/// A communicator provides a cuboid with data of other
/// cuboids.
/** For each cuboid a cuboid neighbourhood is defined. A communicator
 * coordinates the communication of the single cuboids, in which the data
 * is stored locally.
 *
 * There are methods to add needed cells, send, receive and write
 * data of the cell.
 *
 * This class is not intended to be derived from.
 */


template<typename T> class CuboidNeighbourhood3D;
template<typename T> class SuperStructure3D;

template<typename T>
class Communicator3D {

private:
  /// Stores neighbourhoods information for a single cuboid
  std::vector<CuboidNeighbourhood3D<T> > _nh;
  /// Reference to the super structure
  SuperStructure3D<T>& _superStructure;
  /// Number of cuboids in the structure
  int _nC;
  /// Indicates if there was an initialization done
  bool _initDone;
public:
  ~Communicator3D() {};
  /// Constructor
  Communicator3D(SuperStructure3D<T>& superStructure);
  /// Initializes the cuboid neighbourhoods
  void init_nh();

  /// Adds a cell to the cuboid neigbourhood iCloc
  /// that is then able to get data of the cell
  void add_cell(int iCloc, int iX, int iY, int iZ);
  /// Adds all cells of the overlap to its associated cuboid
  /// neigbourhood
  void add_cells(int overlap);
  /// Initializes the communicator
  void init();
  /// Resets all cuboids neighbourhoods
  void reset();

  /// Sends the needed data
  void send();
  /// Receives the needed data
  void receive();
  /// Writes all data to the block lattices of the super structure
  void write();
};

}  // namespace olb

#endif
