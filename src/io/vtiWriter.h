/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007, 2009, 2012, 2015 Mathias J. Krause, Benjamin Förster, Jonas Latt, Tim Dornieden
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
 * A method to write vti data for cuboid geometries
 * (only for uniform grids) -- header file.
 */

#ifndef VTI_WRITER_H
#define VTI_WRITER_H


#include "io/ostreamManager.h"


namespace olb {

template< typename T, typename BaseType> class BlockData3D;
template< typename T, typename BaseType> class SuperData3D;
template< typename T> class LoadBalancer;
template< typename T> class Cuboid3D;
template< typename T> class CuboidGeometry3D;

template<typename T, typename BaseType>
class VTIwriter3D {
public:
  VTIwriter3D();
  /// Write Single Block Data
  static void writeData( std::string const& fName, std::string const& fieldName,
                         BlockData3D<T,BaseType> const& blockData, Cuboid3D<T> const& cuboid);
  /// Write Super Data
  static void writeData( std::string const& fName, std::string const& fieldName,
                         SuperData3D<T,BaseType> const& superData, CuboidGeometry3D<T> const& cGeometry,
                         LoadBalancer<T> const& loadBalancer);
  /// Write Super Data with its own cGeometry and loadBalancer
  static void writeData( std::string const& fName, std::string const& fieldName,
                         SuperData3D<T,BaseType> const& superData);
private:
  static OstreamManager clout;
  /// Write VTK Preamble and PostScript
  static void writePreamble(std::string& fullName, Cuboid3D<T> const& cuboid);
  static void writePreamble(std::string& fullName, int nx, int ny, int nz,
                            T delta, T originX, T originY, T originZ);
  static void writePostScript(std::string& fullName);
  /// Write BlockData3D - Used for single BlockData as well as SuperData
  static void writeBlockData(std::string& fullName, std::string const& fieldName, BlockData3D<T,BaseType> const& blockData,
                             Cuboid3D<T> const& cuboid);
  /// Helper Functions
  static std::string getFullName(std::string const& fName);
};

}  // namespace olb


#endif
