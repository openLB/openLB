/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause
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
 * Dynamics for a generic 3D block -- header file.
 */
#ifndef BLOCK_STRUCTURE_3D_H
#define BLOCK_STRUCTURE_3D_H



namespace olb {

/** An empty hull with left bottom corner at (0,0,0).
 *
 * \param _nx extension in x direction
 * \param _ny extension in y direction
 * \param _nz extension in z direction
 *
 */
class BlockStructure3D {
protected:
  /// Block width
  int _nx;
  /// Block height
  int _ny;
  /// Block lenght
  int _nz;
public:
  BlockStructure3D(int nx, int ny, int nz) : _nx(nx), _ny(ny), _nz(nz) {};
  //template <typename T>
  //BlockStructure3D(Cuboid3D<T> cuboid) : _nx(cuboid.getNx()), _ny(cuboid.getNy()), _nz(cuboid.getNz()) {};
  /// Read only access to block width
  int getNx() const
  {
    return _nx;
  };
  /// Read only access to block height
  int getNy() const
  {
    return _ny;
  };
  /// Read only access to block height
  int getNz() const
  {
    return _nz;
  };
};

}  // namespace olb

#endif
