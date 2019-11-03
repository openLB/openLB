/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause
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
 * A method to write vtk data for block geometries
 * (only for uniform grids) -- header file.
 */

#ifndef BLOCK_VTK_WRITER_3D_H
#define BLOCK_VTK_WRITER_3D_H

#include "io/ostreamManager.h"
#include "functors/lattice/blockBaseF3D.h"

namespace olb {

/** BlockVTKwriter3D writes any BLockF3D to vtk-based output files
 *
 * One can add functors via addFunctors. To write added functors
 * call write(int iT=0).
 *
 * To write a functor without adding him, call
 * call write(BlockF3D<T>& f, int iT).
 *
 * Default output type is binary, to change it have a look on the
 * constructor.
 */
template<typename T>
class BlockVTKwriter3D {
public:
  BlockVTKwriter3D( std::string name, bool binary=true );
  ~BlockVTKwriter3D();
  ///  method calls preamble(), pointData(), data() and coresponding
  ///  closing methods.
  ///  writes given functor
  void write( BlockF3D<T>& f, int iT=0 );
  ///  writes functors stored at pointerVec
  void write( int iT=0 );
  ///  put functor to _pointerVec
  ///  to simplify writing process of several functors
  void addFunctor( BlockF3D<T>& f );
  ///  to clear stored functors
  void clearAddedFunctors();
private:
  ///  writes <VTKFile .... >, <ImageData ... >, <Piece ... > and  <PointData Scalar="..." >
  void preamble( const std::string& fullName, int nx, int ny, int nz,
                 T originX=0, T originY=0, T originZ=0);
  ///  writes </PointData>, </Piece>, </ImageData> and  </VTKFile>
  void closePreamble( const std::string& fullName );
  ///  writes data of given functor f as ASCII
  void writeRawData( const std::string& fullNameVti, BlockF3D<T>& f, int nx, int ny, int nz );
  ///  writes data of given functor f as binary
  void writeRawDataBinary( const std::string& fullNameVti, BlockF3D<T>& f,
                           int nx, int ny, int nz );
private:
  mutable OstreamManager clout;
  ///  determines the name of .vti per iT
  std::string _name;
  ///  default is true, may be changed at constructor
  bool _binary;
  ///  holds added functor, to simplify the use of write function
  std::vector< BlockF3D<T>* > _pointerVec;
};

}  // namespace olb


#endif
