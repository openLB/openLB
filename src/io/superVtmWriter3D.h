/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016-2017 Albert Mink, Maximilian Gaedtke, Markus Morhard Mathias J. Krause
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
 * A method to write vtk data for cuboid geometries
 * (only for uniform grids) -- header file.
 */

#ifndef SUPER_VTM_WRITER_3D_H
#define SUPER_VTM_WRITER_3D_H

#include <sstream>
#include <vector>
#include "io/ostreamManager.h"
#include "functors/lattice/superBaseF3D.h"

namespace olb {

/** SuperVTMwriter3D writes any SuperF3D to vtk-based output files.
 *
 * In .pvd files, there are only links/references to a VTKmultiblock file 'vtm'
 *
 * .pvd file structur
 * the time series is represented by different 'vtm' files.
 *
 * .vtm file
 * This file links cuboids ('vti') and represents the entire data of a single timestep.
 *
 */
template<typename T, typename W=T>
class SuperVTMwriter3D {
public:
  SuperVTMwriter3D( const std::string& name, bool binary=true, bool compress=true );
  ///  writes functors stored in pointerVec
  ///  every thread writes a vti file with data from his cuboids
  ///  the vti files are linked in a pvd file
  void write(int iT=0);
  ///  writes functor instantaneously, same vti-pvd file structure as above
  void write(SuperF3D<T,W>& f, int iT=0);
  void write(std::shared_ptr<SuperF3D<T,W>> ptr_f, int iT=0);
  ///  have to be called before calling write(int iT=0), since it creates
  //   the master pvd file, where all vti are linked!
  void createMasterFile();
  ///  put functor to _pointerVec
  ///  to simplify writing process of several functors
  void addFunctor(SuperF3D<T,W>& f);
  ///  to clear stored functors, not yet used due to lack of necessity
  void clearAddedFunctors();
  /// getter for _name
  std::string getName() const;
private:
  ///  performes <VTKFile ...>, <ImageData ...>, <PieceExtent ...> and <PointData ...>
  void preambleVTI(const std::string& fullName, const Vector<int,3> extent0, const Vector<int,3> extent1,
                   T origin[], T delta);
  ///  performes </ImageData> and </VTKFile>
  void closeVTI(const std::string& fullNamePiece);
  ///  performes <VTKFile ...> and <Collection>
  void preamblePVD(const std::string& fullNamePVD);
  ///  performes </Collection> and </VTKFile>
  void closePVD(const std::string& fullNamePVD);
  ///  performes <VTKFile ...> and <Collection>
  void preambleVTM(const std::string& fullNamePVD);
  ///  performes </Collection> and </VTKFile>
  void closeVTM(const std::string& fullNamePVD);
  ///  performes <DataSet timestep= ... file=namePiece />
  ///  used for linking vti into pvd files
  void dataVTM(int iC, const std::string& fullNamePVD, const std::string& namePiece);
  ///  performes <DataSet timestep= ... file=namePiece />
  ///  *** nasty function ***
  void dataPVDmaster(int iT, const std::string& fullNamePVDMaster,
                     const std::string& namePiece);
  ///  writes given functor f, ascii or base64 or zLib
  void dataArray(const std::string& fullName, SuperF3D<T,W>& f,
                 int iC, const Vector<int,3> extent1);
  ///  performes </PointData> and </Piece>
  void closePiece(const std::string& fullNamePiece);
private:
  mutable OstreamManager clout;
  ///  default is false, call createMasterFile() and it will be true
  bool _createFile;
  ///  determines the name of .vti and .pvd per iT
  std::string const _name;
  ///  holds added functor, to simplify the use of write function
  std::vector< SuperF3D<T,W>* > _pointerVec;
  ///  writing data base64 encoded
  bool _binary;
  ///  writing data zLib compressed
  bool _compress;
};


}  // namespace olb


#endif
