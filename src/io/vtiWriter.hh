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
 * (only for uniform grids) -- generic implementation.
 */

#ifndef VTI_WRITER_HH
#define VTI_WRITER_HH

#include <fstream>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "core/singleton.h"
#include "geometry/cuboid3D.h"
#include "geometry/cuboidGeometry3D.h"
#include "core/blockData3D.h"
#include "core/superData3D.h"
#include "communication/loadBalancer.h"
#include "vtiWriter.h"

namespace olb {

// Initialize Output
template<typename T, typename BaseType>
OstreamManager VTIwriter3D<T,BaseType>::clout("VTIwriter3D");

template<typename T, typename BaseType>
VTIwriter3D<T,BaseType>::VTIwriter3D() {};

template<typename T, typename BaseType>
void VTIwriter3D<T,BaseType>::writeData( std::string const& fName,
    std::string const& fieldName, BlockData3D<T,BaseType> const& blockData,
    Cuboid3D<T> const& cuboid)
{
  // Generate full file path
  std::string fullName = getFullName(fName);
  writePreamble(fullName, cuboid);
  writeBlockData(fullName, fieldName, blockData, cuboid);
  writePostScript(fullName);
}

template<typename T, typename BaseType>
void VTIwriter3D<T,BaseType>::writeData( std::string const& fName, std::string const& fieldName,
    SuperData3D<T,BaseType> const& superData, CuboidGeometry3D<T> const& cGeometry,
    LoadBalancer<T> const& loadBalancer)
{
  // Generate full file path
  std::string fullName = getFullName(fName);

  int rank = 0;
  int size = 1;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
  size = singleton::mpi().getSize();
#endif

  if ( rank == 0 ) {
    // Write XML Preamble with geometry of surrounding cuboid
    writePreamble(fullName, cGeometry.getMotherCuboid());
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().barrier();
#endif

  // Write VTI File in Order of MPI Rank
  for (int iRank = 0; iRank < size; ++iRank) {
    if (rank == iRank) {
      for (int iC = 0; iC < loadBalancer.size(); ++iC) {
        // Write data block by block
        writeBlockData(fullName, fieldName, superData.get(iC), cGeometry.get(loadBalancer.glob(iC)));
      }
    }
#ifdef PARALLEL_MODE_MPI
    // Wait for current writer process
    singleton::mpi().barrier();
#endif
  }
  //
  if (rank == 0) {
    // Close All XML Tags
    writePostScript(fullName);
  }
}

template<typename T, typename BaseType>
void VTIwriter3D<T,BaseType>::writeData( std::string const& fName,
    std::string const& fieldName, SuperData3D<T,BaseType> const& superData)
{
  writeData(fName, fieldName, superData, superData.getCuboidGeometry(), superData.getLoadBalancer());
}

template<typename T, typename BaseType>
void VTIwriter3D<T,BaseType>::writeBlockData(std::string& fullName, std::string const& fieldName,
    BlockData3D<T,BaseType> const& blockData, Cuboid3D<T> const& cuboid)
{
  // Open File
  std::ofstream fout(fullName.c_str(), std::ios::app);
  if (!fout) {
    std::cout << "Error: could not open file " << fullName << std::endl;
  }

  // Calculate extent info from cuboid
  // NOTE: origin is actual coordinate, globX needs to be number of starting node
  Vector<T,3> origin = cuboid.getOrigin();
  T delta = cuboid.getDeltaR();
  int globX = (int) (origin[0] / (T)delta);
  int globY = (int) (origin[1] / (T)delta);
  int globZ = (int) (origin[2] / (T)delta);
  fout << "<Piece Extent=\""
       << globX << " " << cuboid.getNx() - 1 + globX << " "
       << globY << " " << cuboid.getNy() - 1 + globY << " "
       << globZ << " " << cuboid.getNz() - 1 + globZ << "\">\n";

  // Write Point Data
  fout << "<PointData>\n";
  //  fout << "<DataArray type=\"Float32\" Name=\"" << blockData.getName()
  //       << "\" NumberOfComponents=\"" << blockData.getSize() << "\">\n";

  fout << "<DataArray type=\"Float32\" Name=\"" << fieldName
       << "\" NumberOfComponents=\"" << blockData.getSize() << "\">\n";

  // NOTE: Respect ordering of VTI File!
  for (int iZ = 0; iZ < cuboid.getNz(); ++iZ) {
    for (int iY = 0; iY < cuboid.getNy(); ++iY) {
      for (int iX = 0; iX < cuboid.getNx(); ++iX) {
        for (int iSize = 0; iSize < blockData.getSize(); ++iSize) {
          fout << blockData.get(iX, iY, iZ, iSize) << " ";
        }
      }
    }
    fout << "\n";
  }

  // Close All Tags
  fout << "</DataArray>\n";
  fout << "</PointData>\n";
  fout << "</Piece>\n";

  // Close File
  fout.close();
}

template<typename T, typename BaseType>
void VTIwriter3D<T,BaseType>::writePreamble(std::string& fullName, Cuboid3D<T> const& cuboid)
{
  Vector<T,3> origin = cuboid.getOrigin();
  writePreamble(fullName,
                cuboid.getNx() - 1, cuboid.getNy() - 1, cuboid.getNz() - 1,
                cuboid.getDeltaR(),
                origin[0], origin[1], origin[2]);
}

template<typename T, typename BaseType>
void VTIwriter3D<T,BaseType>::writePreamble(std::string& fullName, int nx,
    int ny, int nz, T delta, T originX, T originY, T originZ)
{
  // Open File
  std::ofstream fout(fullName.c_str());
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }

  // Follow VTI standard
  fout << "<?xml version=\"1.0\"?>\n";
  fout << "<VTKFile type=\"ImageData\" version=\"0.1\" "
       << "byte_order=\"LittleEndian\">\n";
  fout << "<ImageData";
  fout << " WholeExtent=\" 0 " << nx << " 0 " << ny << " 0 " << nz << " \"";
  fout << " Origin=\"" << originX << " " << originY << " " << originZ << "\"";
  fout << " Spacing=\"" << delta << " " << delta << " " << delta << "\">\n";

  // Close File
  fout.close();
}

template<typename T, typename BaseType>
void VTIwriter3D<T,BaseType>::writePostScript(std::string& fullName)
{
  std::ofstream fout(fullName.c_str(), std::ios::app );
  if (!fout) {
    clout << "Error: could not open " << fullName << std::endl;
  }

  fout << "</ImageData>\n";
  fout << "</VTKFile>\n";

  fout.close();
}

template<typename T, typename BaseType>
std::string VTIwriter3D<T,BaseType>::getFullName(std::string const& fName)
{
  return singleton::directories().getVtkOutDir() + fName + ".vti";
}


}  // namespace olb

#endif
