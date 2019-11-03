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
 * A method to write vtk data for cuboid geometries
 * (only for uniform grids) -- generic implementation.
 */

#ifndef BLOCK_VTK_WRITER_3D_HH
#define BLOCK_VTK_WRITER_3D_HH

#include <fstream>
#include <iostream>
#include "communication/mpiManager.h"
#include "core/singleton.h"
#include "io/blockVtkWriter3D.h"
#include "io/base64.h"
#include "io/fileName.h"

namespace olb {

template<typename T>
BlockVTKwriter3D<T>::BlockVTKwriter3D( std::string name, bool binary )
  : clout( std::cout,"BlockVTKwriter3D" ), _name(name), _binary(binary)
{}

template<typename T>
BlockVTKwriter3D<T>::~BlockVTKwriter3D ()
{
  clearAddedFunctors();
}

template<typename T>
void BlockVTKwriter3D<T>::write(int iT)
{
  if ( _pointerVec.empty() ) {
    clout << "Error: Please add functor via addFunctor()";
  } else {
    // get first functor
    auto it = _pointerVec.cbegin();

    T originX = 0;
    T originY = 0;
    T originZ = 0;
    int nx = (**it).getBlockStructure().getNx() -1;
    int ny = (**it).getBlockStructure().getNy() -1;
    int nz = (**it).getBlockStructure().getNz() -1;

    std::string fullNameVti = singleton::directories().getVtkOutDir()
                              + createFileName( _name, iT ) + ".vti";

    preamble( fullNameVti, nx,ny,nz, originX,originY,originZ );
    if ( _binary ) {
      // iterate on functors
      for ( auto functor = _pointerVec.cbegin(); functor != _pointerVec.cend(); ++functor) {
        writeRawDataBinary( fullNameVti, **functor, nx, ny, nz);
      }
    } else {
      for ( auto functor = _pointerVec.cbegin(); functor != _pointerVec.cend(); ++functor) {
        writeRawData( fullNameVti, **functor, nx, ny, nz);
      }
    }
    closePreamble( fullNameVti );
  }
}

template<typename T>
void BlockVTKwriter3D<T>::write(BlockF3D<T>& f, int iT)
{
  T originX = 0;
  T originY = 0;
  T originZ = 0;
  int nx = f.getBlockStructure().getNx() -1;
  int ny = f.getBlockStructure().getNy() -1;
  int nz = f.getBlockStructure().getNz() -1;

  std::string fullNameVti = singleton::directories().getVtkOutDir()
                            + createFileName( f.getName(), iT ) + ".vti";

  preamble( fullNameVti, nx,ny,nz, originX,originY,originZ );
  if ( _binary ) {
    writeRawData( fullNameVti, f, nx,ny,nz );
  } else {
    writeRawDataBinary( fullNameVti, f, nx,ny,nz );
  }
  closePreamble( fullNameVti );
}

template<typename T>
void BlockVTKwriter3D<T>::addFunctor(BlockF3D<T>& f)
{
  _pointerVec.push_back(&f);
}

template<typename T>
void BlockVTKwriter3D<T>::clearAddedFunctors()
{
  _pointerVec.clear();
}

template<typename T>
void BlockVTKwriter3D<T>::preamble(const std::string& fullName, int nx, int ny, int nz,
                                   T originX, T originY, T originZ)
{
  if (singleton::mpi().getRank()==0) {
    std::ofstream fout(fullName.c_str());
    if (!fout) {
      clout << "Error: could not open " << fullName << std::endl;
    }
    // spacing is not known for BlockF3D classes
    // prone to error: spacing might correspond to extension in y or z direction
    //                 at the end of the day, is can be fixed by apply a scaling in paraview
    double spacing = 1/double(nx);

    fout << "<?xml version=\"1.0\"?>\n";
    fout << "<VTKFile type=\"ImageData\" version=\"0.1\" "
         << "byte_order=\"LittleEndian\">\n";
    fout << "<ImageData WholeExtent=\""
         << "0" <<" "<< nx <<" "
         << "0" <<" "<< ny <<" "
         << "0" <<" "<< nz
         << "\" Origin=\"" << originX << " " << originY << " " << originZ
         << "\" Spacing=\"" << spacing << " " << spacing << " " << spacing << "\">\n";

    fout << "<Piece Extent=\""
         << 0 <<" "<< nx <<" "
         << 0 <<" "<< ny <<" "
         << 0 <<" "<< nz <<"\">\n";

    fout << "<PointData>\n";
    fout.close();
  }
}

template<typename T>
void BlockVTKwriter3D<T>::closePreamble(const std::string& fullNamePiece)
{
  if (singleton::mpi().getRank()==0) {
    std::ofstream fout(fullNamePiece.c_str(), std::ios::app );
    if (!fout) {
      clout << "Error: could not open " << fullNamePiece << std::endl;
    }
    fout << "</PointData>\n";
    fout << "</Piece>\n";
    fout << "</ImageData>\n";
    fout << "</VTKFile>\n";
    fout.close();
  }
}



template<typename T>
void BlockVTKwriter3D<T>::writeRawData(const std::string& fullNameVti, BlockF3D<T>& f,
                                       int nx, int ny, int nz)
{
  std::ofstream fout(fullNameVti.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNameVti << std::endl;
  }

  if (singleton::mpi().getRank()==0) {
    fout << "<DataArray " ;
    fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
         << "NumberOfComponents=\"" << f.getTargetDim() << "\">\n";
  }

  int i[3] = {int()};
  T evaluated[f.getTargetDim()];
  for (int iDim = 0; iDim < f.getTargetDim(); ++iDim) {
    evaluated[iDim] = T();
  }
  for (i[2] = 0; i[2] < nz+1; ++i[2]) {
    for (i[1] = 0; i[1] < ny+1; ++i[1]) {
      for (i[0] = 0; i[0] < nx+1; ++i[0]) {
        f(evaluated,i);
        for (int iDim = 0; iDim < f.getTargetDim(); ++iDim) {
          if (singleton::mpi().getRank()==0) {
            fout << evaluated[iDim] << " ";
          }
        }
      }
    }
  }
  if (singleton::mpi().getRank()==0) {
    fout << "\n</DataArray>\n";
  }

  fout.close();
}

//  uses base64 encoder to write binary output
//  first number is written by a seperate sizeEncoder
//  this number indicates how many numbers will be stored.
//  then dataEncoder will be called to write output.
//  !!code is fixed to float functor values!!
template<typename T>
void BlockVTKwriter3D<T>::writeRawDataBinary(const std::string& fullNameVti,
    BlockF3D<T>& f, int nx, int ny, int nz)
{
  const char* fileName = fullNameVti.c_str();
  std::ofstream fout(fileName, std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fileName << std::endl;
  }

  if (singleton::mpi().getRank()==0) {
    fout << "<DataArray " ;
    if (f.getTargetDim() == 1) {
      fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
           << "format=\"binary\" encoding=\"base64\">\n";
    } else {
      fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
           << "format=\"binary\" encoding=\"base64\" "
           << "NumberOfComponents=\"" << f.getTargetDim() << "\">\n";
    }
  }
  fout.close();

  std::ofstream ofstr( fileName, std::ios::out | std::ios::app | std::ios::binary );
  if (!ofstr) {
    clout << "Error: could not open " << fileName << std::endl;
  }

  size_t fullSize = f.getTargetDim() * (1 + nx) * (1 + ny) * (1 + nz);
  size_t binarySize = size_t( fullSize * sizeof(float) );
  // writes first number, which have to be the size(byte) of the following data
  Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
  unsigned int uintBinarySize = (unsigned int)binarySize;
  sizeEncoder.encode(&uintBinarySize, 1);
  //  write numbers from functor
  Base64Encoder<float>* dataEncoder = nullptr;
  dataEncoder = new Base64Encoder<float>( ofstr, fullSize );

  int i[3] = {int()};
  T evaluated[f.getTargetDim()];
  for (int iDim = 0; iDim < f.getTargetDim(); ++iDim) {
    evaluated[iDim] = T();
  }
  for (i[2] = 0; i[2] < nz+1; ++i[2]) {
    for (i[1] = 0; i[1] < ny+1; ++i[1]) {
      for (i[0] = 0; i[0] < nx+1; ++i[0]) {
        f(evaluated,i);
        for (int iDim = 0; iDim<f.getTargetDim(); ++iDim) {
          if (singleton::mpi().getRank()==0) {
            const float evaluated2 = float( evaluated[iDim] );
            dataEncoder->encode( &evaluated2, 1 );
          }
        }
      }
    }
  }
  ofstr.close();

  if (singleton::mpi().getRank()==0) {
    std::ofstream foutt(fileName,  std::ios::out | std::ios::app);
    if (!foutt) {
      clout << "Error: could not open " << fileName << std::endl;
    }
    foutt << "\n</DataArray>\n";
    foutt.close();
  }
  delete dataEncoder;
}


}  // namespace olb

#endif
