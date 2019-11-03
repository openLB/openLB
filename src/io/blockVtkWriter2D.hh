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

#ifndef BLOCK_VTK_WRITER_2D_HH
#define BLOCK_VTK_WRITER_2D_HH

#include <fstream>
#include <iostream>
#include "communication/mpiManager.h"
#include "core/singleton.h"
#include "io/blockVtkWriter2D.h"
#include "io/base64.h"
#include "io/fileName.h"

namespace olb {


template<typename T>
BlockVTKwriter2D<T>::BlockVTKwriter2D( std::string name, bool binary )
  : clout( std::cout,"BlockVTKwriter2D" ), _name(name), _binary(binary)
{}

template<typename T>
BlockVTKwriter2D<T>::~BlockVTKwriter2D ()
{
  clearAddedFunctors();
}

//  iteration on _pointerVec is realized by function
//  dataArray() respective dataArrayBinary()
template<typename T>
void BlockVTKwriter2D<T>::write(int iT)
{
  if ( _pointerVec.empty() ) {
    clout << "Error: Please add functor via addFunctor()";
  } else {
    auto it = _pointerVec.cbegin();

    T originX = 0;
    T originY = 0;
    T originZ = 0;
    int nx = (**it).getBlockStructure().getNx() -1;
    int ny = (**it).getBlockStructure().getNy() -1;

    std::string fullNameVti = singleton::directories().getVtkOutDir()
                              + createFileName( _name, iT ) + ".vti";


    preamble( fullNameVti, nx,ny, originX,originY,originZ );
    if ( _binary ) {
      for ( auto functor = _pointerVec.cbegin(); functor != _pointerVec.cend(); ++functor) {
        writeRawDataBinary( fullNameVti, **functor, nx, ny);
      }
    } else {
      for ( auto functor = _pointerVec.cbegin(); functor != _pointerVec.cend(); ++functor) {
        writeRawData( fullNameVti, **functor, nx, ny);
      }
    }
    closePreamble( fullNameVti );
  }
}

template<typename T>
void BlockVTKwriter2D<T>::write(BlockF2D<T>& f, int iT)
{
  T originX = 0;
  T originY = 0;
  T originZ = 0;
  int nx = f.getBlockStructure().getNx() -1;
  int ny = f.getBlockStructure().getNy() -1;

  std::string fullNameVti = singleton::directories().getVtkOutDir()
                            + createFileName( f.getName(), iT ) + ".vti";

  preamble( fullNameVti, nx,ny, originX,originY,originZ );
  if ( _binary ) {
    writeRawDataBinary( fullNameVti, f, nx,ny );
  } else {
    writeRawData( fullNameVti, f, nx,ny );
  }
  closePreamble( fullNameVti );
}

template<typename T>
void BlockVTKwriter2D<T>::addFunctor(BlockF2D<T>& f)
{
  _pointerVec.push_back(&f);
}

template<typename T>
void BlockVTKwriter2D<T>::clearAddedFunctors()
{
  _pointerVec.clear();
}

template<typename T>
void BlockVTKwriter2D<T>::preamble(const std::string& fullName, int nx, int ny,
                                   T originX, T originY, T originZ)
{
  if (singleton::mpi().getRank()==0) {
    std::ofstream fout(fullName.c_str());
    if (!fout) {
      clout << "Error: could not open " << fullName << std::endl;
    }

    // spacing is not known for BlockF2D classes
    // prone to error: spacing might correspond to extension in y direction
    //                 at the end of the day, is can be fixed by apply a scaling in paraview
    double spacing = double(1.0/nx);

    fout << "<?xml version=\"1.0\"?>\n";
    fout << "<VTKFile type=\"ImageData\" version=\"0.1\" "
         << "byte_order=\"LittleEndian\">\n";
    fout << "<ImageData WholeExtent=\""
         << 0 <<" "<< nx <<" "
         << 0 <<" "<< ny <<" "
         << 0 <<" "<< 0
         << "\" Origin=\"" << originX << " " << originY << " " << originZ
         << "\" Spacing=\"" << spacing << " " << spacing << " " << spacing << "\">\n";

    fout << "<Piece Extent=\""
         << 0 <<" "<< nx <<" "
         << 0 <<" "<< ny <<" "
         << 0 <<" "<< 0 <<"\">\n";

    fout << "<PointData>\n";
    fout.close();
  }
}

template<typename T>
void BlockVTKwriter2D<T>::closePreamble(const std::string& fullNamePiece)
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
void BlockVTKwriter2D<T>::writeRawData(const std::string& fullNameVti, BlockF2D<T>& f,
                                       int nx, int ny)
{
  std::ofstream fout(fullNameVti.c_str(), std::ios::app);
  if (!fout) {
    clout << "Error: could not open " << fullNameVti << std::endl;
  }

  if (singleton::mpi().getRank()==0) {
    fout << "<DataArray " ;
    fout << "type=\"Float32\" Name=\"" << f.getName() << "\" "
         << "NumberOfComponents=\"" << f.getTargetDim() <<"\">\n";
  }

  int i[2] = {int()};
  T evaluated[f.getTargetDim()];
  for (int iDim = 0; iDim < f.getTargetDim(); ++iDim) {
    evaluated[iDim] = T();
  }
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

  if (singleton::mpi().getRank()==0) {
    fout << "\n</DataArray>\n";
  }

  fout.close();
}

template<typename T>
void BlockVTKwriter2D<T>::writeRawDataBinary(const std::string& fullNameVti,
    BlockF2D<T>& f, int nx, int ny)
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
           << "NumberOfComponents=\"" << f.getTargetDim() <<"\">\n";
    }
  }
  fout.close();

  std::ofstream ofstr( fileName, std::ios::out | std::ios::app | std::ios::binary );
  if (!ofstr) {
    clout << "Error: could not open " << fileName << std::endl;
  }

  size_t fullSize = f.getTargetDim() * (1 + nx) * (1 + ny) * (1);
  size_t binarySize = size_t( fullSize * sizeof(float) );
  // writes first number, which have to be the size(byte) of the following data
  Base64Encoder<unsigned int> sizeEncoder(ofstr, 1);
  unsigned int uintBinarySize = (unsigned int)binarySize;
  sizeEncoder.encode(&uintBinarySize, 1);
  //  write numbers from functor
  Base64Encoder<float>* dataEncoder = nullptr;
  dataEncoder = new Base64Encoder<float>( ofstr, fullSize );


  int i[2] = {int()};
  T evaluated[f.getTargetDim()];
  for (int iDim = 0; iDim < f.getTargetDim(); ++iDim) {
    evaluated[iDim] = T();
  }
  for (i[1] = 0; i[1] < ny+1; ++i[1]) {
    for (i[0] = 0; i[0] < nx+1; ++i[0]) {
      f(evaluated,i);
      for (int iDim = 0; iDim < f.getTargetDim(); ++iDim) {
        if (singleton::mpi().getRank()==0) {
          const float evaluated2 = float( evaluated[iDim] );
          dataEncoder->encode( &evaluated2, 1 );
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
