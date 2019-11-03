/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Albert Mink, Mathias Krause
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

#ifndef BLOCK_GIF_WRITER_HH
#define BLOCK_GIF_WRITER_HH

#include <fstream>
#include <string>
#include <vector>

#include "io/blockGifWriter.h"
#include "io/fileName.h"
#include "utilities/vectorHelpers.h"
#include "core/singleton.h"
#include "communication/mpiManager.h"

namespace olb {


template< typename T >
BlockGifWriter<T>::BlockGifWriter(std::string const& map)
  : clout(std::cout, "BlockGifWriter"), _colorRange(1024), _numColors(1024),
    _colorMap(graphics::mapGenerators::generateMap<T>(map)), _minValue(0),
    _maxValue(1)
{ }

template< typename T >
void BlockGifWriter<T>::write(BlockF2D<T>& f, T minValue, T maxValue, int iT,
                              std::string const& name)
{
  // [!] exeption image(f) != 1
  if ( f.getTargetDim() != 1 ) {
    clout << "Error: Functor targetDim is not 1. " << std::endl;
    exit(-1);
  } else {
    if ( singleton::mpi().getRank() == 0 ) {
      std::string fullNamePpm;
      if ( name == "emptyName") {
        fullNamePpm = createFileName( singleton::directories().getImageOutDir(),
                                      f.getName(), iT);
      } else {
        fullNamePpm = createFileName( singleton::directories().getImageOutDir(),
                                      name, iT);
      }
      fullNamePpm = fullNamePpm  + ".ppm" ;
      std::ofstream fout( fullNamePpm.c_str() );

      // write header
      fout << "P3\n";
      // dimension of image
      fout << f.getBlockStructure().getNx() << " "
           << f.getBlockStructure().getNy() << "\n";
      // dynamic range
      fout << (_colorRange - 1) << "\n";

      int i[2] = {0,0};
      for (i[1] = f.getBlockStructure().getNy() - 1; i[1] >= 0; --i[1]) {
        for (i[0] = 0; i[0] < f.getBlockStructure().getNx(); ++i[0]) {
          T evaluated[1];
          f(evaluated,i);
          // scales evaluated in [getMinValue(),getMaxValue()] to [0,1]
          evaluated[0] = (minValue - evaluated[0]) / (minValue - maxValue);
          // sets evaluated notin [getMinValue(),getMaxValue()] to 1
          if ( evaluated[0] >= T(1) ) {
            evaluated[0] = 1;
          }
          graphics::rgb<T> color = _colorMap.get(evaluated[0]);
          fout << (int)(color.r * (_colorRange - 1)) << " "
               << (int)(color.g * (_colorRange - 1)) << " "
               << (int)(color.b * (_colorRange - 1)) << "\n";
        }
      }
      fout.close();
    }
  }
}

template< typename T >
void BlockGifWriter<T>::write(BlockF2D<T>& f, int iT, std::string const& name)
{
  // determine min-/maxValue
  int i[2] = {0,0};
  // initialize min-/maxValue
  T minValue[1];
  T maxValue[1];
  f(minValue,i);
  f(maxValue,i);
  for (i[0] = 1; i[0] < f.getBlockStructure().getNx(); i[0]++) {
    for (i[1] = 1; i[1] < f.getBlockStructure().getNy(); i[1]++) {
      T valueTmp[1];
      f(valueTmp,i);
      if (valueTmp[0] < minValue[0]) {
        minValue[0] = valueTmp[0];
      }
      if (valueTmp[0] > maxValue[0]) {
        maxValue[0] = valueTmp[0];
      }
    }
  }
  if (maxValue[0] <= minValue[0]) {
    minValue[0] = T();
    maxValue[0] = T(1);
  }
  // call write() with min-/maxValue
  write(f, minValue[0], maxValue[0], iT, name);
}

//  iteration on _pointerVec is realized by function
//  dataArray() respective dataArrayBinary()
template< typename T >
void BlockGifWriter<T>::write(int iT)
{
  if ( _pointerVec.empty() ) {
    // secure code. doesn't appear on console ??
    clout << "Error: Please add functor via addFunctor()";
  } else {
    int i = 0;
    for ( auto it = _pointerVec.cbegin(); it != _pointerVec.cend(); ++it, ++i) {
      if (_autoScale[i]) {
        write(**it, iT, _name[i]);
      } else {
        write(**it, _minValue[i], _maxValue[i], iT, _name[i]);
      }
    }
  }
}

template< typename T >
void BlockGifWriter<T>::addFunctor(BlockF2D<T>& f, std::string const& name)
{
  _pointerVec.push_back(&f);

  if ( name == "emptyName") {
    _name.push_back(f.getName() );
  }
  _name.push_back(name);

  _autoScale.push_back(true);
  _minValue.push_back(T());
  _maxValue.push_back(T());
}

template< typename T >
void BlockGifWriter<T>::addFunctor(BlockF2D<T>& f, T minValue, T maxValue,
                                   std::string const& name)
{
  _pointerVec.push_back(&f);

  if ( name == "emptyName") {
    _name.push_back(f.getName() );
  }
  _name.push_back(name);

  _autoScale.push_back(false);
  _minValue.push_back(minValue);
  _maxValue.push_back(maxValue);
}




/// platform specific code
//template<typename T>
//void BlockGifWriter<T>::ppm2gif(const std::string& fullNamePpm,
//                                const std::string& fullNameGif)
//{
//  std::string consoleCommand = "convert " + fullNamePpm
//                               + " -resize 600x800^ "
//                               + fullNameGif;
//  system(consoleCommand.c_str());
//  consoleCommand = std::string("rm ") + fullNamePpm;
//  system(consoleCommand.c_str());
//}






}  // namespace olb

#endif
