/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Albert Mink, Mathias J. Krause
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

#ifndef BLOCK_GIF_WRITER_H
#define BLOCK_GIF_WRITER_H

#include <iomanip>
#include "io/colormaps.h"
#include "io/ostreamManager.h"
#include "functors/lattice/blockBaseF2D.h"


namespace olb {

/** BlockGifWriter writes given functor data to image file of format .ppm
 *
 *  There are two different modes
 *  1. mode: image maxValue and minValue are computed every time.
 *           this implies that image data is scaled every time exactly to
 *           interval [0,1]. this may lead to different scaling of the functor data
 *           concerning the time steps
 *  2. mode: image maxValue and minValue are passed by the user once.
 *           this static use, prevents a rescaling effect for each time step.
 *
 *  \param map determines the color of the graphics
 */
template< typename T >
class BlockGifWriter {
public:
  // constructor
  BlockGifWriter( std::string const& map="leeloo" );

  /** writes functor values normed to interval [0,1].
   *  imageValue is computed according to :
   *  if functorValue >= maxValue then the imageValue is clipped to 1
   *  else imageValue = (minValue - functorValue) / (minValue - maxValue)
   */
  void write( BlockF2D<T>& f, T minValue, T maxValue, int iT=0,
              std::string const& name="emptyName" );
  /// writes functor data and determines its min-/maxValue
  void write( BlockF2D<T>& f, int iT=0, std::string const& name="emptyName" );
  ///  writes functors stored at pointerVec
  void write(int iT=0);
  ///  put functor to _pointerVec, to simplify writing process of several functors
  void addFunctor( BlockF2D<T>& f, std::string const& name="emptyName" );
  void addFunctor( BlockF2D<T>& f, T minValue, T maxValue, std::string const& name="emptyName" );

private:
  //  void ppm2gif(const std::string& fullNamePpm, const std::string& fullNameGif);
  mutable OstreamManager clout;
  int _colorRange;    // set to 1024
  int _numColors;     // set to 1024
  graphics::ColorMap<T> _colorMap;

  /// Holds added functor, to simplify the use of write function
  std::vector< BlockF2D<T>* > _pointerVec;
  /// Holds the names of added functors
  std::vector<std::string> _name;
  /// Holds the scale instruction of added functors
  std::vector<bool> _autoScale;
  /// Holds the minValues of added functors
  std::vector<T> _minValue;
  /// Holds the maxValues of added functors
  std::vector<T> _maxValue;
};


}  // namespace olb

#endif
