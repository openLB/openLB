/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Marc Haussmann
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

#ifndef GNUPLOT_HEATMAP_WRITER_H
#define GNUPLOT_HEATMAP_WRITER_H

#include <iomanip>
#include <iostream>
#include <list>

#include "functors/lattice/blockReduction3D2D.h"
#include "functors/lattice/blockReduction2D2D.h"

namespace olb {

namespace heatmap {

template <typename T>
struct plotParam
{
  std::string name;
  bool writeCSV = false;
  int contourlevel = 0;
  std::string colour = std::string("rainbow");
  Vector<T,2u> zoomOrigin;
  Vector<T,2u> zoomExtend;
  T minValue = 0.;
  T maxValue = 0.;

  plotParam():
    zoomOrigin(0.0, 0.0),
    zoomExtend(1.0, 1.0)
  {
  }

};

/** This function is used to plot heat maps as jpeg files.
 * minValue and maxValue set a defined scalar range.
 * contourlevel sets the number of contours in the plot.
 * zoomOrigin and zoomExtend set a zoom scale in the plot:
 * (zoomOrigin(0 to 1) and zoomExtend(0 to 1-zoomOrigin)
 * Available colour schemes are "grey", "pm3d", "blackbody" and "rainbow"
 */
template <typename T>
void write(BlockReduction3D2D<T>& blockReduction, int iT, const plotParam<T> param = {});

/** This function is used to plot heat maps as jpeg files.
 * minValue and maxValue set a defined scalar range.
 * contourlevel sets the number of contours in the plot.
 * zoomOrigin and zoomExtend set a zoom scale in the plot:
 * (zoomOrigin(0 to 1) and zoomExtend(0 to 1-zoomOrigin)
 * Available colour schemes are "grey", "pm3d", "blackbody" and "rainbow"
 */
template <typename T>
void write(BlockReduction2D2D<T>& blockReduction, int iT, const plotParam<T> param = {});




namespace detail {

template <typename T>
void genericHeatMapInterface(const HyperplaneLattice3D<T>& hyperPlane, BlockF2D<T>& blockData, int iT,
                             const plotParam<T>& param);



template <typename T>
struct detailParam
{
  const plotParam<T>* plot;
  BlockF2D<T>* blockData;
  const HyperplaneLattice3D<T>* hyperPlane;
  std::string matrixPath;
  std::string csvPath;
  std::string jpegPath;
  std::string pngPath;
  std::string plotFilePath;
  std::string dir;
  int iT;
  int nx;
  int ny;
  T spacing;
  Vector<T,3> origin;
  Vector<T,3> normal;
  Vector<T,2> zoomMin;
  Vector<T,2> zoomMax;
  std::string quantityname;
};

template <typename T>
void writeHeatMapDataFile(detailParam<T>& param);

template <typename T>
void writeHeatMapPlotFile(detailParam<T>& param);

template< typename T >
void executeGnuplot(detailParam<T>& param);

} // namespace detail

} // namespace heatmap

} // namespace olb

#endif
