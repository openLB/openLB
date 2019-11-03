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

#ifndef GNUPLOT_HEATMAP_WRITER_HH
#define GNUPLOT_HEATMAP_WRITER_HH

#include <fstream>
#include <iostream>
#include <unistd.h>

#include "core/singleton.h"
#include "io/fileName.h"
#include "gnuplotHeatMapWriter.h"
#include "utilities/vectorHelpers.h"

namespace olb {

namespace heatmap {

template <typename T>
void write(BlockReduction3D2D<T>& blockReduction, int iT, const plotParam<T> param) {
  detail::genericHeatMapInterface(blockReduction, blockReduction, iT, param);
}

template <typename T>
void write(BlockReduction2D2D<T>& blockReduction, int iT, const plotParam<T> param) {
  detail::genericHeatMapInterface(blockReduction.getPlaneDiscretizationIn3D(), blockReduction, iT,
                                  param);
}

namespace detail {

template <typename T>
void genericHeatMapInterface(const HyperplaneLattice3D<T>& hyperPlane, BlockF2D<T>& blockData, int iT,
                             const plotParam<T>& plot) {


  if ( blockData.getTargetDim() != 1 ) {
    std::cout << "HeatMap Error: Functor targetDim is not 1. " << std::endl;
    exit(-1);
  }
  else {
    if ( singleton::mpi().getRank() == 0 ) {
      detail::detailParam<T> param;
      param.plot = &plot;
      param.blockData = &blockData;
      param.hyperPlane = &hyperPlane;
      param.dir = singleton::directories().getImageOutDir();
      param.quantityname = param.blockData->getName();
      param.quantityname.erase (param.quantityname.begin(), param.quantityname.begin()+15);
      param.quantityname.pop_back();
      std::string name;
      if (plot.name == "") {
        name = param.quantityname;
      } else {
        name = plot.name;
      }
      param.matrixPath = param.dir + "data/" + name  + ".matrix";
      param.csvPath = createFileName( param.dir + "data/", name, iT) + ".csv";
      param.jpegPath = createFileName( param.dir, name, iT) + ".jpeg";
      param.pngPath = createFileName( param.dir, name, iT) + ".png";
      param.plotFilePath = param.dir + "data/heatMap"+name+".p";

      for (size_t pos = 0; pos < param.plotFilePath.length(); pos++) {
        if (param.plotFilePath.at(pos) == '(' || param.plotFilePath.at(pos) == ')') {
          param.plotFilePath.replace(pos,  1,  "");
        }
      }

      param.nx = param.hyperPlane->getNx();
      param.ny = param.hyperPlane->getNy();
      param.spacing = param.hyperPlane->getPhysSpacing();
      param.origin = param.hyperPlane->getPhysOrigin();
      param.normal = crossProduct3D(param.hyperPlane->getVectorU(),param.hyperPlane->getVectorV());
      param.zoomMin = plot.zoomOrigin;
      param.zoomMax = plot.zoomOrigin + plot.zoomExtend;
      param.iT = iT;

      writeHeatMapDataFile(param);
      writeHeatMapPlotFile(param);
      executeGnuplot(param);
    }
  }
  return;
}

template <typename T>
void writeHeatMapDataFile(detailParam<T>& param) {


  std::ofstream foutMatrix( param.matrixPath.c_str() );

  int i[2] = {0,0};
  for (i[1] = param.ny * param.zoomMin[1]; i[1] < param.ny * param.zoomMax[1]; i[1]++) {
    for (i[0] = param.nx * param.zoomMin[0]; i[0] < param.nx * param.zoomMax[0]; i[0]++) {
      T evaluated[1];
      (*param.blockData)(evaluated,i);
      foutMatrix << evaluated[0] << " ";
      if (i[0] == int(param.nx * param.zoomMax[0]) - 1) {
        foutMatrix  << "\n";
      }
    }
  }
  foutMatrix.close();


  if (param.plot->writeCSV) {
    std::ofstream foutCSV( param.csvPath.c_str() );
    for (i[1] = 0; i[1] < param.ny; i[1]++) {
      for (i[0] = 0; i[0] < param.nx; i[0]++) {
      T evaluated[1];
      Vector<T,3> physPoint = param.hyperPlane->getPhysR(i[0], i[1]);
      (*param.blockData)(evaluated,i);
      foutCSV << physPoint[0] << " " << physPoint[1] << " " << physPoint[2] << " " << evaluated[0] << "\n";
      }
    }
    foutCSV.close();
  }

  return;
}

template <typename T>
void writeHeatMapPlotFile(detailParam<T>& param) {
  std::ofstream fout;

  fout.open(param.plotFilePath.c_str(), std::ios::trunc);

  fout << "if (strstrt(GPVAL_TERMINALS, 'jpeg') > 0) {";
  fout << "set terminal jpeg " << "size " << 1920  << "," << 1080 << "font \",25\"" << "\n";
  fout << "set output '"<< param.jpegPath << "'"<< "\n";
  fout << "} else {";
  fout << "set terminal png " << "size " << 1920  << "," << 1080 << "font \",25\"" << "\n";
  fout << "set output '"<< param.pngPath << "'"<< "\n";
  fout << "}" << "\n";
  fout << "set pm3d map" << "\n";
  fout << "unset key" << "\n";
  fout << "set size ratio -1" << "\n";
  fout << "set size 0.925,1.0" << "\n";
  fout << "set xtics out" << "\n";
  fout << "set ytics out" << "\n";
  fout << "set xtics nomirror" << "\n";
  fout << "set ytics nomirror" << "\n";
  fout << "set pm3d interpolate 0,0" << "\n";
  fout << "set colorbox vertical user origin 0.85,0.1 size 0.025 ,0.8" << "\n";

  if ( util::nearZero(param.normal[0]) && util::nearZero(param.normal[1]) ) {
    fout << "set xlabel \"x-axis in m \"" << "\n"
         << "set ylabel \"y-axis in m \"" << "\n";
  } else if ( util::nearZero(param.normal[0]) && util::nearZero(param.normal[2]) ) {
    fout << "set xlabel \"x-axis in m \"" << "\n"
         << "set ylabel \"z-axis in m \"" << "\n";
    param.origin[1] = param.origin[2];
  } else if ( util::nearZero(param.normal[1]) && util::nearZero(param.normal[2]) ) {
    fout << "set xlabel \"y-axis in m \"" << "\n"
         << "set ylabel \"z-axis in m \"" << "\n";
    param.origin[0] = param.origin[1];
    param.origin[1] = param.origin[2];
  } else {
    fout << "set xlabel \"width in m \"" << "\n"
         << "set ylabel \"height in m \"" << "\n";

  }

  if (param.plot->contourlevel > 0) {
    fout << "set contour base" << "\n";
    fout << "set cntrparam levels " << param.plot->contourlevel << "\n";
    fout << "set cntrparam bspline" << "\n";
    fout << "do for [i=1:"<< param.plot->contourlevel <<"] {" << "\n";
    fout << "set linetype i lc rgb \"black\""<< "\n";
    fout << "}" << "\n";
  }

  fout << "set cblabel offset 0.5 \"" <<param.quantityname;
  if (param.quantityname == "l2(physVelocity)" || param.quantityname == "EuklidNorm(physVelocity)") {
    fout << " in m/s\"" << "\n";
  } else if (param.quantityname == "physPressure") {
    fout << " in Pa\"" << "\n";
  } else {
    fout << "\"\n";
  }

  if (!util::nearZero(param.plot->maxValue - param.plot->minValue)) {
    fout << "set cbrange [" << param.plot->minValue << ":" << param.plot->maxValue <<"]" << "\n";
  }

  fout << "set autoscale fix" << "\n";
  if (param.plot->colour == "grey") {
    fout << "set palette grey" << "\n";
  } else if (param.plot->colour == "pm3d") {

  } else if (param.plot->colour == "blackbody") {
    fout << "set palette defined ( 0 \"black\", 1 \"red\", 2 \"yellow\")" << "\n";
  }
  else {
    fout << "set palette defined ( 0 \"blue\", 1 \"green\", 2 \"yellow\", 3 \"orange\", 4 \"red\" )" << "\n";
  }
  fout << "splot '" << param.matrixPath << "' u ($1*" << param.spacing << "+" << param.origin[0] + int(param.nx * param.zoomMin[0])*param.spacing << "):"
                    << "($2*" << param.spacing << "+" << param.origin[1] + int(param.ny * param.zoomMin[1])*param.spacing <<"):3 matrix with pm3d" << "\n";
  fout.close();
  return;
}

template< typename T >
void executeGnuplot(detailParam<T>& param)
{
#ifdef WIN32
  std::cout << "GNUPLOT WORKS ONLT WITH LINUX" << std::endl;
//  exit (EXIT_FAILURE);
  return;
#endif
#ifndef WIN32
  if (!system(nullptr)) {
    exit (EXIT_FAILURE);
  }
  const std::string command = "gnuplot "+ param.plotFilePath +" > /dev/null &";
  if ( system(command.c_str()) ) {
    std::cout << "Error at GnuplotWriter" << std::endl;
  }
  return;
#endif
}

} // namespace detail

} // namespace heatmap

} // namespace olb

#endif
