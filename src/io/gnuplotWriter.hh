/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Fabian Klemens
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

#ifndef GNUPLOT_WRITER_HH
#define GNUPLOT_WRITER_HH

#include <fstream>
#include <iostream>
#include <unistd.h>
#include "core/singleton.h"
#include "io/fileName.h"
#include "gnuplotWriter.h"
#include "utilities/vectorHelpers.h"

namespace olb {
/// Constructor with name of outputFiles
/// boolean true for real-time plotting //WARNING: experimental!
template< typename T >
Gnuplot<T>::Gnuplot(std::string name, bool liveplot)
  : _name(name),
    _dataFile(singleton::directories().getGnuplotOutDir()+"data/"+_name+".dat"),
    _dir(singleton::directories().getGnuplotOutDir())
{
  _liveplot = liveplot;
  if (singleton::mpi().getRank() == _rank) {
    std::ofstream fout;

    ///add (new) data file
    fout.open(_dataFile.c_str(), std::ios::trunc);
    fout.close();
  }
}

/// writes the data and plot file for two doubles (x and y)
/// plotType indicates whether you want a linegraph 'l' (default) or a scatterplot 'p' (default: 'l')
template< typename T >
void Gnuplot<T>::setData(T xValue, T yValue, std::string name, std::string key, char plotType)
{
  if (_init) {
    _dataSize = 1;
    _key = key;
    _names = {name};
    _plotTypes = {plotType};

    if (_liveplot) {
      writePlotFile("plot");
    }
  }
  writeDataFile(xValue, yValue);

  if (_liveplot && _init) {
    startGnuplot("plot");
  }

  _init = false;
  return;
}

/// writes the data and plot file for two doubles (x and y), where x is increasing integer
template< typename T >
void Gnuplot<T>::setData(bool noXvalue, T yValue, std::string name, std::string key, char plotType)
{
  T xValue = _time;
  setData(xValue, yValue, name, key, {plotType});
  _time++;
}

/// writes the data and plot file for a double and a vector of doubles (x and y1,y2,...)
/// plotType indicates whether you want a linegraph 'l' (default) or a scatterplot 'p': (default: {'l','l'})
/// The position in the vector 'plotType'{'l', 'p'} is linked to the rank of the y-axis (y1, y2) :
/// y1 is plotted in form of a line plot & y2 is plotted in form of a scatterplot
template< typename T >
void Gnuplot<T>::setData(T xValue, std::vector<T> yValues, std::vector<std::string> names, std::string key, std::vector<char> plotType)
{
  if (_init) {
    _dataSize = yValues.size();
    _key = key;
    _names = names;
    _plotTypes = plotType;
    if (_names.size() < _dataSize) {
      for (unsigned int i = _names.size(); i < _dataSize; i++) {
        _names.push_back("");
      }
    }
    if (_plotTypes.size() < _dataSize) {
      for (unsigned int i = _plotTypes.size(); i < _dataSize; i++) {
        _plotTypes.push_back('l');
      }
    }
    if (_liveplot) {
      writePlotFile("plot");
    }
  }
  writeDataFile(xValue,yValues);

  if (_liveplot && _init) {
    startGnuplot("plot");
  }

  _init = false;
  return;
}

/// writes the data and plot file for a double and a vector of doubles (x and y1,y2,...), where x is increasing integer
template< typename T >
void Gnuplot<T>::setData(bool noXvalue, std::vector<T> yValues, std::vector<std::string> names, std::string key, std::vector<char> plotType)
{
  T xValue = _time;
  setData(xValue, yValues, names, key, plotType);
  _time++;
}


/// writes an PDF
template< typename T >
void Gnuplot<T>::writePDF()
{
  if (!_init) {
    writePlotFile("pdf");
    startGnuplot("plotPDF");
  }
  return;
}


/// writes PNGs
/// usage: first argument: numbering of png file
/// second argument: range for the x axis
/// thrid argument: specifies the name of the plot in case the user wants to
/// create more than one plot with the simulation results (default: plotName = "")
/// no arguments: writes consecutive numbers with adaptive xrange
template< typename T >
void Gnuplot<T>::writePNG(int iT, double xRange, std::string plotName)
{
  if (!_init) {
    _iT = iT;
    _xRange = xRange;

    /// initialize the writePlotFile for Gnuplot with the type and the name of the output data
    writePlotFile("png", plotName);
    startGnuplot("plotPNG", plotName);
  }
  return;
}

/// plotName specifies the name of the plot in case the user wants to create more than
/// one plot with the simulation results (default: plotName = "")
template< typename T >
void Gnuplot<T>::writePlotFile(std::string type, std::string plotName)
{
  if (singleton::mpi().getRank() == _rank ) {
    std::ofstream fout;

    std::string plotFile;
    if (_liveplot && type == "plot") {
      plotFile = singleton::directories().getGnuplotOutDir()+"data/plot.p";
    } else if (type == "pdf") {
      plotFile = singleton::directories().getGnuplotOutDir()+"data/plotPDF.p";
    } else if (type == "png") {
      plotFile = singleton::directories().getGnuplotOutDir()+"data/plotPNG"+plotName+".p";
    } else {
      std::cout << "WARNING: invalid Gnuplot type={'', 'plot'; 'pdf', 'png'}" << std::endl;
      exit(-1);
    }

    fout.open(plotFile.c_str(), std::ios::trunc);
    fout << "set key " << _key << "\n";

    if (type=="pdf") {
      fout << "set terminal pdf enhanced" << "\n"
           << "set output '"<<_dir<<_name<<".pdf'" << "\n";
    }
    if (type=="png") {
      if ( !util::nearZero(_xRange+1) ) {
        fout << "set xr[0:"<< _xRange <<"]" << "\n";
      }
      fout << "set terminal png" << "\n"
           << "set output '"<<_dir<<_name;
      if (_iT != -1) {
        fout <<"_"<<_iT;
      }
      fout <<".png'" << "\n";
    }
    /// set the x and y label of the Plot
    fout << "set xlabel '" << _xLabel << "'" << "\n";
    fout << "set ylabel '" << _yLabel << "'" << "\n";

    /// vector which holds the information about the plotType
    /// (e.g. scatterplot 'p' or lineplot 'l': default {'l','l'})
    fout << "plot '"<<_dataFile<<"' u 1:2 w " << _plotTypes[0] << " t '"<< _names[0] << "'";
    for (unsigned int i = 0; i < _dataSize-1; ++i) {
      fout << ", '"<<_dataFile<<"' u 1:" << i+3 << " w " << _plotTypes[i+1] << " t '" << _names[i+1] << "'";
    }
    fout << "\n";
    if (_liveplot && type=="plot") {
      fout << "pause -1" << "\n"
           << "reread" << "\n";
    }
    fout.close();
  }
  return;
}


/// writes the data file for two doubles (x and y)
template< typename T >
void Gnuplot<T>::writeDataFile(T xValue, T yValue)
{
  if (singleton::mpi().getRank() == _rank) {
    std::ofstream fout;
    fout.precision(6);
    fout.open(_dataFile.c_str(), std::ios::app);
    fout << xValue
         << " "
         << yValue
         << std::endl;
    fout.close();
  }
  return;
}

/// set Label of the gnuplotPlot; xLabel and yLabel
template< typename T >
void Gnuplot<T>::setLabel(std::string xLabel, std::string yLabel)
{
  _xLabel = xLabel;
  _yLabel = yLabel;
}


/// writes the data file for one double and a vector of doubles (x and y1,y2,...)
template< typename T >
void Gnuplot<T>::writeDataFile(T xValue, std::vector<T> yValues)
{
  if (singleton::mpi().getRank() == _rank) {
    std::ofstream fout;
    fout.precision(6);
    fout.open(_dataFile.c_str(), std::ios::app);
    fout << xValue;
    for (unsigned int i = 0; i < yValues.size(); i++) {
      fout << " " << yValues[i];
    }
    fout << "\n";
    fout.close();
  }
  return;
}


/// system command to start gnuplot (LINUX ONLY!)
/// plotName indicates the name of the plot in case the user wants to create
/// more than one plot with the simulation results (default: plotName = "")
template< typename T >
void Gnuplot<T>::startGnuplot(std::string plotFile, std::string plotName)
{
#ifdef WIN32
  std::cout << "GNUPLOT WORKS ONLT WITH LINUX" << std::endl;
//  exit (EXIT_FAILURE);
  return;
#endif
#ifndef WIN32
  if (singleton::mpi().getRank() == _rank) {
    if (!system(nullptr)) {
      exit (EXIT_FAILURE);
    }
    const std::string command = "gnuplot -persistent "+_dir+"data/"+plotFile+plotName+".p > /dev/null &";
    if ( system(command.c_str()) ) {
      std::cout << "Error at GnuplotWriter" << std::endl;
    }
  }
  return;
#endif
}
}  // namespace olb

#endif
