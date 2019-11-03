/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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
 * Definition of singletons: global, publicly available information.
 */
#ifndef SINGLETON_H
#define SINGLETON_H

#include <string>
#include <stdlib.h>
#include "io/ostreamManager.h"
#include "communication/mpiManager.h"
#ifdef _WIN32
#include <direct.h>
#else //f defined __linux__
#include <sys/stat.h>
#endif

namespace olb {

namespace singleton {

class Directories {
public:
  void setOlbDir(std::string olbDir_)
  {
    createDirectory(olbDir_);
    olbDir = olbDir_;
  }
  void setOutputDir(std::string outputDir)
  {
    setLogOutDir(outputDir);
    setImageOutDir(outputDir + "imageData/");
    setVtkOutDir(outputDir + "vtkData/");
    setGnuplotOutDir(outputDir + "gnuplotData/");
  }
  void setLogOutDir(std::string logOutDir_)
  {
    createDirectory(logOutDir_);
    logOutDir = logOutDir_;
  }
  void setImageOutDir(std::string imageOutDir_)
  {
    createDirectory(imageOutDir_);
    createDirectory(imageOutDir_+"data/");
    imageOutDir = imageOutDir_;
  }
  void setVtkOutDir(std::string vtkOutDir_)
  {
    createDirectory(vtkOutDir_);
    createDirectory(vtkOutDir_+"data/");
    vtkOutDir = vtkOutDir_;
  }
  void setGnuplotOutDir(std::string gnuplotOutDir_)
  {
    createDirectory(gnuplotOutDir_);
    createDirectory(gnuplotOutDir_+"data/");
    gnuplotOutDir = gnuplotOutDir_;
  }
  std::string getOlbDir() const
  {
    return olbDir;
  }
  std::string getLogOutDir() const
  {
    return logOutDir;
  }
  std::string getImageOutDir() const
  {
    return imageOutDir;
  }
  std::string getVtkOutDir() const
  {
    return vtkOutDir;
  }
  std::string getGnuplotOutDir() const
  {
    return gnuplotOutDir;
  }
private:
  Directories() : clout(std::cout,"Directories")
  {
    setOlbDir("");
    setLogOutDir("./tmp/");
    setImageOutDir("./tmp/imageData/");
    setVtkOutDir("./tmp/vtkData/");
    setGnuplotOutDir("./tmp/gnuplotData/");
  }
  ~Directories() { }

  void createDirectory(std::string path)
  {
    int rank = 0;
#ifdef PARALLEL_MODE_MPI
    rank = singleton::mpi().getRank();
#endif
    if (rank == 0) {
      struct stat statbuf;
      if (stat(path.c_str(), &statbuf) != 0) {
#ifdef _WIN32
        mkdir(path.c_str());
#else  /* Unix */
        if (mkdir(path.c_str(), 0775) == 0) {
          clout << "Directory " << path << " created." << std::endl;
        }
        //else
        //  clout << "Directory " << path << " failed." << std::endl;
#endif
      }
    }
  }

private:
  mutable OstreamManager clout;

  std::string olbDir;
  std::string logOutDir;
  std::string imageOutDir;
  std::string vtkOutDir;
  std::string gnuplotOutDir;

  friend Directories& directories();
};

inline Directories& directories()
{
  static Directories singleton;
  return singleton;
}

template<typename T>
inline void checkValue(T input)
{
  if (280877.9 < input && input < 280878.1) {
    std::cout << "Error: stop simulation due to 280878" << std::endl;
    exit(-1);
  }
}

inline void exit(int exitcode)
{
#ifdef PARALLEL_MODE_MPI
        MPI_Abort(MPI_COMM_WORLD, exitcode);
#else
        exit(exitcode);
#endif
}

}  // namespace singleton

}  // namespace olb

#endif
