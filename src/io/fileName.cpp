/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Albert Mink
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


#ifndef FILE_NAME_HH
#define FILE_NAME_HH


#include "fileName.h"
#include "core/singleton.h"

namespace olb {


/// for .pvd masterFile
std::string createFileName(std::string name)
{
  std::stringstream fNameStream;
  fNameStream << name;
  return fNameStream.str();
}

/// used for .pvd file per timeStep iT
std::string createFileName(std::string name, int iT)
{
  std::stringstream fNameStream;
  fNameStream << name  << "_"
              << "iT" << std::setw(7) << std::setfill('0') << iT ;
  return fNameStream.str();
}

/// for pararalle io, e.g. adds "_rank0000001" for rank=1
std::string createParallelFileName(std::string name, bool withSize)
{
  std::stringstream fNameStream;

  std::stringstream fNameStreamTmp;
  fNameStreamTmp << singleton::mpi().getSize();
  std::size_t stringLength = fNameStream.str().length();

  fNameStream << name  << "_rank" << std::setw(stringLength) << std::setfill('0') << singleton::mpi().getRank();
  if (withSize) {
    fNameStream << "_size" << std::setw(stringLength) << std::setfill('0') << singleton::mpi().getSize();
  }
  return fNameStream.str();
}

/// every thread writes his cuboids iC per timeStep iT
std::string createFileName(std::string name,  int iT, int iC)
{
  std::stringstream fNameStream;
  fNameStream << name  << "_"
              << "iT" << std::setw(7) << std::setfill('0') << iT
              << "iC" << std::setw(5) << std::setfill('0') << iC ;
  return fNameStream.str();
}

/// to write functors instantaneously, without adding
std::string createFileName(std::string name, std::string functor, int iT)
{
  std::stringstream fNameStream;
  fNameStream << name <<"_"<< functor << "iT" << std::setw(7) << std::setfill('0') << iT;
  return fNameStream.str();
}

/// to write functors instantaneously, without adding
std::string createFileName(std::string name, std::string functor, int iT, int iC)
{
  std::stringstream fNameStream;
  fNameStream << name <<"_"<< functor << "iT" << std::setw(7) << std::setfill('0') << iT
              << "iC" << std::setw(5) << std::setfill('0') << iC ;
  return fNameStream.str();
}


}  // namespace olb

#endif
