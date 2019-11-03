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

#ifndef FILE_NAME_H
#define FILE_NAME_H

#include <sstream>
#include <iomanip>


namespace olb {

/** \file
 * These functions help you to create file names. By overloading the
 * the function createFileName() the corresponding version is performed.
 */

/// for .pvd masterFile
std::string createFileName(std::string name);

/// used for .pvd file per timeStep iT
std::string createFileName(std::string name, int iT);

/// for parallel io, e.g. adds "_rank0000001" for rank=1, and optional "_size0000016" if withSize==true
std::string createParallelFileName(std::string name, bool withSize=true);

/// every thread writes his cuboids iC per timeStep iT
std::string createFileName(std::string name,  int iT, int iC);

/// to write functors instantaneously, without adding
std::string createFileName(std::string name, std::string functor, int iT=0);

/// to write functors instantaneously, without adding
std::string createFileName(std::string name, std::string functor, int iT, int iC);


}  // namespace olb

#endif
