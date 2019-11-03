/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2011 Lukas Baron, Mathias Krause
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

#ifndef OSTREAM_MANAGER_HH
#define OSTREAM_MANAGER_HH

#include "ostreamManager.h"

#include "communication/mpiManager.h"
#include "communication/ompManager.h"

namespace olb {

// class OMBuf /////////////////////////////////////////


bool OMBuf::multiOutput = 0;

OMBuf::OMBuf()
  : output(nullptr), text("")
{ }

OMBuf::~OMBuf()
{ }

OMBuf::OMBuf(const OMBuf& rhs)
  : output(rhs.output), text(rhs.text)
{ }

OMBuf& OMBuf::operator=(const OMBuf& rhs)
{
  output = rhs.output;
  text = rhs.text;

  return *this;
}

OMBuf::OMBuf(std::ostream& str, std::string classname)
  : output(&str), text(classname)
{ }

void OMBuf::setMultiOutput(bool b)
{
  (*this).multiOutput = b;
}

int OMBuf::sync()
{
#ifdef PARALLEL_MODE_MPI
  if (multiOutput==true) {
    *output << "["
            << text << ":" << singleton::mpi().getRank()
            << "] " << str();
  } else {  // multiOutput==false
    if (singleton::mpi().getRank()==0) {
      *output << "[" << text << "] " << str();
    }
  }
#elif PARALLEL_MODE_OMP
  if (multiOutput==true) {
    *output << "["
            << text << ":" << omp.get_rank()
            << "] " << str();
  } else {  // multiOutput==false
    if (omp.get_rank()==0) {
      *output << "[" << text << "] " << str();
    }
  }
#else
  *output << "[" << text << "] " << str();
#endif
  str("");
  output->flush();
  return 0;
}

// class OstreamManager /////////////////////////////////

/* this function should not need to be called at all,
 * neither implicitly nor explicitly
OstreamManager::OstreamManager()
    : std::ostream(&buffer), buffer(std::cout, "???")
{ }
*/

OstreamManager::OstreamManager(std::string classname)
  : std::ostream(&buffer), buffer(std::cout, classname)
{ }

OstreamManager::OstreamManager(std::ostream& str, std::string classname)
  : std::ostream(&buffer), buffer(str, classname)
{ }

OstreamManager::OstreamManager(const OstreamManager& rhs)
  : std::ostream(&buffer), buffer(rhs.buffer)
{ }

OstreamManager& OstreamManager::operator=(const OstreamManager& rhs)
{
  buffer = rhs.buffer;
  return *this;
}

OstreamManager::~OstreamManager() { }

void OstreamManager::setMultiOutput(bool b)
{
  buffer.setMultiOutput(b);
}

} // end namespace olb

#endif
