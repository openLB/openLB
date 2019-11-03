/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Mathias Krause
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


#ifdef PARALLEL_MODE_OMP

#include <omp.h>
#include "communication/ompManager.h"
#include "io/ostreamManager.h"

void ompManager::init()
{
  set_dynamic(0);
  size = omp_get_max_threads();
  rank = omp_get_thread_num();
  olb::OstreamManager clout(std::cout,"OmpManager");
  clout << "Sucessfully initialized, numThreads=" << get_size() << std::endl;
}

int ompManager::get_size() const
{
  return size;
}

int ompManager::get_rank() const
{
  return rank;
}

void ompManager::set_dynamic(int dynamicThreads)
{
  omp_set_dynamic(dynamicThreads);
}

ompManager omp;

#endif


