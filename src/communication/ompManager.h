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

#ifndef OMP_MANAGER_H
#define OMP_MANAGER_H

#ifdef PARALLEL_MODE_OMP

class ompManager {

  // This must be "public" since gcc "threadprived"
  // requires that the class ompManager is of pod-type,
  // hence "privat" is not allowed. That is a bug in gcc.
public:
  int size;
  int rank;
public:
  void init();
  int get_size() const;
  int get_rank() const;
  void set_dynamic(int dynamicThreads);
};

extern ompManager omp;
#pragma omp threadprivate (omp)

#endif

#endif
