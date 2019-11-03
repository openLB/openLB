/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Jonas Latt
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
 * LB initialisation routine -- header file.
 */
#ifndef OLB_INIT_H
#define OLB_INIT_H

#include "communication/mpiManager.h"
#include "io/ostreamManager.h"
#include "io/parallelIO.h"
#include "communication/ompManager.h"

namespace olb {

inline void olbInit(int *argc, char ***argv, bool verbose=false)
{

  // create an OstreamManager object in order to enable multi output
  olb::OstreamManager clout(std::cout,"olbInit");
  clout.setMultiOutput(verbose);
  singleton::mpi().init(argc, argv);

#ifdef PARALLEL_MODE_MPI
  /*ParBuf *newCoutBuf = new ParBuf(std::cout.rdbuf());
  ParBuf *newClogBuf = new ParBuf(std::clog.rdbuf());
  ParBuf *newCinBuf  = new ParBuf(std::cin.rdbuf());

  std::cout.rdbuf(newCoutBuf);
  std::clog.rdbuf(newClogBuf);
  std::cin. rdbuf(newCinBuf);*/
#endif

#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel
  omp.init();
#endif

}

}  // namespace olb

#endif
