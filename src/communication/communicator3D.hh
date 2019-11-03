/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias J. Krause
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
 * A communincator provides a cuboids with cells of other
 * cuboids -- generic implementation.
 */


#ifndef COMMUNICATOR_3D_HH
#define COMMUNICATOR_3D_HH

#include "communication/mpiManager.h"
#include <vector>
#include "communication/communicator3D.h"
#include "communication/cuboidNeighbourhood3D.h"
#include "geometry/cuboidGeometry3D.h"
#include "communication/superStructure3D.h"


namespace olb {

/////////////////// Class Communicator3D //////////////////////


template<typename T>
Communicator3D<T>::Communicator3D(SuperStructure3D<T>& superStructure):_superStructure(superStructure)
{
  _initDone = false;
}

template<typename T>
void Communicator3D<T>::init_nh()
{

  _nC = _superStructure.getCuboidGeometry().getNc();

  for (int iC=0; iC<_superStructure.getLoadBalancer().size(); iC++) {
    CuboidNeighbourhood3D<T> nh(_superStructure,_superStructure.getLoadBalancer().glob(iC));
    _nh.push_back(nh);
  }
}

template<typename T>
void Communicator3D<T>::add_cell(int iCloc, int iX, int iY, int iZ)
{

  _nh[iCloc].add_inCell(iX,iY,iZ);
}

template<typename T>
void Communicator3D<T>::add_cells(int overlap)
{

  for (int iC=0; iC<_superStructure.getLoadBalancer().size(); iC++) {
    _nh[iC].add_inCells(overlap);
  }
}

template<typename T>
void Communicator3D<T>::init()
{

  reset();
  for (int iC=0; iC<_superStructure.getLoadBalancer().size(); iC++) {
    _nh[iC].init_inCN();
    for (int i=0; i<_nh[iC].get_inCellsSize(); i++) {
      int ID = _nh[iC].get_inCell(i).latticeR[0];
#ifdef PARALLEL_MODE_MPI
      if ( singleton::mpi().getRank() ==_superStructure.getLoadBalancer().rank(ID) )
#endif
      {
        Cell3D<T> temp;
        temp.physR[0] = _nh[iC].get_inCell(i).physR[0];
        temp.physR[1] = _nh[iC].get_inCell(i).physR[1];
        temp.physR[2] = _nh[iC].get_inCell(i).physR[2];
        _superStructure.getCuboidGeometry().getLatticeR(temp.latticeR, temp.physR);
        temp.latticeR[0]    = _superStructure.getLoadBalancer().glob(iC);
        _nh[_superStructure.getLoadBalancer().loc(ID)].add_outCell(temp);
      }
    }
  }

  for (int iC=0; iC<_superStructure.getLoadBalancer().size(); iC++) {
    _nh[iC].init_outCN();
  }

#ifdef PARALLEL_MODE_MPI
  singleton::MpiNonBlockingHelper helper;
  for (int iC=0; iC<_superStructure.getLoadBalancer().size(); iC++) {
    _nh[iC].finish_comm();
  }
  for (int iC=0; iC<_superStructure.getLoadBalancer().size(); iC++) {
    _nh[iC].bufSend_inCells(helper);
  }
  for (int iC=0; iC<_superStructure.getLoadBalancer().size(); iC++) {
    _nh[iC].recWrite_outCells();
  }
  for (int iC=0; iC<_superStructure.getLoadBalancer().size(); iC++) {
    singleton::mpi().waitAll(helper);
  }
#endif
}

template<typename T>
void Communicator3D<T>::reset()
{

  if (_initDone) {
    for (int iC=0; iC<_superStructure.getLoadBalancer().size(); iC++) {
      _nh[iC].reset();
    }
    _initDone = false;
  }
}

template<typename T>
void Communicator3D<T>::send()
{

  for (int iC=0; iC<_superStructure.getLoadBalancer().size(); iC++) {
    _nh[iC].buffer_outData();
#ifdef PARALLEL_MODE_MPI
    _nh[iC].send_outData();
#endif
  }
}

template<typename T>
void Communicator3D<T>::receive()
{

  for (int iC=0; iC<_superStructure.getLoadBalancer().size(); iC++) {
#ifdef PARALLEL_MODE_MPI
    _nh[iC].receive_inData();
#else
    for (int i=0; i<_nh[iC].get_inCsize(); i++) {
      _nh[iC].get_inData()[_nh[iC].get_inC(i)] =
        _nh[_nh[iC].get_inC(i)].get_outData()[iC];
    }
#endif
  }
#ifdef PARALLEL_MODE_MPI
  for (int iC=0; iC<_superStructure.getLoadBalancer().size(); iC++) {
    _nh[iC].finish_comm();
  }
#endif
}

template<typename T>
void Communicator3D<T>::write()
{

  for (int iC=0; iC<_superStructure.getLoadBalancer().size(); iC++) {
    _nh[iC].write_inData();
  }
}

}

#endif
