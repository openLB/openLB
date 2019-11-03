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
 * The description of a 3D cuboid neighbourhood -- generic implementation.
 */


#ifndef CUBOID_NEIGHBOURHOOD_3D_HH
#define CUBOID_NEIGHBOURHOOD_3D_HH

#include "communication/mpiManager.h"
#include <vector>
#include <string>
#include "communication/cuboidNeighbourhood3D.h"
#include "geometry/cuboidGeometry3D.h"
#include "dynamics/dynamics.h"
#include "core/cell.h"
#include "communication/superStructure3D.h"


namespace olb {


//////////////// Class CuboidNeighbourhood3D //////////////////

template<typename T>
CuboidNeighbourhood3D<T>::CuboidNeighbourhood3D(
  SuperStructure3D<T>& superStructure, int iC):_superStructure(superStructure)
{
  _iCglob        = iC;
  _nC            = _superStructure.getCuboidGeometry().getNc();
  _deltaC        = _superStructure.getCuboidGeometry().get(iC).getDeltaR();
  _nData         = _superStructure.getDataSize();
  _nDataType     = _superStructure.getDataTypeSize();
  _initInCNdone  = false;
  _initOutCNdone = false;
}

template<typename T>
CuboidNeighbourhood3D<T>::CuboidNeighbourhood3D (
  CuboidNeighbourhood3D<T> const& rhs ):_superStructure(rhs._superStructure)
{
  _iCglob        = rhs._iCglob;
  _nC            = rhs._nC;
  _deltaC        = rhs._deltaC;
  _inCells       = rhs._inCells;
  _outCells      = rhs._outCells;
  _inC           = rhs._inC;
  _inN           = rhs._inN;
  _outC          = rhs._outC;
  _outN          = rhs._outN;
  _nData         = rhs._nData;
  _nDataType     = rhs._nDataType;
  _initInCNdone  = false;
  _initOutCNdone = false;
}

template<typename T>
CuboidNeighbourhood3D<T> CuboidNeighbourhood3D<T>::operator= (
  CuboidNeighbourhood3D<T> rhs )
{
  CuboidNeighbourhood3D<T> tmp(rhs);
  return tmp;
}

template<typename T>
CuboidNeighbourhood3D<T>::~CuboidNeighbourhood3D<T>()
{
  reset();
}


template<typename T>
Cell3D<T> const& CuboidNeighbourhood3D<T>::get_inCell(int i) const
{
  return _inCells[i];
}

template<typename T>
int CuboidNeighbourhood3D<T>::get_inCellsSize() const
{
  return _inCells.size();
}

template<typename T>
int const& CuboidNeighbourhood3D<T>::get_inC(int i) const
{
  return _inC[i];
}

template<typename T>
int CuboidNeighbourhood3D<T>::get_inCsize() const
{
  return _inC.size();
}

template<typename T>
bool** CuboidNeighbourhood3D<T>::get_inData()
{
  return _inData;
}

template<typename T>
bool** CuboidNeighbourhood3D<T>::get_outData()
{
  return _outData;
}


template<typename T>
void CuboidNeighbourhood3D<T>::add_inCell(Cell3D<T> cell)
{
  _inCells.push_back(cell);
}

template<typename T>
void CuboidNeighbourhood3D<T>::add_outCell(Cell3D<T> cell)
{
  _outCells.push_back(cell);
}

template<typename T>
void CuboidNeighbourhood3D<T>::add_inCell(int iX, int iY, int iZ)
{

  Cell3D<T> found;
  found.latticeR[0] = _iCglob;
  found.latticeR[1] = iX;
  found.latticeR[2] = iY;
  found.latticeR[3] = iZ;

  _superStructure.getCuboidGeometry().getPhysR(found.physR, found.latticeR);
  std::vector<T> tmp(found.physR,found.physR + 3);
  if (_superStructure.getCuboidGeometry().getC(tmp, found.latticeR[0]) ) {
    for (unsigned i=0; i<_inCells.size(); i++) {
      if (_inCells[i]==found) {
        return;
      }
    }
    _inCells.push_back(found);
  }
}

template<typename T>
void CuboidNeighbourhood3D<T>::add_inCells(int overlap)
{

  int nX  = _superStructure.getCuboidGeometry().get(_iCglob).getNx();
  int nY  = _superStructure.getCuboidGeometry().get(_iCglob).getNy();
  int nZ  = _superStructure.getCuboidGeometry().get(_iCglob).getNz();

  for (int iX=0; iX<nX+2*overlap; iX++) {
    for (int iY=0; iY<nY+2*overlap; iY++) {
      for (int iZ=0; iZ<nZ+2*overlap; iZ++) {
        if (iX < overlap || iX > nX + overlap - 1 ||
            iY < overlap || iY > nY + overlap - 1 ||
            iZ < overlap || iZ > nZ + overlap - 1 ) {
          Cell3D<T> found;
          found.latticeR[0] = _iCglob;
          found.latticeR[1] = iX - overlap;
          found.latticeR[2] = iY - overlap;
          found.latticeR[3] = iZ - overlap;

          _superStructure.getCuboidGeometry().getPhysR(found.physR, found.latticeR);
          std::vector<T> tmp(found.physR,found.physR + 3);
          if (_superStructure.getCuboidGeometry().getC(tmp, found.latticeR[0]) ) {
            _inCells.push_back(found);
          }
        }
      }
    }
  }
}

template<typename T>
void CuboidNeighbourhood3D<T>::init_inCN()
{

  _inC.clear();
  _inN.clear();

  _inData = new bool* [_nC];
  _inDataCoordinates = new T* [_nC];
  _tempInCN = new int [_nC];
  for (int i=0; i<_nC; i++) {
    _tempInCN[i]=0;
  }

  for (unsigned i=0; i<_inCells.size(); i++) {
    _tempInCN[_inCells[i].latticeR[0]]++;
  }
  for (int i=0; i<_nC; i++) {
    if (_tempInCN[i]!=0) {
      _inC.push_back(i);
      _inN.push_back(_tempInCN[i]);
#ifdef PARALLEL_MODE_MPI
      _inData[i] = new bool [_tempInCN[i]*_nData*_nDataType];
      _inDataCoordinates[i] = new T [_tempInCN[i]*3];
#endif
    }
#ifdef PARALLEL_MODE_MPI
    else {
      _inData[i] = NULL;
      _inDataCoordinates[i] = NULL;
    }
#endif
  }

#ifdef PARALLEL_MODE_MPI
  int counter=0;
  for (int i=0; i<_nC; i++) {
    // int dRank = _superStructure.getLoadBalancer().rank(i);
    counter++;
  }
  _mpiNbHelper.allocate(counter);
  counter=0;
  for (int i=0; i<_nC; i++) {
    int dRank = _superStructure.getLoadBalancer().rank(i);
    singleton::mpi().iSend(&_tempInCN[i] , 1, dRank, _mpiNbHelper.get_mpiRequest(counter), _iCglob);
    counter++;
  }
#endif

  _initInCNdone = true;
}

template<typename T>
void CuboidNeighbourhood3D<T>::init_outCN()
{

  _outC.clear();
  _outN.clear();
  _outData = new bool* [_nC];
  _outDataCoordinates = new T* [_nC];

  std::vector<int> temp(_nC,0);

  for (unsigned i=0; i<_outCells.size(); i++) {
    temp[_outCells[i].latticeR[0]]++;
  }

  for (int i=0; i<_nC; i++) {
#ifdef PARALLEL_MODE_MPI
    int sRank = _superStructure.getLoadBalancer().rank(i);
    singleton::mpi().receive(&temp[i], 1, sRank, i);
#endif
    if (temp[i]!=0) {
      _outC.push_back(i);
      _outN.push_back(temp[i]);
    }
    _outData[i] = new bool [temp[i]*_nData*_nDataType];
    _outDataCoordinates[i] = new T [temp[i]*3];
  }

  _initOutCNdone = true;
}

#ifdef PARALLEL_MODE_MPI
template<typename T>
void CuboidNeighbourhood3D<T>::bufSend_inCells(singleton::MpiNonBlockingHelper& helper)
{

  helper.free();

  std::vector<int> temp(_nC,0);
  for (unsigned i=0; i<_inCells.size(); i++) {
    int iC = _inCells[i].latticeR[0];
    if (singleton::mpi().getRank() != _superStructure.getLoadBalancer().rank(iC)) {
      _inDataCoordinates[iC][3*temp[iC]] = _inCells[i].physR[0];
      _inDataCoordinates[iC][3*temp[iC]+1] = _inCells[i].physR[1];
      _inDataCoordinates[iC][3*temp[iC]+2] = _inCells[i].physR[2];
      temp[iC]++;
    }
  }

  int counter=0;
  for (unsigned iC=0; iC<_inC.size(); iC++) {
    //int dRank = _superStructure.get_load().rank(_inC[iC]);
    //if ( singleton::mpi().getRank() != dRank )
    counter++;
  }

  helper.allocate(counter);
  counter=0;
  for (unsigned iC=0; iC<_inC.size(); iC++) {
    int dRank = _superStructure.getLoadBalancer().rank(_inC[iC]);
    //if ( singleton::mpi().getRank() != dRank ) {
    singleton::mpi().iSend( _inDataCoordinates[_inC[iC]],
                            _inN[iC]*3, dRank, helper.get_mpiRequest(counter), _iCglob);
    counter++;
    //}
  }
}
#endif

template<typename T>
void CuboidNeighbourhood3D<T>::recWrite_outCells()
{

#ifdef PARALLEL_MODE_MPI
  for (unsigned iC=0; iC<_outC.size(); iC++) {
    int sRank = _superStructure.getLoadBalancer().rank(_outC[iC]);
    singleton::mpi().receive(_outDataCoordinates[_outC[iC]], _outN[iC]*3, sRank,_outC[iC]);
    if ( singleton::mpi().getRank() != sRank ) {
      //singleton::mpi().receive(_outDataCoordinates[_outC[iC]], _outN[iC]*3, sRank, _outC[iC]);
      Cell3D<T> found;
      for (int i=0; i<_outN[iC]; i++) {
        found.physR[0] = _outDataCoordinates[_outC[iC]][3*i];
        found.physR[1] = _outDataCoordinates[_outC[iC]][3*i+1];
        found.physR[2] = _outDataCoordinates[_outC[iC]][3*i+2];
        _superStructure.getCuboidGeometry().getLatticeR(found.latticeR, found.physR);
        found.latticeR[0] = _outC[iC];
        _outCells.push_back(found);
      }
    }
  }
#endif
}

template<typename T>
void CuboidNeighbourhood3D<T>::finish_comm()
{

#ifdef PARALLEL_MODE_MPI
  singleton::mpi().waitAll(_mpiNbHelper);
#endif

}

template<typename T>
void CuboidNeighbourhood3D<T>::buffer_outData()
{

  std::vector<int> temp(_nC,0);
  int iCloc = _superStructure.getLoadBalancer().loc(_iCglob);
  for (unsigned i=0; i<_outCells.size(); i++) {
    int iC = _outCells[i].latticeR[0];
    // WARNING: Here is interpolation needed if globX, globY
    // are not integers. This needs to be fixed if one will
    // use unstructured grids.
    //for (int iData=0; iData<_nData; iData++) {
    //  memcpy(_outData[iC] + (temp[iC]*_nData + iData)*_nDataType, _superStructure(iCloc,iX+overlap,iY+overlap,iZ+overlap,iData), _nDataType);
    //}
    memcpy(_outData[iC] + temp[iC]*_nData*_nDataType, _superStructure(iCloc,_outCells[i].latticeR[1],_outCells[i].latticeR[2],_outCells[i].latticeR[3],0), _nDataType*_nData);
    temp[iC]++;
  }
}

template<typename T>
void CuboidNeighbourhood3D<T>::send_outData()
{
#ifdef PARALLEL_MODE_MPI
  for (unsigned iC=0; iC<_outC.size(); iC++) {
    int dRank = _superStructure.getLoadBalancer().rank(_outC[iC]);
    singleton::mpi().iSend( _outData[_outC[iC]],
                            _outN[iC]*_nData*_nDataType, dRank, _mpiNbHelper.get_mpiRequest(iC), _iCglob);
  }
#endif
}

template<typename T>
void CuboidNeighbourhood3D<T>::receive_inData()
{
#ifdef PARALLEL_MODE_MPI
  for (unsigned iC=0; iC<_inC.size(); iC++) {
    int sRank = _superStructure.getLoadBalancer().rank(_inC[iC]);
    singleton::mpi().receive(_inData[_inC[iC]], _inN[iC]*_nData*_nDataType, sRank,_inC[iC]);
  }
#endif
}

template<typename T>
void CuboidNeighbourhood3D<T>::write_inData()
{

  int iCloc = _superStructure.getLoadBalancer().loc(_iCglob);
  std::vector<int> temp(_nC,0);
  for (unsigned i=0; i<_inCells.size(); i++) {
    int iC = _inCells[i].latticeR[0];
    //for (int iData=0; iData<_nData; iData++) {
    //  memcpy(_superStructure(iCloc,iX+overlap,iY+overlap,iZ+overlap,iData), _inData[iC] + (temp[iC]*_nData + iData)*_nDataType, _nDataType);
    //}
    //memcpy(_superStructure(iCloc,iX+overlap,iY+overlap,iZ+overlap,0), _inData[iC] + temp[iC]*_nData*_nDataType, _nData*_nDataType);
    //temp[iC]++;
    memcpy(_superStructure(iCloc,_inCells[i].latticeR[1],_inCells[i].latticeR[2],_inCells[i].latticeR[3],0), _inData[iC] + temp[iC]*_nData*_nDataType, _nData*_nDataType);
    temp[iC]++;
  }
}

template<typename T>
void CuboidNeighbourhood3D<T>::reset()
{

  if (_initInCNdone) {
#ifdef PARALLEL_MODE_MPI
    for (int iC=0; iC<_nC; iC++) {
      delete[] _inData[iC];
      delete[] _inDataCoordinates[iC];
    }
#endif
    delete[] _inData;
    delete[] _inDataCoordinates;
    delete[] _tempInCN;
    _initInCNdone = false;
  }
  if (_initOutCNdone) {
    for (int iC=0; iC<_nC; iC++) {
      delete[] _outData[iC];
      delete[] _outDataCoordinates[iC];
    }
    delete[] _outData;
    delete[] _outDataCoordinates;
#ifdef PARALLEL_MODE_MPI
    _mpiNbHelper.free();
#endif
    _initOutCNdone = false;
  }
  _inCells.clear();
  _outCells.clear();
  _inC.clear();
  _outC.clear();
  _inN.clear();
  _outN.clear();
}

}  // namespace olb

#endif
