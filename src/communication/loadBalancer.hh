/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias Krause, Peter Weisbrod
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


#ifndef LOAD_BALANCER_HH
#define LOAD_BALANCER_HH

#include "communication/loadBalancer.h"

namespace olb {


template<typename T>
LoadBalancer<T>::LoadBalancer(int size) : _size(size)
{}

template<typename T>
LoadBalancer<T>::LoadBalancer(int size, std::map<int,int>& loc, std::vector<int>& glob, std::map<int,int>& rank)
  : _size(size), _loc(loc), _glob(glob), _rank(rank)
{}

template<typename T>
LoadBalancer<T>::~LoadBalancer()
{}

template<typename T>
void LoadBalancer<T>::swap(LoadBalancer<T>& loadBalancer)
{
  std::swap(_size, loadBalancer._size);
  _loc.swap(loadBalancer._loc);
  _glob.swap(loadBalancer._glob);
  _rank.swap(loadBalancer._rank);
}

template<typename T>
bool LoadBalancer<T>::isLocal(const int& glob)
{
  return rank(glob) == singleton::mpi().getRank();
}

template<typename T>
int LoadBalancer<T>::loc(const int& glob)
{
  return _loc[glob];
}

template<typename T>
int LoadBalancer<T>::loc(int glob) const
{
  std::map<int,int>::const_iterator iter = _loc.find(glob);
  return iter->second;
}

template<typename T>
int LoadBalancer<T>::glob(int loc) const
{
  return _glob[loc];
}

template<typename T>
int LoadBalancer<T>::rank(const int& glob)
{
  return _rank[glob];
}

template<typename T>
int LoadBalancer<T>::rank(int glob) const
{
  std::map<int,int>::const_iterator iter = _rank.find(glob);
  return iter->second;
}

template<typename T>
int LoadBalancer<T>::size() const
{
  return _size;
}

template<typename T>
bool LoadBalancer<T>::operator==(const LoadBalancer<T>& rhs) const
{
  return _size == rhs._size &&
         _loc == rhs._loc &&
         _glob == rhs._glob &&
         _rank == rhs._rank;
}


template<typename T>
size_t LoadBalancer<T>::getNblock() const
{
  return   4     // _size, plus vector length of _loc, _glob and _rank
           + _loc.size()
           + _rank.size()
           + _glob.size();
}


template<typename T>
size_t LoadBalancer<T>::getSerializableSize() const
{
  return   sizeof(int)         // _size
           + 3 * sizeof(size_t)  // vector length of _loc, _glob and _rank
           + _loc.size() * sizeof(std::pair<int, int>)
           + _rank.size() * sizeof(std::pair<int, int>)
           + _glob.size() * sizeof(int);
}


template<typename T>
bool* LoadBalancer<T>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  size_t sizeBufferIndex = 0;
  bool* dataPtr = nullptr;

  this->registerVar<int>            (iBlock, sizeBlock, currentBlock, dataPtr, _size);
  this->registerMap<int, int>       (iBlock, sizeBlock, currentBlock, sizeBufferIndex, dataPtr, _loc, loadingMode);
  this->registerStdVectorOfVars<int>(iBlock, sizeBlock, currentBlock, sizeBufferIndex, dataPtr, _glob, loadingMode);
  this->registerMap<int, int>       (iBlock, sizeBlock, currentBlock, sizeBufferIndex, dataPtr, _rank, loadingMode);

  return dataPtr;
}


template<typename T>
void LoadBalancer<T>::print(bool multiOutput) const
{
  OstreamManager clout(std::cout,"LoadBalancer");
  clout.setMultiOutput(multiOutput);
  for (unsigned i = 0; i < this->_glob.size(); i++) {
    clout << "glob[" << i << "]=" << this->_glob[i] << std::endl;
  }
  for (auto it = this->_loc.cbegin(); it != this->_loc.cend();
       it++) {
    clout << "loc[" << (*it).first << "]=" << (*it).second << std::endl;
  }
  for (auto it = this->_rank.cbegin(); it != this->_rank.cend();
       it++) {
    clout << "rank[" << (*it).first << "]=" << (*it).second << std::endl;
  }
  clout.setMultiOutput(false);
}


}  // namespace olb
#endif
