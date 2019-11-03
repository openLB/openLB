/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Mathias J. Krause
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
 * Representation of the 2D block geometry view -- generic implementation.
 */

#ifndef BLOCK_GEOMETRY_VIEW_2D_HH
#define BLOCK_GEOMETRY_VIEW_2D_HH


#include <vector>
#include "geometry/blockGeometryView2D.h"


namespace olb {

template<typename T>
BlockGeometryView2D<T>::BlockGeometryView2D
(BlockGeometryStructure2D<T>& originalBlockGeometry, int x0, int x1, int y0,
 int y1)
  : BlockGeometryStructure2D<T>(originalBlockGeometry.getIcGlob()),
    BlockStructure2D(x1-x0+1, y1-y0+1),
    _originalBlockGeometry(&originalBlockGeometry), _x0(x0), _y0(y0),
    _nx(x1-x0+1), _ny(y1-y0+1)
{
  this->_statistics = BlockGeometryStatistics2D<T>(this);
  addToStatisticsList( &(this->_statistics.getStatisticsStatus()) );
}

template<typename T>
BlockGeometryView2D<T>::BlockGeometryView2D(BlockGeometryView2D const& rhs)
  : BlockGeometryStructure2D<T>(rhs),
    BlockStructure2D(0,0)
{
  _originalBlockGeometry = rhs._originalBlockGeometry;
  _x0 = rhs._x0;
  _y0 = rhs._y0;
  _nx = rhs._nx;
  _ny = rhs._ny;
  this->_iCglob = rhs._iCglob;
  this->_statistics = BlockGeometryStatistics2D<T>(this);
  addToStatisticsList( &(this->_statistics.getStatisticsStatus()) );
}

template<typename T>
BlockGeometryView2D<T>& BlockGeometryView2D<T>::operator=(BlockGeometryView2D const& rhs)
{
  _originalBlockGeometry = rhs._originalBlockGeometry;
  _x0 = rhs._x0;
  _y0 = rhs._y0;
  _nx = rhs._nx;
  _ny = rhs._ny;
  this->_iCglob = rhs._iCglob;
  this->_statistics = BlockGeometryStatistics2D<T>(this);
  addToStatisticsList( &(this->_statistics.getStatisticsStatus()) );
  return *this;
}

template<typename T>
BlockGeometryView2D<T>::~BlockGeometryView2D()
{
  removeFromStatisticsList( &(this->_statistics.getStatisticsStatus()) );
}

template<typename T>
BlockStructure2D& BlockGeometryView2D<T>::getBlockStructure()
{
  return *this;
}

template<typename T>
BlockGeometryStatistics2D<T>& BlockGeometryView2D<T>::getStatistics(bool verbose)
{
  return this->_statistics;
}

template<typename T>
BlockGeometryStatistics2D<T> const& BlockGeometryView2D<T>::getStatistics(bool verbose) const
{
  return this->_statistics;
}

template<typename T>
Vector<T,2> BlockGeometryView2D<T>::getOrigin() const
{
  Vector<T,2> origin;
  origin[0] = _originalBlockGeometry->getOrigin()[0] + _x0*getDeltaR();
  origin[1] = _originalBlockGeometry->getOrigin()[1] + _y0*getDeltaR();
  return origin;
}

template<typename T>
const T BlockGeometryView2D<T>::getDeltaR() const
{
  return _originalBlockGeometry->getDeltaR();
}

template<typename T>
int BlockGeometryView2D<T>::getNx() const
{
  return _nx;
}

template<typename T>
int BlockGeometryView2D<T>::getNy() const
{
  return _ny;
}


template<typename T>
int& BlockGeometryView2D<T>::get(int iX, int iY)
{
  return _originalBlockGeometry->get(_x0+iX, _y0+iY);
}

template<typename T>
const int& BlockGeometryView2D<T>::get(int iX, int iY) const
{
  return _originalBlockGeometry->get(_x0+iX, _y0+iY);
}

template<typename T>
int BlockGeometryView2D<T>::getMaterial(int iX, int iY) const
{
  return _originalBlockGeometry->getMaterial(_x0+iX, _y0+iY);
}

template<typename T>
void BlockGeometryView2D<T>::getPhysR(T physR[2], const int& iX, const int& iY) const
{
  _originalBlockGeometry->getPhysR(physR, _x0 + iX, _y0 + iY);
  return;
}

template<typename T>
void BlockGeometryView2D<T>::addToStatisticsList(bool* statisticStatus)
{
  _originalBlockGeometry->addToStatisticsList(statisticStatus);
}

template<typename T>
void BlockGeometryView2D<T>::removeFromStatisticsList(bool* statisticStatus)
{
  _originalBlockGeometry->removeFromStatisticsList(statisticStatus);
}

} // namespace olb

#endif
