/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007, 2008 Jonas Latt
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

#ifndef POST_PROCESSING_HH
#define POST_PROCESSING_HH

#include "blockLattice2D.h"
#include "blockLattice3D.h"
//#include <cmath>
//#include <numeric>
//#include <limits>
//#include "util.h"

namespace olb {

////////////////////// Class PostProcessorGenerator2D /////////////////

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>::PostProcessorGenerator2D (
  int x0_, int x1_, int y0_, int y1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_)
{ }

template<typename T, typename DESCRIPTOR>
void PostProcessorGenerator2D<T,DESCRIPTOR>::shift(int deltaX, int deltaY)
{
  x0 += deltaX;
  x1 += deltaX;
  y0 += deltaY;
  y1 += deltaY;
}

template<typename T, typename DESCRIPTOR>
bool PostProcessorGenerator2D<T,DESCRIPTOR>::
extract(int x0_, int x1_, int y0_, int y1_)
{
  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {
    x0 = newX0;
    x1 = newX1;
    y0 = newY0;
    y1 = newY1;
    return true;
  } else {
    return false;
  }
}


////////////////////// Class LatticeCouplingGenerator2D /////////////////

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>::LatticeCouplingGenerator2D (
  int x0_, int x1_, int y0_, int y1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_)
{ }

template<typename T, typename DESCRIPTOR>
void LatticeCouplingGenerator2D<T,DESCRIPTOR>::shift(int deltaX, int deltaY)
{
  x0 += deltaX;
  x1 += deltaX;
  y0 += deltaY;
  y1 += deltaY;
}

template<typename T, typename DESCRIPTOR>
bool LatticeCouplingGenerator2D<T,DESCRIPTOR>::extract(int x0_, int x1_, int y0_, int y1_)
{
  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {
    x0 = newX0;
    x1 = newX1;
    y0 = newY0;
    y1 = newY1;
    return true;
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
void LatticeCouplingGenerator2D<T,DESCRIPTOR>::reset(int x0_, int x1_, int y0_, int y1_)
{
  x0 = x0_;
  x1 = x1_;
  y0 = y0_;
  y1 = y1_;
}


////////////////////// Class StatisticsPostProcessor2D //////////////

template<typename T, typename DESCRIPTOR>
StatisticsPostProcessor2D<T,DESCRIPTOR>::StatisticsPostProcessor2D()
{ }

#ifndef PARALLEL_MODE_OMP
template<typename T, typename DESCRIPTOR>
void StatisticsPostProcessor2D<T,DESCRIPTOR>::process (
  BlockLattice2D<T,DESCRIPTOR>& blockLattice )
{
  blockLattice.getStatistics().reset();
}
#endif


#ifdef PARALLEL_MODE_OMP
template<typename T, typename DESCRIPTOR>
void StatisticsPostProcessor2D<T,DESCRIPTOR>::process (
  BlockLattice2D<T,DESCRIPTOR>& blockLattice )
{
  #pragma omp parallel
  blockLattice.getStatistics().reset();


  int numCells     = 0;
  T avRho    = T();
  T avEnergy = T();
  T maxU     = T();

  #pragma omp parallel
  {
    #pragma omp critical
    {
      numCells       += blockLattice.getStatistics().getNumCells();
      avRho          += blockLattice.getStatistics().getAverageRho()
      *blockLattice.getStatistics().getNumCells();
      avEnergy       += blockLattice.getStatistics().getAverageEnergy()
      *blockLattice.getStatistics().getNumCells();
      if (maxU<blockLattice.getStatistics().getMaxU() )
      {
        maxU        = blockLattice.getStatistics().getMaxU();
      }
    }
  }
  if (numCells==0) {
    // avoid division by zero
    avRho = T();
    avEnergy = T();
    maxU = T();
    numCells = 0;
  } else {
    avRho    = avRho / numCells;
    avEnergy = avEnergy / numCells;
  }
  #pragma omp parallel
  blockLattice.getStatistics().reset(avRho,avEnergy, maxU, numCells);
}
#endif


////////////////////// Class StatPPGenerator2D //////////////

template<typename T, typename DESCRIPTOR>
StatPPGenerator2D<T,DESCRIPTOR>::StatPPGenerator2D()
  : PostProcessorGenerator2D<T,DESCRIPTOR>(-1,-1,-1,-1)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>* StatPPGenerator2D<T,DESCRIPTOR>::generate() const
{
  return new StatisticsPostProcessor2D<T,DESCRIPTOR>;
}


template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
StatPPGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new StatPPGenerator2D;
}


////////////////////// Class PostProcessorGenerator3D /////////////////

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>::PostProcessorGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{ }

template<typename T, typename DESCRIPTOR>
void PostProcessorGenerator3D<T,DESCRIPTOR>::shift (
  int deltaX, int deltaY, int deltaZ )
{
  x0 += deltaX;
  x1 += deltaX;
  y0 += deltaY;
  y1 += deltaY;
  z0 += deltaZ;
  z1 += deltaZ;
}

template<typename T, typename DESCRIPTOR>
bool PostProcessorGenerator3D<T,DESCRIPTOR>::
extract(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {
    x0 = newX0;
    x1 = newX1;
    y0 = newY0;
    y1 = newY1;
    z0 = newZ0;
    z1 = newZ1;
    return true;
  } else {
    return false;
  }
}

////////////////////// Class LatticeCouplingGenerator3D /////////////////

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator3D<T,DESCRIPTOR>::LatticeCouplingGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{ }

template<typename T, typename DESCRIPTOR>
void LatticeCouplingGenerator3D<T,DESCRIPTOR>::shift (
  int deltaX, int deltaY, int deltaZ, int iC_)
{
  x0 += deltaX;
  x1 += deltaX;
  y0 += deltaY;
  y1 += deltaY;
  z0 += deltaZ;
  z1 += deltaZ;
  iC = iC_;
}

template<typename T, typename DESCRIPTOR>
bool LatticeCouplingGenerator3D<T,DESCRIPTOR>::
extract(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {
    x0 = newX0;
    x1 = newX1;
    y0 = newY0;
    y1 = newY1;
    z0 = newZ0;
    z1 = newZ1;
    return true;
  } else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
void LatticeCouplingGenerator3D<T,DESCRIPTOR>::
reset(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  x0 = x0_;
  x1 = x1_;
  y0 = y0_;
  y1 = y1_;
  z0 = z0_;
  z1 = z1_;
}

////////////////////// Class StatisticsPostProcessor3D //////////////

template<typename T, typename DESCRIPTOR>
StatisticsPostProcessor3D<T,DESCRIPTOR>::StatisticsPostProcessor3D()
{ }

#ifndef PARALLEL_MODE_OMP
template<typename T, typename DESCRIPTOR>
void StatisticsPostProcessor3D<T,DESCRIPTOR>::process (
  BlockLattice3D<T,DESCRIPTOR>& blockLattice )
{
  blockLattice.getStatistics().reset();
}
#endif
#ifdef PARALLEL_MODE_OMP
template<typename T, typename DESCRIPTOR>
void StatisticsPostProcessor3D<T,DESCRIPTOR>::process (
  BlockLattice3D<T,DESCRIPTOR>& blockLattice )
{
  #pragma omp parallel
  blockLattice.getStatistics().reset();


  int numCells     = 0;
  T avRho    = T();
  T avEnergy = T();
  T maxU     = T();

  #pragma omp parallel
  {
    #pragma omp critical
    {
      numCells       += blockLattice.getStatistics().getNumCells();
      avRho          += blockLattice.getStatistics().getAverageRho()
      *blockLattice.getStatistics().getNumCells();
      avEnergy       += blockLattice.getStatistics().getAverageEnergy()
      *blockLattice.getStatistics().getNumCells();
      if (maxU<blockLattice.getStatistics().getMaxU() )
      {
        maxU        = blockLattice.getStatistics().getMaxU();
      }
    }
  }
  if (numCells==0) {
    // avoid division by zero
    avRho = T();
    avEnergy = T();
    maxU = T();
    numCells = 0;
  } else {
    avRho    = avRho / numCells;
    avEnergy = avEnergy / numCells;
  }
  #pragma omp parallel
  blockLattice.getStatistics().reset(avRho,avEnergy, maxU, numCells);
}
#endif

////////////////////// Class StatPPGenerator3D //////////////

template<typename T, typename DESCRIPTOR>
StatPPGenerator3D<T,DESCRIPTOR>::StatPPGenerator3D()
  : PostProcessorGenerator3D<T,DESCRIPTOR>(-1,-1,-1,-1,-1,-1)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* StatPPGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new StatisticsPostProcessor3D<T,DESCRIPTOR>;
}


template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>* StatPPGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new StatPPGenerator3D;
}




}  // namespace olb

#endif
