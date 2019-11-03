/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2008 Jonas Latt
 *  OMP parallel code by Mathias Krause, Copyright (C) 2007
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
 * The dynamics of a 2D block lattice -- generic implementation.
 */
#ifndef BLOCK_LATTICE_2D_HH
#define BLOCK_LATTICE_2D_HH

#include <algorithm>
#include "blockLattice2D.h"
#include "dynamics/dynamics.h"
#include "dynamics/lbHelpers.h"
#include "util.h"
#include "communication/loadBalancer.h"
#include "communication/blockLoadBalancer.h"
#include "functors/lattice/indicator/blockIndicatorF2D.h"
#include "communication/ompManager.h"

namespace olb {

////////////////////// Class BlockLattice2D /////////////////////////

/** \param nx lattice width (first index)
 *  \param ny lattice height (second index)
 */
template<typename T, typename DESCRIPTOR>
BlockLattice2D<T,DESCRIPTOR>::BlockLattice2D(int nx, int ny)
  : BlockLatticeStructure2D<T,DESCRIPTOR>(nx,ny)
{
  allocateMemory();
  resetPostProcessors();
#ifdef PARALLEL_MODE_OMP
  statistics = new LatticeStatistics<T>* [3*omp.get_size()];
  #pragma omp parallel
  {
    statistics[omp.get_rank() + omp.get_size()]
      = new LatticeStatistics<T>;
    statistics[omp.get_rank()] = new LatticeStatistics<T>;
    statistics[omp.get_rank() + 2*omp.get_size()]
      = new LatticeStatistics<T>;
  }
#else
  statistics = new LatticeStatistics<T>;
#endif
}

/** During destruction, the memory for the lattice and the contained
 * cells is released. However, the dynamics objects pointed to by
 * the cells must be deleted manually by the user.
 */
template<typename T, typename DESCRIPTOR>
BlockLattice2D<T,DESCRIPTOR>::~BlockLattice2D()
{
  releaseMemory();
  clearPostProcessors();
  clearLatticeCouplings();
#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel
  {
    delete statistics[omp.get_rank()];
  }
  delete statistics;
#else
  delete statistics;
#endif
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::initialize()
{
  postProcess();
}

template<typename T, typename DESCRIPTOR>
Dynamics<T,DESCRIPTOR>* BlockLattice2D<T,DESCRIPTOR>::getDynamics (int iX, int iY)
{
  return grid[iX][iY].getDynamics();
}

/** The dynamics object is not duplicated: all cells of the rectangular
 * domain point to the same dynamics.
 *
 * The dynamics object is not owned by the BlockLattice2D object, its
 * memory management is in charge of the user.
 */
template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::defineDynamics (
  int x0, int x1, int y0, int y1, Dynamics<T,DESCRIPTOR>* dynamics )
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      grid[iX][iY].defineDynamics(dynamics);
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::defineDynamics (
  int iX, int iY, Dynamics<T,DESCRIPTOR>* dynamics )
{
  OLB_PRECONDITION(iX>=0 && iX<this->_nx);
  OLB_PRECONDITION(iY>=0 && iY<this->_ny);

  grid[iX][iY].defineDynamics(dynamics);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::defineDynamics (
  BlockIndicatorF2D<T>& indicator, Dynamics<T,DESCRIPTOR>* dynamics)
{
  int latticeR[2];
  for (latticeR[0] = 0; latticeR[0] < this->_nx; ++latticeR[0]) {
    for (latticeR[1] = 0; latticeR[1] < this->_ny; ++latticeR[1]) {
      if (indicator(latticeR)) {
        get(latticeR).defineDynamics(dynamics);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::defineDynamics(
  BlockGeometryStructure2D<T>& blockGeometry, int material, Dynamics<T,DESCRIPTOR>* dynamics)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometry, std::vector<int>(1, material));
  defineDynamics(indicator, dynamics);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::collide(int x0, int x1, int y0, int y1)
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);

  int iX;
#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for schedule(dynamic,1)
#endif
  for (iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      grid[iX][iY].collide(getStatistics());
      grid[iX][iY].revert();
    }
  }
}

/** \sa collide(int,int,int,int) */
template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::collide()
{
  collide(0, this->_nx-1, 0, this->_ny-1);
}

/** The distribution functions never leave the rectangular domain. On the
 * domain boundaries, the (outgoing) distribution functions that should
 * be streamed outside are simply left untouched.
 * The post-processing steps are not automatically invoked by this method,
 * as they are in the method stream(). If you want them to be executed, you
 * must explicitly call the method postProcess().
 * \sa stream()
 */
template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::stream(int x0, int x1, int y0, int y1)
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);

  static const int vicinity = descriptors::vicinity<DESCRIPTOR>();

  bulkStream(x0+vicinity,x1-vicinity,y0+vicinity,y1-vicinity);

  boundaryStream(x0,x1,y0,y1, x0,x0+vicinity-1, y0,y1);
  boundaryStream(x0,x1,y0,y1, x1-vicinity+1,x1, y0,y1);
  boundaryStream(x0,x1,y0,y1, x0+vicinity,x1-vicinity, y0,y0+vicinity-1);
  boundaryStream(x0,x1,y0,y1, x0+vicinity,x1-vicinity, y1-vicinity+1,y1);
}

/** At the end of this method, the post-processing steps are automatically
 * invoked.
 * \sa stream(int,int,int,int)
 */
template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::stream(bool periodic)
{
  stream(0, this->_nx-1, 0, this->_ny-1);

  if (periodic) {
    makePeriodic();
  }

  postProcess();
  getStatistics().incrementTime();
}

/** This operation is more efficient than a successive application of
 * collide(int,int,int,int) and stream(int,int,int,int), because memory
 * is traversed only once instead of twice.
 * The post-processing steps are not automatically invoked by this method,
 * as they are in the method stream(). If you want them to be executed, you
 * must explicitly call the method postProcess().
 * \sa collideAndStream()
 */
template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::collideAndStream(int x0, int x1, int y0, int y1)
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);

  static const int vicinity = descriptors::vicinity<DESCRIPTOR>();

  // First, do the collision on cells within a boundary envelope of width
  // equal to the range of the lattice vectors (e.g. 1 for D2Q9)
  collide(x0,x0+vicinity-1, y0,y1);
  collide(x1-vicinity+1,x1, y0,y1);
  collide(x0+vicinity,x1-vicinity, y0,y0+vicinity-1);
  collide(x0+vicinity,x1-vicinity, y1-vicinity+1,y1);

  // Then, do the efficient collideAndStream algorithm in the bulk,
  // excluding the envelope (this is efficient because there is no
  // if-then-else statement within the loop, given that the boundary
  // region is excluded)
  bulkCollideAndStream(x0+vicinity,x1-vicinity,y0+vicinity,y1-vicinity);

  // Finally, do streaming in the boundary envelope to conclude the
  // collision-stream cycle
  boundaryStream(x0,x1,y0,y1, x0,x0+vicinity-1,y0,y1);
  boundaryStream(x0,x1,y0,y1, x1-vicinity+1,x1,y0,y1);
  boundaryStream(x0,x1,y0,y1, x0+vicinity,x1-vicinity, y0,y0+vicinity-1);
  boundaryStream(x0,x1,y0,y1, x0+vicinity,x1-vicinity, y1-vicinity+1,y1);
}

/** At the end of this method, the post-processing steps are automatically
 * invoked.
 * \sa collideAndStream(int,int,int,int) */
template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::collideAndStream(bool periodic)
{
  collideAndStream(0, this->_nx-1, 0, this->_ny-1);

  if (periodic) {
    makePeriodic();
  }

  postProcess();
  getStatistics().incrementTime();
}

template<typename T, typename DESCRIPTOR>
T BlockLattice2D<T,DESCRIPTOR>::computeAverageDensity ( int x0, int x1, int y0, int y1) const
{
  T sumRho = T();
  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      T rho, u[DESCRIPTOR::d];
      get(iX,iY).computeRhoU(rho, u);
      sumRho += rho;
    }
  }
  return sumRho / (T)(x1-x0+1) / (T)(y1-y0+1);
}

template<typename T, typename DESCRIPTOR>
T BlockLattice2D<T,DESCRIPTOR>::computeAverageDensity() const
{
  return computeAverageDensity(0, this->_nx-1, 0, this->_ny-1);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::computeStress(int iX, int iY, T pi[util::TensorVal<DESCRIPTOR>::n])
{
  grid[iX][iY].computeStress(pi);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::stripeOffDensityOffset ( int x0, int x1, int y0, int y1, T offset )
{
  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
        get(iX,iY)[iPop] -= descriptors::t<T,DESCRIPTOR>(iPop) * offset;
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::stripeOffDensityOffset(T offset)
{
  stripeOffDensityOffset(0, this->_nx-1, 0, this->_ny-1, offset);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::forAll (
  int x0, int x1, int y0, int y1, WriteCellFunctional<T,DESCRIPTOR> const& application )
{
  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      int pos[] = {iX, iY};
      application.apply( get(iX,iY), pos );
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::forAll(WriteCellFunctional<T,DESCRIPTOR> const& application)
{
  forAll(0, this->_nx-1, 0, this->_ny-1, application);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::addPostProcessor (
  PostProcessorGenerator2D<T,DESCRIPTOR> const& ppGen )
{
  postProcessors.push_back(ppGen.generate());
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::resetPostProcessors()
{
  clearPostProcessors();
  StatPPGenerator2D<T,DESCRIPTOR> statPPGenerator;
  addPostProcessor(statPPGenerator);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::clearPostProcessors()
{
  typename PostProcVector::iterator ppIt = postProcessors.begin();
  for (; ppIt != postProcessors.end(); ++ppIt) {
    delete *ppIt;
  }
  postProcessors.clear();
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::postProcess()
{
  for (unsigned iPr=0; iPr<postProcessors.size(); ++iPr) {
    postProcessors[iPr] -> process(*this);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::postProcess(int x0_, int x1_, int y0_, int y1_)
{
  for (unsigned iPr=0; iPr<postProcessors.size(); ++iPr) {
    postProcessors[iPr] -> processSubDomain(*this, x0_, x1_, y0_, y1_);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::addLatticeCoupling (
  LatticeCouplingGenerator2D<T,DESCRIPTOR> const& lcGen,
  std::vector<SpatiallyExtendedObject2D*> partners )
{
  latticeCouplings.push_back(lcGen.generate(partners));
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::executeCoupling()
{
  for (unsigned iPr=0; iPr<latticeCouplings.size(); ++iPr) {
    latticeCouplings[iPr] -> process(*this);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::executeCoupling(int x0_, int x1_, int y0_, int y1_)
{
  for (unsigned iPr=0; iPr<latticeCouplings.size(); ++iPr) {
    latticeCouplings[iPr] -> processSubDomain(*this, x0_, x1_, y0_, y1_);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::clearLatticeCouplings()
{
  typename PostProcVector::iterator ppIt = latticeCouplings.begin();
  for (; ppIt != latticeCouplings.end(); ++ppIt) {
    delete *ppIt;
  }
  latticeCouplings.clear();
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T>& BlockLattice2D<T,DESCRIPTOR>::getStatistics()
{
#ifdef PARALLEL_MODE_OMP
  return *statistics[omp.get_rank()];
#else
  return *statistics;
#endif
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T> const& BlockLattice2D<T,DESCRIPTOR>::getStatistics() const
{
#ifdef PARALLEL_MODE_OMP
  return *statistics[omp.get_rank()];
#else
  return *statistics;
#endif
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::allocateMemory()
{
  // The conversions to size_t ensure 64-bit compatibility. Note that
  //   nx and ny are of type int, which might by 32-bit types, even on
  //   64-bit platforms. Therefore, nx*ny may lead to a type overflow.
  rawData = new Cell<T,DESCRIPTOR> [(size_t)(this->_nx)*(size_t)(this->_ny)];
  grid    = new Cell<T,DESCRIPTOR>* [(size_t)(this->_nx)];
  for (int iX=0; iX<this->_nx; ++iX) {
    grid[iX] = rawData + (size_t)iX*(size_t)(this->_ny);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::releaseMemory()
{
  delete [] rawData;
  delete [] grid;
}

/** This method is slower than bulkStream(int,int,int,int), because it must
 * be verified which distribution functions are to be kept from leaving
 * the domain.
 * \sa stream(int,int,int,int)
 * \sa stream()
 */
template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::boundaryStream (
  int lim_x0, int lim_x1, int lim_y0, int lim_y1,
  int x0, int x1, int y0, int y1 )
{
  OLB_PRECONDITION(lim_x0>=0 && lim_x1<this->_nx);
  OLB_PRECONDITION(lim_x1>=lim_x0);
  OLB_PRECONDITION(lim_y0>=0 && lim_y1<this->_ny);
  OLB_PRECONDITION(lim_y1>=lim_y0);

  OLB_PRECONDITION(x0>=lim_x0 && x1<=lim_x1);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=lim_y0 && y1<=lim_y1);
  OLB_PRECONDITION(y1>=y0);

  int iX;

#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for
#endif
  for (iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iPop=1; iPop<=DESCRIPTOR::q/2; ++iPop) {
        int nextX = iX + descriptors::c<DESCRIPTOR>(iPop,0);
        int nextY = iY + descriptors::c<DESCRIPTOR>(iPop,1);
        if (nextX>=lim_x0 && nextX<=lim_x1 && nextY>=lim_y0 && nextY<=lim_y1) {
          std::swap(grid[iX][iY][iPop+DESCRIPTOR::q/2],
                    grid[nextX][nextY][iPop]);
        }
      }
    }
  }
}

/** This method is faster than boundaryStream(int,int,int,int), but it
 * is erroneous when applied to boundary cells.
 * \sa stream(int,int,int,int)
 * \sa stream()
 */
template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::bulkStream (
  int x0, int x1, int y0, int y1 )
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);

  int iX;
#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for
#endif
  for (iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iPop=1; iPop<=DESCRIPTOR::q/2; ++iPop) {
        int nextX = iX + descriptors::c<DESCRIPTOR>(iPop,0);
        int nextY = iY + descriptors::c<DESCRIPTOR>(iPop,1);
        std::swap(grid[iX][iY][iPop+DESCRIPTOR::q/2],
                  grid[nextX][nextY][iPop]);
      }
    }
  }
}

#ifndef PARALLEL_MODE_OMP  // OpenMP parallel version is at the
// end of this file
/** This method is fast, but it is erroneous when applied to boundary
 * cells.
 * \sa collideAndStream(int,int,int,int)
 * \sa collideAndStream()
 */
template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::bulkCollideAndStream (
  int x0, int x1, int y0, int y1 )
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      grid[iX][iY].collide(getStatistics());
      lbHelpers<T,DESCRIPTOR>::swapAndStream2D(grid, iX, iY);
    }
  }
}
#endif // not defined PARALLEL_MODE_OMP

template<typename T, typename DESCRIPTOR>
std::size_t BlockLattice2D<T,DESCRIPTOR>::getNblock() const
{
  return 2 + rawData[0].getNblock() * this->_nx * this->_ny;
}

template<typename T, typename DESCRIPTOR>
std::size_t BlockLattice2D<T,DESCRIPTOR>::getSerializableSize() const
{
  return 2 * sizeof(int) + rawData[0].getSerializableSize() * this->_nx * this->_ny;
}

template<typename T, typename DESCRIPTOR>
bool* BlockLattice2D<T,DESCRIPTOR>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar                      (iBlock, sizeBlock, currentBlock, dataPtr, this->_nx);
  registerVar                      (iBlock, sizeBlock, currentBlock, dataPtr, this->_ny);
  registerSerializablesOfConstSize (iBlock, sizeBlock, currentBlock, dataPtr, rawData,
                                    (size_t) this->_nx * this->_ny, loadingMode);
  return dataPtr;
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::periodicEdge(int x0, int x1, int y0, int y1)
{
  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iPop=1; iPop<=DESCRIPTOR::q/2; ++iPop) {
        int nextX = iX + descriptors::c<DESCRIPTOR>(iPop,0);
        int nextY = iY + descriptors::c<DESCRIPTOR>(iPop,1);
        if ( nextX<0 || nextX>=this->_nx ||
             nextY<0 || nextY>=this->_ny ) {
          nextX = (nextX+this->_nx)%this->_nx;
          nextY = (nextY+this->_ny)%this->_ny;
          std::swap (
            grid[iX][iY][iPop+DESCRIPTOR::q/2],
            grid[nextX][nextY][iPop] );
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::makePeriodic()
{
  static const int vicinity = descriptors::vicinity<DESCRIPTOR>();
  int maxX = this->_nx-1;
  int maxY = this->_ny-1;
  periodicEdge(0,vicinity-1, 0,maxY);
  periodicEdge(maxX-vicinity+1,maxX, 0,maxY);
  periodicEdge(vicinity,maxX-vicinity, 0,vicinity-1);
  periodicEdge(vicinity,maxX-vicinity, maxY-vicinity+1,maxY);
}

//// OpenMP implementation of the method bulkCollideAndStream,
//   by Mathias Krause                                         ////

#ifdef PARALLEL_MODE_OMP
template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::bulkCollideAndStream (
  int x0, int x1, int y0, int y1 )
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);

  if (omp.get_size() <= x1-x0+1) {
    #pragma omp parallel
    {
      BlockLoadBalancer<T> loadbalance(omp.get_rank(), omp.get_size(), x1-x0+1, x0);
      int iX, iY, iPop;

      iX=loadbalance.firstGlobNum();
      for (int iY=y0; iY<=y1; ++iY)
      {
        grid[iX][iY].collide(getStatistics());
        grid[iX][iY].revert();
      }

      for (iX=loadbalance.firstGlobNum()+1; iX<=loadbalance.lastGlobNum(); ++iX)
      {
        for (iY=y0; iY<=y1; ++iY) {
          grid[iX][iY].collide(getStatistics());
          /** The method beneath doesnt work with Intel compiler 9.1044 and 9.1046 for Itanium prozessors
           *    lbHelpers<T,DESCRIPTOR>::swapAndStream2D(grid, iX, iY);
           *  Therefore we use:
           */
          int half = DESCRIPTOR::q/2;
          for (int iPop=1; iPop<=half; ++iPop) {
            int nextX = iX + descriptors::c<DESCRIPTOR>(iPop,0);
            int nextY = iY + descriptors::c<DESCRIPTOR>(iPop,1);
            T fTmp                   = grid[iX][iY][iPop];
            grid[iX][iY][iPop]       = grid[iX][iY][iPop+half];
            grid[iX][iY][iPop+half]  = grid[nextX][nextY][iPop];
            grid[nextX][nextY][iPop] = fTmp;
          }
        }
      }

      #pragma omp barrier

      iX=loadbalance.firstGlobNum();
      for (iY=y0; iY<=y1; ++iY)
      {
        for (iPop=1; iPop<=DESCRIPTOR::q/2; ++iPop) {
          int nextX = iX + descriptors::c<DESCRIPTOR>(iPop,0);
          int nextY = iY + descriptors::c<DESCRIPTOR>(iPop,1);
          std::swap(grid[iX][iY][iPop+DESCRIPTOR::q/2],
                    grid[nextX][nextY][iPop]);
        }
      }
    }
  }
  else {
    for (int iX=x0; iX<=x1; ++iX) {
      for (int iY=y0; iY<=y1; ++iY) {
        grid[iX][iY].collide(getStatistics());
        lbHelpers<T,DESCRIPTOR>::swapAndStream2D(grid, iX, iY);
      }
    }
  }
}
#endif // defined PARALLEL_MODE_OMP

}  // namespace olb

#endif
