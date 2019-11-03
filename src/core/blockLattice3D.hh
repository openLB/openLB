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
 * The dynamics of a 3D block lattice -- generic implementation.
 */
#ifndef BLOCK_LATTICE_3D_HH
#define BLOCK_LATTICE_3D_HH

#include <algorithm>

#include "util.h"
#include "blockLattice3D.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "dynamics/dynamics.h"
#include "dynamics/lbHelpers.h"
#include "communication/loadBalancer.h"
#include "communication/blockLoadBalancer.h"
#include "communication/ompManager.h"

namespace olb {

////////////////////// Class BlockLattice3D /////////////////////////

/** \param nx_ lattice width (first index)
 *  \param ny_ lattice height (second index)
 *  \param nz_ lattice depth (third index)
 */
template<typename T, typename DESCRIPTOR>
BlockLattice3D<T,DESCRIPTOR>::BlockLattice3D(int nx, int ny, int nz, BlockGeometry3D<T>& geometry)
  : BlockLatticeStructure3D<T,DESCRIPTOR>(nx,ny,nz),
    geometry_(geometry)
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
  statistics->initialize();
#endif
}

/** During destruction, the memory for the lattice and the contained
 * cells is released. However, the dynamics objects pointed to by
 * the cells must be deleted manually by the user.
 */
template<typename T, typename DESCRIPTOR>
BlockLattice3D<T,DESCRIPTOR>::~BlockLattice3D()
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
void BlockLattice3D<T,DESCRIPTOR>::initialize()
{
  postProcess();
}

template<typename T, typename DESCRIPTOR>
Dynamics<T,DESCRIPTOR>* BlockLattice3D<T,DESCRIPTOR>::getDynamics (
  int iX, int iY, int iZ)
{
  return get(iX, iY, iZ).getDynamics();
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::defineDynamics (
  int iX, int iY, int iZ, Dynamics<T,DESCRIPTOR>* dynamics )
{
  OLB_PRECONDITION(iX>=0 && iX<this->_nx);
  OLB_PRECONDITION(iY>=0 && iY<this->_ny);
  OLB_PRECONDITION(iZ>=0 && iZ<this->_nz);

  get(iX, iY, iZ).defineDynamics(dynamics);
}

/** The dynamics object is not duplicated: all cells of the rectangular
 * domain point to the same dynamics.
 *
 * The dynamics object is not owned by the BlockLattice3D object, its
 * memory management must be taken care of by the user.
 */
template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::defineDynamics (
  int x0, int x1, int y0, int y1, int z0, int z1,
  Dynamics<T,DESCRIPTOR>* dynamics )
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);
  OLB_PRECONDITION(z0>=0 && z1<this->_nz);
  OLB_PRECONDITION(z1>=z0);

  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        get(iX, iY, iZ).defineDynamics(dynamics);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::defineDynamics (
  BlockIndicatorF3D<T>& indicator, Dynamics<T,DESCRIPTOR>* dynamics)
{
  for (int iX = 0; iX < this->_nx; ++iX) {
    for (int iY = 0; iY < this->_ny; ++iY) {
      for (int iZ = 0; iZ < this->_nz; ++iZ) {
        if (indicator(iX, iY, iZ)) {
          get(iX, iY, iZ).defineDynamics(dynamics);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::defineDynamics (
  BlockGeometryStructure3D<T>& blockGeometry, int material, Dynamics<T,DESCRIPTOR>* dynamics)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometry, std::vector<int>(1, material));
  defineDynamics(indicator, dynamics);
}

/**
 * This method is automatically parallelized if your compiler understands
 * OpenMP
 */
template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::collide (
  int x0, int x1, int y0, int y1, int z0, int z1)
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);
  OLB_PRECONDITION(z0>=0 && z1<this->_nz);
  OLB_PRECONDITION(z1>=z0);

  int iX;
#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for schedule(dynamic,1)
#endif
  for (iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        grid[iX][iY][iZ].collide(getStatistics());
        grid[iX][iY][iZ].revert();
      }
    }
  }
}

/** \sa collide(int,int,int,int) */
template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::collide()
{
  collide(0, this->_nx-1, 0, this->_ny-1, 0, this->_nz-1);
}

/** The distribution function never leave the rectangular domain. On the
 * domain boundaries, the (outgoing) distribution functions that should
 * be streamed outside are simply left untouched.
 */
template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::stream(int x0, int x1, int y0, int y1, int z0, int z1)
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);
  OLB_PRECONDITION(z0>=0 && z1<this->_nz);
  OLB_PRECONDITION(z1>=z0);

  static const int vicinity = descriptors::vicinity<DESCRIPTOR>();

  bulkStream(x0+vicinity,x1-vicinity, y0+vicinity,y1-vicinity, z0+vicinity,z1-vicinity);

  boundaryStream(x0,x1,y0,y1,z0,z1, x0,x0+vicinity-1, y0,y1, z0,z1);
  boundaryStream(x0,x1,y0,y1,z0,z1, x1-vicinity+1,x1, y0,y1, z0,z1);
  boundaryStream(x0,x1,y0,y1,z0,z1, x0+vicinity,x1-vicinity, y0,y0+vicinity-1, z0,z1);
  boundaryStream(x0,x1,y0,y1,z0,z1, x0+vicinity,x1-vicinity, y1-vicinity+1,y1, z0,z1);
  boundaryStream(x0,x1,y0,y1,z0,z1, x0+vicinity,x1-vicinity, y0+vicinity,y1-vicinity, z0,z0+vicinity-1);
  boundaryStream(x0,x1,y0,y1,z0,z1, x0+vicinity,x1-vicinity, y0+vicinity,y1-vicinity, z1-vicinity+1,z1);
}

/** Post-processing steps are called at the end of this method.
 * \sa stream(int,int,int,int,int,int) */
template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::stream(bool periodic)
{
  stream(0, this->_nx-1, 0, this->_ny-1, 0, this->_nz-1);

  if (periodic) {
    makePeriodic();
  }

  postProcess();
  getStatistics().incrementTime();
}

/** This operation is more efficient than a successive application of
 * collide(int,int,int,int,int,int) and stream(int,int,int,int,int,int),
 * because memory is traversed only once instead of twice.
 */
template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::collideAndStream(int x0, int x1, int y0, int y1, int z0, int z1)
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);
  OLB_PRECONDITION(z0>=0 && z1<this->_nz);
  OLB_PRECONDITION(z1>=z0);

  static const int vicinity = descriptors::vicinity<DESCRIPTOR>();

  collide(x0,x0+vicinity-1, y0,y1, z0,z1);
  collide(x1-vicinity+1,x1, y0,y1, z0,z1);
  collide(x0+vicinity,x1-vicinity, y0,y0+vicinity-1, z0,z1);
  collide(x0+vicinity,x1-vicinity, y1-vicinity+1,y1, z0,z1);
  collide(x0+vicinity,x1-vicinity, y0+vicinity,y1-vicinity, z0,z0+vicinity-1);
  collide(x0+vicinity,x1-vicinity, y0+vicinity,y1-vicinity, z1-vicinity+1,z1);

  bulkCollideAndStream(x0+vicinity,x1-vicinity, y0+vicinity,y1-vicinity, z0+vicinity,z1-vicinity);

  boundaryStream(x0,x1,y0,y1,z0,z1, x0,x0+vicinity-1, y0,y1, z0,z1);
  boundaryStream(x0,x1,y0,y1,z0,z1, x1-vicinity+1,x1, y0,y1, z0,z1);
  boundaryStream(x0,x1,y0,y1,z0,z1, x0+vicinity,x1-vicinity, y0,y0+vicinity-1, z0,z1);
  boundaryStream(x0,x1,y0,y1,z0,z1, x0+vicinity,x1-vicinity, y1-vicinity+1,y1, z0,z1);
  boundaryStream(x0,x1,y0,y1,z0,z1, x0+vicinity,x1-vicinity, y0+vicinity,y1-vicinity, z0,z0+vicinity-1);
  boundaryStream(x0,x1,y0,y1,z0,z1, x0+vicinity,x1-vicinity, y0+vicinity,y1-vicinity, z1-vicinity+1,z1);
}

/** Post-processing steps are called at the end of this method.
 * \sa collideAndStream(int,int,int,int,int,int) */
template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::collideAndStream(bool periodic)
{
  collideAndStream(0, this->_nx-1, 0, this->_ny-1, 0, this->_nz-1);

  if (periodic) {
    makePeriodic();
  }

  postProcess();
  getStatistics().incrementTime();
}

template<typename T, typename DESCRIPTOR>
T BlockLattice3D<T,DESCRIPTOR>::computeAverageDensity (
  int x0, int x1, int y0, int y1, int z0, int z1) const
{
  T sumRho = T();
  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        T rho, u[DESCRIPTOR::d];
        get(iX,iY,iZ).computeRhoU(rho, u);
        sumRho += rho;
      }
    }
  }
  return sumRho / (T)(x1-x0+1) / (T)(y1-y0+1) / (T)(z1-z0+1);
}

template<typename T, typename DESCRIPTOR>
T BlockLattice3D<T,DESCRIPTOR>::computeAverageDensity() const
{
  return computeAverageDensity(0, this->_nx-1, 0, this->_ny-1, 0, this->_nz-1);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::computeStress(int iX, int iY, int iZ,
    T pi[util::TensorVal<DESCRIPTOR >::n])
{
  grid[iX][iY][iZ].computeStress(pi);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::stripeOffDensityOffset (
  int x0, int x1, int y0, int y1, int z0, int z1, T offset )
{
  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        //if (offset<-42000.) {
        //T rho = get(iX,iY,iZ).computeRho();
        // if (rho<0) {
        //for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
        //  if (get(iX,iY,iZ)[iPop] + descriptors::t<T,DESCRIPTOR>(iPop) < T() ) {
        //    get(iX,iY,iZ)[iPop] = -descriptors::t<T,DESCRIPTOR>(iPop)+0.0000001;
        //  }
        //  else if(rho>1.)
        // get(iX,iY,iZ)[iPop] -= descriptors::t<T,DESCRIPTOR>(iPop) * (rho-1.);
        //}
        //}
        //}
        //else {
        // only stripe off if rho stays positive
        //if (get(iX,iY,iZ).computeRho()>offset) {
        for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
          get(iX,iY,iZ)[iPop] -= descriptors::t<T,DESCRIPTOR>(iPop) * offset;
        }
        // }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::stripeOffDensityOffset(T offset)
{
  stripeOffDensityOffset(0, this->_nx-1, 0, this->_ny-1, 0, this->_nz-1, offset);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::forAll (
  int x0, int x1, int y0, int y1, int z0, int z1,
  WriteCellFunctional<T,DESCRIPTOR> const& application )
{
  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        int pos[] = {iX, iY, iZ};
        application.apply( get(iX,iY,iZ), pos );
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::forAll(WriteCellFunctional<T,DESCRIPTOR> const& application)
{
  forAll(0, this->_nx-1, 0, this->_ny-1, 0, this->_nz-1, application);
}


template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::addPostProcessor (
  PostProcessorGenerator3D<T,DESCRIPTOR> const& ppGen )
{
  postProcessors.push_back(ppGen.generate());
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::resetPostProcessors()
{
  clearPostProcessors();
  StatPPGenerator3D<T,DESCRIPTOR> statPPGenerator;
  addPostProcessor(statPPGenerator);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::clearPostProcessors()
{
  typename std::vector<PostProcessor3D<T,DESCRIPTOR>*>::iterator ppIt = postProcessors.begin();
  for (; ppIt != postProcessors.end(); ++ppIt) {
    delete *ppIt;
  }
  postProcessors.clear();
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::postProcess()
{
  for (unsigned iPr=0; iPr<postProcessors.size(); ++iPr) {
    postProcessors[iPr] -> process(*this);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::postProcess (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  for (unsigned iPr=0; iPr<postProcessors.size(); ++iPr) {
    postProcessors[iPr] -> processSubDomain(*this, x0_, x1_, y0_, y1_, z0_, z1_);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::addLatticeCoupling (
  LatticeCouplingGenerator3D<T,DESCRIPTOR> const& lcGen,
  std::vector<SpatiallyExtendedObject3D*> partners )
{
  latticeCouplings.push_back(lcGen.generate(partners));
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::executeCoupling()
{
  for (unsigned iPr=0; iPr<latticeCouplings.size(); ++iPr) {
    latticeCouplings[iPr] -> process(*this);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::executeCoupling (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  for (unsigned iPr=0; iPr<latticeCouplings.size(); ++iPr) {
    latticeCouplings[iPr] -> processSubDomain(*this, x0_, x1_, y0_, y1_, z0_, z1_);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::clearLatticeCouplings()
{
  typename std::vector<PostProcessor3D<T,DESCRIPTOR>*>::iterator ppIt = latticeCouplings.begin();
  for (; ppIt != latticeCouplings.end(); ++ppIt) {
    delete *ppIt;
  }
  latticeCouplings.clear();
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T>& BlockLattice3D<T,DESCRIPTOR>::getStatistics()
{
#ifdef PARALLEL_MODE_OMP
  return *statistics[omp.get_rank()];
#else
  return *statistics;
#endif
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T> const&
BlockLattice3D<T,DESCRIPTOR>::getStatistics() const
{
#ifdef PARALLEL_MODE_OMP
  return *statistics[omp.get_rank()];
#else
  return *statistics;
#endif
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::allocateMemory()
{
  // The conversions to size_t ensure 64-bit compatibility. Note that
  //   nx, ny and nz are of type int, which might by 32-bit types, even on
  //   64-bit platforms. Therefore, nx*ny*nz may lead to a type overflow.
  rawData = new Cell<T,DESCRIPTOR> [(size_t)this->_nx*(size_t)this->_ny*(size_t)this->_nz];
  grid    = new Cell<T,DESCRIPTOR>** [(size_t)this->_nx];
  for (int iX=0; iX<this->_nx; ++iX) {
    grid[iX] = new Cell<T,DESCRIPTOR>* [(size_t)this->_ny];
    for (int iY=0; iY<this->_ny; ++iY) {
      grid[iX][iY] = rawData + (size_t)this->_nz*((size_t)iY+(size_t)this->_ny*(size_t)iX);
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::releaseMemory()
{
  delete [] rawData;
  for (int iX=0; iX<this->_nx; ++iX) {
    delete [] grid[iX];
  }
  delete [] grid;
}

/** This method is slower than bulkStream(int,int,int,int), because it must
 * be verified which distribution functions are to be kept from leaving
 * the domain.
 * \sa stream(int,int,int,int)
 * \sa stream()
 */
template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::boundaryStream (
  int lim_x0, int lim_x1, int lim_y0, int lim_y1, int lim_z0, int lim_z1,
  int x0, int x1, int y0, int y1, int z0, int z1 )
{
  OLB_PRECONDITION(lim_x0>=0 && lim_x1<this->_nx);
  OLB_PRECONDITION(lim_x1>=lim_x0);
  OLB_PRECONDITION(lim_y0>=0 && lim_y1<this->_ny);
  OLB_PRECONDITION(lim_y1>=lim_y0);
  OLB_PRECONDITION(lim_z0>=0 && lim_z1<this->_nz);
  OLB_PRECONDITION(lim_z1>=lim_z0);

  OLB_PRECONDITION(x0>=lim_x0 && x1<=lim_x1);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=lim_y0 && y1<=lim_y1);
  OLB_PRECONDITION(y1>=y0);
  OLB_PRECONDITION(z0>=lim_z0 && z1<=lim_z1);
  OLB_PRECONDITION(z1>=z0);

  int iX;

#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for
#endif
  for (iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        for (int iPop=1; iPop<=DESCRIPTOR::q/2; ++iPop) {
          int nextX = iX + descriptors::c<DESCRIPTOR>(iPop,0);
          int nextY = iY + descriptors::c<DESCRIPTOR>(iPop,1);
          int nextZ = iZ + descriptors::c<DESCRIPTOR>(iPop,2);
          if ( nextX>=lim_x0 && nextX<=lim_x1 &&
               nextY>=lim_y0 && nextY<=lim_y1 &&
               nextZ>=lim_z0 && nextZ<=lim_z1 ) {
            std::swap(grid[iX][iY][iZ][iPop+DESCRIPTOR::q/2],
                      grid[nextX][nextY][nextZ][iPop]);
          }
        }
      }
    }
  }
}

/** This method is faster than boundaryStream(int,int,int,int,int,int), but it
 * is erroneous when applied to boundary cells.
 * \sa stream(int,int,int,int,int,int)
 * \sa stream()
 */
template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::bulkStream (
  int x0, int x1, int y0, int y1, int z0, int z1 )
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);
  OLB_PRECONDITION(z0>=0 && z1<this->_nz);
  OLB_PRECONDITION(z1>=z0);

  int iX;
#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for
#endif
  for (iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        for (int iPop=1; iPop<=DESCRIPTOR::q/2; ++iPop) {
          int nextX = iX + descriptors::c<DESCRIPTOR>(iPop,0);
          int nextY = iY + descriptors::c<DESCRIPTOR>(iPop,1);
          int nextZ = iZ + descriptors::c<DESCRIPTOR>(iPop,2);
          std::swap(grid[iX][iY][iZ][iPop+DESCRIPTOR::q/2],
                    grid[nextX][nextY][nextZ][iPop]);
        }
      }
    }
  }
}

#ifndef PARALLEL_MODE_OMP // OpenMP parallel version is at the end
// of this file
/** This method is fast, but it is erroneous when applied to boundary
 * cells.
 * \sa collideAndStream(int,int,int,int,int,int)
 * \sa collideAndStream()
 */
template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::bulkCollideAndStream (
  int x0, int x1, int y0, int y1, int z0, int z1 )
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);
  OLB_PRECONDITION(z0>=0 && z1<this->_nz);
  OLB_PRECONDITION(z1>=z0);

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        grid[iX][iY][iZ].collide(getStatistics());
        lbHelpers<T,DESCRIPTOR>::swapAndStream3D(grid, iX, iY, iZ);
      }
    }
  }
}
#endif  // not defined PARALLEL_MODE_OMP

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::makePeriodic()
{
  static const int vicinity = descriptors::vicinity<DESCRIPTOR>();
  int maxX = this->getNx()-1;
  int maxY = this->getNy()-1;
  int maxZ = this->getNz()-1;
  periodicSurface(0,             vicinity-1,    0,    maxY,              0,       maxZ);
  periodicSurface(maxX-vicinity+1,     maxX,    0,    maxY,              0,       maxZ);
  periodicSurface(vicinity,      maxX-vicinity, 0,    vicinity-1,        0,       maxZ);
  periodicSurface(vicinity,      maxX-vicinity, maxY-vicinity+1, maxY,   0,       maxZ);
  periodicSurface(vicinity,      maxX-vicinity, vicinity, maxY-vicinity, 0, vicinity-1);
  periodicSurface(vicinity,      maxX-vicinity, vicinity, maxY-vicinity, maxZ-vicinity+1, maxZ);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::periodicSurface (
  int x0, int x1, int y0, int y1, int z0, int z1)
{
  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        for (int iPop=1; iPop<=DESCRIPTOR::q/2; ++iPop) {
          int nextX = iX + descriptors::c<DESCRIPTOR>(iPop,0);
          int nextY = iY + descriptors::c<DESCRIPTOR>(iPop,1);
          int nextZ = iZ + descriptors::c<DESCRIPTOR>(iPop,2);
          if ( nextX<0 || nextX>=this->getNx() ||
               nextY<0 || nextY>=this->getNy() ||
               nextZ<0 || nextZ>=this->getNz() ) {
            nextX = (nextX+this->getNx())%this->getNx();
            nextY = (nextY+this->getNy())%this->getNy();
            nextZ = (nextZ+this->getNz())%this->getNz();
            std::swap (
              grid[iX][iY][iZ]         [iPop+DESCRIPTOR::q/2],
              grid[nextX][nextY][nextZ][iPop] );
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
std::size_t BlockLattice3D<T,DESCRIPTOR>::getNblock() const
{
  return 3 + rawData[0].getNblock() * this->_nx * this->_ny * this->_nz;
}


template<typename T, typename DESCRIPTOR>
std::size_t BlockLattice3D<T,DESCRIPTOR>::getSerializableSize() const
{
  return 3 * sizeof(int) + rawData[0].getSerializableSize() * this->_nx * this->_ny * this->_nz;
}


template<typename T, typename DESCRIPTOR>
bool* BlockLattice3D<T,DESCRIPTOR>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar                      (iBlock, sizeBlock, currentBlock, dataPtr, this->_nx);
  registerVar                      (iBlock, sizeBlock, currentBlock, dataPtr, this->_ny);
  registerVar                      (iBlock, sizeBlock, currentBlock, dataPtr, this->_nz);
  registerSerializablesOfConstSize (iBlock, sizeBlock, currentBlock, dataPtr, rawData,
                                    (size_t) this->_nx * this->_ny * this->_nz, loadingMode);

  return dataPtr;
}

template<typename T, typename DESCRIPTOR>
int BlockLattice3D<T,DESCRIPTOR>::getMaterial(int iX, int iY, int iZ)
{
  return geometry_.getMaterial(iX, iY, iZ);
}

//// OpenMP implementation of the method bulkCollideAndStream,
//   by Mathias Krause                                         ////
#ifdef PARALLEL_MODE_OMP
template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::bulkCollideAndStream (
  int x0, int x1, int y0, int y1, int z0, int z1 )
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);
  OLB_PRECONDITION(z0>=0 && z1<this->_nz);
  OLB_PRECONDITION(z1>=z0);

  if (omp.get_size() <= x1-x0+1) {
    #pragma omp parallel
    {
      BlockLoadBalancer<T> loadbalance(omp.get_rank(), omp.get_size(), x1-x0+1, x0);
      int iX, iY, iZ, iPop;

      iX=loadbalance.firstGlobNum();
      for (int iY=y0; iY<=y1; ++iY)
      {
        for (int iZ=z0; iZ<=z1; ++iZ) {
          grid[iX][iY][iZ].collide(getStatistics());
          grid[iX][iY][iZ].revert();
        }
      }

      for (iX=loadbalance.firstGlobNum()+1; iX<=loadbalance.lastGlobNum(); ++iX)
      {
        for (iY=y0; iY<=y1; ++iY) {
          for (iZ=z0; iZ<=z1; ++iZ) {
            grid[iX][iY][iZ].collide(getStatistics());
            /** The method beneath doesnt work with Intel
             *  compiler 9.1044 and 9.1046 for Itanium prozessors
             *    lbHelpers<T,DESCRIPTOR>::swapAndStream3D(grid, iX, iY, iZ);
             *  Therefore we use:
             */
            int half = DESCRIPTOR::q/2;
            for (int iPop=1; iPop<=half; ++iPop) {
              int nextX = iX + descriptors::c<DESCRIPTOR>(iPop,0);
              int nextY = iY + descriptors::c<DESCRIPTOR>(iPop,1);
              int nextZ = iZ + descriptors::c<DESCRIPTOR>(iPop,2);
              T fTmp                          = grid[iX][iY][iZ][iPop];
              grid[iX][iY][iZ][iPop]          = grid[iX][iY][iZ][iPop+half];
              grid[iX][iY][iZ][iPop+half]     = grid[nextX][nextY][nextZ][iPop];
              grid[nextX][nextY][nextZ][iPop] = fTmp;
            }
          }
        }
      }

      #pragma omp barrier
      iX=loadbalance.firstGlobNum();
      for (iY=y0; iY<=y1; ++iY)
      {
        for (iZ=z0; iZ<=z1; ++iZ) {
          for (iPop=1; iPop<=DESCRIPTOR::q/2; ++iPop) {
            int nextX = iX + descriptors::c<DESCRIPTOR>(iPop,0);
            int nextY = iY + descriptors::c<DESCRIPTOR>(iPop,1);
            int nextZ = iZ + descriptors::c<DESCRIPTOR>(iPop,2);
            std::swap(grid[iX][iY][iZ][iPop+DESCRIPTOR::q/2],
                      grid[nextX][nextY][nextZ][iPop]);
          }
        }
      }
    }
  }
  else {
    for (int iX=x0; iX<=x1; ++iX) {
      for (int iY=y0; iY<=y1; ++iY) {
        for (int iZ=z0; iZ<=z1; ++iZ) {
          grid[iX][iY][iZ].collide(getStatistics());
          lbHelpers<T,DESCRIPTOR>::swapAndStream3D(grid, iX, iY, iZ);
        }
      }
    }
  }
}

#endif // defined PARALLEL_MODE_OMP

}  // namespace olb

#endif
