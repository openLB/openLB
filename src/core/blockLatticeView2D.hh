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

/** \file
 * Dynamics for a generic 2D block lattice view -- generic implementation.
 */
#ifndef BLOCK_LATTICE_VIEW_2D_HH
#define BLOCK_LATTICE_VIEW_2D_HH

#include "blockLatticeView2D.h"
#include "cell.h"

namespace olb {

////////////////////// Class BlockLatticeView2D /////////////////////////

template<typename T, typename DESCRIPTOR>
BlockLatticeView2D<T,DESCRIPTOR>::BlockLatticeView2D (
  BlockLatticeStructure2D<T,DESCRIPTOR>& originalLattice_ )
  : BlockLatticeStructure2D<T,DESCRIPTOR>(originalLattice_.getNx(),originalLattice_.getNy()), originalLattice(&originalLattice_),
    x0(0), y0(0)
{ }

template<typename T, typename DESCRIPTOR>
BlockLatticeView2D<T,DESCRIPTOR>::BlockLatticeView2D (
  BlockLatticeStructure2D<T,DESCRIPTOR>& originalLattice_,
  int x0_, int x1_, int y0_, int y1_ )
  : BlockLatticeStructure2D<T,DESCRIPTOR>(x1_-x0_+1,y1_-y0_+1), originalLattice(&originalLattice_),
    x0(x0_), y0(y0_)
{
  OLB_PRECONDITION(x1_ < originalLattice->getNx());
  OLB_PRECONDITION(x0_ <= x1_);
  OLB_PRECONDITION(y1_ < originalLattice->getNy());
  OLB_PRECONDITION(y0_ <= y1_);
}

template<typename T, typename DESCRIPTOR>
BlockLatticeView2D<T,DESCRIPTOR>::~BlockLatticeView2D()
{
}


template<typename T, typename DESCRIPTOR>
BlockLatticeView2D<T,DESCRIPTOR>::BlockLatticeView2D(BlockLatticeView2D<T,DESCRIPTOR> const& rhs)
  : BlockLatticeStructure2D<T,DESCRIPTOR>(rhs._nx,rhs._ny), originalLattice(rhs.originalLattice),
    x0(rhs.x0), y0(rhs.y0)
{ }

template<typename T, typename DESCRIPTOR>
BlockLatticeView2D<T,DESCRIPTOR>& BlockLatticeView2D<T,DESCRIPTOR>::operator= (
  BlockLatticeView2D<T,DESCRIPTOR> const& rhs )
{
  BlockLatticeView2D<T,DESCRIPTOR> tmp(rhs);
  swap(tmp);
  return *this;
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::swap (
  BlockLatticeView2D<T,DESCRIPTOR>& rhs)
{
  std::swap(this->_nx, rhs._nx);
  std::swap(this->_ny, rhs._ny);
  std::swap(x0, rhs.x0);
  std::swap(y0, rhs.y0);
  std::swap(originalLattice, rhs.originalLattice);
}

template<typename T, typename DESCRIPTOR>
Cell<T,DESCRIPTOR>& BlockLatticeView2D<T,DESCRIPTOR>::get(int iX, int iY)
{
  return originalLattice->get(iX+x0, iY+y0);
}

template<typename T, typename DESCRIPTOR>
Cell<T,DESCRIPTOR>& BlockLatticeView2D<T,DESCRIPTOR>::get(int latticeR[])
{
  return get(latticeR[0], latticeR[1]);
}

template<typename T, typename DESCRIPTOR>
Cell<T,DESCRIPTOR> const& BlockLatticeView2D<T,DESCRIPTOR>::get (
  int iX, int iY ) const
{
  return originalLattice->get(iX+x0, iY+y0);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::initialize()
{
  originalLattice->initialize();
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::defineDynamics (
  int x0_, int x1_, int y0_, int y1_,
  Dynamics<T,DESCRIPTOR>* dynamics )
{
  originalLattice->defineDynamics( x0_+x0, x1_+x0,
                                   y0_+y0, y1_+y0,
                                   dynamics );

}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::defineDynamics (
  int iX, int iY, Dynamics<T,DESCRIPTOR>* dynamics )
{
  originalLattice->defineDynamics( iX+x0, iY+y0, dynamics );
}

template<typename T, typename DESCRIPTOR>
Dynamics<T,DESCRIPTOR>* BlockLatticeView2D<T,DESCRIPTOR>::getDynamics (
  int iX, int iY)
{
  return originalLattice->getDynamics( iX+x0, iY+y0 );
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::collide (
  int x0_, int x1_, int y0_, int y1_ )
{
  originalLattice->collide( x0_+x0, x1_+x0,
                            y0_+y0, y1_+y0 );
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::collide()
{
  originalLattice->collide( x0, x0+this->_nx-1, y0, y0+this->_ny-1);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::stream(int x0_, int x1_, int y0_, int y1_)
{
  originalLattice->stream(x0_+x0, x1_+x0, y0_+y0, y1_+y0);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::stream(bool periodic)
{
  OLB_PRECONDITION(!periodic);
  originalLattice->stream( x0, x0+this->_nx-1, y0, y0+this->_ny-1);
  postProcess();
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::collideAndStream(int x0_, int x1_, int y0_, int y1_)
{
  originalLattice->collideAndStream(x0_+x0, x1_+x0, y0_+y0, y1_+y0);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::collideAndStream(bool periodic)
{
  OLB_PRECONDITION(!periodic);
  originalLattice->collideAndStream( x0, x0+this->_nx-1, y0, y0+this->_ny-1);
  postProcess();
}

template<typename T, typename DESCRIPTOR>
T BlockLatticeView2D<T,DESCRIPTOR>::computeAverageDensity (
  int x0_, int x1_, int y0_, int y1_ ) const
{
  return originalLattice->computeAverageDensity( x0_+x0, x1_+x0,
         y0_+y0, y1_+y0 );
}

template<typename T, typename DESCRIPTOR>
T BlockLatticeView2D<T,DESCRIPTOR>::computeAverageDensity() const
{
  return originalLattice->computeAverageDensity( x0, x0+this->_nx-1, y0, y0+this->_ny-1);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::computeStress(int iX, int iY, T pi[util::TensorVal<DESCRIPTOR >::n])
{
    originalLattice->computeStress( iX + x0, iY + y0, pi);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::stripeOffDensityOffset (
  int x0_, int x1_, int y0_, int y1_, T offset )
{
  originalLattice->stripeOffDensityOffset( x0_+x0, x1_+x0,
      y0_+y0, y1_+y0, offset );
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::stripeOffDensityOffset(T offset)
{
  originalLattice->stripeOffDensityOffset(x0, x0+this->_nx-1, y0, y0+this->_ny-1, offset);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::forAll (
  int x0_, int x1_, int y0_, int y1_, WriteCellFunctional<T,DESCRIPTOR> const& application )
{
  originalLattice->forAll( x0_+x0, x1_+x0, y0_+y0, y1_+y0, application );
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::forAll(WriteCellFunctional<T,DESCRIPTOR> const& application)
{
  originalLattice->forAll(x0, x0+this->_nx-1, y0, y0+this->_ny-1, application);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::addPostProcessor (
  PostProcessorGenerator2D<T,DESCRIPTOR> const& ppGen)
{
  PostProcessorGenerator2D<T,DESCRIPTOR>* shiftedGen = ppGen.clone();
  shiftedGen->shift(x0,y0);
  originalLattice->addPostProcessor(*shiftedGen);
  delete shiftedGen;
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::resetPostProcessors()
{
  originalLattice->resetPostProcessors();
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::postProcess()
{
  originalLattice -> postProcess(x0, x0+this->_nx-1, y0, y0+this->_ny-1);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::postProcess (
  int x0_, int x1_, int y0_, int y1_ )
{
  originalLattice -> postProcess( x0_+x0, x1_+x0, y0_+y0, y1_+y0 );
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::addLatticeCoupling (
  LatticeCouplingGenerator2D<T,DESCRIPTOR> const& lcGen,
  std::vector<SpatiallyExtendedObject2D*> partners )
{
  LatticeCouplingGenerator2D<T,DESCRIPTOR>* shiftedGen = lcGen.clone();
  shiftedGen->shift(x0,y0);
  originalLattice->addLatticeCoupling(*shiftedGen, partners);
  delete shiftedGen;
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::executeCoupling()
{
  originalLattice -> executeCoupling(x0, x0+this->_nx-1, y0, y0+this->_ny-1);
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeView2D<T,DESCRIPTOR>::executeCoupling (
  int x0_, int x1_, int y0_, int y1_ )
{
  originalLattice -> executeCoupling( x0_+x0, x1_+x0, y0_+y0, y1_+y0 );
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T>& BlockLatticeView2D<T,DESCRIPTOR>::getStatistics()
{
  return originalLattice->getStatistics();
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T> const& BlockLatticeView2D<T,DESCRIPTOR>::getStatistics() const
{
  return originalLattice->getStatistics();
}

}  // namespace olb

#endif
