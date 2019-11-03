/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

#ifndef FD_BOUNDARIES_2D_HH
#define FD_BOUNDARIES_2D_HH

#include "boundaryPostProcessors2D.h"
#include "core/finiteDifference2D.h"
#include "core/blockLattice2D.h"
#include "core/util.h"
#include "dynamics/lbHelpers.h"
#include "dynamics/firstOrderLbHelpers.h"

namespace olb {

///////////  StraightFdBoundaryProcessor2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, int direction, int orientation>
StraightFdBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>::
StraightFdBoundaryProcessor2D(int x0_, int x1_, int y0_, int y1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);
}

template<typename T, typename DESCRIPTOR, int direction,int orientation>
void StraightFdBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  using namespace olb::util::tensorIndices2D;

  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    int iX;

#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      T dx_u[DESCRIPTOR::d], dy_u[DESCRIPTOR::d];
      for (int iY=newY0; iY<=newY1; ++iY) {
        Cell<T,DESCRIPTOR>& cell = blockLattice.get(iX,iY);
        Dynamics<T,DESCRIPTOR>* dynamics = blockLattice.getDynamics(iX, iY);

        T rho, u[DESCRIPTOR::d];
        cell.computeRhoU(rho,u);

        interpolateGradients<0>(blockLattice, dx_u, iX, iY);
        interpolateGradients<1>(blockLattice, dy_u, iX, iY);
        T dx_ux = dx_u[0];
        T dy_ux = dy_u[0];
        T dx_uy = dx_u[1];
        T dy_uy = dy_u[1];
        T omega = dynamics->getOmega();
        T sToPi = - rho / descriptors::invCs2<T,DESCRIPTOR>() / omega;
        T pi[util::TensorVal<DESCRIPTOR >::n];
        pi[xx] = (T)2 * dx_ux * sToPi;
        pi[yy] = (T)2 * dy_uy * sToPi;
        pi[xy] = (dx_uy + dy_ux) * sToPi;

        // Computation of the particle distribution functions
        // according to the regularized formula

        T uSqr = util::normSqr<T,2>(u);
        for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
          cell[iPop] = dynamics -> computeEquilibrium(iPop,rho,u,uSqr) +
                       firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, int direction,int orientation>
void StraightFdBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

template<typename T, typename DESCRIPTOR, int direction,int orientation>
template<int deriveDirection>
void StraightFdBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>::
interpolateGradients(BlockLattice2D<T,DESCRIPTOR> const& blockLattice,
                     T velDeriv[DESCRIPTOR::d], int iX, int iY) const
{
  fd::DirectedGradients2D<T,DESCRIPTOR,direction,orientation,direction==deriveDirection>::
  interpolateVector(velDeriv, blockLattice, iX, iY);
}

////////  StraightFdBoundaryProcessorGenerator2D ////////////////////////////////

template<typename T, typename DESCRIPTOR, int direction,int orientation>
StraightFdBoundaryProcessorGenerator2D<T,DESCRIPTOR, direction,orientation>::
StraightFdBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_)
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_)
{ }

template<typename T, typename DESCRIPTOR, int direction,int orientation>
PostProcessor2D<T,DESCRIPTOR>*
StraightFdBoundaryProcessorGenerator2D<T,DESCRIPTOR,direction,orientation>::generate() const
{
  return new StraightFdBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>
         ( this->x0, this->x1, this->y0, this->y1);
}

template<typename T, typename DESCRIPTOR, int direction,int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
StraightFdBoundaryProcessorGenerator2D<T,DESCRIPTOR,direction,orientation>::clone() const
{
  return new StraightFdBoundaryProcessorGenerator2D<T,DESCRIPTOR,direction,orientation>
         (this->x0, this->x1, this->y0, this->y1);
}


////////  StraightConvectionBoundaryProcessor2D ////////////////////////////////

template<typename T, typename DESCRIPTOR, int direction, int orientation>
StraightConvectionBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>::
StraightConvectionBoundaryProcessor2D(int x0_, int x1_, int y0_, int y1_, T* uAv_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), uAv(uAv_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);

  saveCell = new T** [(size_t)(x1_-x0_+1)];
  for (int iX=0; iX<=x1_-x0_; ++iX) {
    saveCell[iX] = new T* [(size_t)(y1_-y0_+1)];
    for (int iY=0; iY<=y1_-y0_; ++iY) {
      saveCell[iX][iY] = new T [(size_t)(DESCRIPTOR::q)];
      for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
        // default set to -1 in order to avoid wrong results at first call
        saveCell[iX][iY][iPop] = T(-1);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
StraightConvectionBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>::
~StraightConvectionBoundaryProcessor2D()
{
  for (int iX=0; iX<=x1-x0; ++iX) {
    for (int iY=0; iY<=y1-y0; ++iY) {
      delete [] saveCell[iX][iY];
    }
    delete [] saveCell[iX];
  }
  delete [] saveCell;
}

template<typename T, typename DESCRIPTOR, int direction,int orientation>
void StraightConvectionBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  using namespace olb::util::tensorIndices2D;

  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    int iX;
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        Cell<T,DESCRIPTOR>& cell = blockLattice.get(iX,iY);
        for (int iPop = 0; iPop < DESCRIPTOR::q ; ++iPop) {
          if (descriptors::c<DESCRIPTOR>(iPop,direction)==-orientation) {
            // using default -1 avoids wrong first call
            if (!util::nearZero(1 + saveCell[iX-newX0][iY-newY0][iPop]) ){
              cell[iPop] = saveCell[iX-newX0][iY-newY0][iPop];
            }
          }
        }

        T rho0, u0[2];
        T rho1, u1[2];
        T rho2, u2[2];
        if (direction==0) {
          blockLattice.get(iX,iY).computeRhoU(rho0,u0);
          blockLattice.get(iX-orientation,iY).computeRhoU(rho1,u1);
          blockLattice.get(iX-orientation*2,iY).computeRhoU(rho2,u2);
        } else {
          blockLattice.get(iX,iY).computeRhoU(rho0,u0);
          blockLattice.get(iX,iY-orientation).computeRhoU(rho1,u1);
          blockLattice.get(iX,iY-orientation*2).computeRhoU(rho2,u2);
        }

        // rho0 = T(1); rho1 = T(1); rho2 = T(1);

        T uDelta[2];
        T uAverage = rho0*u0[direction];
        if (uAv!=nullptr) {
          uAverage = *uAv;
        }
        uDelta[0]=-uAverage*0.5*(3*rho0*u0[0]-4*rho1*u1[0]+rho2*u2[0]);
        uDelta[1]=-uAverage*0.5*(3*rho0*u0[1]-4*rho1*u1[1]+rho2*u2[1]);

        for (int iPop = 0; iPop < DESCRIPTOR::q ; ++iPop) {
          if (descriptors::c<DESCRIPTOR>(iPop,direction) == -orientation) {
            saveCell[iX-newX0][iY-newY0][iPop] = cell[iPop] + descriptors::invCs2<T,DESCRIPTOR>()*descriptors::t<T,DESCRIPTOR>(iPop)*(uDelta[0]*descriptors::c<DESCRIPTOR>(iPop,0)+uDelta[1]*descriptors::c<DESCRIPTOR>(iPop,1));
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, int direction,int orientation>
void StraightConvectionBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}


////////  StraightConvectionBoundaryProcessorGenerator2D ////////////////////////////////

template<typename T, typename DESCRIPTOR, int direction,int orientation>
StraightConvectionBoundaryProcessorGenerator2D<T,DESCRIPTOR, direction,orientation>::
StraightConvectionBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_, T* uAv_)
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_), uAv(uAv_)
{ }

template<typename T, typename DESCRIPTOR, int direction,int orientation>
PostProcessor2D<T,DESCRIPTOR>*
StraightConvectionBoundaryProcessorGenerator2D<T,DESCRIPTOR,direction,orientation>::generate() const
{
  return new StraightConvectionBoundaryProcessor2D<T,DESCRIPTOR,direction,orientation>
         ( this->x0, this->x1, this->y0, this->y1, uAv);
}

template<typename T, typename DESCRIPTOR, int direction,int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
StraightConvectionBoundaryProcessorGenerator2D<T,DESCRIPTOR,direction,orientation>::clone() const
{
  return new StraightConvectionBoundaryProcessorGenerator2D<T,DESCRIPTOR,direction,orientation>
         (this->x0, this->x1, this->y0, this->y1, uAv);
}


////////  SlipBoundaryProcessor2D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
SlipBoundaryProcessor2D<T,DESCRIPTOR>::
SlipBoundaryProcessor2D(int x0_, int x1_, int y0_, int y1_, int discreteNormalX, int discreteNormalY)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);
  reflectionPop[0] = 0;
  for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
    reflectionPop[iPop] = 0;
    // iPop are the directions which ointing into the fluid, discreteNormal is pointing outwarts
    if (descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY < 0) {
      // std::cout << "-----" <<s td::endl;
      int mirrorDirection0;
      int mirrorDirection1;
      int mult = 1;
      if (discreteNormalX*discreteNormalY==0) {
        mult = 2;
      }
      mirrorDirection0 = (descriptors::c<DESCRIPTOR>(iPop,0) - mult*(descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY )*discreteNormalX);
      mirrorDirection1 = (descriptors::c<DESCRIPTOR>(iPop,1) - mult*(descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY )*discreteNormalY);

      // computes mirror jPop
      for (reflectionPop[iPop] = 1; reflectionPop[iPop] < DESCRIPTOR::q ; reflectionPop[iPop]++) {
        if (descriptors::c<DESCRIPTOR>(reflectionPop[iPop],0)==mirrorDirection0 && descriptors::c<DESCRIPTOR>(reflectionPop[iPop],1)==mirrorDirection1 ) {
          break;
        }
      }
      //std::cout <<iPop << " to "<< jPop <<" for discreteNormal= "<< discreteNormalX << "/"<<discreteNormalY <<std::endl;
    }
  }
}

template<typename T, typename DESCRIPTOR>
void SlipBoundaryProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    int iX;
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
          if (reflectionPop[iPop]!=0) {
            //do reflection
            blockLattice.get(iX,iY)[iPop] = blockLattice.get(iX,iY)[reflectionPop[iPop]];
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void SlipBoundaryProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

////////  SlipBoundaryProcessorGenerator2D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
SlipBoundaryProcessorGenerator2D<T,DESCRIPTOR>::
SlipBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_, int discreteNormalX_, int discreteNormalY_)
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_), discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>*
SlipBoundaryProcessorGenerator2D<T,DESCRIPTOR>::generate() const
{
  return new SlipBoundaryProcessor2D<T,DESCRIPTOR>
         ( this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
SlipBoundaryProcessorGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new SlipBoundaryProcessorGenerator2D<T,DESCRIPTOR>
         (this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY);
}

////////  PartialSlipBoundaryProcessor2D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
PartialSlipBoundaryProcessor2D<T,DESCRIPTOR>::
PartialSlipBoundaryProcessor2D(T tuner_, int x0_, int x1_, int y0_, int y1_, int discreteNormalX, int discreteNormalY)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), tuner(tuner_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);
  reflectionPop[0] = 0;
  for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
    reflectionPop[iPop] = 0;
    // iPop are the directions which ointing into the fluid, discreteNormal is pointing outwarts
    if (descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY < 0) {
      //std::cout << "-----" <<std::endl;
      int mirrorDirection0;
      int mirrorDirection1;
      int mult = 1;
      if (discreteNormalX*discreteNormalY==0) {
        mult = 2;
      }
      mirrorDirection0 = (descriptors::c<DESCRIPTOR>(iPop,0) - mult*(descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY )*discreteNormalX);
      mirrorDirection1 = (descriptors::c<DESCRIPTOR>(iPop,1) - mult*(descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY )*discreteNormalY);

      // computes mirror jPop
      for (reflectionPop[iPop] = 1; reflectionPop[iPop] < DESCRIPTOR::q ; reflectionPop[iPop]++) {
        if (descriptors::c<DESCRIPTOR>(reflectionPop[iPop],0)==mirrorDirection0 && descriptors::c<DESCRIPTOR>(reflectionPop[iPop],1)==mirrorDirection1 ) {
          break;
        }
      }
      //std::cout <<iPop << " to "<< jPop <<" for discreteNormal= "<< discreteNormalX << "/"<<discreteNormalY <<std::endl;
    }
  }
}

template<typename T, typename DESCRIPTOR>
void PartialSlipBoundaryProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    int iX;
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
          if (reflectionPop[iPop]!=0) {
            //do reflection
            blockLattice.get(iX,iY)[iPop] = tuner*blockLattice.get(iX,iY)[reflectionPop[iPop]];
          }
        }
        for (int iPop = 1; iPop < DESCRIPTOR::q/2 ; ++iPop) {
          T provv = blockLattice.get(iX,iY)[descriptors::opposite<DESCRIPTOR>(iPop)];
          blockLattice.get(iX,iY)[descriptors::opposite<DESCRIPTOR>(iPop)] +=
                                         (1.-tuner)*blockLattice.get(iX,iY)[iPop];
          blockLattice.get(iX,iY)[iPop] += (1.-tuner)*provv;
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void PartialSlipBoundaryProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

////////  PartialSlipBoundaryProcessorGenerator2D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
PartialSlipBoundaryProcessorGenerator2D<T,DESCRIPTOR>::
PartialSlipBoundaryProcessorGenerator2D(T tuner_, int x0_, int x1_, int y0_, int y1_, int discreteNormalX_, int discreteNormalY_)
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_), discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_), tuner(tuner_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>*
PartialSlipBoundaryProcessorGenerator2D<T,DESCRIPTOR>::generate() const
{
  return new PartialSlipBoundaryProcessor2D<T,DESCRIPTOR>
         (tuner, this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
PartialSlipBoundaryProcessorGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new PartialSlipBoundaryProcessorGenerator2D<T,DESCRIPTOR>
         (tuner, this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY);
}

/////////// OuterVelocityCornerProcessor2D /////////////////////////////////////

template<typename T, typename DESCRIPTOR, int xNormal,int yNormal>
OuterVelocityCornerProcessor2D<T, DESCRIPTOR, xNormal, yNormal>::
OuterVelocityCornerProcessor2D(int x_, int y_)
  : x(x_), y(y_)
{ }

template<typename T, typename DESCRIPTOR, int xNormal,int yNormal>
void OuterVelocityCornerProcessor2D<T, DESCRIPTOR, xNormal, yNormal>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  using namespace olb::util::tensorIndices2D;

  T rho10 = blockLattice.get(x-1*xNormal, y-0*yNormal).computeRho();
  T rho01 = blockLattice.get(x-0*xNormal, y-1*yNormal).computeRho();

  T rho20 = blockLattice.get(x-2*xNormal, y-0*yNormal).computeRho();
  T rho02 = blockLattice.get(x-0*xNormal, y-2*yNormal).computeRho();

  T rho = (T)2/(T)3*(rho01+rho10) - (T)1/(T)6*(rho02+rho20);

  T dx_u[DESCRIPTOR::d], dy_u[DESCRIPTOR::d];
  fd::DirectedGradients2D<T, DESCRIPTOR, 0, xNormal, true>::interpolateVector(dx_u, blockLattice, x,y);
  fd::DirectedGradients2D<T, DESCRIPTOR, 1, yNormal, true>::interpolateVector(dy_u, blockLattice, x,y);
  T dx_ux = dx_u[0];
  T dy_ux = dy_u[0];
  T dx_uy = dx_u[1];
  T dy_uy = dy_u[1];

  Cell<T,DESCRIPTOR>& cell = blockLattice.get(x,y);
  Dynamics<T,DESCRIPTOR>* dynamics = blockLattice.getDynamics(x, y);
  T omega = dynamics -> getOmega();

  T sToPi = - rho / descriptors::invCs2<T,DESCRIPTOR>() / omega;
  T pi[util::TensorVal<DESCRIPTOR >::n];
  pi[xx] = (T)2 * dx_ux * sToPi;
  pi[yy] = (T)2 * dy_uy * sToPi;
  pi[xy] = (dx_uy + dy_ux) * sToPi;

  // Computation of the particle distribution functions
  // according to the regularized formula
  T u[DESCRIPTOR::d];
  blockLattice.get(x,y).computeU(u);

  T uSqr = util::normSqr<T,2>(u);
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] =
      dynamics -> computeEquilibrium(iPop,rho,u,uSqr) +
      firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi);
  }
}

template<typename T, typename DESCRIPTOR, int xNormal,int yNormal>
void OuterVelocityCornerProcessor2D<T, DESCRIPTOR, xNormal, yNormal>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_ )
{
  if (util::contained(x, y, x0_, x1_, y0_, y1_)) {
    process(blockLattice);
  }
}


////////  OuterVelocityCornerProcessorGenerator2D ////////////////////////////

template<typename T, typename DESCRIPTOR, int xNormal,int yNormal>
OuterVelocityCornerProcessorGenerator2D<T, DESCRIPTOR, xNormal, yNormal>::
OuterVelocityCornerProcessorGenerator2D(int x_, int y_)
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x_, x_, y_, y_)
{ }

template<typename T, typename DESCRIPTOR, int xNormal,int yNormal>
PostProcessor2D<T,DESCRIPTOR>*
OuterVelocityCornerProcessorGenerator2D<T, DESCRIPTOR, xNormal, yNormal>::generate() const
{
  return new OuterVelocityCornerProcessor2D<T, DESCRIPTOR, xNormal, yNormal>
         ( this->x0, this->y0);
}

template<typename T, typename DESCRIPTOR, int xNormal,int yNormal>
PostProcessorGenerator2D<T,DESCRIPTOR>*
OuterVelocityCornerProcessorGenerator2D<T, DESCRIPTOR, xNormal, yNormal>::
clone() const
{
  return new OuterVelocityCornerProcessorGenerator2D<T, DESCRIPTOR, xNormal, yNormal>
         ( this->x0, this->y0);
}


template<typename T, typename DESCRIPTOR>
FreeEnergyWallProcessor2D<T,DESCRIPTOR>::
FreeEnergyWallProcessor2D(int x0_, int x1_, int y0_, int y1_, int discreteNormalX_, int discreteNormalY_, T addend_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_),
    discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_),
    addend(addend_)
{ }

template<typename T, typename DESCRIPTOR>
void FreeEnergyWallProcessor2D<T,DESCRIPTOR>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        T rhoAvg = blockLattice.get(iX-discreteNormalX, iY-discreteNormalY).computeRho();
        T rhoTmp = 0.;
        for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
            rhoTmp += blockLattice.get(iX,iY)[iPop];
        }
        T rho = rhoAvg + addend;
        rho -= rhoTmp;
        blockLattice.get(iX,iY)[0] = rho-1.;
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyWallProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
processSubDomain(blockLattice, x0, x1, y0, y1);
}

template<typename T, typename DESCRIPTOR>
FreeEnergyWallProcessorGenerator2D<T,DESCRIPTOR>::
FreeEnergyWallProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_, int discreteNormalX_,
    int discreteNormalY_, T addend_)
: PostProcessorGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_), discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_), addend(addend_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>*
FreeEnergyWallProcessorGenerator2D<T,DESCRIPTOR>::generate() const
{
return new FreeEnergyWallProcessor2D<T,DESCRIPTOR>
     (this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY, addend);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
FreeEnergyWallProcessorGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new FreeEnergyWallProcessorGenerator2D<T,DESCRIPTOR>
         (this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY, addend);
}


////////  FreeEnergyChemPotBoundaryProcessor2D ////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyChemPotBoundaryProcessor2D<T,DESCRIPTOR>::
FreeEnergyChemPotBoundaryProcessor2D(
  int x0_, int x1_, int y0_, int y1_, int discreteNormalX_, int discreteNormalY_, int latticeNumber_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), discreteNormalX(discreteNormalX_),
    discreteNormalY(discreteNormalY_), latticeNumber(latticeNumber_)
{ }

template<typename T, typename DESCRIPTOR>
void FreeEnergyChemPotBoundaryProcessor2D<T,DESCRIPTOR>::
processSubDomain(
  BlockLattice2D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        blockLattice.get(iX,iY).template setField<descriptors::CHEM_POTENTIAL>(
          blockLattice.get(iX-discreteNormalX, iY-discreteNormalY).template getField<descriptors::CHEM_POTENTIAL>()
        );
        if (latticeNumber == 1) {
          T rho0 = blockLattice.get(iX, iY).computeRho();
          T rho1 = blockLattice.get(iX-discreteNormalX, iY-discreteNormalY).computeRho();
          *(blockLattice.get(iX,iY).template getFieldPointer<descriptors::CHEM_POTENTIAL>()) += (rho1 / rho0 - 1) / descriptors::invCs2<T,DESCRIPTOR>();
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyChemPotBoundaryProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
processSubDomain(blockLattice, x0, x1, y0, y1);
}

////////  FreeEnergyChemPotBoundaryProcessorGenerator2D ////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyChemPotBoundaryProcessorGenerator2D<T,DESCRIPTOR>::
FreeEnergyChemPotBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_,
  int discreteNormalX_, int discreteNormalY_, int latticeNumber_)
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_), discreteNormalX(discreteNormalX_),
    discreteNormalY(discreteNormalY_), latticeNumber(latticeNumber_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>*
FreeEnergyChemPotBoundaryProcessorGenerator2D<T,DESCRIPTOR>::generate() const
{
return new FreeEnergyChemPotBoundaryProcessor2D<T,DESCRIPTOR>
     (this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY, latticeNumber);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
FreeEnergyChemPotBoundaryProcessorGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new FreeEnergyChemPotBoundaryProcessorGenerator2D<T,DESCRIPTOR>
         (this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY, latticeNumber);
}


////////  FreeEnergyConvectiveProcessor2D ////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyConvectiveProcessor2D<T,DESCRIPTOR>::
FreeEnergyConvectiveProcessor2D(
  int x0_, int x1_, int y0_, int y1_, int discreteNormalX_, int discreteNormalY_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_),
    discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_)
{ }

template<typename T, typename DESCRIPTOR>
void FreeEnergyConvectiveProcessor2D<T,DESCRIPTOR>::
processSubDomain(
  BlockLattice2D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        T rho, rho0, rho1, u[2];
        rho0 = blockLattice.get(iX,iY).computeRho();
        blockLattice.get(iX-discreteNormalX,iY-discreteNormalY).computeRhoU(rho1, u);
        T uPerp = 0;
        if (discreteNormalX == 0) {
          uPerp = discreteNormalY * u[1];
        } else if (discreteNormalY == 0) {
          uPerp = discreteNormalX * u[0];
        }
        rho = (rho0 + uPerp * rho1) / (1. + uPerp);
        blockLattice.get(iX,iY).defineRho(rho);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyConvectiveProcessor2D<T,DESCRIPTOR>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
processSubDomain(blockLattice, x0, x1, y0, y1);
}

////////  FreeEnergyConvectiveProcessorGenerator2D ////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyConvectiveProcessorGenerator2D<T,DESCRIPTOR>::
FreeEnergyConvectiveProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_,
  int discreteNormalX_, int discreteNormalY_)
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_),
    discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor2D<T,DESCRIPTOR>*
FreeEnergyConvectiveProcessorGenerator2D<T,DESCRIPTOR>::generate() const
{
return new FreeEnergyConvectiveProcessor2D<T,DESCRIPTOR>
     (this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>*
FreeEnergyConvectiveProcessorGenerator2D<T,DESCRIPTOR>::clone() const
{
  return new FreeEnergyConvectiveProcessorGenerator2D<T,DESCRIPTOR>
         (this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY);
}

}  // namespace olb

#endif
