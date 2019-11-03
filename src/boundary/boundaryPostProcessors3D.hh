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

#ifndef BOUNDARY_POST_PROCESSORS_3D_HH
#define BOUNDARY_POST_PROCESSORS_3D_HH

#include "boundaryPostProcessors3D.h"
#include "core/finiteDifference3D.h"
#include "core/blockLattice3D.h"
#include "dynamics/firstOrderLbHelpers.h"
#include "core/util.h"
#include "utilities/vectorHelpers.h"

namespace olb {

////////  PlaneFdBoundaryProcessor3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, int direction, int orientation>
PlaneFdBoundaryProcessor3D<T,DESCRIPTOR,direction,orientation>::
PlaneFdBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void PlaneFdBoundaryProcessor3D<T,DESCRIPTOR,direction,orientation>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  using namespace olb::util::tensorIndices3D;

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {
    int iX;

#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      T dx_u[DESCRIPTOR::d], dy_u[DESCRIPTOR::d], dz_u[DESCRIPTOR::d];
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          Cell<T,DESCRIPTOR>& cell = blockLattice.get(iX,iY,iZ);
          Dynamics<T,DESCRIPTOR>* dynamics = blockLattice.getDynamics(iX, iY, iZ);
          T rho, u[DESCRIPTOR::d];
          cell.computeRhoU(rho,u);

          interpolateGradients<0> ( blockLattice, dx_u, iX, iY, iZ );
          interpolateGradients<1> ( blockLattice, dy_u, iX, iY, iZ );
          interpolateGradients<2> ( blockLattice, dz_u, iX, iY, iZ );
          T dx_ux = dx_u[0];
          T dy_ux = dy_u[0];
          T dz_ux = dz_u[0];
          T dx_uy = dx_u[1];
          T dy_uy = dy_u[1];
          T dz_uy = dz_u[1];
          T dx_uz = dx_u[2];
          T dy_uz = dy_u[2];
          T dz_uz = dz_u[2];
          T omega = dynamics->getOmega();
          T sToPi = - rho / descriptors::invCs2<T,DESCRIPTOR>() / omega;
          T pi[util::TensorVal<DESCRIPTOR >::n];
          pi[xx] = (T)2 * dx_ux * sToPi;
          pi[yy] = (T)2 * dy_uy * sToPi;
          pi[zz] = (T)2 * dz_uz * sToPi;
          pi[xy] = (dx_uy + dy_ux) * sToPi;
          pi[xz] = (dx_uz + dz_ux) * sToPi;
          pi[yz] = (dy_uz + dz_uy) * sToPi;

          // Computation of the particle distribution functions
          // according to the regularized formula
          T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
          for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop)
            cell[iPop] = dynamics -> computeEquilibrium(iPop,rho,u,uSqr) +
                         firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void PlaneFdBoundaryProcessor3D<T,DESCRIPTOR,direction,orientation>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}


template<typename T, typename DESCRIPTOR, int direction, int orientation>
template<int deriveDirection>
void PlaneFdBoundaryProcessor3D<T,DESCRIPTOR,direction,orientation>::
interpolateGradients(BlockLattice3D<T,DESCRIPTOR> const& blockLattice, T velDeriv[DESCRIPTOR::d],
                     int iX, int iY, int iZ) const
{
  fd::DirectedGradients3D<T, DESCRIPTOR, direction, orientation, deriveDirection, direction==deriveDirection>::
  interpolateVector(velDeriv, blockLattice, iX, iY, iZ);
}


////////  PlaneFdBoundaryProcessorGenerator3D ///////////////////////////////

template<typename T, typename DESCRIPTOR, int direction, int orientation>
PlaneFdBoundaryProcessorGenerator3D<T,DESCRIPTOR,direction,orientation>::
PlaneFdBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_)
{ }

template<typename T, typename DESCRIPTOR, int direction, int orientation>
PostProcessor3D<T,DESCRIPTOR>* PlaneFdBoundaryProcessorGenerator3D<T,DESCRIPTOR,direction,orientation>::
generate() const
{
  return new PlaneFdBoundaryProcessor3D<T,DESCRIPTOR, direction,orientation>
         ( this->x0, this->x1, this->y0, this->y1, this->z0, this->z1 );
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
PostProcessorGenerator3D<T,DESCRIPTOR>*
PlaneFdBoundaryProcessorGenerator3D<T,DESCRIPTOR,direction,orientation>::clone() const
{
  return new PlaneFdBoundaryProcessorGenerator3D<T,DESCRIPTOR,direction,orientation>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1);
}


////////  StraightConvectionBoundaryProcessorGenerator3D ////////////////////////////////

template<typename T, typename DESCRIPTOR, int direction, int orientation>
StraightConvectionBoundaryProcessor3D<T,DESCRIPTOR,direction,orientation>::
StraightConvectionBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T* uAv_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), uAv(uAv_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);

  saveCell = new T*** [(size_t)(x1_-x0_+1)];
  for (int iX=0; iX<=x1_-x0_; ++iX) {
    saveCell[iX] = new T** [(size_t)(y1_-y0_+1)];
    for (int iY=0; iY<=y1_-y0_; ++iY) {
      saveCell[iX][iY] = new T* [(size_t)(z1_-z0_+1)];
      for (int iZ=0; iZ<=z1_-z0_; ++iZ) {
        saveCell[iX][iY][iZ] = new T [(size_t)(DESCRIPTOR::q)];
        for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
          // default set to -1 in order to avoid wrong results at first call
          saveCell[iX][iY][iZ][iPop] = T(-1);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
StraightConvectionBoundaryProcessor3D<T,DESCRIPTOR,direction,orientation>::
~StraightConvectionBoundaryProcessor3D()
{
  for (int iX=0; iX<=x1-x0; ++iX) {
    for (int iY=0; iY<=y1-y0; ++iY) {
      for (int iZ=0; iZ<=z1-z0; ++iZ) {
        delete [] saveCell[iX][iY][iZ];
      }
      delete [] saveCell[iX][iY];
    }
    delete [] saveCell[iX];
  }
  delete [] saveCell;
}

template<typename T, typename DESCRIPTOR, int direction,int orientation>
void StraightConvectionBoundaryProcessor3D<T,DESCRIPTOR,direction,orientation>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    int iX;
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          Cell<T,DESCRIPTOR>& cell = blockLattice.get(iX,iY,iZ);
          for (int iPop = 0; iPop < DESCRIPTOR::q ; ++iPop) {
            if (descriptors::c<DESCRIPTOR>(iPop,direction)==-orientation) {
              // using default -1 avoids wrong first call
              if (!util::nearZero(1 + saveCell[iX-newX0][iY-newY0][iZ-newZ0][iPop]) ) {
                cell[iPop] = saveCell[iX-newX0][iY-newY0][iZ-newZ0][iPop];
              }
            }
          }

          T rho0, u0[3];
          T rho1, u1[3];
          T rho2, u2[3];
          if (direction==0) {
            blockLattice.get(iX,iY,iZ).computeRhoU(rho0,u0);
            blockLattice.get(iX-orientation,iY,iZ).computeRhoU(rho1,u1);
            blockLattice.get(iX-orientation*2,iY,iZ).computeRhoU(rho2,u2);
          }
          else if (direction==1) {
            blockLattice.get(iX,iY,iZ).computeRhoU(rho0,u0);
            blockLattice.get(iX,iY-orientation,iZ).computeRhoU(rho1,u1);
            blockLattice.get(iX,iY-orientation*2,iZ).computeRhoU(rho2,u2);
          }
          else {
            blockLattice.get(iX,iY,iZ).computeRhoU(rho0,u0);
            blockLattice.get(iX,iY,iZ-orientation).computeRhoU(rho1,u1);
            blockLattice.get(iX,iY,iZ-orientation*2).computeRhoU(rho2,u2);
          }

          // rho0 = T(1); rho1 = T(1); rho2 = T(1);

          T uDelta[3];
          T uAverage = rho0*u0[direction];
          if (uAv!=nullptr) {
            uAverage = *uAv * rho0;
          }
          uDelta[0]=-uAverage*0.5*(3*rho0*u0[0]-4*rho1*u1[0]+rho2*u2[0]);
          uDelta[1]=-uAverage*0.5*(3*rho0*u0[1]-4*rho1*u1[1]+rho2*u2[1]);
          uDelta[2]=-uAverage*0.5*(3*rho0*u0[2]-4*rho1*u1[2]+rho2*u2[2]);

          for (int iPop = 0; iPop < DESCRIPTOR::q ; ++iPop) {
            if (descriptors::c<DESCRIPTOR>(iPop,direction) == -orientation) {
              saveCell[iX-newX0][iY-newY0][iZ-newZ0][iPop] = cell[iPop] + descriptors::invCs2<T,DESCRIPTOR>()*descriptors::t<T,DESCRIPTOR>(iPop)*(uDelta[0]*descriptors::c<DESCRIPTOR>(iPop,0)+uDelta[1]*descriptors::c<DESCRIPTOR>(iPop,1)+uDelta[2]*descriptors::c<DESCRIPTOR>(iPop,2));
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, int direction,int orientation>
void StraightConvectionBoundaryProcessor3D<T,DESCRIPTOR,direction,orientation>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}


////////  StraightConvectionBoundaryProcessorGenerator2D ////////////////////////////////

template<typename T, typename DESCRIPTOR, int direction,int orientation>
StraightConvectionBoundaryProcessorGenerator3D<T,DESCRIPTOR, direction,orientation>::
StraightConvectionBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T* uAv_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_), uAv(uAv_)
{ }

template<typename T, typename DESCRIPTOR, int direction,int orientation>
PostProcessor3D<T,DESCRIPTOR>*
StraightConvectionBoundaryProcessorGenerator3D<T,DESCRIPTOR,direction,orientation>::generate() const
{
  return new StraightConvectionBoundaryProcessor3D<T,DESCRIPTOR,direction,orientation>
         ( this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, uAv);
}

template<typename T, typename DESCRIPTOR, int direction,int orientation>
PostProcessorGenerator3D<T,DESCRIPTOR>*
StraightConvectionBoundaryProcessorGenerator3D<T,DESCRIPTOR,direction,orientation>::clone() const
{
  return new StraightConvectionBoundaryProcessorGenerator3D<T,DESCRIPTOR,direction,orientation>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, uAv);
}


////////  OuterVelocityEdgeProcessor3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
OuterVelocityEdgeProcessor3D<T,DESCRIPTOR, plane,normal1,normal2>::
OuterVelocityEdgeProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{
  OLB_PRECONDITION (
    (plane==2 && x0==x1 && y0==y1) ||
    (plane==1 && x0==x1 && z0==z1) ||
    (plane==0 && y0==y1 && z0==z1)     );

}

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
void OuterVelocityEdgeProcessor3D<T,DESCRIPTOR, plane,normal1,normal2>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  using namespace olb::util::tensorIndices3D;

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect ( x0, x1, y0, y1, z0, z1,
                         x0_, x1_, y0_, y1_, z0_, z1_,
                         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {
    int iX;

#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          Cell<T,DESCRIPTOR>& cell = blockLattice.get(iX,iY,iZ);
          Dynamics<T,DESCRIPTOR>* dynamics = blockLattice.getDynamics(iX, iY, iZ);

          T rho10 = getNeighborRho(iX,iY,iZ,1,0, blockLattice);
          T rho01 = getNeighborRho(iX,iY,iZ,0,1, blockLattice);
          T rho20 = getNeighborRho(iX,iY,iZ,2,0, blockLattice);
          T rho02 = getNeighborRho(iX,iY,iZ,0,2, blockLattice);
          T rho = (T)2/(T)3*(rho01+rho10)-(T)1/(T)6*(rho02+rho20);

          T dA_uB_[3][3];
          interpolateGradients<plane,0>            ( blockLattice, dA_uB_[0], iX, iY, iZ );
          interpolateGradients<direction1,normal1> ( blockLattice, dA_uB_[1], iX, iY, iZ );
          interpolateGradients<direction2,normal2> ( blockLattice, dA_uB_[2], iX, iY, iZ );
          T dA_uB[3][3];
          for (int iBeta=0; iBeta<3; ++iBeta) {
            dA_uB[plane][iBeta]      = dA_uB_[0][iBeta];
            dA_uB[direction1][iBeta] = dA_uB_[1][iBeta];
            dA_uB[direction2][iBeta] = dA_uB_[2][iBeta];
          }
          T omega = dynamics -> getOmega();
          T sToPi = - rho / descriptors::invCs2<T,DESCRIPTOR>() / omega;
          T pi[util::TensorVal<DESCRIPTOR >::n];
          pi[xx] = (T)2 * dA_uB[0][0] * sToPi;
          pi[yy] = (T)2 * dA_uB[1][1] * sToPi;
          pi[zz] = (T)2 * dA_uB[2][2] * sToPi;
          pi[xy] = (dA_uB[0][1]+dA_uB[1][0]) * sToPi;
          pi[xz] = (dA_uB[0][2]+dA_uB[2][0]) * sToPi;
          pi[yz] = (dA_uB[1][2]+dA_uB[2][1]) * sToPi;

          // Computation of the particle distribution functions
          // according to the regularized formula
          T u[DESCRIPTOR::d];
          cell.computeU(u);
          T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);

          for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
            cell[iPop] = dynamics -> computeEquilibrium(iPop,rho,u,uSqr) +
                         firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi);
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
void OuterVelocityEdgeProcessor3D<T,DESCRIPTOR, plane,normal1,normal2>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
T OuterVelocityEdgeProcessor3D<T,DESCRIPTOR, plane,normal1,normal2>::
getNeighborRho(int x, int y, int z, int step1, int step2, BlockLattice3D<T,DESCRIPTOR> const& blockLattice)
{
  int coords[3] = {x, y, z};
  coords[direction1] += -normal1*step1;
  coords[direction2] += -normal2*step2;
  return blockLattice.get(coords[0], coords[1], coords[2]).computeRho();
}

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
template<int deriveDirection, int orientation>
void OuterVelocityEdgeProcessor3D<T,DESCRIPTOR, plane,normal1,normal2>::
interpolateGradients(BlockLattice3D<T,DESCRIPTOR> const& blockLattice,
                     T velDeriv[DESCRIPTOR::d],
                     int iX, int iY, int iZ) const
{
  fd::DirectedGradients3D<T,DESCRIPTOR,deriveDirection,orientation,deriveDirection,deriveDirection!=plane>::
  interpolateVector(velDeriv, blockLattice, iX, iY, iZ);
}

////////  OuterVelocityEdgeProcessorGenerator3D ///////////////////////////////

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, plane,normal1,normal2>::
OuterVelocityEdgeProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_)
{ }

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
PostProcessor3D<T,DESCRIPTOR>*
OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, plane,normal1,normal2>::
generate() const
{
  return new OuterVelocityEdgeProcessor3D < T,DESCRIPTOR, plane,normal1,normal2 >
         ( this->x0, this->x1, this->y0, this->y1, this->z0, this->z1);
}

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,DESCRIPTOR>*
OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, plane,normal1,normal2>::clone() const
{
  return new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, plane,normal1,normal2 >
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1);
}

/////////// OuterVelocityCornerProcessor3D /////////////////////////////////////

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
OuterVelocityCornerProcessor3D<T, DESCRIPTOR, xNormal, yNormal, zNormal>::
OuterVelocityCornerProcessor3D ( int x_, int y_, int z_ )
  : x(x_), y(y_), z(z_)
{ }

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
void OuterVelocityCornerProcessor3D<T, DESCRIPTOR, xNormal, yNormal, zNormal>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  using namespace olb::util::tensorIndices3D;
  Cell<T,DESCRIPTOR>& cell = blockLattice.get(x,y,z);
  Dynamics<T,DESCRIPTOR>* dynamics = blockLattice.getDynamics(x, y, z);

  T rho100 = blockLattice.get(x - 1*xNormal, y - 0*yNormal, z - 0*zNormal).computeRho();
  T rho010 = blockLattice.get(x - 0*xNormal, y - 1*yNormal, z - 0*zNormal).computeRho();
  T rho001 = blockLattice.get(x - 0*xNormal, y - 0*yNormal, z - 1*zNormal).computeRho();
  T rho200 = blockLattice.get(x - 2*xNormal, y - 0*yNormal, z - 0*zNormal).computeRho();
  T rho020 = blockLattice.get(x - 0*xNormal, y - 2*yNormal, z - 0*zNormal).computeRho();
  T rho002 = blockLattice.get(x - 0*xNormal, y - 0*yNormal, z - 2*zNormal).computeRho();
  T rho = (T)4/(T)9 * (rho001 + rho010 + rho100) - (T)1/(T)9 * (rho002 + rho020 + rho200);

  T dx_u[DESCRIPTOR::d], dy_u[DESCRIPTOR::d], dz_u[DESCRIPTOR::d];
  fd::DirectedGradients3D<T, DESCRIPTOR, 0, xNormal, 0, true>::interpolateVector(dx_u, blockLattice, x,y,z);
  fd::DirectedGradients3D<T, DESCRIPTOR, 1, yNormal, 0, true>::interpolateVector(dy_u, blockLattice, x,y,z);
  fd::DirectedGradients3D<T, DESCRIPTOR, 2, zNormal, 0, true>::interpolateVector(dz_u, blockLattice, x,y,z);

  T dx_ux = dx_u[0];
  T dy_ux = dy_u[0];
  T dz_ux = dz_u[0];
  T dx_uy = dx_u[1];
  T dy_uy = dy_u[1];
  T dz_uy = dz_u[1];
  T dx_uz = dx_u[2];
  T dy_uz = dy_u[2];
  T dz_uz = dz_u[2];
  T omega = dynamics -> getOmega();
  T sToPi = - rho / descriptors::invCs2<T,DESCRIPTOR>() / omega;
  T pi[util::TensorVal<DESCRIPTOR >::n];
  pi[xx] = (T)2 * dx_ux * sToPi;
  pi[yy] = (T)2 * dy_uy * sToPi;
  pi[zz] = (T)2 * dz_uz * sToPi;
  pi[xy] = (dx_uy + dy_ux) * sToPi;
  pi[xz] = (dx_uz + dz_ux) * sToPi;
  pi[yz] = (dy_uz + dz_uy) * sToPi;

  // Computation of the particle distribution functions
  // according to the regularized formula
  T u[DESCRIPTOR::d];
  cell.computeU(u);
  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);

  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = dynamics -> computeEquilibrium(iPop,rho,u,uSqr) +
                 firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi);
  }

}

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
void OuterVelocityCornerProcessor3D<T, DESCRIPTOR, xNormal, yNormal, zNormal>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  if (util::contained(x, y, z, x0_, x1_, y0_, y1_, z0_, z1_)) {
    process(blockLattice);
  }
}

////////  OuterVelocityCornerProcessorGenerator3D ///////////////////////////////

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
OuterVelocityCornerProcessorGenerator3D<T,DESCRIPTOR, xNormal,yNormal,zNormal>::
OuterVelocityCornerProcessorGenerator3D(int x_, int y_, int z_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x_,x_, y_,y_, z_,z_)
{ }

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
PostProcessor3D<T,DESCRIPTOR>*
OuterVelocityCornerProcessorGenerator3D<T,DESCRIPTOR, xNormal,yNormal,zNormal>::
generate() const
{
  return new OuterVelocityCornerProcessor3D<T,DESCRIPTOR, xNormal,yNormal,zNormal>
         ( this->x0, this->y0, this->z0 );
}

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,DESCRIPTOR>*
OuterVelocityCornerProcessorGenerator3D<T,DESCRIPTOR, xNormal,yNormal,zNormal>::clone() const
{
  return new OuterVelocityCornerProcessorGenerator3D<T,DESCRIPTOR, xNormal, yNormal, zNormal>
         (this->x0, this->y0, this->z0);
}


////////  SlipBoundaryProcessor3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
SlipBoundaryProcessor3D<T,DESCRIPTOR>::
SlipBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX, int discreteNormalY, int discreteNormalZ)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);
  int mirrorDirection0;
  int mirrorDirection1;
  int mirrorDirection2;
  int mult = 2 / (discreteNormalX*discreteNormalX + discreteNormalY*discreteNormalY + discreteNormalZ*discreteNormalZ);
  reflectionPop[0] = 0;
  for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
    reflectionPop[iPop] = 0;
    // iPop are the directions which pointing into the fluid, discreteNormal is pointing outwarts
    int scalarProduct = descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY + descriptors::c<DESCRIPTOR>(iPop,2)*discreteNormalZ;
    if ( scalarProduct < 0) {
      // bounce back for the case discreteNormalX = discreteNormalY = discreteNormalZ = 1, that is mult=0
      if (mult == 0) {
        mirrorDirection0 = -descriptors::c<DESCRIPTOR>(iPop,0);
        mirrorDirection1 = -descriptors::c<DESCRIPTOR>(iPop,1);
        mirrorDirection2 = -descriptors::c<DESCRIPTOR>(iPop,2);
      }
      else {
        mirrorDirection0 = descriptors::c<DESCRIPTOR>(iPop,0) - mult*scalarProduct*discreteNormalX;
        mirrorDirection1 = descriptors::c<DESCRIPTOR>(iPop,1) - mult*scalarProduct*discreteNormalY;
        mirrorDirection2 = descriptors::c<DESCRIPTOR>(iPop,2) - mult*scalarProduct*discreteNormalZ;
      }

      // run through all lattice directions and look for match of direction
      for (int i = 1; i < DESCRIPTOR::q; i++) {
        if (descriptors::c<DESCRIPTOR>(i,0)==mirrorDirection0
            && descriptors::c<DESCRIPTOR>(i,1)==mirrorDirection1
            && descriptors::c<DESCRIPTOR>(i,2)==mirrorDirection2) {
          reflectionPop[iPop] = i;
          break;
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void SlipBoundaryProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    int iX;
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
            if (reflectionPop[iPop]!=0) {
              //do reflection
              blockLattice.get(iX,iY,iZ)[iPop] = blockLattice.get(iX,iY,iZ)[reflectionPop[iPop]];
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void SlipBoundaryProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

////////  SlipBoundaryProcessorGenerator3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
SlipBoundaryProcessorGenerator3D<T,DESCRIPTOR>::
SlipBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_), discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_), discreteNormalZ(discreteNormalZ_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
SlipBoundaryProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new SlipBoundaryProcessor3D<T,DESCRIPTOR>
         ( this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, discreteNormalX, discreteNormalY, discreteNormalZ);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
SlipBoundaryProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new SlipBoundaryProcessorGenerator3D<T,DESCRIPTOR>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, discreteNormalX, discreteNormalY, discreteNormalZ);
}

////////  PartialSlipBoundaryProcessor3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
PartialSlipBoundaryProcessor3D<T,DESCRIPTOR>::
PartialSlipBoundaryProcessor3D(T tuner_, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX, int discreteNormalY, int discreteNormalZ)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), tuner(tuner_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);
  reflectionPop[0] = 0;
  for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
    reflectionPop[iPop] = 0;
    // iPop are the directions which pointing into the fluid, discreteNormal is pointing outwarts
    if (descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY + descriptors::c<DESCRIPTOR>(iPop,2)*discreteNormalZ < 0) {
      //std::cout << "-----" <<std::endl;
      int mirrorDirection0;
      int mirrorDirection1;
      int mirrorDirection2;
      int mult = 2 / (discreteNormalX*discreteNormalX + discreteNormalY*discreteNormalY + discreteNormalZ*discreteNormalZ);

      mirrorDirection0 = (descriptors::c<DESCRIPTOR>(iPop,0) - mult*(descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY + descriptors::c<DESCRIPTOR>(iPop,2)*discreteNormalZ)*discreteNormalX);
      mirrorDirection1 = (descriptors::c<DESCRIPTOR>(iPop,1) - mult*(descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY + descriptors::c<DESCRIPTOR>(iPop,2)*discreteNormalZ)*discreteNormalY);
      mirrorDirection2 = (descriptors::c<DESCRIPTOR>(iPop,2) - mult*(descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY + descriptors::c<DESCRIPTOR>(iPop,2)*discreteNormalZ)*discreteNormalZ);

      // bounce back for the case discreteNormalX = discreteNormalY = discreteNormalZ = 1, that is mult=0
      if (mult == 0) {
        mirrorDirection0 = -descriptors::c<DESCRIPTOR>(iPop,0);
        mirrorDirection1 = -descriptors::c<DESCRIPTOR>(iPop,1);
        mirrorDirection2 = -descriptors::c<DESCRIPTOR>(iPop,2);
      }

      // computes mirror jPop
      for (reflectionPop[iPop] = 1; reflectionPop[iPop] < DESCRIPTOR::q ; reflectionPop[iPop]++) {
        if (descriptors::c<DESCRIPTOR>(reflectionPop[iPop],0)==mirrorDirection0 && descriptors::c<DESCRIPTOR>(reflectionPop[iPop],1)==mirrorDirection1 && descriptors::c<DESCRIPTOR>(reflectionPop[iPop],2)==mirrorDirection2) {
          break;
        }
      }
      //std::cout <<iPop << " to "<< jPop <<" for discreteNormal= "<< discreteNormalX << "/"<<discreteNormalY <<std::endl;
    }
  }
}

template<typename T, typename DESCRIPTOR>
void PartialSlipBoundaryProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    int iX;
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
            if (reflectionPop[iPop]!=0) {
              //do reflection
              blockLattice.get(iX,iY,iZ)[iPop] = tuner*blockLattice.get(iX,iY,iZ)[reflectionPop[iPop]];
            }
          }
          for (int iPop = 1; iPop < DESCRIPTOR::q/2 ; ++iPop) {
            T provv = blockLattice.get(iX,iY,iZ)[descriptors::opposite<DESCRIPTOR>(iPop)];
            blockLattice.get(iX,iY,iZ)[descriptors::opposite<DESCRIPTOR>(iPop)] +=
              (1.-tuner)*blockLattice.get(iX,iY,iZ)[iPop];
            blockLattice.get(iX,iY,iZ)[iPop] += (1.-tuner)*provv;
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void PartialSlipBoundaryProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

////////  PartialSlipBoundaryProcessorGenerator3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
PartialSlipBoundaryProcessorGenerator3D<T,DESCRIPTOR>::
PartialSlipBoundaryProcessorGenerator3D(T tuner_, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_), discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_), discreteNormalZ(discreteNormalZ_), tuner(tuner_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
PartialSlipBoundaryProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new PartialSlipBoundaryProcessor3D<T,DESCRIPTOR>
         (tuner, this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, discreteNormalX, discreteNormalY, discreteNormalZ);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
PartialSlipBoundaryProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new PartialSlipBoundaryProcessorGenerator3D<T,DESCRIPTOR>
         (tuner, this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, discreteNormalX, discreteNormalY, discreteNormalZ);
}

////////  FreeEnergyWallProcessor3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyWallProcessor3D<T,DESCRIPTOR>::FreeEnergyWallProcessor3D(
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
  int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_, T addend_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), discreteNormalX(discreteNormalX_),
    discreteNormalY(discreteNormalY_), discreteNormalZ(discreteNormalZ_), addend(addend_)
{ }

template<typename T, typename DESCRIPTOR>
void FreeEnergyWallProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          T rhoBulk = blockLattice.get(iX-discreteNormalX, iY-discreteNormalY, iZ-discreteNormalZ).computeRho();
          T rhoTmp = 0;
          for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
            rhoTmp += blockLattice.get(iX,iY,iZ)[iPop];
          }
          T rho = rhoBulk + addend;
          blockLattice.get(iX,iY,iZ)[0] = rho - rhoTmp - 1;
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyWallProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

////////  FreeEnergyWallProcessorGenerator3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyWallProcessorGenerator3D<T,DESCRIPTOR>::FreeEnergyWallProcessorGenerator3D(
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
  int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_, T addend_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_), discreteNormalX(discreteNormalX_),
    discreteNormalY(discreteNormalY_), discreteNormalZ(discreteNormalZ_), addend(addend_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
FreeEnergyWallProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new FreeEnergyWallProcessor3D<T,DESCRIPTOR>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1,
          discreteNormalX, discreteNormalY, discreteNormalZ, addend);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
FreeEnergyWallProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new FreeEnergyWallProcessorGenerator3D<T,DESCRIPTOR>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1,
          discreteNormalX, discreteNormalY, discreteNormalZ, addend);
}


////////  FreeEnergyChemPotBoundaryProcessor3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyChemPotBoundaryProcessor3D<T,DESCRIPTOR>::
FreeEnergyChemPotBoundaryProcessor3D(
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX_,
  int discreteNormalY_, int discreteNormalZ_, int latticeNumber_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), discreteNormalX(discreteNormalX_),
    discreteNormalY(discreteNormalY_), discreteNormalZ(discreteNormalZ_), latticeNumber(latticeNumber_)
{ }

template<typename T, typename DESCRIPTOR>
void FreeEnergyChemPotBoundaryProcessor3D<T,DESCRIPTOR>::
processSubDomain(
  BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          blockLattice.get(iX,iY,iZ).template setField<descriptors::CHEM_POTENTIAL>(
            blockLattice.get(iX-discreteNormalX, iY-discreteNormalY, iZ-discreteNormalZ).template getField<descriptors::CHEM_POTENTIAL>()
          );
          if (latticeNumber == 1) {
            T rho0 = blockLattice.get(iX,iY,iZ).computeRho();
            T rho1 = blockLattice.get(iX-discreteNormalX, iY-discreteNormalY, iZ-discreteNormalZ).computeRho();
            *(blockLattice.get(iX,iY,iZ).template getFieldPointer<descriptors::CHEM_POTENTIAL>()) += (rho1 / rho0 - 1) / descriptors::invCs2<T,DESCRIPTOR>();
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyChemPotBoundaryProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

////////  FreeEnergyChemPotBoundaryProcessorGenerator3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyChemPotBoundaryProcessorGenerator3D<T,DESCRIPTOR>::
FreeEnergyChemPotBoundaryProcessorGenerator3D(
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
  int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_, int latticeNumber_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_), discreteNormalX(discreteNormalX_),
    discreteNormalY(discreteNormalY_), discreteNormalZ(discreteNormalZ_), latticeNumber(latticeNumber_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
FreeEnergyChemPotBoundaryProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new FreeEnergyChemPotBoundaryProcessor3D<T,DESCRIPTOR>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1,
          discreteNormalX, discreteNormalY, discreteNormalZ, latticeNumber);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
FreeEnergyChemPotBoundaryProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new FreeEnergyChemPotBoundaryProcessorGenerator3D<T,DESCRIPTOR>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1,
          discreteNormalX, discreteNormalY, discreteNormalZ, latticeNumber);
}


////////  FreeEnergyConvectiveProcessor3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyConvectiveProcessor3D<T,DESCRIPTOR>::
FreeEnergyConvectiveProcessor3D(
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
  int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), discreteNormalX(discreteNormalX_),
    discreteNormalY(discreteNormalY_), discreteNormalZ(discreteNormalZ_)
{ }

template<typename T, typename DESCRIPTOR>
void FreeEnergyConvectiveProcessor3D<T,DESCRIPTOR>::
processSubDomain(
  BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          T rho, uPerp, rho0, rho1, u[3];
          rho0 = blockLattice.get(iX,iY,iZ).computeRho();
          blockLattice.get(iX-discreteNormalX,iY-discreteNormalY,iZ-discreteNormalZ).computeRhoU(rho1, u);

          if (discreteNormalZ == 0) {
            if (discreteNormalY == 0) {
              if (discreteNormalX < 0) {
                uPerp = u[0];
              }
              else {
                uPerp = -u[0];
              }
            }
            else if (discreteNormalX == 0) {
              if (discreteNormalY < 0) {
                uPerp = u[1];
              }
              else {
                uPerp = -u[1];
              }
            }
            else {
              uPerp = sqrt(u[0] * u[0] + u[1] * u[1]);
            }
          }
          else if (discreteNormalY == 0) {
            if (discreteNormalX == 0) {
              if (discreteNormalZ < 0) {
                uPerp = u[2];
              }
              else {
                uPerp = -u[2];
              }
            }
            else {
              uPerp = sqrt(u[0] * u[0] + u[2] * u[2]);
            }
          }
          else if (discreteNormalX == 0) {
            uPerp = sqrt(u[1] * u[1] + u[2] * u[2]);
          }
          else {
            uPerp = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
          }

          rho = (rho0 + uPerp * rho1) / (1. + uPerp);
          blockLattice.get(iX,iY,iZ).defineRho(rho);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void FreeEnergyConvectiveProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

////////  FreeEnergyConvectiveProcessorGenerator3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
FreeEnergyConvectiveProcessorGenerator3D<T,DESCRIPTOR>::
FreeEnergyConvectiveProcessorGenerator3D(
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
  int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_), discreteNormalX(discreteNormalX_),
    discreteNormalY(discreteNormalY_), discreteNormalZ(discreteNormalZ_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
FreeEnergyConvectiveProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new FreeEnergyConvectiveProcessor3D<T,DESCRIPTOR>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1,
          discreteNormalX, discreteNormalY, discreteNormalZ);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
FreeEnergyConvectiveProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new FreeEnergyConvectiveProcessorGenerator3D<T,DESCRIPTOR>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1,
          discreteNormalX, discreteNormalY, discreteNormalZ);
}

}  // namespace olb

#endif
