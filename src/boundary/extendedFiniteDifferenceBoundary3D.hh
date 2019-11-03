/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Orestis Malaspinas, Jonas Latt
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

#ifndef EXTENDED_FINITE_DIFFERENCE_BOUNDARY_3D_HH
#define EXTENDED_FINITE_DIFFERENCE_BOUNDARY_3D_HH

#include "extendedFiniteDifferenceBoundary3D.h"
#include "core/finiteDifference3D.h"
#include "core/blockLattice3D.h"
#include "core/util.h"
#include "dynamics/lbHelpers.h"
#include "dynamics/firstOrderLbHelpers.h"
#include "boundaryInstantiator3D.h"
#include "momentaOnBoundaries3D.h"


namespace olb {

////////  ExtendedFdPlaneBoundaryPostProcessor3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, int direction, int orientation>
ExtendedFdPlaneBoundaryPostProcessor3D <T,DESCRIPTOR,direction,orientation>::
ExtendedFdPlaneBoundaryPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void ExtendedFdPlaneBoundaryPostProcessor3D<T,DESCRIPTOR,direction,orientation>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  typedef DESCRIPTOR L;
  using namespace olb::util::tensorIndices3D;
  typedef lbHelpers<T,DESCRIPTOR> lbH;
  enum {x,y,z};

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          Cell<T,DESCRIPTOR>& cell = blockLattice.get(iX,iY,iZ);
          T rho, u[L::d];
          cell.computeRhoU(rho,u);
          T dx_U[DESCRIPTOR::d], dy_U[DESCRIPTOR::d], dz_U[DESCRIPTOR::d];
          interpolateGradients<0>(blockLattice, dx_U, iX, iY, iZ);
          interpolateGradients<1>(blockLattice, dy_U, iX, iY, iZ);
          interpolateGradients<2>(blockLattice, dz_U, iX, iY, iZ);

          T rhoGradU[L::d][L::d];
          rhoGradU[x][x] = rho *dx_U[x];
          rhoGradU[x][y] = rho *dx_U[y];
          rhoGradU[x][z] = rho *dx_U[z];
          rhoGradU[y][x] = rho *dy_U[x];
          rhoGradU[y][y] = rho *dy_U[y];
          rhoGradU[y][z] = rho *dy_U[z];
          rhoGradU[z][x] = rho *dz_U[x];
          rhoGradU[z][y] = rho *dz_U[y];
          rhoGradU[z][z] = rho *dz_U[z];

          T omega = blockLattice.getDynamics(iX, iY, iZ) -> getOmega();
          T sToPi = - (T)1 / descriptors::invCs2<T,DESCRIPTOR>() / omega;
          T pi[util::TensorVal<DESCRIPTOR >::n];

          pi[xx] = (T)2 * rhoGradU[x][x] * sToPi;
          pi[yy] = (T)2 * rhoGradU[y][y] * sToPi;
          pi[zz] = (T)2 * rhoGradU[z][z] * sToPi;
          pi[xy] = (rhoGradU[x][y] + rhoGradU[y][x]) * sToPi;
          pi[xz] = (rhoGradU[x][z] + rhoGradU[z][x]) * sToPi;
          pi[yz] = (rhoGradU[y][z] + rhoGradU[z][y]) * sToPi;

          // here ends the "regular" fdBoudaryCondition
          // implemented in OpenLB

          T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);

          // first we compute the term
          // (c_{i\alpha} \nabla_\beta)(rho*u_\alpha*u_\beta)
          T dx_rho, dy_rho, dz_rho;
          interpolateGradients<0>(blockLattice, dx_rho, iX, iY, iZ);
          interpolateGradients<1>(blockLattice, dy_rho, iX, iY, iZ);
          interpolateGradients<2>(blockLattice, dz_rho, iX, iY, iZ);
          for (int iPop = 0; iPop < L::q; ++iPop) {
            T cGradRhoUU = T();
            for (int iAlpha=0; iAlpha < L::d; ++iAlpha) {
              cGradRhoUU += descriptors::c<L>(iPop,iAlpha) * (
                              dx_rho*u[iAlpha]*u[x] +
                              dx_U[iAlpha]*rho*u[x] +
                              dx_U[x]*rho*u[iAlpha] + //end of dx derivative
                              dy_rho*u[iAlpha]*u[y] +
                              dy_U[iAlpha]*rho*u[y] +
                              dy_U[y]*rho*u[iAlpha] +//end of dy derivative
                              dz_rho*u[iAlpha]*u[z] +
                              dz_U[iAlpha]*rho*u[z] +
                              dz_U[z]*rho*u[iAlpha]);
            }

            // then we compute the term
            // c_{i\gamma}\nabla_{\gamma}(\rho*u_\alpha * u_\beta)
            T cDivRhoUU[L::d][L::d]; //first step towards QcdivRhoUU
            for (int iAlpha = 0; iAlpha < L::d; ++iAlpha) {
              for (int iBeta = 0; iBeta < L::d; ++iBeta) {
                cDivRhoUU[iAlpha][iBeta] = descriptors::c<L>(iPop,x)*
                                           (dx_rho*u[iAlpha]*u[iBeta] +
                                            dx_U[iAlpha]*rho*u[iBeta] +
                                            dx_U[iBeta]*rho*u[iAlpha])
                                           +descriptors::c<L>(iPop,y)*
                                           (dy_rho*u[iAlpha]*u[iBeta] +
                                            dy_U[iAlpha]*rho*u[iBeta] +
                                            dy_U[iBeta]*rho*u[iAlpha])
                                           +descriptors::c<L>(iPop,z)*
                                           (dz_rho*u[iAlpha]*u[iBeta] +
                                            dz_U[iAlpha]*rho*u[iBeta] +
                                            dz_U[iBeta]*rho*u[iAlpha]);
              }
            }

            //Finally we can compute
            // Q_{i\alpha\beta}c_{i\gamma}\nabla_{\gamma}(\rho*u_\alpha * u_\beta)
            // and Q_{i\alpha\beta}\rho\nabla_{\alpha}u_\beta
            T qCdivRhoUU = T();
            T qRhoGradU = T();
            for (int iAlpha = 0; iAlpha < L::d; ++iAlpha) {
              for (int iBeta = 0; iBeta < L::d; ++iBeta) {
                int ci_ci = descriptors::c<L>(iPop,iAlpha)*descriptors::c<L>(iPop,iBeta);
                qCdivRhoUU  += ci_ci * cDivRhoUU[iAlpha][iBeta];
                qRhoGradU += ci_ci * rhoGradU[iAlpha][iBeta];
                if (iAlpha == iBeta) {
                  qCdivRhoUU -= cDivRhoUU[iAlpha][iBeta]/descriptors::invCs2<T,L>();
                  qRhoGradU  -= rhoGradU[iAlpha][iBeta]/descriptors::invCs2<T,L>();
                }
              }
            }

            // we then can reconstruct the value of the populations
            // according to the complete C-E expansion term
            cell[iPop] = lbH::equilibrium(iPop,rho,u,uSqr)
                         - descriptors::t<T,L>(iPop) * descriptors::invCs2<T,L>() / omega
                         * (qRhoGradU - cGradRhoUU + 0.5*descriptors::invCs2<T,L>()*qCdivRhoUU);
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void ExtendedFdPlaneBoundaryPostProcessor3D<T,DESCRIPTOR,direction,orientation>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
template<int deriveDirection>
void ExtendedFdPlaneBoundaryPostProcessor3D<T,DESCRIPTOR,direction,orientation>::
interpolateGradients(BlockLattice3D<T,DESCRIPTOR> const& blockLattice,
                     T velDeriv[DESCRIPTOR::d],
                     int iX, int iY, int iZ) const
{
  fd::DirectedGradients3D<T, DESCRIPTOR, direction, orientation, deriveDirection, direction==deriveDirection>::
  interpolateVector(velDeriv, blockLattice, iX, iY, iZ);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
template<int deriveDirection>
void ExtendedFdPlaneBoundaryPostProcessor3D<T,DESCRIPTOR,direction,orientation>::
interpolateGradients(BlockLattice3D<T,DESCRIPTOR> const& blockLattice,
                     T& rhoDeriv, int iX, int iY, int iZ) const
{
  fd::DirectedGradients3D<T, DESCRIPTOR, direction, orientation, deriveDirection, direction==deriveDirection>::
  interpolateScalar(rhoDeriv, blockLattice, iX, iY, iZ);
}


////////  ExtendedFdPlaneBoundaryProcessorGenertor3D ///////////////////////////////

template<typename T, typename DESCRIPTOR, int direction, int orientation>
ExtendedFdPlaneBoundaryProcessorGenerator3D<T,DESCRIPTOR,direction,orientation>::
ExtendedFdPlaneBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_)
{ }

template<typename T, typename DESCRIPTOR, int direction, int orientation>
PostProcessor3D<T,DESCRIPTOR>*
ExtendedFdPlaneBoundaryProcessorGenerator3D<T,DESCRIPTOR,direction,orientation>::generate() const
{
  return new ExtendedFdPlaneBoundaryPostProcessor3D<T,DESCRIPTOR,direction,orientation>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ExtendedFdPlaneBoundaryProcessorGenerator3D<T,DESCRIPTOR,direction,orientation>::clone() const
{
  return new ExtendedFdPlaneBoundaryProcessorGenerator3D<T,DESCRIPTOR,direction,orientation>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1);
}


////////// Class ExtendedFdBoundaryManager3D /////////////////////////////////////////

template<typename T, typename DESCRIPTOR, class MixinDynamics>
class ExtendedFdBoundaryManager3D {
public:
  template<int direction, int orientation> static Momenta<T,DESCRIPTOR>*
  getVelocityBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,DESCRIPTOR>*
  getVelocityBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int direction, int orientation> static Momenta<T,DESCRIPTOR>*
  getPressureBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,DESCRIPTOR>*
  getPressureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getPressureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int direction, int orientation> static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1, T* uAv=NULL);

  template<int plane, int normal1, int normal2> static Momenta<T,DESCRIPTOR>*
  getExternalVelocityEdgeMomenta();
  template<int plane, int normal1, int normal2> static Dynamics<T,DESCRIPTOR>*
  getExternalVelocityEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int plane, int normal1, int normal2> static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getExternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int plane, int normal1, int normal2> static Momenta<T,DESCRIPTOR>*
  getInternalVelocityEdgeMomenta();
  template<int plane, int normal1, int normal2> static Dynamics<T,DESCRIPTOR>*
  getInternalVelocityEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int plane, int normal1, int normal2> static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getInternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int xNormal, int yNormal, int zNormal> static Momenta<T,DESCRIPTOR>*
  getExternalVelocityCornerMomenta();
  template<int xNormal, int yNormal, int zNormal> static Dynamics<T,DESCRIPTOR>*
  getExternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int xNormal, int yNormal, int zNormal> static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getExternalVelocityCornerProcessor(int x, int y, int z);

  template<int xNormal, int yNormal, int zNormal> static Momenta<T,DESCRIPTOR>*
  getInternalVelocityCornerMomenta();
  template<int xNormal, int yNormal, int zNormal> static Dynamics<T,DESCRIPTOR>*
  getInternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int xNormal, int yNormal, int zNormal> static PostProcessorGenerator3D<T,DESCRIPTOR>*
  getInternalVelocityCornerProcessor(int x, int y, int z);
};


template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getVelocityBoundaryMomenta()
{
  return new BasicDirichletBM<T,DESCRIPTOR,VelocityBM, direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getVelocityBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return new ExtendedFdPlaneBoundaryProcessorGenerator3D
         <T,DESCRIPTOR, direction,orientation>(x0,x1, y0,y1, z0,z1);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getPressureBoundaryMomenta()
{
  return new BasicDirichletBM<T,DESCRIPTOR,PressureBM, direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getPressureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getPressureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return new ExtendedFdPlaneBoundaryProcessorGenerator3D
         <T,DESCRIPTOR, direction,orientation>(x0,x1, y0,y1, z0,z1);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1, T* uAv)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Momenta<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getExternalVelocityEdgeMomenta()
{
  return new FixedVelocityBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Dynamics<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return new OuterVelocityEdgeProcessorGenerator3D<T,DESCRIPTOR, plane,normal1,normal2>(x0,x1, y0,y1, z0,z1);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Momenta<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getInternalVelocityEdgeMomenta()
{
  return new InnerEdgeVelBM3D<T,DESCRIPTOR, plane,normal1,normal2>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
Dynamics<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityEdgeDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getExternalVelocityCornerMomenta()
{
  return new FixedVelocityBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityCornerProcessor(int x, int y, int z)
{
  return new OuterVelocityCornerProcessorGenerator3D<T,DESCRIPTOR, xNormal,yNormal,zNormal> (x,y,z);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::getInternalVelocityCornerMomenta()
{
  return new InnerCornerVelBM3D<T,DESCRIPTOR, xNormal,yNormal,zNormal>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ExtendedFdBoundaryManager3D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityCornerProcessor(int x, int y, int z)
{
  return nullptr;
}


////////// Factory function //////////////////////////////////////////////////

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
OnLatticeBoundaryCondition3D<T,DESCRIPTOR>*
createExtendedFdBoundaryCondition3D(BlockLatticeStructure3D<T,DESCRIPTOR>& block)
{
  return new BoundaryConditionInstantiator3D <
         T, DESCRIPTOR,
         ExtendedFdBoundaryManager3D<T,DESCRIPTOR, MixinDynamics> > (block);
}


}  // namespace olb

#endif
