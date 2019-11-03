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

#ifndef EXTENDED_FINITE_DIFFERENCE_BOUNDARY_2D_HH
#define EXTENDED_FINITE_DIFFERENCE_BOUNDARY_2D_HH

#include "extendedFiniteDifferenceBoundary2D.h"
#include "core/finiteDifference2D.h"
#include "core/blockLattice2D.h"
#include "core/util.h"
#include "dynamics/lbHelpers.h"
#include "dynamics/firstOrderLbHelpers.h"
#include "boundaryInstantiator2D.h"


namespace olb {


///////////  ExtendedStraightFdBoundaryPostProcessor2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, int direction, int orientation>
ExtendedStraightFdBoundaryPostProcessor2D<T,DESCRIPTOR,direction,orientation>::
ExtendedStraightFdBoundaryPostProcessor2D(int x0_, int x1_, int y0_, int y1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void ExtendedStraightFdBoundaryPostProcessor2D<T,DESCRIPTOR,direction,orientation>::
processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  using namespace olb::util::tensorIndices2D;
  typedef lbHelpers<T,DESCRIPTOR> lbH;
  typedef DESCRIPTOR L;
  enum {x,y};

  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {
    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        Cell<T,DESCRIPTOR>& cell = blockLattice.get(iX,iY);
        T rho, u[L::d];
        cell.computeRhoU(rho,u);

        T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);

        T dx_U[L::d], dy_U[L::d];
        interpolateGradients<0>(blockLattice, dx_U, iX, iY);
        interpolateGradients<1>(blockLattice, dy_U, iX, iY);

        T rhoGradU[L::d][L::d];
        rhoGradU[x][x] = rho * dx_U[x];
        rhoGradU[x][y] = rho * dx_U[y];
        rhoGradU[y][x] = rho * dy_U[x];
        rhoGradU[y][y] = rho * dy_U[y];

        T omega = blockLattice.getDynamics(iX, iY) -> getOmega();
        T sToPi = - (T)1 / descriptors::invCs2<T,DESCRIPTOR>() / omega;

        T pi[util::TensorVal<DESCRIPTOR >::n];
        pi[xx] = (T)2 * rhoGradU[x][x] * sToPi;
        pi[yy] = (T)2 * rhoGradU[y][y] * sToPi;
        pi[xy] = (rhoGradU[x][y] + rhoGradU[y][x]) * sToPi;
        // here ends the "regular" fdBoudaryCondition
        // implemented in OpenLB

        // first we compute the term
        // (c_{i\alpha} \nabla_\beta)(rho*u_\alpha*u_\beta)
        T dx_rho, dy_rho;
        interpolateGradients<0>(blockLattice, dx_rho, iX, iY);
        interpolateGradients<1>(blockLattice, dy_rho, iX, iY);
        for (int iPop = 0; iPop < L::q; ++iPop) {
          T cGradRhoUU = T();
          for (int iAlpha=0; iAlpha < L::d; ++iAlpha) {
            cGradRhoUU += descriptors::c<L>(iPop,iAlpha) * (
                            dx_rho*u[iAlpha]*u[x] +
                            dx_U[iAlpha]*rho*u[x] +
                            dx_U[x]*rho*u[iAlpha] + //end of dx derivatice
                            dy_rho*u[iAlpha]*u[y] +
                            dy_U[iAlpha]*rho*u[y] +
                            dy_U[y]*rho*u[iAlpha]);
          }

          // then we compute the term
          // c_{i\gamma}\nabla_{\gamma}(\rho*u_\alpha * u_\beta)
          T cDivRhoUU[L::d][L::d]; //first step towards QcdivRhoUU
          for (int iAlpha = 0; iAlpha < L::d; ++iAlpha) {
            for (int iBeta = 0; iBeta < L::d; ++iBeta) {
              cDivRhoUU[iAlpha][iBeta] = descriptors::c<L>(iPop,x) *
                                         (dx_rho*u[iAlpha]*u[iBeta] +
                                          dx_U[iAlpha]*rho*u[iBeta] +
                                          dx_U[iBeta]*rho*u[iAlpha])
                                         + descriptors::c<L>(iPop,y) *
                                         (dy_rho*u[iAlpha]*u[iBeta] +
                                          dy_U[iAlpha]*rho*u[iBeta] +
                                          dy_U[iBeta]*rho*u[iAlpha]);
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
              qRhoGradU   += ci_ci * rhoGradU[iAlpha][iBeta];
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

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void ExtendedStraightFdBoundaryPostProcessor2D<T,DESCRIPTOR,direction,orientation>::
process(BlockLattice2D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
template<int deriveDirection>
void ExtendedStraightFdBoundaryPostProcessor2D<T,DESCRIPTOR,direction,orientation>::
interpolateGradients(BlockLattice2D<T,DESCRIPTOR> const& blockLattice,
                     T velDeriv[DESCRIPTOR::d], int iX, int iY) const
{
  fd::DirectedGradients2D<T, DESCRIPTOR, direction, orientation, direction==deriveDirection>::
  interpolateVector(velDeriv, blockLattice, iX, iY);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
template<int deriveDirection>
void ExtendedStraightFdBoundaryPostProcessor2D<T,DESCRIPTOR,direction,orientation>::
interpolateGradients(BlockLattice2D<T,DESCRIPTOR> const& blockLattice, T& rhoDeriv, int iX, int iY) const
{
  fd::DirectedGradients2D<T, DESCRIPTOR, direction, orientation, direction==deriveDirection>::
  interpolateScalar(rhoDeriv, blockLattice, iX, iY);
}


////////  ExtendedStraightFdBoundaryPostProcessorGenerator ////////////////////////////////

template<typename T, typename DESCRIPTOR, int direction, int orientation>
ExtendedStraightFdBoundaryProcessorGenerator2D<T,DESCRIPTOR,direction,orientation>::
ExtendedStraightFdBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_)
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_)
{ }

template<typename T, typename DESCRIPTOR, int direction, int orientation>
PostProcessor2D<T,DESCRIPTOR>*
ExtendedStraightFdBoundaryProcessorGenerator2D<T,DESCRIPTOR,direction,orientation>::generate() const
{
  return new ExtendedStraightFdBoundaryPostProcessor2D<T,DESCRIPTOR,direction,orientation>
         ( this->x0, this->x1, this->y0, this->y1 );
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
ExtendedStraightFdBoundaryProcessorGenerator2D<T,DESCRIPTOR,direction,orientation>::clone() const
{
  return new ExtendedStraightFdBoundaryProcessorGenerator2D<T,DESCRIPTOR,direction,orientation>
         (this->x0, this->x1, this->y0, this->y1);
}


////////// Class ExtendedFdBoundaryManager2D /////////////////////////////////////////

template<typename T, typename DESCRIPTOR, class MixinDynamics>
class ExtendedFdBoundaryManager2D {
public:
  template<int direction, int orientation> static Momenta<T,DESCRIPTOR>*
  getVelocityBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,DESCRIPTOR>*
  getVelocityBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1);

  template<int direction, int orientation> static Momenta<T,DESCRIPTOR>*
  getPressureBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,DESCRIPTOR>*
  getPressureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getPressureBoundaryProcessor(int x0, int x1, int y0, int y1);

  template<int direction, int orientation> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, T* uAv=NULL);

  template<int xNormal, int yNormal> static Momenta<T,DESCRIPTOR>*
  getExternalVelocityCornerMomenta();
  template<int xNormal, int yNormal> static Dynamics<T,DESCRIPTOR>*
  getExternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int xNormal, int yNormal> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getExternalVelocityCornerProcessor(int x, int y);

  template<int xNormal, int yNormal> static Momenta<T,DESCRIPTOR>*
  getInternalVelocityCornerMomenta();
  template<int xNormal, int yNormal> static Dynamics<T,DESCRIPTOR>*
  getInternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  template<int xNormal, int yNormal> static PostProcessorGenerator2D<T,DESCRIPTOR>*
  getInternalVelocityCornerProcessor(int x, int y);
};


template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>* ExtendedFdBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getVelocityBoundaryMomenta()
{
  return new BasicDirichletBM<T,DESCRIPTOR,VelocityBM,direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* ExtendedFdBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getVelocityBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>* ExtendedFdBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1)
{
  return new ExtendedStraightFdBoundaryProcessorGenerator2D<T,DESCRIPTOR,direction,orientation>
         (x0, x1, y0, y1);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,DESCRIPTOR>* ExtendedFdBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getPressureBoundaryMomenta()
{
  return new BasicDirichletBM<T,DESCRIPTOR,PressureBM,direction,orientation>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,DESCRIPTOR>* ExtendedFdBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getPressureBoundaryDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
ExtendedFdBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getPressureBoundaryProcessor(int x0, int x1, int y0, int y1)
{
  return new ExtendedStraightFdBoundaryProcessorGenerator2D<T,DESCRIPTOR,direction,orientation>
         (x0, x1, y0, y1);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,DESCRIPTOR>*
ExtendedFdBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, T* uAv)
{
  return nullptr;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Momenta<T,DESCRIPTOR>*
ExtendedFdBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getExternalVelocityCornerMomenta()
{
  return new FixedVelocityBM<T,DESCRIPTOR>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Dynamics<T,DESCRIPTOR>* ExtendedFdBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getExternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
PostProcessorGenerator2D<T,DESCRIPTOR>*
ExtendedFdBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getExternalVelocityCornerProcessor(int x, int y)
{
  return new OuterVelocityCornerProcessorGenerator2D<T,DESCRIPTOR, xNormal,yNormal> (x,y);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Momenta<T,DESCRIPTOR>*
ExtendedFdBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getInternalVelocityCornerMomenta()
{
  return new InnerCornerVelBM2D<T,DESCRIPTOR, xNormal,yNormal>;
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
Dynamics<T,DESCRIPTOR>* ExtendedFdBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::
getInternalVelocityCornerDynamics(T omega, Momenta<T,DESCRIPTOR>& momenta)
{
  return new CombinedRLBdynamics<T,DESCRIPTOR, MixinDynamics>(omega, momenta);
}

template<typename T, typename DESCRIPTOR, class MixinDynamics>
template<int xNormal, int yNormal>
PostProcessorGenerator2D<T,DESCRIPTOR>*
ExtendedFdBoundaryManager2D<T,DESCRIPTOR,MixinDynamics>::getInternalVelocityCornerProcessor (int x, int y)
{
  return nullptr;
}


////////// Factory functions //////////////////////////////////////////////////

template<typename T, typename DESCRIPTOR, typename MixinDynamics>
OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
createExtendedFdBoundaryCondition2D(BlockLatticeStructure2D<T,DESCRIPTOR>& block)
{
  return new BoundaryConditionInstantiator2D <
         T, DESCRIPTOR,
         ExtendedFdBoundaryManager2D<T,DESCRIPTOR, MixinDynamics> > (block);
}


}  // namespace olb

#endif
