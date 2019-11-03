/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2016 Robin Trunk
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

#include "advectionDiffusionBoundaryPostProcessor3D.h"
#include "core/blockLattice3D.h"
#include "core/util.h"
#include "dynamics/lbHelpers.h"
#include "dynamics/firstOrderLbHelpers.h"

namespace olb {

////////  ConvectionBoundaryProcessor3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
ConvectionBoundaryProcessor3D<T,DESCRIPTOR>::
ConvectionBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                              int discreteNormalX, int discreteNormalY, int discreteNormalZ)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);

  interpolationPop[0] = 0;
  for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
    interpolationPop[iPop] = 0;
    // find incoming iPop from material 0
    if (descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY + descriptors::c<DESCRIPTOR>(iPop,2)*discreteNormalZ > 0) {
      // check for material number of neighbours has to be one level higher
      interpolationPop[iPop] = 1;
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ConvectionBoundaryProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_,
                 int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
            if (interpolationPop[iPop]!=0) {
              //do reflection
              blockLattice.get(iX,iY,iZ)[iPop] =
                ( blockLattice.get(iX+descriptors::c<DESCRIPTOR>(iPop,0),iY+descriptors::c<DESCRIPTOR>(iPop,1),iZ+descriptors::c<DESCRIPTOR>(iPop,2))[iPop]
                  + blockLattice.get(iX+2*descriptors::c<DESCRIPTOR>(iPop,0),iY+2*descriptors::c<DESCRIPTOR>(iPop,1),iZ+2*descriptors::c<DESCRIPTOR>(iPop,2))[iPop]
                ) * 0.5;
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ConvectionBoundaryProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

////////  ZeroDistributionBoundaryProcessor3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
ZeroDistributionBoundaryProcessor3D<T,DESCRIPTOR>::
ZeroDistributionBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_,
                                    int z1_, int discreteNormalX, int discreteNormalY, int discreteNormalZ)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);

  resetPop[0] = 0;
  for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
    resetPop[iPop] = 0;
    // find incoming iPop from material 0
    if (descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY + descriptors::c<DESCRIPTOR>(iPop,2)*discreteNormalZ > 0) {
      resetPop[iPop] = 1;
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ZeroDistributionBoundaryProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_,
                 int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
            if (resetPop[iPop]!=0) {
              blockLattice.get(iX,iY,iZ)[iPop] = -descriptors::t<T,DESCRIPTOR>(iPop);
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ZeroDistributionBoundaryProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

////////  ExtFieldBoundaryProcessor3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
ExtFieldBoundaryProcessor3D<T,DESCRIPTOR>::
ExtFieldBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                            int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_, int offset_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_),
    discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_),
    discreteNormalZ(discreteNormalZ_), offset(offset_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);
  par = true;
}

template<typename T, typename DESCRIPTOR>
void ExtFieldBoundaryProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_,
                 int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;

  int off = (par) ? 3 : 0;
  int off2 = (par) ? 0 : 3;

  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          // TODO: Update to use descriptor fields
          T* velNeighbour = &blockLattice.get(iX+discreteNormalX,iY+discreteNormalY,iZ+discreteNormalZ)[offset+off];
          T* velCell = &blockLattice.get(iX,iY,iZ)[offset+off2];
          for (unsigned iD = 0; iD < DESCRIPTOR::d; ++iD) {
            velCell[iD] = velNeighbour[iD];
          }
        }
      }
    }
  }
  par = !par;
}

template<typename T, typename DESCRIPTOR>
void ExtFieldBoundaryProcessor3D<T,DESCRIPTOR>::
process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

////////  ConvectionBoundaryProcessorGenerator3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
ConvectionBoundaryProcessorGenerator3D<T,DESCRIPTOR>::
ConvectionBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_,
                                       int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_),
    discreteNormalZ(discreteNormalZ_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
ConvectionBoundaryProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new ConvectionBoundaryProcessor3D<T,DESCRIPTOR>(this->x0, this->x1, this->y0,
         this->y1, this->z0, this->z1,
         discreteNormalX,
         discreteNormalY,
         discreteNormalZ);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ConvectionBoundaryProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new ConvectionBoundaryProcessorGenerator3D<T,DESCRIPTOR>(this->x0, this->x1,
         this->y0, this->y1, this->z0, this->z1,
         discreteNormalX, discreteNormalY, discreteNormalZ);
}

////////  ZeroDistributionBoundaryProcessorGenerator3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
ZeroDistributionBoundaryProcessorGenerator3D<T,DESCRIPTOR>::
ZeroDistributionBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_,
    int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_),
    discreteNormalZ(discreteNormalZ_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
ZeroDistributionBoundaryProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new ZeroDistributionBoundaryProcessor3D<T,DESCRIPTOR>(this->x0, this->x1,
         this->y0, this->y1, this->z0, this->z1,
         discreteNormalX, discreteNormalY, discreteNormalZ);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ZeroDistributionBoundaryProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new ZeroDistributionBoundaryProcessorGenerator3D<T,DESCRIPTOR>(this->x0,
         this->x1, this->y0, this->y1, this->z0, this->z1,
         discreteNormalX, discreteNormalY, discreteNormalZ);
}

////////  ExtFieldBoundaryProcessorGenerator3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
ExtFieldBoundaryProcessorGenerator3D<T,DESCRIPTOR>::
ExtFieldBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_,
                                     int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_, int offset_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_),
    discreteNormalZ(discreteNormalZ_), offset(offset_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
ExtFieldBoundaryProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new ExtFieldBoundaryProcessor3D<T,DESCRIPTOR>(this->x0, this->x1, this->y0,
         this->y1, this->z0, this->z1, discreteNormalX, discreteNormalY,
         discreteNormalZ, offset);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ExtFieldBoundaryProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new ExtFieldBoundaryProcessorGenerator3D<T,DESCRIPTOR>(this->x0, this->x1,
         this->y0, this->y1, this->z0, this->z1,
         discreteNormalX, discreteNormalY, discreteNormalZ, offset);
}

}  // namespace olb
