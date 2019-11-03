/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2016 Robin Trunk
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  Generic version of the collision, which modifies the particle
 *  distribution functions, by Orestis Malaspinas.
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

#include "core/postProcessing.h"
#include "core/blockLattice3D.h"

namespace olb {

/**
* This class interpolates missing f_i from values
* near the boundary to get a more stable outflow
* condition for the density. It is assumed that the
* next two cells in direction of this f_i have
* viable values.
*/
template<typename T, typename DESCRIPTOR>
class ConvectionBoundaryProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  ConvectionBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_,
                                int z1_, int discreteNormalX_,
                                int discreteNormalY_, int discreteNormalZ_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain ( BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_ , int z0_, int z1_) override;
private:
  int interpolationPop[DESCRIPTOR::q];
  int x0, x1, y0, y1, z0, z1;
};

template<typename T, typename DESCRIPTOR>
class ConvectionBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  ConvectionBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_,
                                         int z0_, int z1_, int discreteNormalX_,
                                         int discreteNormalY_, int discreteNormalZ_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  int discreteNormalX;
  int discreteNormalY;
  int discreteNormalZ;
};

/**
* This class copies missing values in the
* external field from the neighbour in normal direction.
* Therefore it is assumed this neighbour is a fluid cell.
*/
template<typename T, typename DESCRIPTOR>
class ExtFieldBoundaryProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  ExtFieldBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_,
                              int z1_, int discreteNormalX_, int discreteNormalY_,
                              int discreteNormalZ_, int offset_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain ( BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_ , int z0_, int z1_) override;
private:
  int x0, x1, y0, y1, z0, z1;
  int discreteNormalX, discreteNormalY, discreteNormalZ;
  int offset;
  bool par;
};

template<typename T, typename DESCRIPTOR>
class ExtFieldBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  ExtFieldBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_,
                                       int z1_, int discreteNormalX_, int discreteNormalY_,
                                       int discreteNormalZ_, int offset_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  int discreteNormalX, discreteNormalY, discreteNormalZ;
  int offset;
};

/**
* This class resets some values of the distribution
* on the boundary that can have arbitrary values
* to be zero and thus ensures a correct computation
* of the density that is about to leave the domain.
*/
template<typename T, typename DESCRIPTOR>
class ZeroDistributionBoundaryProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  ZeroDistributionBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_,
                                      int z1_, int discreteNormalX_, int discreteNormalY_,
                                      int discreteNormalZ_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain ( BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_ , int z0_, int z1_) override;
private:
  int resetPop[DESCRIPTOR::q];
  int x0, x1, y0, y1, z0, z1;
};

template<typename T, typename DESCRIPTOR>
class ZeroDistributionBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  ZeroDistributionBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_,
      int z0_, int z1_, int discreteNormalX_,
      int discreteNormalY_, int discreteNormalZ_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  int discreteNormalX;
  int discreteNormalY;
  int discreteNormalZ;
};
}

