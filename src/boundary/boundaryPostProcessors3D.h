/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  Generic collision, which modifies the particle distribution
 *  functions, implemented by Orestis Malaspinas, 2007
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

#ifndef BOUNDARY_POST_PROCESSORS_3D_H
#define BOUNDARY_POST_PROCESSORS_3D_H

#include "core/postProcessing.h"
#include "momentaOnBoundaries.h"
#include "core/blockLattice3D.h"

namespace olb {

/**
* This class computes the skordos BC
* on a plane wall in 3D but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template<typename T, typename DESCRIPTOR, int direction, int orientation>
class PlaneFdBoundaryProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  PlaneFdBoundaryProcessor3D (int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) override;
private:
  template<int deriveDirection>
  void interpolateGradients (
    BlockLattice3D<T,DESCRIPTOR> const& blockLattice,
    T velDeriv[DESCRIPTOR::d], int iX, int iY, int iZ ) const;
private:
  int x0, x1, y0, y1, z0, z1;
};

template<typename T, typename DESCRIPTOR, int direction, int orientation>
class PlaneFdBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  PlaneFdBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
};


/**
* This class computes a convection BC on a flat wall in 2D
*/
template<typename T, typename DESCRIPTOR, int direction, int orientation>
class StraightConvectionBoundaryProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  StraightConvectionBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T* uAv_ = NULL);
  ~StraightConvectionBoundaryProcessor3D() override;
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain ( BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_ , int z0_, int z1_) override;
private:
  int x0, x1, y0, y1, z0, z1;
  T**** saveCell;
  T* uAv;
};

template<typename T, typename DESCRIPTOR, int direction, int orientation>
class StraightConvectionBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  StraightConvectionBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T* uAv_ = NULL);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  T* uAv;
};

/**
* This class computes the skordos BC
* on a convex edge wall in 3D but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template<typename T, typename DESCRIPTOR,
         int plane, int normal1, int normal2>
class OuterVelocityEdgeProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  enum { direction1 = (plane+1)%3, direction2 = (plane+2)%3 };
public:
  OuterVelocityEdgeProcessor3D (
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ );
  int extent() const override
  {
    return 2;
  }
  int extent(int whichDirection) const override
  {
    return 2;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_,
                                int z0_, int z1_ ) override;
private:
  T getNeighborRho(int x, int y, int z, int step1, int step2,
                   BlockLattice3D<T,DESCRIPTOR> const& blockLattice);
  template<int deriveDirection, int orientation>
  void interpolateGradients (
    BlockLattice3D<T,DESCRIPTOR> const& blockLattice,
    T velDeriv[DESCRIPTOR::d], int iX, int iY, int iZ ) const;
private:
  int x0, x1, y0, y1, z0, z1;
};

template<typename T, typename DESCRIPTOR,
         int plane, int normal1, int normal2>
class OuterVelocityEdgeProcessorGenerator3D
  : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  OuterVelocityEdgeProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_,
                                        int z0_, int z1_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
};


template<typename T, typename DESCRIPTOR,
         int xNormal, int yNormal, int zNormal>
class OuterVelocityCornerProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  OuterVelocityCornerProcessor3D(int x_, int y_, int z_);
  int extent() const override
  {
    return 2;
  }
  int extent(int whichDirection) const override
  {
    return 2;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_,
                                int z0_, int z1_ ) override;
private:
  int x,y,z;
};

template<typename T, typename DESCRIPTOR,
         int xNormal, int yNormal, int zNormal>
class OuterVelocityCornerProcessorGenerator3D
  : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  OuterVelocityCornerProcessorGenerator3D(int x_, int y_, int z_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
};

/**
* This class computes a slip BC in 3D
*/

template<typename T, typename DESCRIPTOR>
class SlipBoundaryProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  SlipBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_);
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
                                  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) override;
private:
  int reflectionPop[DESCRIPTOR::q];
  int x0, x1, y0, y1, z0, z1;
};


template<typename T, typename DESCRIPTOR>
class SlipBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  SlipBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  int discreteNormalX;
  int discreteNormalY;
  int discreteNormalZ;
};

/**
* This class computes a partial slip BC in 3D
*/

template<typename T, typename DESCRIPTOR>
class PartialSlipBoundaryProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  PartialSlipBoundaryProcessor3D(T tuner_, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_);
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
                                  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) override;
private:
  int reflectionPop[DESCRIPTOR::q];
  int x0, x1, y0, y1, z0, z1;
  T tuner;
};


template<typename T, typename DESCRIPTOR>
class PartialSlipBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  PartialSlipBoundaryProcessorGenerator3D(T tuner_, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  int discreteNormalX;
  int discreteNormalY;
  int discreteNormalZ;
  T tuner;
};

/// PostProcessor for the wetting boundary condition in the free energy model. This is 
/// required to set rho on the boundary (using the denisty of the neighbouring cell in 
/// direction of inwards facing normal at the boundary), as the coupling between the 
/// lattices requires the calculation of a density gradient.
template<typename T, typename DESCRIPTOR>
class FreeEnergyWallProcessor3D : public LocalPostProcessor3D<T, DESCRIPTOR> {
public:
  FreeEnergyWallProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
          int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_, T addend_);
  int extent() const override
  {
    return 2;
  }
  int extent(int whichDirection) const override
  {
    return 2;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) override;
private:
  int x0, x1, y0, y1, z0, z1;
  int discreteNormalX, discreteNormalY, discreteNormalZ;
  T addend;
};

/// Generator class for the FreeEnergyWall PostProcessor handling the wetting boundary condition.
template<typename T, typename DESCRIPTOR>
class FreeEnergyWallProcessorGenerator3D : public PostProcessorGenerator3D<T, DESCRIPTOR> {
public:
  FreeEnergyWallProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
          int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_, T addend_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  int discreteNormalX;
  int discreteNormalY;
  int discreteNormalZ;
  T addend;
};


/// PostProcessor for the chemical potential boundary condition in the free energy model.
/// The chemical potentials on the boundary are set equal to the chemical potential on the
/// fluid cell normal to the boundary. This is necessary because the coupling between the 
/// lattices requires the calculation of the gradient of the chemical potential.
///
/// It would be preferable if this were implemented as a lattice coupling that ran
/// between the chemical potential and force lattice couplings. However there is no
/// access to the discrete normals in lattice couplings.
template<typename T, typename DESCRIPTOR>
class FreeEnergyChemPotBoundaryProcessor3D : public LocalPostProcessor3D<T, DESCRIPTOR> {
public:
  FreeEnergyChemPotBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
          int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_, int latticeNumber_);
  int extent() const override
  {
    return 2;
  }
  int extent(int whichDirection) const override
  {
    return 2;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) override;
private:
  int x0, x1, y0, y1, z0, z1;
  int discreteNormalX, discreteNormalY, discreteNormalZ;
  int latticeNumber;
};

/// Generator class for the FreeEnergyChemPotBoundary PostProcessor.
template<typename T, typename DESCRIPTOR>
class FreeEnergyChemPotBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T, DESCRIPTOR> {
public:
  FreeEnergyChemPotBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
          int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_, int latticeNumber_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  int discreteNormalX;
  int discreteNormalY;
  int discreteNormalZ;
  int latticeNumber;
};


/// PostProcessor for the density / velocity outflow boundaries in the free energy model.
/// The density / order parameters are prescribed to the outflow nodes such that they obey
/// the local-velocity convective boundary condition given in Lou, Gou, Shi (2013).
template<typename T, typename DESCRIPTOR>
class FreeEnergyConvectiveProcessor3D : public LocalPostProcessor3D<T, DESCRIPTOR> {
public:
  FreeEnergyConvectiveProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
          int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_);
  int extent() const override
  {
    return 2;
  }
  int extent(int whichDirection) const override
  {
    return 2;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) override;
private:
  int x0, x1, y0, y1, z0, z1;
  int discreteNormalX, discreteNormalY, discreteNormalZ;
};

/// Generator class for the FreeEnergyConvective post processor.
template<typename T, typename DESCRIPTOR>
class FreeEnergyConvectiveProcessorGenerator3D : public PostProcessorGenerator3D<T, DESCRIPTOR> {
public:
  FreeEnergyConvectiveProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
          int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  int discreteNormalX;
  int discreteNormalY;
  int discreteNormalZ;
};

}

#endif
