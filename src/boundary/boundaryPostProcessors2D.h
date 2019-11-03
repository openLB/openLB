/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

#ifndef FD_BOUNDARIES_2D_H
#define FD_BOUNDARIES_2D_H

#include "core/postProcessing.h"
#include "momentaOnBoundaries.h"
#include "core/blockLattice2D.h"

namespace olb {

/**
* This class computes the skordos BC
* on a flat wall in 2D but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template<typename T, typename DESCRIPTOR, int direction, int orientation>
class StraightFdBoundaryProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  StraightFdBoundaryProcessor2D(int x0_, int x1_, int y0_, int y1_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain ( BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_ ) override;
private:
  template<int deriveDirection>
  void interpolateGradients (
    BlockLattice2D<T,DESCRIPTOR> const& blockLattice,
    T velDeriv[DESCRIPTOR::d], int iX, int iY ) const;
private:
  int x0, x1, y0, y1;
};

template<typename T, typename DESCRIPTOR, int direction, int orientation>
class StraightFdBoundaryProcessorGenerator2D : public PostProcessorGenerator2D<T,DESCRIPTOR> {
public:
  StraightFdBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_);
  PostProcessor2D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator2D<T,DESCRIPTOR>*  clone() const override;
};

/**
* This class computes a convection BC on a flat wall in 2D
*/
template<typename T, typename DESCRIPTOR, int direction, int orientation>
class StraightConvectionBoundaryProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  StraightConvectionBoundaryProcessor2D(int x0_, int x1_, int y0_, int y1_, T* uAv_ = NULL);
  ~StraightConvectionBoundaryProcessor2D() override;
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain ( BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_ ) override;
private:
  int x0, x1, y0, y1;
  T*** saveCell;
  T* uAv;
};

template<typename T, typename DESCRIPTOR, int direction, int orientation>
class StraightConvectionBoundaryProcessorGenerator2D : public PostProcessorGenerator2D<T,DESCRIPTOR> {
public:
  StraightConvectionBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_, T* uAv_ = NULL);
  PostProcessor2D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator2D<T,DESCRIPTOR>*  clone() const override;
private:
  T* uAv;
};

/**
* This class computes a slip BC in 2D
*/

template<typename T, typename DESCRIPTOR>
class SlipBoundaryProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  SlipBoundaryProcessor2D(int x0_, int x1_, int y0_, int y1_, int discreteNormalX_, int discreteNormalY_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain ( BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_ ) override;
private:
  int reflectionPop[DESCRIPTOR::q];
  int x0, x1, y0, y1;
};


template<typename T, typename DESCRIPTOR>
class SlipBoundaryProcessorGenerator2D : public PostProcessorGenerator2D<T,DESCRIPTOR> {
public:
  SlipBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_, int discreteNormalX_, int discreteNormalY_);
  PostProcessor2D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator2D<T,DESCRIPTOR>*  clone() const override;
private:
  int discreteNormalX;
  int discreteNormalY;
};

/**
* This class computes a partial slip BC in 2D
*/

template<typename T, typename DESCRIPTOR>
class PartialSlipBoundaryProcessor2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  PartialSlipBoundaryProcessor2D(T tuner_, int x0_, int x1_, int y0_, int y1_, int discreteNormalX_, int discreteNormalY_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain ( BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_ ) override;
private:
  int reflectionPop[DESCRIPTOR::q];
  int x0, x1, y0, y1;
  T tuner;
};


template<typename T, typename DESCRIPTOR>
class PartialSlipBoundaryProcessorGenerator2D : public PostProcessorGenerator2D<T,DESCRIPTOR> {
public:
  PartialSlipBoundaryProcessorGenerator2D(T tuner_, int x0_, int x1_, int y0_, int y1_, int discreteNormalX_, int discreteNormalY_);
  PostProcessor2D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator2D<T,DESCRIPTOR>*  clone() const override;
private:
  int discreteNormalX;
  int discreteNormalY;
  T tuner;
};

/**
* This class computes the skordos BC in 2D on a convex
* corner but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template<typename T, typename DESCRIPTOR, int xNormal,int yNormal>
class OuterVelocityCornerProcessor2D : public LocalPostProcessor2D<T, DESCRIPTOR> {
public:
  OuterVelocityCornerProcessor2D(int x_, int y_);
  int extent() const override
  {
    return 2;
  }
  int extent(int whichDirection) const override
  {
    return 2;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                int x0_,int x1_,int y0_,int y1_ ) override;
private:
  int x, y;
};

template<typename T, typename DESCRIPTOR, int xNormal,int yNormal>
class OuterVelocityCornerProcessorGenerator2D : public PostProcessorGenerator2D<T, DESCRIPTOR> {
public:
  OuterVelocityCornerProcessorGenerator2D(int x_, int y_);
  PostProcessor2D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator2D<T,DESCRIPTOR>*  clone() const override;
};


/// PostProcessor for the wetting boundary condition in the free energy model. This is 
/// required to set rho on the boundary (using the denisty of the neighbouring cell in 
/// direction of inwards facing normal at the boundary), as the coupling between the 
/// lattices requires the calculation of a density gradient.
template<typename T, typename DESCRIPTOR>
class FreeEnergyWallProcessor2D : public LocalPostProcessor2D<T, DESCRIPTOR> {
public:
  FreeEnergyWallProcessor2D(int x0_, int x1_, int y0_, int y1_,
          int discreteNormalX_, int discreteNormalY_, T addend_);
  int extent() const override
  {
    return 2;
  }
  int extent(int whichDirection) const override
  {
    return 2;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                int x0_,int x1_,int y0_,int y1_ ) override;
private:
  int x0, x1, y0, y1;
  int discreteNormalX, discreteNormalY;
  T addend;
};

/// Generator class for the FreeEnergyWall PostProcessor handling the wetting boundary condition.
template<typename T, typename DESCRIPTOR>
class FreeEnergyWallProcessorGenerator2D : public PostProcessorGenerator2D<T, DESCRIPTOR> {
public:
  FreeEnergyWallProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_,
          int discreteNormalX_, int discreteNormalY_, T addend_);
  PostProcessor2D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator2D<T,DESCRIPTOR>*  clone() const override;
private:
  int discreteNormalX;
  int discreteNormalY;
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
class FreeEnergyChemPotBoundaryProcessor2D : public LocalPostProcessor2D<T, DESCRIPTOR> {
public:
  FreeEnergyChemPotBoundaryProcessor2D(int x0_, int x1_, int y0_, int y1_,
          int discreteNormalX_, int discreteNormalY_, int latticeNumber_);
  int extent() const override
  {
    return 2;
  }
  int extent(int whichDirection) const override
  {
    return 2;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                int x0_,int x1_,int y0_,int y1_ ) override;
private:
  int x0, x1, y0, y1;
  int discreteNormalX, discreteNormalY;
  int latticeNumber;
};

/// Generator class for the FreeEnergyChemPotBoundary PostProcessor.
template<typename T,typename DESCRIPTOR>
class FreeEnergyChemPotBoundaryProcessorGenerator2D : public PostProcessorGenerator2D<T, DESCRIPTOR> {
public:
  FreeEnergyChemPotBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_,
          int discreteNormalX_, int discreteNormalY_, int latticeNumber_);
  PostProcessor2D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator2D<T,DESCRIPTOR>*  clone() const override;
private:
  int discreteNormalX;
  int discreteNormalY;
  int latticeNumber;
};


/// PostProcessor for pressure / velocity outflow boundaries in the free energy model. 
/// The density / order parameters are prescribed to the outflow nodes such that they
/// obey the local-velocity convective boundary condition given in Lou, Gou, Shi (2013).
template<typename T, typename DESCRIPTOR>
class FreeEnergyConvectiveProcessor2D : public LocalPostProcessor2D<T, DESCRIPTOR> {
public:
  FreeEnergyConvectiveProcessor2D(int x0_, int x1_, int y0_, int y1_,
          int discreteNormalX_, int discreteNormalY_);
  int extent() const override
  {
    return 2;
  }
  int extent(int whichDirection) const override
  {
    return 2;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                int x0_,int x1_,int y0_,int y1_ ) override;
private:
  int x0, x1, y0, y1;
  int discreteNormalX, discreteNormalY;
};

/// Generator class for the FreeEnergyConvective post processor.
template<typename T, typename DESCRIPTOR>
class FreeEnergyConvectiveProcessorGenerator2D : public PostProcessorGenerator2D<T, DESCRIPTOR> {
public:
  FreeEnergyConvectiveProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_,
          int discreteNormalX_, int discreteNormalY_);
  PostProcessor2D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator2D<T,DESCRIPTOR>*  clone() const override;
private:
  int discreteNormalX;
  int discreteNormalY;
};

}

#endif
