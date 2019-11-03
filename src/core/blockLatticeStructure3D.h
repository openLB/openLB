/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2018 Jonas Latt, Adrian Kummerlaender
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
 * Dynamics for a generic 3D block structure -- header file.
 */
#ifndef BLOCK_LATTICE_STRUCTURE_3D_H
#define BLOCK_LATTICE_STRUCTURE_3D_H

#include <vector>
#include "cell.h"
#include "blockData3D.h"
#include "blockStructure3D.h"
#include "postProcessing.h"
#include "serializer.h"
#include "spatiallyExtendedObject3D.h"
#include "geometry/blockGeometryStructure3D.h"
#include "latticeStatistics.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"


namespace olb {
template<typename T, typename S> class AnalyticalF3D;

template<typename T, typename DESCRIPTOR> struct Dynamics;
template<typename T, typename DESCRIPTOR> class Cell;
template<typename T, typename DESCRIPTOR> struct WriteCellFunctional;
template<typename T> class BlockIndicatorF3D;
template<typename T> class IndicatorSphere3D;


/** BlockLatticeStructure3D is a interface class for defining dynamics on a BlockStructure3D.
 *  The pure virtual methods are wrapped by generic methods.
 */
template<typename T, typename DESCRIPTOR>
class BlockLatticeStructure3D : public BlockStructure3D, public SpatiallyExtendedObject3D {
public:
  BlockLatticeStructure3D(int nx, int ny, int nz) : BlockStructure3D(nx,ny,nz) {};
  ~BlockLatticeStructure3D() override { }
public:
  /// Define rho on a domain described by an indicator
  /**
   * \param indicator Block indicator describing the target domain
   * \param rho       Analytical functor (global)
   **/
  virtual void defineRho(BlockIndicatorF3D<T>& indicator,
                         AnalyticalF3D<T,T>& rho);
  /// Define rho on a domain with a particular material number
  virtual void defineRho(BlockGeometryStructure3D<T>& blockGeometry, int material,
                         AnalyticalF3D<T,T>& rho);

  /// Define u on a domain described by an indicator
  /**
   * \param indicator Block indicator describing the target domain
   * \param u         Analytical functor (global)
   **/
  virtual void defineU(BlockIndicatorF3D<T>& indicator,
                       AnalyticalF3D<T,T>& u);
  /// Define u on a domain with a particular material number
  virtual void defineU(BlockGeometryStructure3D<T>& blockGeometry, int material,
                       AnalyticalF3D<T,T>& u);

  /// Define rho and u on a domain described by an indicator
  /**
   * \param indicator Block indicator describing the target domain
   * \param rho       Analytical functor (global)
   * \param u         Analytical functor (global)
   **/
  virtual void defineRhoU(BlockIndicatorF3D<T>& indicator,
                          AnalyticalF3D<T,T>& rho, AnalyticalF3D<T,T>& u);
  /// Define rho and u on a domain with a particular material number
  virtual void defineRhoU(BlockGeometryStructure3D<T>& blockGeometry, int material,
                          AnalyticalF3D<T,T>& rho, AnalyticalF3D<T,T>& u);

  /// Define a population on a domain described by an indicator
  /**
   * \param indicator Block indicator describing the target domain
   * \param Pop       Analytical functor (global), target dimension DESCRIPTOR::q
   **/
  virtual void definePopulations(BlockIndicatorF3D<T>& indicator,
                                 AnalyticalF3D<T,T>& Pop);
  /// Define a population on a domain with a particular material number
  virtual void definePopulations(BlockGeometryStructure3D<T>& blockGeometry, int material,
                                 AnalyticalF3D<T,T>& Pop);
  /**
   * \param indicator Block indicator describing the target domain
   * \param Pop       Block functor, target dimension DESCRIPTOR::q
   **/
  virtual void definePopulations(BlockIndicatorF3D<T>& indicator,
                                 BlockF3D<T>& Pop);
  /// Define a population on a domain with a particular material number
  virtual void definePopulations(BlockGeometryStructure3D<T>& blockGeometry, int material,
                                 BlockF3D<T>& Pop);

  /// Define a field on a domain described by an indicator
  /**
   * \param indicator     Block indicator describing the target domain
   * \param field         Analytical functor (global)
   **/
  template <typename FIELD>
  void defineField(BlockIndicatorF3D<T>& indicator,
                   AnalyticalF3D<T,T>& field);
  /// Define a field on a domain with a particular material number
  template <typename FIELD>
  void defineField(BlockGeometryStructure3D<T>& blockGeometry, int material,
                   AnalyticalF3D<T,T>& field);
  /// Define a field on a domain described by an analytical indicator
  /**
   * \param indicatorF Domain indicator to be reduced to BlockIndicatorFfromIndicatorF3D
   **/
  template <typename FIELD>
  void defineField(BlockGeometryStructure3D<T>& blockGeometry,
                   IndicatorF3D<T>& indicatorF,
                   AnalyticalF3D<T,T>& field);

  /// Initialize by equilibrium on a domain described by an indicator
  /**
   * \param indicator Block indicator describing the target domain
   * \param rho       Analytical functor (global)
   * \param u         Analytical functor (global)
   **/
  virtual void iniEquilibrium(BlockIndicatorF3D<T>& indicator,
                              AnalyticalF3D<T,T>& rho, AnalyticalF3D<T,T>& u);
  /// Initialize by equilibrium on a domain with a particular material number
  virtual void iniEquilibrium(BlockGeometryStructure3D<T>& blockGeometry, int material,
                              AnalyticalF3D<T,T>& rho, AnalyticalF3D<T,T>& u);

  // pure virtual member functions
  virtual Cell<T,DESCRIPTOR>& get(int iX, int iY, int iZ) =0;
  virtual Cell<T,DESCRIPTOR>& get(const int latticeR[]) =0;
  virtual Cell<T,DESCRIPTOR> const& get(int iX, int iY, int iZ) const =0;

  virtual void initialize() =0;

  /// Define the dynamics on a lattice site
  virtual void defineDynamics(int iX, int iY, int iZ, Dynamics<T,DESCRIPTOR>* dynamics) = 0;
  /// Define the dynamics on a 3D sub-box
  virtual void defineDynamics(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                              Dynamics<T,DESCRIPTOR>* dynamics) = 0;
  /// Define the dynamics on a domain described by an indicator
  virtual void defineDynamics(BlockIndicatorF3D<T>& indicator,
                              Dynamics<T,DESCRIPTOR>* dynamics) = 0;
  /// Define the dynamics on a domain with a particular material number
  virtual void defineDynamics(BlockGeometryStructure3D<T>& blockGeometry, int material,
                              Dynamics<T,DESCRIPTOR>* dynamics) = 0;

  virtual Dynamics<T,DESCRIPTOR>* getDynamics(int iX, int iY, int iZ) = 0;

  virtual void collide(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) =0;
  virtual void collide() =0;

  virtual void stream(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) =0;
  virtual void stream(bool periodic=false) =0;

  virtual void collideAndStream(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) =0;
  virtual void collideAndStream(bool periodic=false) =0;

  virtual T computeAverageDensity(int x0_, int x1_, int y0_, int y1_, int z0_,
                                  int z1_) const =0;
  virtual T computeAverageDensity() const =0;
  virtual void computeStress(int iX, int iY, int iZ, T pi[util::TensorVal<DESCRIPTOR >::n]) = 0;

  virtual void stripeOffDensityOffset(int x0_, int x1_, int y0_, int y1_, int z0_,
                                      int z1_, T offset) =0;
  virtual void stripeOffDensityOffset(T offset) =0;

  virtual void forAll(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                      WriteCellFunctional<T,DESCRIPTOR> const& application) =0;
  virtual void forAll(WriteCellFunctional<T,DESCRIPTOR> const& application) =0;

  virtual void addPostProcessor(PostProcessorGenerator3D<T,DESCRIPTOR> const& ppGen) =0;
  virtual void resetPostProcessors() =0;
  virtual void postProcess(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) =0;
  virtual void postProcess() =0;

  virtual void addLatticeCoupling(LatticeCouplingGenerator3D<T,DESCRIPTOR> const& lcGen,
                                  std::vector<SpatiallyExtendedObject3D*> partners ) =0;
  virtual void executeCoupling(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) =0;
  virtual void executeCoupling() =0;

  virtual LatticeStatistics<T>& getStatistics() =0;
  virtual LatticeStatistics<T> const& getStatistics() const =0;
};


////////// FREE FUNCTIONS //////////

template <typename T, typename DESCRIPTOR>
void setBlockExternalParticleField( BlockGeometryStructure3D<T>& blockGeometry, AnalyticalF3D<T,T>& velocity,
                                    SmoothIndicatorF3D<T,T,true>& sIndicator,
                                    BlockLattice3D<T,DESCRIPTOR>& extendedBlockLattice );

}  // namespace olb

#endif
