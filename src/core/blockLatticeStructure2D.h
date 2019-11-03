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
 * Dynamics for a generic 2D block structure -- header file.
 */
#ifndef BLOCK_LATTICE_STRUCTURE_2D_H
#define BLOCK_LATTICE_STRUCTURE_2D_H

#include <vector>
#include "cell.h"
#include "blockStructure2D.h"
#include "postProcessing.h"
#include "serializer.h"
#include "spatiallyExtendedObject2D.h"
#include "geometry/blockGeometryStructure2D.h"
#include "latticeStatistics.h"
#include "functors/analytical/analyticalF.h"


namespace olb {

template<typename T, typename DESCRIPTOR> struct Dynamics;
template<typename T, typename DESCRIPTOR> class Cell;
template<typename T, typename DESCRIPTOR> struct WriteCellFunctional;
template<typename T> class BlockIndicatorF2D;


/// An interface to all the variants of (more or less) regular lattices.
template<typename T, typename DESCRIPTOR>
class BlockLatticeStructure2D : public BlockStructure2D, public SpatiallyExtendedObject2D {
public:
  BlockLatticeStructure2D(int nx, int ny) : BlockStructure2D(nx,ny) {};
  ~BlockLatticeStructure2D() override { }
public:
  virtual void defineRho(BlockIndicatorF2D<T>& indicator,
                         AnalyticalF2D<T,T>& rho);
  virtual void defineRho(BlockGeometryStructure2D<T>& blockGeometry, int material,
                         AnalyticalF2D<T,T>& rho);
  virtual void defineU(BlockIndicatorF2D<T>& indicator,
                       AnalyticalF2D<T,T>& u);
  virtual void defineU(BlockGeometryStructure2D<T>& blockGeometry, int material,
                       AnalyticalF2D<T,T>& u);
  virtual void defineRhoU(BlockIndicatorF2D<T>& indicator,
                          AnalyticalF2D<T,T>& rho, AnalyticalF2D<T,T>& u);
  virtual void defineRhoU(BlockGeometryStructure2D<T>& blockGeometry, int material,
                          AnalyticalF2D<T,T>& rho, AnalyticalF2D<T,T>& u);
  virtual void definePopulations(BlockIndicatorF2D<T>& indicator,
                                 AnalyticalF2D<T,T>& Pop);
  virtual void definePopulations(BlockGeometryStructure2D<T>& blockGeometry, int material,
                                 AnalyticalF2D<T,T>& Pop);

  template <typename FIELD>
  void defineField(BlockIndicatorF2D<T>& indicator,
                   AnalyticalF2D<T,T>& field);
  template <typename FIELD>
  void defineField(BlockGeometryStructure2D<T>& blockGeometry, int material,
                   AnalyticalF2D<T,T>& field);
  template <typename FIELD>
  void defineField(BlockGeometryStructure2D<T>& blockGeometry,
                   IndicatorF2D<T>& indicator,
                   AnalyticalF2D<T,T>& field);
  template <typename FIELD>
  void addField(BlockGeometryStructure2D<T>& blockGeometry,
                IndicatorF2D<T>& indicator,
                AnalyticalF2D<T,T>& field);
  template <typename FIELD>
  void addField(BlockGeometryStructure2D<T>& blockGeometry,
                IndicatorF2D<T>& indicator,
                AnalyticalF2D<T,T>& field, AnalyticalF2D<T,T>& porous);
  template <typename FIELD>
  void multiplyField(BlockGeometryStructure2D<T>& blockGeometry,
                     IndicatorF2D<T>& indicator,
                     AnalyticalF2D<T,T>& field);

  virtual void iniEquilibrium(BlockIndicatorF2D<T>& indicator,
                              AnalyticalF2D<T,T>& rho, AnalyticalF2D<T,T>& u);
  virtual void iniEquilibrium(BlockGeometryStructure2D<T>& blockGeometry, int material,
                              AnalyticalF2D<T,T>& rho, AnalyticalF2D<T,T>& u);

  // pure virtual member functions
  virtual Cell<T,DESCRIPTOR>& get(int iX, int iY) =0;
  virtual Cell<T,DESCRIPTOR>& get(int latticeR[]) =0;
  virtual Cell<T,DESCRIPTOR> const& get(int iX, int iY) const =0;
  virtual void initialize() =0;
  virtual void defineDynamics(int x0_, int x1_, int y0_, int y1_,
                              Dynamics<T,DESCRIPTOR>* dynamics ) =0;
  virtual void defineDynamics(int iX, int iY, Dynamics<T,DESCRIPTOR>* dynamics ) =0;
  virtual Dynamics<T,DESCRIPTOR>* getDynamics(int iX, int iY) = 0;
  virtual void collide(int x0_, int x1_, int y0_, int y1_) =0;
  virtual void collide() =0;
  virtual void stream(int x0_, int x1_, int y0_, int y1_) =0;
  virtual void stream(bool periodic=false) =0;
  virtual void collideAndStream(int x0_, int x1_, int y0_, int y1_) =0;
  virtual void collideAndStream(bool periodic=false) =0;
  virtual T computeAverageDensity(int x0_, int x1_, int y0_, int y1_) const =0;
  virtual T computeAverageDensity() const =0;
  virtual void computeStress(int iX, int iY, T pi[util::TensorVal<DESCRIPTOR >::n]) = 0;
  virtual void stripeOffDensityOffset(int x0_, int x1_, int y0_, int y1_,
                                      T offset ) =0;
  virtual void stripeOffDensityOffset(T offset) =0;
  virtual void forAll(int x0_, int x1_, int y0_, int y1_,
                      WriteCellFunctional<T,DESCRIPTOR> const& application) =0;
  virtual void forAll(WriteCellFunctional<T,DESCRIPTOR> const& application) =0;
  virtual void addPostProcessor(PostProcessorGenerator2D<T,DESCRIPTOR> const& ppGen) =0;
  virtual void resetPostProcessors() =0;
  virtual void postProcess(int x0_, int x1_, int y0_, int y1_) =0;
  virtual void postProcess() =0;
  virtual void addLatticeCoupling(LatticeCouplingGenerator2D<T,DESCRIPTOR> const& lcGen,
                                  std::vector<SpatiallyExtendedObject2D*> partners ) =0;
  virtual void executeCoupling(int x0_, int x1_, int y0_, int y1_) =0;
  virtual void executeCoupling() =0;
  virtual LatticeStatistics<T>& getStatistics() =0;
  virtual LatticeStatistics<T> const& getStatistics() const =0;

};

////////// FREE FUNCTIONS //////////

template <typename T, typename DESCRIPTOR>
void setBlockExternalParticleField( BlockGeometryStructure2D<T>& blockGeometry, AnalyticalF2D<T,T>& velocity,
                                    SmoothIndicatorF2D<T,T,true>& sIndicator,
                                    BlockLattice2D<T,DESCRIPTOR>& extendedBlockLattice );


}  // namespace olb

#endif
