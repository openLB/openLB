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

/** \file
 * Dynamics for a generic 2D block lattice view -- header file.
 */
#ifndef BLOCK_LATTICE_VIEW_2D_H
#define BLOCK_LATTICE_VIEW_2D_H

#include <vector>
#include "blockLatticeStructure2D.h"
#include "geometry/blockGeometry2D.h"
#include "geometry/blockGeometryStatistics2D.h"

namespace olb {

/// A rectangular extract from a given BlockLatticeStructure.
/** Attention:
 * - This class can only be applied to already existing BlockLattices.
 * - The postProcessors of the original BlockLattice are not called any
 *   more. New, appropriate PostProcessors are generated instead.
 */
template<typename T, typename DESCRIPTOR>
class BlockLatticeView2D : public BlockLatticeStructure2D<T,DESCRIPTOR> {
public:
  BlockLatticeView2D(BlockLatticeStructure2D<T,DESCRIPTOR>& originalLattice_);
  BlockLatticeView2D(BlockLatticeStructure2D<T,DESCRIPTOR>& originalLattice_,
                     int x0_, int x1_, int y0_, int y1_);
  ~BlockLatticeView2D() override;
  BlockLatticeView2D(BlockLatticeView2D const& rhs);
  BlockLatticeView2D<T,DESCRIPTOR>& operator=
  (BlockLatticeView2D<T,DESCRIPTOR> const& rhs);
  void swap(BlockLatticeView2D<T,DESCRIPTOR>& rhs);

  Cell<T,DESCRIPTOR>& get(int iX, int iY) override;
  Cell<T,DESCRIPTOR>& get(int latticeR[]) override;
  Cell<T,DESCRIPTOR> const& get(int iX, int iY) const override;
  void initialize() override;
  void defineDynamics (
    int x0_, int x1_, int y0_, int y1_,
    Dynamics<T,DESCRIPTOR>* dynamics ) override;
  void defineDynamics(int iX, int iY, Dynamics<T,DESCRIPTOR>* dynamics ) override;
  Dynamics<T,DESCRIPTOR>* getDynamics(int iX, int iY) override;
  void collide(int x0_, int x1_, int y0_, int y1_) override;
  void collide() override;
  void stream(int x0_, int x1_, int y0_, int y1_) override;
  void collideAndStream(int x0_, int x1_, int y0_, int y1_) override;
  void stream(bool periodic=false) override;
  void collideAndStream(bool periodic=false) override;
  T computeAverageDensity(int x0_, int x1_, int y0_, int y1_) const override;
  T computeAverageDensity() const override;
  void computeStress(int iX, int iY, T pi[util::TensorVal<DESCRIPTOR >::n]) override;
  void stripeOffDensityOffset (
    int x0_, int x1_, int y0_, int y1_, T offset ) override;
  void stripeOffDensityOffset(T offset) override;
  void forAll(int x0_, int x1_, int y0_, int y1_,
                      WriteCellFunctional<T,DESCRIPTOR> const& application) override;
  void forAll(WriteCellFunctional<T,DESCRIPTOR> const& application) override;
  void addPostProcessor (
    PostProcessorGenerator2D<T,DESCRIPTOR> const& ppGen) override;
  void resetPostProcessors() override;
  void postProcess(int x0_, int x1_, int y0_, int y1_) override;
  void postProcess() override;
  void addLatticeCoupling (
    LatticeCouplingGenerator2D<T,DESCRIPTOR> const& lcGen,
    std::vector<SpatiallyExtendedObject2D*> partners ) override;
  void executeCoupling(int x0_, int x1_, int y0_, int y1_) override;
  void executeCoupling() override;
  LatticeStatistics<T>& getStatistics() override;
  LatticeStatistics<T> const& getStatistics() const override;
private:
  BlockLatticeStructure2D<T,DESCRIPTOR>  *originalLattice;
  int                          x0, y0;
};

}  // namespace olb

#endif
