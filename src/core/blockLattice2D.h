/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2008 Jonas Latt
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
 * The dynamics of a 2D block lattice -- header file.
 */
#ifndef BLOCK_LATTICE_2D_H
#define BLOCK_LATTICE_2D_H

#include <vector>
#include "olbDebug.h"
#include "postProcessing.h"
#include "blockLatticeStructure2D.h"
#include "core/cell.h"
#include "latticeStatistics.h"
#include "serializer.h"

/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T> class BlockGeometryStructure2D;
template<typename T, typename DESCRIPTOR> struct Dynamics;
template<typename T> class BlockIndicatorF2D;

/** A regular lattice for highly efficient 2D LB dynamics.
 * A block lattice contains a regular array of Cell objects and
 * some useful methods to execute the LB dynamics on the lattice.
 *
 * This class is not intended to be derived from.
 */
template<typename T, typename DESCRIPTOR>
class BlockLattice2D final : public BlockLatticeStructure2D<T,DESCRIPTOR>, public Serializable {
public:
  typedef std::vector<PostProcessor2D<T,DESCRIPTOR>*> PostProcVector;
private:
  /// Actual data array
  Cell<T,DESCRIPTOR>      *rawData;
  /// 3D-Array pointing to rawData; grid[iX] points to beginning of y-array in rawData
  Cell<T,DESCRIPTOR>      **grid;
  PostProcVector       postProcessors, latticeCouplings;
#ifdef PARALLEL_MODE_OMP
  LatticeStatistics<T> **statistics;
#else
  LatticeStatistics<T> *statistics;
#endif

public:
  /// Construction of an nx_ by ny_ lattice
  BlockLattice2D(int nx, int ny);
  /// Destruction of the lattice
  ~BlockLattice2D() override;
  /// Copy construction
  BlockLattice2D(BlockLattice2D<T,DESCRIPTOR> const& rhs) = delete;
  // Move constructor
  BlockLattice2D(BlockLattice2D&&) = default;
  /// Copy assignment
  BlockLattice2D& operator=(BlockLattice2D<T,DESCRIPTOR> const& rhs) = delete;
public:
  /// Read/write access to lattice cells
  Cell<T,DESCRIPTOR>& get(int iX, int iY) override
  {
    OLB_PRECONDITION(iX<this->_nx);
    OLB_PRECONDITION(iY<this->_ny);
    return grid[iX][iY];
  }
  /// Read/write access to lattice cells
  Cell<T,DESCRIPTOR>& get(int latticeR[]) override
  {
    return get(latticeR[0], latticeR[1]);
  }
  /// Read only access to lattice cells
  Cell<T,DESCRIPTOR> const& get(int iX, int iY) const override
  {
    OLB_PRECONDITION(iX<this->_nx);
    OLB_PRECONDITION(iY<this->_ny);
    return grid[iX][iY];
  }
  /// Initialize the lattice cells to become ready for simulation
  void initialize() override;
  /// Get the dynamics for the specific point
  Dynamics<T,DESCRIPTOR>* getDynamics(int iX, int iY) override;
  /// Define the dynamics on a rectangular domain
  void defineDynamics (int x0, int x1, int y0, int y1,
                       Dynamics<T,DESCRIPTOR>* dynamics ) override;
  /// Define the dynamics on a lattice site
  void defineDynamics(int iX, int iY, Dynamics<T,DESCRIPTOR>* dynamics) override;
  /// Define the dynamics on a domain described by an indicator
  void defineDynamics(BlockIndicatorF2D<T>& indicator,
                      Dynamics<T,DESCRIPTOR>* dynamics);
  /// Define the dynamics on a domain described by a material number
  void defineDynamics(BlockGeometryStructure2D<T>& blockGeometry, int material,
                      Dynamics<T,DESCRIPTOR>* dynamics);

  /// Apply collision step to a rectangular domain
  void collide(int x0, int x1, int y0, int y1) override;
  /// Apply collision step to the whole domain
  void collide() override;
  /// Apply streaming step to a rectangular domain
  void stream(int x0, int x1, int y0, int y1) override;
  /// Apply streaming step to the whole domain
  void stream(bool periodic=false) override;
  /// Apply first collision, then streaming step to a rectangular domain
  void collideAndStream(int x0, int x1, int y0, int y1) override;
  /// Apply first collision, then streaming step to the whole domain
  void collideAndStream(bool periodic=false) override;
  /// Compute the average density within a rectangular domain
  T computeAverageDensity(int x0, int x1, int y0, int y1) const override;
  /// Compute the average density within the whole domain
  T computeAverageDensity() const override;
  /// Compute components of the stress tensor on the cell.
  void computeStress(int iX, int iY, T pi[util::TensorVal<DESCRIPTOR >::n]) override;
  /// Subtract a constant offset from the density within the whole domain
  void stripeOffDensityOffset (
    int x0, int x1, int y0, int y1, T offset ) override;
  /// Subtract a constant offset from the density within a rect. domain
  void stripeOffDensityOffset(T offset) override;
  /// Apply an operation to all cells of a sub-domain
  void forAll(int x0_, int x1_, int y0_, int y1_,
              WriteCellFunctional<T,DESCRIPTOR> const& application) override;
  /// Apply an operation to all cells
  void forAll(WriteCellFunctional<T,DESCRIPTOR> const& application) override;
  /// Add a non-local post-processing step
  void addPostProcessor (    PostProcessorGenerator2D<T,DESCRIPTOR> const& ppGen ) override;
  /// Clean up all non-local post-processing steps
  void resetPostProcessors() override;
  /// Execute post-processing on a sub-lattice
  void postProcess(int x0_, int x1_, int y0_, int y1_) override;
  /// Execute post-processing steps
  void postProcess() override;
  /// Add a non-local post-processing step which couples together lattices
  void addLatticeCoupling( LatticeCouplingGenerator2D<T,DESCRIPTOR> const& lcGen,
                           std::vector<SpatiallyExtendedObject2D*> partners ) override;
  /// Execute couplings on a sub-lattice
  void executeCoupling(int x0_, int x1_, int y0_, int y1_) override;
  /// Execute couplings
  void executeCoupling() override;
  /// Return a handle to the LatticeStatistics object
  LatticeStatistics<T>& getStatistics() override;
  /// Return a constant handle to the LatticeStatistics object
  LatticeStatistics<T> const& getStatistics() const override;
public:
  /// Apply streaming step to bulk (non-boundary) cells
  void bulkStream(int x0, int x1, int y0, int y1);
  /// Apply streaming step to boundary cells
  void boundaryStream ( int lim_x0, int lim_x1, int lim_y0, int lim_y1, int x0,
                        int x1, int y0, int y1 );
  /// Apply collision and streaming step to bulk (non-boundary) cells
  void bulkCollideAndStream(int x0, int x1, int y0, int y1);


  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

private:
  /// Helper method for memory allocation
  void allocateMemory();
  /// Helper method for memory de-allocation
  void releaseMemory();
  /// Release memory for post processors
  void clearPostProcessors();
  /// Release memory for lattice couplings
  void clearLatticeCouplings();
  void periodicEdge(int x0, int x1, int y0, int y1);
  void makePeriodic();
};

}  // namespace olb

#endif
