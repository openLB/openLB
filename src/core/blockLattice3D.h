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
 * The dynamics of a 3D block lattice -- header file.
 */
#ifndef BLOCK_LATTICE_3D_H
#define BLOCK_LATTICE_3D_H

#include <vector>
#include "olbDebug.h"
#include "core/cell.h"
#include "postProcessing.h"
#include "blockLatticeStructure3D.h"
#include "geometry/blockGeometry3D.h"
#include "latticeStatistics.h"
#include "serializer.h"


/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T> class BlockGeometryStructure3D;
template<typename T> class BlockIndicatorF3D;
template<typename T, typename DESCRIPTOR> struct Dynamics;


/** BlockLattice3D is a regular lattice for highly efficient 3D LB dynamics.
 * A block lattice contains an array of Cell objects and
 * some useful methods to execute the LB dynamics on the lattice.
 *
 * This class is not intended to be derived from.
 */
template<typename T, typename DESCRIPTOR>
class BlockLattice3D final : public BlockLatticeStructure3D<T,DESCRIPTOR>, public Serializable {
public:
  typedef std::vector<PostProcessor3D<T,DESCRIPTOR>*> PostProcVector;

private:
  /// Actual data array
  Cell<T,DESCRIPTOR>      *rawData;
  /// 3D-Array pointing to rawData; grid[iX][iY] points to beginning of z-array in rawData
  Cell<T,DESCRIPTOR>      ***grid;
  PostProcVector       postProcessors, latticeCouplings;
#ifdef PARALLEL_MODE_OMP
  LatticeStatistics<T> **statistics;
#else
  LatticeStatistics<T> *statistics;
#endif
  BlockGeometry3D<T>& geometry_;

public:
  /// Construction of an nx_ by ny_ by nz_ lattice
  BlockLattice3D(int nx_, int ny_, int nz_, BlockGeometry3D<T>& geometry);
  /// Destruction of the lattice
  ~BlockLattice3D() override;
  /// Copy construction
  BlockLattice3D(BlockLattice3D<T,DESCRIPTOR> const& rhs) = delete;
  /// Move constructor
  BlockLattice3D(BlockLattice3D<T,DESCRIPTOR>&&) = default;
  /// Copy assignment
  BlockLattice3D& operator=(BlockLattice3D<T,DESCRIPTOR> const& rhs) = delete;

  /// Read/write access to lattice cells
  Cell<T,DESCRIPTOR>& get(int iX, int iY, int iZ) override
  {
    OLB_PRECONDITION(iX<this->_nx);
    OLB_PRECONDITION(iY<this->_ny);
    OLB_PRECONDITION(iZ<this->_nz);
    return grid[iX][iY][iZ];
  }
  /// Read/write access to lattice cells
  Cell<T,DESCRIPTOR>& get(const int latticeR[]) override
  {
    return get(latticeR[0], latticeR[1], latticeR[2]);
  }
  /// Read only access to lattice cells
  Cell<T,DESCRIPTOR> const& get(int iX, int iY, int iZ) const override
  {
    OLB_PRECONDITION(iX<this->_nx);
    OLB_PRECONDITION(iY<this->_ny);
    OLB_PRECONDITION(iZ<this->_nz);
    return grid[iX][iY][iZ];
  }

  /// Initialize the lattice cells to become ready for simulation
  void initialize() override;

  /// Get the dynamics on a lattice site
  Dynamics<T,DESCRIPTOR>* getDynamics(int iX, int iY, int iZ) override;

  /// Define the dynamics on a lattice site
  void defineDynamics(int iX, int iY, int iZ, Dynamics<T,DESCRIPTOR>* dynamics) override;
  /// Define the dynamics on a 3D sub-box
  void defineDynamics(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                      Dynamics<T,DESCRIPTOR>* dynamics ) override;
  /// Define the dynamics on a domain described by an indicator
  void defineDynamics(BlockIndicatorF3D<T>& indicator,
                      Dynamics<T,DESCRIPTOR>* dynamics) override;
  /// Define the dynamics on a domain with a particular material number
  void defineDynamics(BlockGeometryStructure3D<T>& blockGeometry, int material,
                      Dynamics<T,DESCRIPTOR>* dynamics) override;

  /// Apply collision step to a 3D sub-box
  void collide(int x0_, int x1_, int y0_, int y1_,
               int z0_, int z1_) override;
  /// Apply collision step to the whole domain
  void collide() override;
  /// Apply streaming step to a 3D sub-box
  void stream(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
  /// Apply streaming step to the whole domain
  void stream(bool periodic=false) override;
  /// Apply first collision, then streaming step to a 3D sub-box
  void collideAndStream(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
  /// Apply first collision, then streaming step to the whole domain
  void collideAndStream(bool periodic=false) override;
  /// Compute the average density within a rectangular domain
  T computeAverageDensity(int x0_, int x1_, int y0_, int y1_,
                          int z0_, int z1_ ) const override;
  /// Compute the average density within the whole domain
  T computeAverageDensity() const override;
  /// Compute components of the stress tensor on the cell.
  void computeStress(int iX, int iY, int iZ, T pi[util::TensorVal<DESCRIPTOR >::n]) override;
  /// Subtract a constant offset from the density within the whole domain
  void stripeOffDensityOffset (
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T offset ) override;
  /// Subtract a constant offset from the density within a rect. domain
  void stripeOffDensityOffset(T offset) override;
  /// Apply an operation to all cells of a sub-domain
  void forAll(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
              WriteCellFunctional<T,DESCRIPTOR> const& application) override;
  /// Apply an operation to all cells
  void forAll(WriteCellFunctional<T,DESCRIPTOR> const& application) override;
  /// Add a non-local post-processing step
  void addPostProcessor (
    PostProcessorGenerator3D<T,DESCRIPTOR> const& ppGen ) override;
  /// Clean up all non-local post-processing steps
  void resetPostProcessors() override;
  /// Execute post-processing on a sub-lattice
  void postProcess(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
  /// Execute post-processing steps
  void postProcess() override;
  /// Add a non-local post-processing step
  void addLatticeCoupling (
    LatticeCouplingGenerator3D<T,DESCRIPTOR> const& lcGen,
    std::vector<SpatiallyExtendedObject3D*> partners ) override;
  /// Execute couplings on a sub-lattice
  void executeCoupling(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
  /// Execute couplings steps
  void executeCoupling() override;
  /// Return a handle to the LatticeStatistics object
  LatticeStatistics<T>& getStatistics() override;
  /// Return a constant handle to the LatticeStatistics object
  LatticeStatistics<T> const& getStatistics() const override;

  /// Apply streaming step to bulk (non-boundary) cells
  void bulkStream(int x0, int x1, int y0, int y1, int z0, int z1);
  /// Apply streaming step to boundary cells
  void boundaryStream (
    int lim_x0, int lim_x1, int lim_y0, int lim_y1,
    int lim_z0, int lim_z1,
    int x0, int x1, int y0, int y1, int z0, int z1 );
  /// Apply collision and streaming step to bulk (non-boundary) cells
  void bulkCollideAndStream(int x0, int x1, int y0, int y1, int z0, int z1);

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

  /// Get material number of cell
  int getMaterial(int iX, int iY, int iZ);

private:
  /// Helper method for memory allocation
  void allocateMemory();
  /// Helper method for memory de-allocation
  void releaseMemory();
  /// Release memory for post processors
  void clearPostProcessors();
  /// Release memory for post processors
  void clearLatticeCouplings();
  /// Make the lattice periodic in all directions
  void makePeriodic();

  void periodicSurface(int x0, int x1, int y0, int y1, int z0, int z1);

};

}  // namespace olb

#endif
