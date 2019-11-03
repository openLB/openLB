/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2014 Mathias J. Krause, Peter Weisbrod
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
 * Representation of a parallel 3D geometry -- header file.
 */

/// A super geometry represents a parallel voxel mesh
/** A super geometry consists of a number of block geometries,
 * where the material numbers are stored. It is constructed
 * from a cuboid geometry. All coboids of the cuboid geometry
 * are assigned to block geometries which are extended by an
 * overlap in order to enable efficient parallelisation.
 *
 * By the class access is provided to the material numbers of
 * the mesh. Methods for renaming materials are provided as
 * well as a statistic class.
 *
 * This class is not intended to be derived from.
 */


#ifndef SUPER_GEOMETRY_3D_H
#define SUPER_GEOMETRY_3D_H

#include <vector>
#include <iostream>
#include <string>

#include "geometry/cuboidGeometry3D.h"
#include "geometry/superGeometryStatistics3D.h"
#include "geometry/blockGeometry3D.h"
#include "geometry/blockGeometryView3D.h"
#include "communication/superStructure3D.h"
#include "communication/loadBalancer.h"
#include "functors/analytical/indicator/indicatorF3D.h"
#include "io/ostreamManager.h"
#include "utilities/functorPtr.h"


/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T> class CuboidGeometry3D;
template<typename T> class BlockGeometry3D;
template<typename T> class BlockGeometryView3D;
template<typename T> class IndicatorF3D;
template<typename T> class SuperIndicatorF3D;
template<typename T> class SuperStructure3D;
template<typename T> class SuperGeometryStatistics3D;

template<typename T>
class SuperGeometry3D : public SuperStructure3D<T> {

private:
  /// Vector of block geometries with overlap
  std::vector<BlockGeometry3D<T> > _extendedBlockGeometries;
  /// Vector of block geometries without overlap
  std::vector<BlockGeometryView3D<T> > _blockGeometries;
  /// Statistic class
  SuperGeometryStatistics3D<T> _statistics;
  /// class specific output stream
  mutable OstreamManager clout;

public:
  /// Constructor
  SuperGeometry3D(CuboidGeometry3D<T>& cuboidGeometry,
                  LoadBalancer<T>& lb, int overlap = 2);
  /// Copy constructor
  SuperGeometry3D(SuperGeometry3D const& rhs);
  /// Copy assignment
  SuperGeometry3D& operator=(SuperGeometry3D const& rhs);

  /// Interface for the communicator class: Write access to the memory of the data of the super structure
  bool* operator() (int iCloc, int iX, int iY, int iZ, int iData) override;
  /// Interface for the communicator class: Read only access to the dim of the data of the super structure
  int getDataSize() const override;
  /// Interface for the communicator class: Read only access to the data type dim of the data of the super structure
  int getDataTypeSize() const override;

  /// Write access to the material numbers, error handling: stops the program if data is not available
  int& set(int iCglob, int iXloc, int iYloc, int iZloc); //TODO to be removed set->get, problem: with get calling wrong function
  /// Read only access to the material numbers, error handling: returns 0 if data is not available
  int const&  get(int iCglob, int iXloc, int iYloc, int iZloc) const;
  /// Read only access to the material numbers with global communication to all ranks
  int getAndCommunicate(int iCglob, int iXloc, int iYloc, int iZloc) const;
  /// Write access to the material numbers, error handling: stops the program if data is not available
  int& set(std::vector<int> latticeR); //TODO to be removed set->get, problem: with get calling wrong function
  /// Read only access to the material numbers, error handling: returns 0 if data is not available
  int const& get(std::vector<int> latticeR) const;
  int const& get(const int latticeR[4]) const;

  /// Read only access to the material numbers with global communication to all ranks
  int getAndCommunicate(std::vector<int> latticeR) const;

  /// Transforms a lattice to physical position (SI unites)
  std::vector<T> getPhysR(int iCglob, int iX, int iY, int iZ) const;
  /// Transforms a lattice to physical position (SI unites)
  std::vector<T> getPhysR(std::vector<int> latticeR) const;
  /// Returns the physical position to the given lattice position respecting periodicity for the overlap nodes which are not in the mother cuboid for the case the flag periodicityOn[iDim]=true if the   physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry
  void getPhysR(T physR[3], const int& iCglob,  const int& iX, const int& iY, const int& iZ) const;
  /// Returns the physical position to the given lattice position respecting periodicity for the overlap nodes which are not in the mother cuboid for the case the flag periodicityOn[iDim]=true
  void getPhysR(T physR[3], const int latticeR[4]) const;

  /// Read and write access to a single extended block geometry
  BlockGeometry3D<T>& getExtendedBlockGeometry(int locIC);
  /// Read only access to a single extended block geometry
  BlockGeometry3D<T> const& getExtendedBlockGeometry(int locIC) const;
  /// Read and write access to a single block geometry
  BlockGeometryView3D<T>& getBlockGeometry(int locIC);
  /// Read only access to a single block geometry
  BlockGeometryView3D<T> const& getBlockGeometry(int locIC) const;

  /// Returns the statistics object
  SuperGeometryStatistics3D<T>& getStatistics();
  /// Read and write access to the statistic status flag, update needed = true
  bool& getStatisticsStatus();
  /// Read only access to the statistic status flag, update needed = true
  bool const& getStatisticsStatus() const;
  /// Updates the super geometry at the boundaries if needed and afterwards the statisics if needed
  void updateStatistics(bool verbose=true);

  /// Executes an outer cleaning
  int clean(bool verbose=true);
  int clean(int material, bool verbose=true);
  /// Removes not needed material fluids from the outer domain
  int outerClean(bool verbose=true);
  /// inner cleaning for all boundary types
  int innerClean(bool verbose=true);
  /// inner cleaning for specific boundary types
  int innerClean(int material, bool verbose=true);
  /// check for errors (searches for all outer voxels (=0) with an inner voxel (=1) as a direct neighbour)
  bool checkForErrors(bool verbose=true);

  /// replace one material with another
  void rename(int fromM, int toM);
  /// replace one material that fulfills an indicator functor condition with another
  void rename(int fromM, int toM, FunctorPtr<IndicatorF3D<T>>&& condition);
  /// replace one material with another respecting an offset (overlap)
  void rename(int fromM, int toM, unsigned offsetX, unsigned offsetY, unsigned offsetZ);
  /// renames all voxels of material fromM to toM if the number of voxels given by testDirection is of material testM
  void rename(int fromM, int toM, int testM, std::vector<int> testDirection);
  /// renames all boundary voxels of material fromBcMat to toBcMat if two neighbour voxel in the direction of the discrete normal are fluid voxel with material fluidM in the region where the indicator function is fulfilled
  void rename(int fromBcMat, int toBcMat, int fluidMat, FunctorPtr<IndicatorF3D<T>>&& condition);
  /// Copy a layer of material numbers inside an indicator in a discrete normal direction
  void copyMaterialLayer(IndicatorF3D<T>& condition, int discreteNormal[3], int numberOfLayers=3);

  /// Prints some information about the super geometry
  void print();

  template<template<typename> class Indicator, typename... Args>
  std::unique_ptr<SuperIndicatorF3D<T>> getIndicator(Args&&... args)
  {
    static_assert(std::is_base_of<SuperIndicatorF3D<T>, Indicator<T>>::value,
                  "Indicator to be constructed is SuperIndicatorF3D implementation");

    return std::unique_ptr<SuperIndicatorF3D<T>>(
             new Indicator<T>(*this, std::forward<Args>(args)...)
           );
  }

  /**
   * Returns a material indicator using the given vector of materials
   *
   * \param  materials Materials to be indicated
   * \returns          Unique ownership of the constructed indicator.
   *                   May be stored or passed directly to e.g. defineDynamics
   **/
  std::unique_ptr<SuperIndicatorF3D<T>> getMaterialIndicator(std::vector<int>&& materials);
  /**
   * Returns a material indicator using a single material number
   *
   * \param material Material to be indicated
   * \returns        Unique ownership of the constructed indicator.
   *                 May be stored or passed directly to e.g. defineDynamics
   **/
  std::unique_ptr<SuperIndicatorF3D<T>> getMaterialIndicator(int material);

};

} // namespace olb

#endif
