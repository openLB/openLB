/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Mathias J. Krause
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
 * Representation of a 2d block geometry structure -- header file.
 */

#ifndef BLOCK_GEOMETRY_STRUCTURE_2D_H
#define BLOCK_GEOMETRY_STRUCTURE_2D_H

#include <vector>
#include <list>

#include "functors/analytical/indicator/indicatorF2D.h"
#include "geometry/blockGeometryStatistics2D.h"
#include "io/ostreamManager.h"


/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T> class BlockGeometryStatistics2D;
template<typename T> class BlockIndicatorF2D;

/// Representation of a block geometry structure
/** This pure virtual class provides an interface for the classes
 * block geometry and block geometry view. It presents a volume
 * of voxels. Different types are given my material numbers
 * which is imporant e.g. to work with different boundaries
 * (like for inflow/output regions). Several methods are
 * provided to rename material numbers of a type at several
 * nodes at once. Further clearing methods are provided to
 * rename material in order to obtain a robust mesh as well
 * as a method which checks the mesh for possible errors. The
 * member variable statistics provides integral information
 * about the mesh and their materials.
 *
 * This class is intended to be derived from.
 */
template<typename T>
class BlockGeometryStructure2D {

protected:
  /// Number of the cuboid, default=-1
  int _iCglob;
  /// Statistic class
  BlockGeometryStatistics2D<T> _statistics;

  /// class specific output stream
  mutable OstreamManager clout;

public:
  /// Constructor
  BlockGeometryStructure2D(int iCglob=-1);
  /// Destructor
  virtual ~BlockGeometryStructure2D() {};

  /// Returns the underlying block structure
  virtual BlockStructure2D& getBlockStructure() = 0;

  /// Read only access to the global iC number which is given !=-1 if the block geometries are part of a super geometry
  virtual int const& getIcGlob() const;
  /// Read and write access to the statistic object
  virtual BlockGeometryStatistics2D<T>& getStatistics(bool verbose=true) = 0;
  /// Read only access to the statistic object
  virtual BlockGeometryStatistics2D<T> const& getStatistics(bool verbose=true) const = 0;

  /// Returns the position of the block origin which is the node (iX=0/iY=0) in physical units (meter)
  virtual Vector<T,2> getOrigin() const = 0;
  /// Returns the extend in x direction of the block in lattice units
  virtual int getNx() const = 0;
  /// Returns the extend in y direction of the block in lattice units
  virtual int getNy() const = 0;
  /// Returns the extend of the block in lattice units
  virtual Vector<int,2> const getExtend() const;
  /// Returns the spacing in physical units (meter)
  virtual const T getDeltaR() const = 0;

  /// Transforms lattice to physical coordinates (wrapped from cuboid geometry)
  virtual void getPhysR(T physR[2], const int& iX, const int& iY) const = 0;
  /// Transforms lattice to physical coordinates (wrapped from cuboid geometry)
  virtual void getPhysR(T physR[2], const int latticeR[2]) const;

  // TODO to be removed old once
  /// returns the (iX,iY) entry in the 2D scalar field
  virtual int getMaterial(int iX, int iY) const = 0;

  /// Write access to a material number
  virtual int& get(int iX, int iY) = 0;
  /// Read only access to a material number
  virtual int const& get(int iX, int iY) const = 0;
  /// Write access to a material number
  virtual int& get(std::vector<int> latticeR);
  /// Read only access to a material number
  virtual int const& get(std::vector<int> latticeR) const;

  /// Changes all materials which are not 0 or 1 to 0 if there is no neighbour with material 1
  virtual int clean(bool verbose=true);
  /// Changes all materials with material 1 to 0 if there is a neighbour with material 0
  virtual int outerClean(bool verbose=true);
  /// Changes all materials which are not 0 or 1 to 1 if there is a non robust constiallation
  virtual int innerClean(bool verbose=true);
  /// Changes all materials with material fromM to 1 if there is a non robust constiallation
  virtual int innerClean(int fromM, bool verbose=true);

  /// Resets all cell materials inside of a domain to 0
  virtual void reset(IndicatorF2D<T>& domain);

  /// Returns true if there is at least one stream direction found for a boundary voxel with coordinates (iX,iY) and if its neighbor is a bulk voxel.
  /**
   * It further fills the vector streamDirection with true if this direction is pointing to a bulk voxel.
   **/
  template<typename V, typename DESCRIPTOR >
  bool findStreamDirections(int iX, int iY, BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator, bool streamDirections[]);
  /// Returns true if there is at least one stream diection found for a voxel with coordinates (iX,iY) of material and if the neighbouring material is one of the bulkMaterials in the list. It further fills the vector streamDirection with true if this direction is pointing to a bulkMaterial.
  template<typename V, typename DESCRIPTOR >
  bool findStreamDirections(int iX, int iY, int material, std::list<int> bulkMaterials, bool streamDirections[]);
  /// Returns the coordinates (iX,iY) of a voxel with a given material number (material) if there exists an neighbourhood of size (offsetX,offsetY) only with voxels of the  given material number
  virtual bool find(int material, unsigned offsetX, unsigned offsetY, int& iX, int& iY);
  /// Returns true if at position (iX,iY) and in a neighbourhood of size (offsetX,offsetY) only voxels with a given material number (material) are there
  virtual bool check(int material, int iX, int iY, unsigned offsetX, unsigned offsetY);
  /// Checks for errors (searches for all outer voxels (=0) with an inner voxel (=1) as a direct neighbour)
  virtual bool checkForErrors(bool verbose=true) const;

  /// Replaces all material numbers (fromM) to another (toM)
  virtual void rename(int fromM, int toM);
  /// Replaces all material numbers (fromM) to another (toM) if an indicator functor condition is fulfilled
  virtual void rename(int fromM, int toM, IndicatorF2D<T>& condition);
  /// Replaces all material numbers (fromM) to another (toM) if all materials in the neighbourhood (iX-offsetX,..,iX,..,ix+offsetX), .. are of the original material number (fromM)
  virtual void rename(int fromM, int toM, unsigned offsetX, unsigned offsetY);
  /// Replaces all material numbers (fromM) to another (toM) if all materials in the neighbourhood (iX+1,iX+2,..,ix+testDirection[0]), .. are of another material number (testM)
  virtual void rename(int fromM, int toM, int testM, std::vector<int> testDirection);
  /// Replaces all material numbers (fromM) to another (toM) if all materials in the neighbourhood (iX+discreteNormal[0],iX+2*discreteNormal[0]), .. are of another material number (testM) and if an indicator functor condition is fulfilled
  virtual void rename(int fromM, int toM, int fluidM, IndicatorF2D<T>& condition, std::vector<int> discreteNormal);
  /// Replaces all material numbers (fromM) to another (toM) if all materials in the neighbourhood (iX+discreteNormal[0],iX+2*discreteNormal[0]), .. are of another material number (fluidM) and if an indicator functor condition is fulfilled, the discreteNormal is computed from all fromM which fulfill the indicator functor condition
  virtual void rename(int fromM, int toM, int fluidM, IndicatorF2D<T>& condition);

  /// Replaces all material numbers (fromM) to another (toM) using a seed point and max. directions indicated by offsetX,Y != 0
  virtual void regionGrowing(int fromM, int toM, int seedX, int seedY, int offsetX, int offsetY, std::map<std::vector<int>, int >* tmp=nullptr);

public:

  /// Adds a pointer to the flag statistic status which is depended on the data of an block geometry
  virtual void addToStatisticsList(bool* statisticStatus) = 0;
  /// Removes a pointer if it exists from the flag statistic status which is depended on the data of an block geometry
  virtual void removeFromStatisticsList(bool* statisticStatus) = 0;
};

} // namespace olb

#endif
