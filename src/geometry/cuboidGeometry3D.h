/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias J. Krause
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
 * The description of a vector of 3D cuboid -- header file.
 */


#ifndef CUBOID_GEOMETRY_3D_H
#define CUBOID_GEOMETRY_3D_H

#include <vector>
#include <fstream>
#include "core/singleton.h"
#include "core/serializer.h"
#include "io/ostreamManager.h"
#include "io/xmlReader.h"
#include "core/vector.h"
#include "geometry/cuboid3D.h"

/// All OpenLB code is contained in this namespace.
namespace olb {


template< typename T> class Cuboid3D;
template <typename T> class IndicatorF3D;
template <typename T> class LoadBalancer;


/// A cuboid geometry represents a voxel mesh
/** A cuboid geometry is given by a number of cuboids. To represent
 * a connected domain it is required that the distance between two
 * neighbouring cuboids is less than the smallest delta of them.
 *
 * By the class, access is provided to the cuboids. Core methods of
 * the class are transforming lattice to physical positions in the
 * corresponding unit systems.
 *
 * WARNING:
 * At the moment there are only cuboids with a constant delta possible
 * and the distance between two neighbouring cuboids must be delta
 * since an interpolation operator in time and space is missing in
 * cuboidNeighbourhood and superLattice.
 *
 * This class is not intended to be derived from.
 */
template<typename T>
class CuboidGeometry3D : public BufferSerializable {

private:
  /// Vector of the cuboids
  std::vector<Cuboid3D<T> > _cuboids;
  /// Cuboid which contains all other cuboids
  Cuboid3D<T> _motherCuboid;
  /// Periodicity flag
  Vector<bool,3> _periodicityOn;
  /// class specific ostream
  mutable OstreamManager clout;

public:
  /// Constructs empty Geometry
  CuboidGeometry3D();
  /// Constructs a cuboid geometry with a cubic shape of size nX times nY times nZ with origin originR=(originX, originY, originZ) that consists of nC cuboids
  CuboidGeometry3D(T originX, T originY, T originZ, T deltaR, int nX, int nY, int nZ, int nC=1);
  /// Constructor with vectors
  CuboidGeometry3D(std::vector<T> origin, T deltaR, std::vector<int> extent, int nC=1);
  /// Constructs a cuboid structure with a uniform spacing of voxelSize which consists of  nC cuboids, the cuboids not needed are removed and too big ones are shrinked
  CuboidGeometry3D(IndicatorF3D<T>& indicatorF, T voxelSize, int nC=1);
  CuboidGeometry3D(std::shared_ptr<IndicatorF3D<T>> indicator_sharedPtrF, T voxelSize, int nC=1);
  /// Constructs a cuboid structure with a uniform spacing of voxelSize which consists of nC cuboids, the cuboids not needed are removed and too big ones are shrinked. Uses an iterative method: Halves the largest child cuboid until the set number of cuboids is reached. Largest cuboid is determined by either volume or weight as choosen in minimizeBy.
  CuboidGeometry3D(IndicatorF3D<T>& indicatorF, T voxelSize, int nC, std::string minimizeBy);
  CuboidGeometry3D(std::shared_ptr<IndicatorF3D<T>> indicator_sharedPtrF, T voxelSize, int nC, std::string minimizeBy);
  /// Destructs CuboidGeometry
  virtual ~CuboidGeometry3D();


  /// Read and write access to a single cuboid
  Cuboid3D<T>& get(int iC);
  /// Read access to a single cuboid
  Cuboid3D<T> const& get(int iC) const;
  /// Returns the smallest cuboid that includes all cuboids of the structure
  Cuboid3D<T> getMotherCuboid();
  Cuboid3D<T> const& getMotherCuboid() const;
  /// Set flag to enable/disable periodicity depending of direction. Be aware that not all directions are true to ensure boundary conditions like for velocity are not disturbed.
  void setPeriodicity(bool periodicityX, bool periodicityY, bool periodicityZ);

  /// Gives for a given point (globX/globY/globZ) the related cuboidID
  /// and _p if the point is not in any of the cuboid _childrenQ
  int get_iC(T globX, T globY, T globZ, int offset = 0) const; //TODO old ones
  int get_iC(Vector<T,3>, int offset = 0) const;
  /// This function checks if the points (globX/globY/globZ) and
  /// (globX + orientationX/delta/globY + orientationY/delta/
  /// globZ + orientationZ/delta) is in a cuboid.
  /// It gives the related cuboidID and _p if the points are
  /// not in any of the cuboids.
  /// abs(orientationX) = abs(orientationY) = abs(orientationY) = 1
  /// must be satisfied
  int get_iC(T globX, T globY, T globZ, int orientationX, int orientationY, int orientationZ) const; //TODO old ones
  /// Returns true and the cuboid number of the nearest lattice position to the given physical position if the physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry
  bool getC(T physR[3], int& iC) const;
  /// Returns true and the cuboid number of the nearest lattice position to the given physical position if the physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry
  bool getC(std::vector<T> physR, int& iC) const; //TODO new one
  /// Returns true and the cuboid number of the nearest lattice position to the given physical position if the physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry
  bool getC(const Vector<T,3>& physR, int& iC) const;
  /// Returns true and the nearest lattice position to the given physical position if the physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry
  bool getLatticeR(int latticeR[4], const T physR[3]) const;
  /// Returns true and the floor lattice position to the given physical position if the physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry
  bool getFloorLatticeR(const std::vector<T>& physR, std::vector<int>& latticeR) const;
  /// Returns true and the floor lattice position to the given physical position if the physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry
  bool getFloorLatticeR(const Vector<T,3>& physR, Vector<int,4>& latticeR) const;
  /// Returns the physical position to the given lattice position respecting periodicity for the overlap nodes which are not in the mother cuboid for the case the flag periodicityOn[iDim]=true if the   physical position is within any of the cuboids with an overlap of 1/2*delta belonging to the cuboid geometry
  void getPhysR(T physR[3], const int& iCglob,  const int& iX, const int& iY, const int& iZ) const;
  /// Returns the physical position to the given lattice position respecting periodicity for the overlap nodes which are not in the mother cuboid for the case the flag periodicityOn[iDim]=true
  void getPhysR(T physR[3], const int latticeR[4]) const;

  /// Stores the iC of the neighbouring cuboids in a vector;
  void getNeighbourhood(int cuboid, std::vector<int>& neighbours, int offset = 0);
  /// Returns the number of cuboids in the structure
  int getNc() const;
  /// Returns the minimum of the ratio nX/NY in the structure
  T getMinRatio() const;
  /// Returns the maximum of the ratio nX/NY in the structure
  T getMaxRatio() const;
  /// Returns the minimum coordinate in the structure
  Vector<T,3> getMinPhysR() const;
  /// Returns the maximum coordinate in the structure
  Vector<T,3> getMaxPhysR() const;
  /// Returns the minimum volume in the structure
  T getMinPhysVolume() const;
  /// Returns the maximum volume in the structure
  T getMaxPhysVolume() const;
  /// Returns the minimum number of nodes in the structure
  size_t getMinLatticeVolume() const;
  /// Returns the maximum number of nodes in the structure
  size_t getMaxLatticeVolume() const;
  /// Returns the minimum number of nodes in the structure inside the indicator
  size_t getMinLatticeWeight() const;
  /// Returns the maximum number of nodes in the structure inside the indicator
  size_t getMaxLatticeWeight() const;
  /// Returns the minimum delta in the structure
  T getMinDeltaR() const;
  /// Returns the maximum delta in the structure
  T getMaxDeltaR() const;

  /// Compares two CuboidGeometries
  bool operator==(CuboidGeometry3D<T>& rhs);

  /// Swaps data from input into object
  void swap(CuboidGeometry3D<T>& rhs);
  /// Swaps the vector of cuboids
  void swapCuboids(std::vector< Cuboid3D<T> >& cuboids);
  /// Replace the vector of cuboids
  void replaceCuboids(std::vector< Cuboid3D<T> >& cuboids);

  /// Sets the number of full cells of each cuboid
  void setWeights(IndicatorF3D<T>& indicatorF);
  /// Resets the cuboid array
  void clearCuboids()
  {
    _cuboids.clear();
  }
  /// Adds a cuboid
  void add(Cuboid3D<T> cuboid);
  /// Splits cuboid iC, removes it and adds p cuboids of same volume
  void split(int iC, int p);
  /// Splits cuboid iC, removes it and adds p cuboids of same weight
  void splitByWeight(int iC, int p, IndicatorF3D<T>& indicatorF);
  /// Removes the cuboid iC
  void remove(int iC);
  /// Removes all cuboids where indicatorF = 0
  void remove(IndicatorF3D<T>& indicatorF);
  /// Shrink cuboid iC so that no empty planes are left
  void shrink(int iC, IndicatorF3D<T>& indicatorF);
  /// Shrink all cuboids so that no empty planes are left
  void shrink(IndicatorF3D<T>& indicatorF);

  /// Number of data blocks for the serializer interface
  size_t getNblock() const override;
  /// Binary size for the serializer interface
  size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

  /// Prints cuboid geometry details
  void print() const;
  /// Prints cuboid geometry details plus details of all cuboids
  void printExtended();

  /// Save CuboidGeometry into an existing XML File
  void writeToExistingFile(std::string completeFileName, LoadBalancer<T>& loadBalancer);
  /// Save CuboidGeometry into XML File
  void writeToFile(std::string fileName, LoadBalancer<T>& loadBalancer);

private:
  /// Helper Function to create cuboid parameters for XML tag
  std::string _cuboidParameters(Cuboid3D<T> const& cub);
};






/// Helper Function to retrieve nData-dimensional std::vector of type S from space separated tag
template<typename S>
std::vector<S> getDataFromTag(XMLreader const& reader, std::string attrName, int nData)
{
  std::vector<S> values(nData, S());
  std::stringstream extstr(reader.getAttribute(attrName));
  for (auto& valueI: values) {
    extstr >> valueI;
  }
  return values;
}


/// Load CuboidGeometry from XML File
template<typename T>
CuboidGeometry3D<T>* createCuboidGeometry(std::string fileName)
{
  OstreamManager clout("saveCuboidGeometry3D");
  std::string fname = singleton::directories().getLogOutDir() + fileName + ".xml";

  XMLreader reader(fname);

  std::vector<T> origin = getDataFromTag<T>(reader["CuboidGeometry"], "origin", 3);
  std::vector<int> extent = getDataFromTag<int>(reader["CuboidGeometry"], "extent", 3);
  T deltaR = getDataFromTag<T>(reader["CuboidGeometry"], "deltaR", 1)[0];
  size_t weight = getDataFromTag<size_t>(reader["CuboidGeometry"], "weight", 1)[0];
  int refinementLevel = getDataFromTag<int>(reader["CuboidGeometry"], "refinementLevel", 1)[0];

  CuboidGeometry3D<T>* cGeo = new CuboidGeometry3D<T> (origin, deltaR, extent);
  cGeo->getMotherCuboid().setRefinementLevel(refinementLevel);
  cGeo->getMotherCuboid().setWeight(weight);
  cGeo->clearCuboids();

  for ( XMLreader* cub: reader["CuboidGeometry"] ) {
    origin = getDataFromTag<T>(*cub, "origin", 3);
    extent = getDataFromTag<int>(*cub, "extent", 3);
    deltaR = getDataFromTag<T>(*cub, "deltaR", 1)[0];
    weight = getDataFromTag<int>(*cub, "weight", 1)[0];
    refinementLevel = getDataFromTag<int>(*cub, "refinementLevel", 1)[0];

    cGeo->add( Cuboid3D<T>(origin, deltaR, extent, refinementLevel) );
    cGeo->get(cGeo->getNc() - 1).setWeight(weight);
  }

  return cGeo;
}

}  // namespace olb

#endif
