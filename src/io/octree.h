/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Thomas Henn
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
 * Octree
 */

#ifndef OCTREE_H
#define OCTREE_H

#include <iostream>
#include <set>

#include "core/singleton.h"
#include "core/vector.h"


using namespace olb::util;
/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T>
class STLmesh;

template<typename T>
struct STLtriangle;

template <typename T>
class Octree {
public:
  /*
   * Constructs Octree containing triangles of an STLmesh.
   * \param center Centerpoint
   * \param rad Radius
   * \param mesh STLmesh
   * \param maxDepth Maximal depth of tree
   * \param overlap Triangles within rad+overlap are added to this Octree
   */
  Octree(Vector<T,3> center, T rad, STLmesh<T>* mesh, short maxDepth, T overlap = 0., Octree<T>* parent=nullptr);
  /// Destructor destructs
  ~Octree();
  /// Find the node containing the first param with remaining maxDepth
  Octree<T>* find(const Vector<T,3>&,const int& maxDepth = 0);
  /// Write Octree
  void write(const Vector<T,3>& pt, const std::string no);
  /// Write Octree
  void write(const int, const std::string);
  /// Write Octree
  void write(const std::string);
  /// Test intersection of ray with all triangles in Octree
  /// returns number of intersections
  int testIntersection(const Vector<T,3>& pt,const Vector<T,3>& dir, bool print = false);
  /// Test intersection of ray with all triangles in Octree
  /// q contains point of closest intersection to pt in direction direction
  bool closestIntersection(const Vector<T,3>& pt, const Vector<T,3>& direction, Vector<T,3>& q, T& a);
  /// Test intersection of ray with all triangles in Octree
  /// q contains point of closest intersection to pt in direction direction
  /// tri contains triangle with closest intersection
  bool closestIntersection(const Vector<T,3>& pt, const Vector<T,3>& direction, Vector<T,3>& q, T& a, STLtriangle<T>& tri, const T& rad = 0., bool print = false);
  /// Test intersection of sphere moving along ray with radius rad
  /// q contains point of closest intersection to pt in direction direction
  /// tri contains triangle with closest intersection
  bool closestIntersectionSphere(const Vector<T,3>& pt, const T& rad, const Vector<T,3>& direction, Vector<T,3>& q, T& a, STLtriangle<T>& tri);
  /// It's complicated. Computes intersections of a ray with triangles inside this Octree. Sets _inside depending on value of rayInside and changes rayInside depending on the number of intersections. Also takes into account if the intersections happen before or after the center.
  void checkRay(const Vector<T,3>& pt,const Vector<T,3>& dir, unsigned short& rayInside);
  /// Computes intersection of ray with Octree boundaries
  void intersectRayNode(const Vector<T,3>& pt, const Vector<T,3>& dir, Vector<T,3>& s);
  /// Computes all centerpoints of Octree
  void getCenterpoints(std::vector<std::vector<T> >& pts);
  /// Collectes all leafs
  void getLeafs(std::vector<Octree<T>* >& pts);
  /// Return status of _isLeaf;
  bool isLeaf();
  /// Sets Inside
  inline void setInside(bool ins)
  {
    _inside = ins;
  };
  /// Gets Inside
  inline bool getInside()
  {
    return _inside;
  };
  /// Gets _boundarNode
  inline bool getBoundaryNode()
  {
    return _boundaryNode;
  };
  /// Gets Maxdepth
  inline int getMaxdepth() const
  {
    return _maxDepth;
  };
  /// Gets numbers of triangles contained by this Octree
  inline const std::vector<unsigned int>& getTriangles() const
  {
    return _triangles;
  };
  /// Gets centerpoint
  inline const Vector<T,3>& getCenter() const
  {
    return _center;
  };
  /// Gets radius
  inline const T getRadius() const
  {
    return _radius;
  };
  /// Prints console output
  void print();
  /// Returns set of indices of all triangles in nodes containing a line.
  void trianglesOnLine(const Vector<T,3>& pt1, const Vector<T,3>& pt2, std::set<unsigned int>& tris);
  /// Returns reference to _mesh
  inline STLmesh<T>* getMesh()
  {
    return _mesh;
  }

protected:
  ///_vector _triangles contains number of triangles
  std::vector<unsigned int> _triangles;
  Vector<T,3> _center;
  T _radius;
  STLmesh<T>* _mesh;
  short _maxDepth;
  bool _isLeaf;
  bool _boundaryNode;
  bool _inside;
  Octree<T> *_parent;
  Octree<T> **_child;
  void findTriangles(T overlap = 0.);
  bool AABBTri(const STLtriangle<T>& tri, T overlap = 0.);
};


}
#endif
