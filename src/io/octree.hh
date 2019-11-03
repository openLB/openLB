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
 * Octree - adapted from http://www.flipcode.com/archives/Octree_Implementation.shtml
 */

#ifndef OCTREE_HH
#define OCTREE_HH

#include <iostream>

#include "core/singleton.h"

using namespace olb::util;
/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T>
class STLmesh;

template<typename T>
struct STLtriangle;

template <typename T>
class Octree;

template <typename T>
Octree<T>::Octree(Vector<T,3> center, T rad, STLmesh<T>* mesh, short maxDepth, T overlap, Octree<T>* parent)
  : _center(center), _radius(rad), _mesh(mesh), _maxDepth(maxDepth),_isLeaf(false), _boundaryNode(false), _inside(false), _parent(parent), _child(nullptr)
{


  findTriangles(overlap);
  //  cout << _triangles.size() << std::endl;
  if (_triangles.size() > 0 && 0 < _maxDepth) {

    _child = new Octree<T>*[8];

    Vector<T,3> tmpCenter = _center;
    T tmpRad = _radius/2.;
    tmpCenter[0] = _center[0] - tmpRad;
    tmpCenter[1] = _center[1] - tmpRad;
    tmpCenter[2] = _center[2] + tmpRad;
    _child[0] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth-1, overlap, this);

    tmpCenter[0] = _center[0] + tmpRad;
    tmpCenter[1] = _center[1] - tmpRad;
    tmpCenter[2] = _center[2] + tmpRad;
    _child[1] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth-1, overlap, this);

    tmpCenter[0] = _center[0] - tmpRad;
    tmpCenter[1] = _center[1] - tmpRad;
    tmpCenter[2] = _center[2] - tmpRad;
    _child[2] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth-1, overlap, this);

    tmpCenter[0] = _center[0] + tmpRad;
    tmpCenter[1] = _center[1] - tmpRad;
    tmpCenter[2] = _center[2] - tmpRad;
    _child[3] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth-1, overlap, this);

    tmpCenter[0] = _center[0] - tmpRad;
    tmpCenter[1] = _center[1] + tmpRad;
    tmpCenter[2] = _center[2] + tmpRad;
    _child[4] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth-1, overlap, this);

    tmpCenter[0] = _center[0] + tmpRad;
    tmpCenter[1] = _center[1] + tmpRad;
    tmpCenter[2] = _center[2] + tmpRad;
    _child[5] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth-1, overlap, this);

    tmpCenter[0] = _center[0] - tmpRad;
    tmpCenter[1] = _center[1] + tmpRad;
    tmpCenter[2] = _center[2] - tmpRad;
    _child[6] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth-1, overlap, this);

    tmpCenter[0] = _center[0] + tmpRad;
    tmpCenter[1] = _center[1] + tmpRad;
    tmpCenter[2] = _center[2] - tmpRad;
    _child[7] = new Octree<T>(tmpCenter, tmpRad, _mesh, _maxDepth-1, overlap, this);

  } else {
    _isLeaf = true;
    if (_triangles.size() > 0 ) {
      _boundaryNode = true;
    }
  }
}

template <typename T>
Octree<T>::~Octree()
{
  if (_maxDepth!=0 && !_isLeaf) {
    for (int i=0; i<8; i++) {
      delete _child[i];
    }
    delete[] _child;
  }
}

template <typename T>
void Octree<T>::findTriangles(T overlap)
{
  if (_parent == nullptr) {
    _triangles.reserve(_mesh->triangleSize());
    for (unsigned int i=0; i<_mesh->triangleSize(); ++i) {
      if (AABBTri(_mesh->getTri(i))) {
        _triangles.push_back(i);
      }
    }
  } else {
    std::vector<unsigned int>::iterator it;
    for (it = _parent->_triangles.begin(); it!=_parent->_triangles.end(); ++it) {
      if (AABBTri(_mesh->getTri(*it), overlap)) {
        _triangles.push_back(*it);
      }
    }
  }
}

template <typename T>
bool Octree<T>::AABBTri(const STLtriangle<T>& tri, T overlap)
{
  std::vector<T> v0(3,T()), v1(3,T()), v2(3,T()), f0(3,T()), f1(3,T()), f2(3,T()), e(3, T());

  /* Test intersection cuboids - triangle
  * Intersection test after Christer Ericson - Real time Collision Detection p.
  * TestTriangleAABB p.171 */
  Vector<T,3> c(_center);
  T eps = 2e-16;

  for (int j=0; j<3; j++) {
    v0[j] = tri.point[0].r[j]-_center[j];
    v1[j] = tri.point[1].r[j]-_center[j];
    v2[j] = tri.point[2].r[j]-_center[j];
    e[j] = _radius*1.01 + overlap; // + std::numeric_limits<T>::epsilon(); // *1.01;
  }
  for (int j=0; j<3; j++) {
    f0[j] = v1[j] - v0[j];
    f1[j] = v2[j] - v1[j];
    f2[j] = v0[j] - v2[j];
  }
  T p0=T(), p1=T(), r=T();
  //test a00
  p0 = v0[2]*v1[1]-v0[1]*v1[2];
  p1 = v2[2]*v1[1]-v2[2]*v0[1]+v0[2]*v2[1]-v1[2]*v2[1];
  r = e[1] * std::fabs(f0[2]) + e[2]*std::fabs(f0[1]);
  T mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r+eps) {
    return false;
  }

  // test a01
  p0 = v0[1]*v1[2]-v0[1]*v2[2]-v0[2]*v1[1]+v0[2]*v2[1];
  p1 = -v1[1]*v2[2]+v1[2]*v2[1];
  r = e[1] * std::fabs(f1[2]) + e[2]*std::fabs(f1[1]);
  mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r+eps) {
    return false;
  }

  // test a02
  p0 = v0[1]*v2[2]-v0[2]*v2[1];
  p1 = v0[1]*v1[2]-v0[2]*v1[1]+v1[1]*v2[2]-v1[2]*v2[1];
  r = e[1]*std::fabs(f2[2]) + e[2]*std::fabs(f2[1]);
  mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r+eps) {
    return false;
  }

  // test a10
  p0 = v0[0]*v1[2]-v0[2]*v1[0];
  p1 = v0[0]*v2[2]-v0[2]*v2[0]-v1[0]*v2[2]+v1[2]*v2[0];
  r = e[0]*std::fabs(f0[2]) + e[2]*std::fabs(f0[0]);
  mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r+eps) {
    return false;
  }

  // test a11
  p0 = -v0[0]*v1[2]+v0[0]*v2[2]+v0[2]*v1[0]-v0[2]*v2[0];
  p1 = v1[0]*v2[2]-v1[2]*v2[0];
  r =  (T)(e[0]*std::fabs(f1[2])+e[2]*std::fabs(f1[0]));
  mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r+eps) {
    return false;
  }

  // test a12
  p0 = -v0[0]*v2[2]+v0[2]*v2[0];
  p1 = -v0[0]*v1[2]+v0[2]*v1[0]-v1[0]*v2[2]+v1[2]*v2[0];
  r = e[0]*std::fabs(f2[2])+e[2]*std::fabs(f2[0]);
  mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r+eps) {
    return false;
  }

  // test a20
  p0 = -v0[0]*v1[1]+v0[1]*v1[0];
  p1 = -v0[0]*v2[1]+v0[1]*v2[0]+v1[0]*v2[1]-v1[1]*v2[0];
  r = e[0]*std::fabs(f0[1])+e[1]*std::fabs(f0[0]);
  mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r+eps) {
    return false;
  }

  // test a21
  p0 = v0[0]*v1[1]-v0[0]*v2[1]-v0[1]*v1[0]+v0[1]*v2[0];
  p1 = -v1[0]*v2[1]+v1[1]*v2[0];
  r = e[0]*std::fabs(f1[1])+e[1]*std::fabs(f1[0]);
  mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r+eps) {
    return false;
  }

  // test a22
  p0 = v0[0]*v2[1]-v0[1]*v2[0];
  p1 = v0[0]*v1[1]-v0[1]*v1[0]+v1[0]*v2[1]-v1[1]*v2[0];
  r = e[0]*std::fabs(f2[1])+e[1]*std::fabs(f2[0]);
  mmm = std::max<T>(-std::max<T>(p0, p1), std::min<T>(p0, p1));
  if (mmm > r+eps) {
    return false;
  }

  if (std::max(std::max(v0[0], v1[0]), v2[0]) < -e[0] || std::min(std::min(v0[0], v1[0]), v2[0]) > e[0]) {
    return false;
  }
  if (std::max(std::max(v0[1], v1[1]), v2[1]) < -e[1] || std::min(std::min(v0[1], v1[1]), v2[1]) > e[1]) {
    return false;
  }
  if (std::max(std::max(v0[2], v1[2]), v2[2]) < -e[2] || std::min(std::min(v0[2], v1[2]), v2[2]) > e[2]) {
    return false;
  }

  /* Test intersection cuboids - triangle plane*/
  r = e[0]*std::fabs(tri.normal[0]) + e[1]*std::fabs(tri.normal[1]) + e[2]*std::fabs(tri.normal[2]);
  T s =  tri.normal[0]*c[0] + tri.normal[1]*c[1] + tri.normal[2]*c[2] - tri.d;
  return (fabs(s) <= r);
}

template<typename T>
Octree<T>* Octree<T>::find(const Vector<T,3>& pt,const int& maxDepth)
{
  //  clout << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;

  if (_isLeaf || maxDepth == _maxDepth) {
    if (std::abs(_center[0] - pt[0]) < _radius + std::numeric_limits<T>::epsilon() &&
        std::abs(_center[1] - pt[1]) < _radius + std::numeric_limits<T>::epsilon() &&
        std::abs(_center[2] - pt[2]) < _radius + std::numeric_limits<T>::epsilon()) {
      //       clout << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
      return this;
    } else {
      OstreamManager clout(std::cout, "Octree");
      clout << "Point: " << std::setprecision(10) << pt[0]<< " " <<pt[1]<< " " <<pt[2]<< " " <<std::endl;
      clout << "Center: " << std::setprecision(10) << _center[0] << " " << _center[1] << " " << _center[2] << " " << std::endl;
      clout << "Radius: "  << std::setprecision(10)<< _radius << std::endl;
      //throw std::runtime_error("[Octree->find] Point outside of geometry.");
      return nullptr;
    }
  } else {
    if (pt[0] < _center[0]) {
      if (pt[1] < _center[1]) {
        if (pt[2] < _center[2]) {
          return _child[2]->find(pt, maxDepth);
        } else {
          return _child[0]->find(pt, maxDepth);
        }
      } else {
        if (pt[2] < _center[2]) {
          return _child[6]->find(pt, maxDepth);
        } else {
          return _child[4]->find(pt, maxDepth);
        }
      }
    } else {
      if (pt[1] < _center[1]) {
        if (pt[2] < _center[2]) {
          return _child[3]->find(pt, maxDepth);
        } else {
          return _child[1]->find(pt, maxDepth);
        }
      } else {
        if (pt[2] < _center[2]) {
          return _child[7]->find(pt, maxDepth);
        } else {
          return _child[5]->find(pt, maxDepth);
        }
      }
    }
  }
}

template<typename T>
int Octree<T>::testIntersection(const Vector<T,3>& pt,const Vector<T,3>& dir, bool print)
{
  int intersections = 0;
  Vector<T,3> q;
  std::vector<Vector<T,3> > qs;
  T a;
#ifdef OLB_DEBUG
  if (print) {
    OstreamManager clout(std::cout, "Octree");
    clout << "Center: " << _center[0] << " " << _center[1] << " "<< _center[2] << " Dir: " << dir[0] << " " << dir[1] << " "<< dir[2] << std::endl;
    clout << "Tris: ";
    for (unsigned k=0; k<_triangles.size(); ++k) {
      clout << _triangles[k] << " ";
    }
    clout << std::endl;
  }
#endif
  for (unsigned k=0; k<_triangles.size(); ++k) {
    if (_mesh->getTri(_triangles[k]).testRayIntersect(pt, dir, q, a)) {
      if (std::fabs(_center[0]-q[0]) <= _radius + std::numeric_limits<T>::epsilon() + 1/1000. * _radius && std::fabs(_center[1]-q[1]) <= _radius + std::numeric_limits<T>::epsilon() + 1/1000. * _radius && std::fabs(_center[2]-q[2]) <= _radius + std::numeric_limits<T>::epsilon() + 1/1000. * _radius) {
        bool newpoint=true;
        for (unsigned i=0; i<qs.size(); i++) {
          newpoint = ( !util::nearZero(q[0]-qs[i][0]) || !util::nearZero(q[1]-qs[i][1]) || !util::nearZero(q[2]-qs[i][2]) );
        }
        if (newpoint) {
          qs.push_back(q);
          intersections++;
        }
#ifdef OLB_DEBUG
        if (print) {
          OstreamManager clout(std::cout, "Octree");
          clout << "Q: " << q[0] << " " <<q[1] << " "<<q[2] << " inside "<< intersections << std::endl;
        }
#endif
      }
#ifdef OLB_DEBUG
      else if (print) {
        OstreamManager clout(std::cout, "Octree");
        clout << "Q: " << q[0] << " " <<q[1] << " "<<q[2] << " outside"<<std::endl;
      }
#endif
    }
  }
#ifdef OLB_DEBUG
  if (intersections == 0 && print) {
    OstreamManager clout(std::cout, "Octree");
    clout << "No Intersection found!" << std::endl;
  }
#endif
  return intersections;
}

template<typename T>
void Octree<T>::checkRay(const Vector<T,3>& pt,const Vector<T,3>& dir, unsigned short& rayInside)
{
  unsigned short left=0, right=0;
  Vector<T,3> dirNormed(dir);
  dirNormed.normalize();
  dirNormed*= _radius * 2.;
  Vector<T,3> q;
  std::vector<Vector<T,3> > qs;
  T a = 1.;

  for (unsigned int k=0; k<_triangles.size(); ++k) {
    if (_mesh->getTri(_triangles[k]).testRayIntersect(pt, dirNormed, q, a, 0.) && a<1.) {
      bool newpoint=true;
      for (unsigned int i=0; i<qs.size(); i++) {
        newpoint &= ( !util::nearZero(q[0]-qs[i][0]) || !util::nearZero(q[1]-qs[i][1]) || !util::nearZero(q[2]-qs[i][2]) );
      }
      if (newpoint) {
        qs.push_back(q);
        if (a < .5) {
          left++;
        } else {
          right++;
        }
      }
    }
  }
  rayInside += left;
  rayInside %=2;
  setInside(rayInside);
  rayInside += right;
  rayInside %=2;
}

template<typename T>
void Octree<T>::getCenterpoints(std::vector<std::vector<T> >& pts)
{
  if (_isLeaf) {
    pts.push_back(_center);
  } else {
    for (int i=0; i<8; i++) {
      _child[i]->getCenterpoints(pts);
    }
  }
}

template<typename T>
void Octree<T>::getLeafs(std::vector<Octree<T>* >& pts)
{
  if (_isLeaf) {
    pts.push_back(this);
  } else {
    for (int i=0; i<8; i++) {
      _child[i]->getLeafs(pts);
    }
  }
}

template<typename T>
bool Octree<T>::isLeaf()
{
  return _isLeaf;
}

template<typename T>
void Octree<T>::write(const Vector<T,3>& pt,const std::string no)
{
  if (_triangles.size()>0 && (std::fabs(pt[0]-_center[0]) < _radius && std::fabs(pt[1]-_center[1]) < _radius && std::fabs(pt[2]-_center[2]) < _radius)) {
    std::string fullName = singleton::directories().getVtkOutDir() + "Octree_" + no + ".stl";
    std::ofstream f(fullName.c_str());
    if (!f) {
      std::cerr << "[Octree] could not open file: " << fullName << std::endl;
    }
    f << "solid ascii" << std::endl << std::flush;
    std::vector<unsigned int>::iterator it = _triangles.begin();
    for (; it != _triangles.end(); ++it) {
      f << "facet normal" << _mesh->getTri(*it).normal[0] << " "  << _mesh->getTri(*it).normal[1] << " "  << _mesh->getTri(*it).normal[2] << " " <<std::endl;
      f << "    outer loop\n";
      f << "        vertex " << _mesh->getTri(*it).point[0].r[0] << " " << _mesh->getTri(*it).point[0].r[1] << " " << _mesh->getTri(*it).point[0].r[2] << "\n";
      f << "        vertex " << _mesh->getTri(*it).point[1].r[0] << " " << _mesh->getTri(*it).point[1].r[1] << " " << _mesh->getTri(*it).point[1].r[2] << "\n";
      f << "        vertex " << _mesh->getTri(*it).point[2].r[0] << " " << _mesh->getTri(*it).point[2].r[1] << " " << _mesh->getTri(*it).point[2].r[2] << "\n";
      f << "    endloop\n";
      f << "endfacet\n";
    }
    f.close();
  }
  if (!_isLeaf) {
    for (int i=0; i<8; i++) {
      std::stringstream istr;
      istr << i;
      _child[i]->write(pt, no+istr.str());
    }
  }
}

template<typename T>
void Octree<T>::write(const int depth,const std::string no)
{
  if (_triangles.size()>0 && _maxDepth == depth) {
    std::string fullName = singleton::directories().getVtkOutDir() + "Octree_" + no + ".stl";
    std::ofstream f(fullName.c_str());
    if (!f) {
      std::cerr << "[Octree] could not open file: " << fullName << std::endl;
    }
    f << "solid ascii" << std::endl << std::flush;
    std::vector<unsigned int>::iterator it = _triangles.begin();
    for (; it != _triangles.end(); ++it) {
      f << "facet normal" << _mesh->getTri(*it).normal[0] << " "  << _mesh->getTri(*it).normal[1] << " "  << _mesh->getTri(*it).normal[2] << " " <<std::endl;
      f << "    outer loop\n";
      f << "        vertex " << _mesh->getTri(*it).point[0].r[0] << " " << _mesh->getTri(*it).point[0].r[1] << " " << _mesh->getTri(*it).point[0].r[2] << "\n";
      f << "        vertex " << _mesh->getTri(*it).point[1].r[0] << " " << _mesh->getTri(*it).point[1].r[1] << " " << _mesh->getTri(*it).point[1].r[2] << "\n";
      f << "        vertex " << _mesh->getTri(*it).point[2].r[0] << " " << _mesh->getTri(*it).point[2].r[1] << " " << _mesh->getTri(*it).point[2].r[2] << "\n";
      f << "    endloop\n";
      f << "endfacet\n";
    }
    f.close();
  }
  if (!_isLeaf) {
    for (int i=0; i<8; i++) {
      std::stringstream istr;
      istr << i;
      _child[i]->write(depth, no+istr.str());
    }
  }
}


template<typename T>
void Octree<T>::write(const std::string fName)
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif
  if (rank==0) {
    std::vector<Octree<T>* > leafs;
    getLeafs(leafs);
    typename std::vector<Octree<T>* >::iterator it = leafs.begin();

    std::string fullName = singleton::directories().getVtkOutDir() + fName + ".vtk";
    std::ofstream f(fullName.c_str());
    if (!f) {
      std::cerr << "[Particles3D] could not open file: " << fullName << std::endl;
    }
    f << "# vtk DataFile Version 2.0" << std::endl << std::flush;
    f << "Octree"<< std::endl << std::flush;
    f << "ASCII"<< std::endl << std::flush;
    f << "DATASET UNSTRUCTURED_GRID"<< std::endl << std::flush;
    std::stringstream points;
    std::stringstream cells;
    std::stringstream cell_types;
    // std::stringstream point_data;
    std::stringstream cell_data;
    std::stringstream cell_leaf;
    std::stringstream cell_boundary;

    points << "POINTS " << leafs.size()*8 << " float" << std::endl;
    cells << "CELLS " << leafs.size() << " " << leafs.size()*9 << std::endl;
    cell_types << "CELL_TYPES " << leafs.size() << std::endl;
    cell_data << "CELL_DATA " << leafs.size() << std::endl;
    cell_data << "SCALARS insideout int" << std::endl;
    cell_data << "LOOKUP_TABLE default" << std::endl;
    cell_leaf << "SCALARS leaf int" << std::endl;
    cell_leaf << "LOOKUP_TABLE default" << std::endl;
    cell_boundary << "SCALARS boundary int" << std::endl;
    cell_boundary << "LOOKUP_TABLE default" << std::endl;

    Vector<T,3> center;
    Vector<T,3> pt;

    T rad;
    int i=0;
    for (; it != leafs.end(); ++it) {
      center  = (*it)->getCenter();
      rad = (*it)->getRadius();

      pt[0] = center[0] - rad;
      pt[1] = center[1] - rad;
      pt[2] = center[2] - rad;
      points << pt[0]<< " " << pt[1]<< " " << pt[2] << " ";

      pt[0] = center[0] + rad;
      pt[1] = center[1] - rad;
      pt[2] = center[2] - rad;
      points << pt[0]<< " " << pt[1]<< " " << pt[2] << " ";

      pt[0] = center[0] - rad;
      pt[1] = center[1] + rad;
      pt[2] = center[2] - rad;
      points << pt[0]<< " " << pt[1]<< " " << pt[2] << " ";

      pt[0] = center[0] + rad;
      pt[1] = center[1] + rad;
      pt[2] = center[2] - rad;
      points << pt[0]<< " " << pt[1]<< " " << pt[2] << " ";

      pt[0] = center[0] - rad;
      pt[1] = center[1] - rad;
      pt[2] = center[2] + rad;
      points << pt[0]<< " " << pt[1]<< " " << pt[2] << " ";

      pt[0] = center[0] + rad;
      pt[1] = center[1] - rad;
      pt[2] = center[2] + rad;
      points << pt[0]<< " " << pt[1]<< " " << pt[2] << " ";

      pt[0] = center[0] - rad;
      pt[1] = center[1] + rad;
      pt[2] = center[2] + rad;
      points << pt[0]<< " " << pt[1]<< " " << pt[2] << " ";

      pt[0] = center[0] + rad;
      pt[1] = center[1] + rad;
      pt[2] = center[2] + rad;
      points << pt[0]<< " " << pt[1]<< " " << pt[2] << " "  << std::endl;

      cells << "8 ";
      for (int j=0; j<8; j++) {
        cells << i+j << " ";
      }
      i+=8;
      cells << std::endl;

      cell_types << 11 << std::endl;

      cell_data << (*it)->getInside() << " "<< std::endl;
      cell_leaf << (*it)->getMaxdepth() << " "<< std::endl;
      cell_boundary << (*it)->getBoundaryNode() << " " << std::endl;
    }

    f << points.str() << cells.str() << cell_types.str() << cell_data.str() << cell_leaf.str() << cell_boundary.str();

    // f << "POINT_DATA 0\nCELL_DATA 0\n" << std::endl;
    f.close();
  }
  OstreamManager clout(std::cout, "Octree");
  /*if (_verbose)*/ clout << "Write ... OK" << std::endl;
}

template<typename T>
bool Octree<T>::closestIntersectionSphere(const Vector<T,3>& pt, const T& rad, const Vector<T,3>& direction, Vector<T,3>& q, T& a, STLtriangle<T>& tri)
{
  a = std::numeric_limits<T>::infinity();
  T alpha = T();
  std::vector<T> qtmp(3, T());
  bool found = false;
  for (unsigned int k=0; k<_triangles.size(); ++k) {
    if (_mesh->getTri(_triangles[k]).testMovingSphereIntersect(pt, rad, direction, qtmp, alpha)) {
      if (alpha < a) {
        a = alpha;
        q = qtmp;
        found = true;
        tri = _mesh->getTri(_triangles[k]);
      }
    }
  }
  return found;
}

template<typename T>
bool Octree<T>::closestIntersection(const Vector<T,3>& pt, const Vector<T,3>& direction, Vector<T,3>& q, T& a, STLtriangle<T>& tri, const T& rad, bool print)
{
  a = std::numeric_limits<T>::infinity();
  T alpha = T();
  Vector<T,3> qtmp;
  bool found = false;
#ifdef OLB_DEBUG
  if (print) {
    std::cout << "Tri Size: " << _triangles.size() << std::endl;
  }
#endif

  for (unsigned int k=0; k<_triangles.size(); ++k) {
    if (_mesh->getTri(_triangles[k]).testRayIntersect(pt, direction, qtmp, alpha, rad)) {
      if (print) {
        std::cout << "Found intersection!" << std::endl;
      }
      if (alpha < a) {
        a = alpha;
        q = qtmp;
        found = true;
        tri = _mesh->getTri(_triangles[k]);
      }
    }
  }
  //  std::cout << a << std::endl;
  return found;
}

template<typename T>
bool Octree<T>::closestIntersection(const Vector<T,3>& pt, const Vector<T,3>& direction, Vector<T,3>& q, T& a)
{
  STLtriangle<T> tri;
  return closestIntersection(pt, direction, q, a, tri, 0.);
}

template<typename T>
void Octree<T>::intersectRayNode(const Vector<T,3>& pt, const Vector<T,3>& dir, Vector<T,3>& s)
{
  T t,d;
  s*=T();
  //Plane Normals outside

  if (dir[0] > 0.) {
    // n = {1, 0, 0}
    d = _center[0] + _radius;
    t = (d - pt[0])/dir[0];
    s = pt + t*dir;
    if (std::fabs(s[1]-_center[1]) < _radius && std::fabs(s[2]-_center[2]) < _radius) {
      return;
    }
  } else if (dir[0] < 0.) {
    // n = {-1, 0, 0}
    d = _center[0] - _radius;
    t = (d - pt[0])/dir[0];
    s = pt + t*dir;
    if (std::fabs(s[1]-_center[1]) < _radius && std::fabs(s[2]-_center[2]) < _radius) {
      return;
    }
  }

  if (dir[1] > 0.) {
    d = _center[1] + _radius;
    t = (d - pt[1])/dir[1];
    s = pt + t*dir;
    if (std::fabs(s[0]-_center[0]) < _radius && std::fabs(s[2]-_center[2]) < _radius) {
      return;
    }
  } else if (dir[1] < 0.) {
    // n = {0, 0, -1}
    d = _center[1] - _radius;
    t = (d - pt[1])/dir[1];
    s = pt + t*dir;
    if (std::fabs(s[0]-_center[0]) < _radius && std::fabs(s[2]-_center[2]) < _radius) {
      return;
    }
  }

  if (dir[2] > 0.) {
    // n = {0, 0, 1}
    d = _center[2] + _radius;
    t = (d - pt[2])/dir[2];
    s = pt + t*dir;
    if (std::fabs(s[0]-_center[0]) < _radius && std::fabs(s[1]-_center[1]) < _radius) {
      return;
    }
  } else if (dir[2] < 0.) {
    // n = {0, 0, -1}
    d = _center[2] - _radius;
    t = (d - pt[2])/dir[2];
    s = pt + t*dir;
    if (std::fabs(s[0]-_center[0]) < _radius && std::fabs(s[1]-_center[1]) < _radius) {
      return;
    }
  }
}

template <typename T>
void Octree<T>::print()
{
  OstreamManager clout(std::cout, "Octree");
  clout << "radius=" << _radius << "; center=(" << _center[0] << "," << _center[1] << "," << _center[2] << ")" << std::endl;
}

template <typename T>
void Octree<T>::trianglesOnLine(const Vector<T,3>& pt1, const Vector<T,3>& pt2, std::set<unsigned int>& tris)
{
  tris.clear();
  std::vector<T> line = pt2-pt1;
  std::vector<T> s = pt1;
  T lineNorm2 = line[0] * line[0] + line[1] * line[1] + line[2] * line[2];
  T dist2 = T();
  Octree<T>* node = NULL;
  int it = 0;
  while (dist2 < lineNorm2 && it < 50) {
    node = find(s);
    tris.insert(node->_triangles.begin(), node->_triangles.end());
    node->intersectRayNode(s, line, s);
    for (int i=0; i<3; i++) {
      s[i] = s[i] + line[i]*_radius*0.001 /* *node->getRadius()*/;
    }
    it++;
    dist2 = (pt1[0]-s[0])*(pt1[0]-s[0]) + (pt1[1]-s[1])*(pt1[1]-s[1]) +  (pt1[2]-s[2])*(pt1[2]-s[2]);
  }


}

}

#endif
