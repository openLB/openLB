/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2010-2015 Thomas Henn, Mathias J. Krause
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
 * Input in STL format -- header file. - nope
 */

#ifndef STL_READER_HH
#define STL_READER_HH


#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include "core/singleton.h"
#include "communication/mpiManager.h"
#include "octree.hh"
#include "stlReader.h"

using namespace olb::util;

/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T>
void STLtriangle<T>::init()
{
  Vector<T,3> A = point[0].r;
  Vector<T,3> B = point[1].r;
  Vector<T,3> C = point[2].r;
  Vector<T,3> b, c;
  T bb = 0., bc = 0., cc = 0.;

  for (int i = 0; i < 3; i++) {
    b[i] = B[i] - A[i];
    c[i] = C[i] - A[i];
    bb += b[i] * b[i];
    bc += b[i] * c[i];
    cc += c[i] * c[i];
  }

  normal[0] = b[1] * c[2] - b[2] * c[1];
  normal[1] = b[2] * c[0] - b[0] * c[2];
  normal[2] = b[0] * c[1] - b[1] * c[0];

  T norm = sqrt(
             std::pow(normal[0], 2) + std::pow(normal[1], 2) + std::pow(normal[2], 2));
  normal[0] /= norm;
  normal[1] /= norm;
  normal[2] /= norm;

  T D = 1.0 / (cc * bb - bc * bc);
  T bbD = bb * D;
  T bcD = bc * D;
  T ccD = cc * D;

  kBeta = 0.;
  kGamma = 0.;
  d = 0.;

  for (int i = 0; i < 3; i++) {
    uBeta[i] = b[i] * ccD - c[i] * bcD;
    uGamma[i] = c[i] * bbD - b[i] * bcD;
    kBeta -= A[i] * uBeta[i];
    kGamma -= A[i] * uGamma[i];
    d += A[i] * normal[i];
  }
}

template<typename T>
std::vector<T> STLtriangle<T>::getE0()
{
  Vector<T,3> vec;
  vec[0] = point[0].r[0] - point[1].r[0];
  vec[1] = point[0].r[1] - point[1].r[1];
  vec[2] = point[0].r[2] - point[1].r[2];
  return vec;
}

template<typename T>
std::vector<T> STLtriangle<T>::getE1()
{
  Vector<T,3> vec;
  vec[0] = point[0].r[0] - point[2].r[0];
  vec[1] = point[0].r[1] - point[2].r[1];
  vec[2] = point[0].r[2] - point[2].r[2];
  return vec;
}

/* Schnitttest nach
 * http://www.uninformativ.de/bin/RaytracingSchnitttests-76a577a-CC-BY.pdf
 *
 * Creative Commons Namensnennung 3.0 Deutschland
 * http://creativecommons.org/licenses/by/3.0/de/
 *
 * P. Hofmann, 22. August 2010
 *
 */
template<typename T>
bool STLtriangle<T>::testRayIntersect(const Vector<T,3>& pt,
                                      const Vector<T,3>& dir,
                                      Vector<T,3>& q, T& alpha, const T& rad,
                                      bool print)
{
  T rn = 0.;
  Vector<T,3> testPt = pt + rad * normal;
  Vector<T,3> help;

  for (int i = 0; i < 3; i++) {
    rn += dir[i] * normal[i];
  }
#ifdef OLB_DEBUG

  if (print) {
    std::cout << "Pt: " << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;
  }
  if (print)
    std::cout << "testPt: " << testPt[0] << " " << testPt[1] << " " << testPt[2]
              << std::endl;
  if (print)
    std::cout << "PosNeg: "
              << normal[0] * testPt[0] + normal[1] * testPt[1] + normal[2] * testPt[2]
              - d << std::endl;
  if (print)
    std::cout << "Normal: " << normal[0] << " " << normal[1] << " " << normal[2]
              << std::endl;
#endif

  // Schnitttest Flugrichtung -> Ebene
  if (fabs(rn) < std::numeric_limits<T>::epsilon()) {
#ifdef OLB_DEBUG
    if (print) {
      std::cout << "FALSE 1" << std::endl;
    }
#endif
    return false;
  }
  alpha = d - testPt[0] * normal[0] - testPt[1] * normal[1] - testPt[2] * normal[2];
  //  alpha -= testPt[i] * normal[i];
  alpha /= rn;

  // Abstand Partikel Ebene
  if (alpha < -std::numeric_limits<T>::epsilon()) {
#ifdef OLB_DEBUG
    if (print) {
      std::cout << "FALSE 2" << std::endl;
    }
#endif
    return false;
  }
  for (int i = 0; i < 3; i++) {
    q[i] = testPt[i] + alpha * dir[i];
  }
  double beta = kBeta;
  for (int i = 0; i < 3; i++) {
    beta += uBeta[i] * q[i];
  }
#ifdef OLB_DEBUG
  T dist = sqrt(
             pow(q[0] - testPt[0], 2) + pow(q[1] - testPt[1], 2)
             + pow(q[2] - testPt[2], 2));
#endif

  // Schnittpunkt q in der Ebene?
  if (beta < -std::numeric_limits<T>::epsilon()) {
#ifdef OLB_DEBUG

    if (print) {
      std::cout << "FALSE 3 BETA " << beta << " DIST " << dist << std::endl;
    }
#endif
    return false;
  }
  double gamma = kGamma;
  for (int i = 0; i < 3; i++) {
    gamma += uGamma[i] * q[i];
  }
  if (gamma < -std::numeric_limits<T>::epsilon()) {
#ifdef OLB_DEBUG
    if (print) {
      std::cout << "FALSE 4 GAMMA " << gamma << " DIST " << dist << std::endl;
    }
#endif
    return false;
  }
  if (1. - beta - gamma < -std::numeric_limits<T>::epsilon()) {
#ifdef OLB_DEBUG
    if (print)
      std::cout << "FALSE 5 VAL " << 1 - beta - gamma << " DIST " << dist
                << std::endl;
#endif
    return false;
  }
#ifdef OLB_DEBUG
  if (print) {
    std::cout << "TRUE" << " GAMMA " << gamma << " BETA " << beta << std::endl;
  }
#endif
  return true;
}

/**
 * computes closest Point in a triangle to another point.
 * source: Real-Time Collision Detection. Christer Ericson. ISBN-10: 1558607323
 */
template<typename T>
Vector<T,3> STLtriangle<T>::closestPtPointTriangle(
  const Vector<T,3>& pt) const
{

  const T nEps = -std::numeric_limits<T>::epsilon();
  const T Eps = std::numeric_limits<T>::epsilon();

  Vector<T,3> ab = point[1].r - point[0].r;
  Vector<T,3> ac = point[2].r - point[0].r;
  Vector<T,3> bc = point[2].r - point[1].r;

  T snom = (pt - point[0].r)*ab;
  T sdenom = (pt - point[1].r)*(point[0].r - point[1].r);

  T tnom = (pt - point[0].r)*ac;
  T tdenom = (pt - point[2].r)*(point[0].r - point[2].r);

  if (snom < nEps && tnom < nEps) {
    return point[0].r;
  }

  T unom = (pt - point[1].r)*bc;
  T udenom = (pt - point[2].r)*(point[1].r - point[2].r);

  if (sdenom < nEps && unom < nEps) {
    return point[1].r;
  }
  if (tdenom < nEps && udenom < nEps) {
    return point[2].r;
  }

  T vc = normal*crossProduct3D(point[0].r - pt, point[1].r - pt);

  if (vc < nEps && snom > Eps && sdenom > Eps) {
    return point[0].r + snom / (snom + sdenom) * ab;
  }

  T va = normal*crossProduct3D(point[1].r - pt, point[2].r - pt);

  if (va < nEps && unom > Eps && udenom > Eps) {
    return point[1].r + unom / (unom + udenom) * bc;
  }

  T vb = normal*crossProduct3D(point[2].r - pt, point[0].r - pt);

  if (vb < nEps && tnom > Eps && tdenom > Eps) {
    return point[0].r + tnom / (tnom + tdenom) * ac;
  }

  T u = va / (va + vb + vc);
  T v = vb / (va + vb + vc);
  T w = 1. - u - v;

  return u * point[0].r + v * point[1].r + w * point[2].r;
}

template<typename T>
STLmesh<T>::STLmesh(std::string fName, T stlSize)
  : _fName(fName),
    _min(T()),
    _max(T()),
    _maxDist2(0),
    clout(std::cout, "STLmesh")
{
  std::ifstream f(fName.c_str(), std::ios::in);
  _triangles.reserve(10000);
  if (!f.good()) {
    throw std::runtime_error("STL File not valid.");
  }
  char buf[6];
  buf[5] = 0;
  f.read(buf, 5);
  const std::string asciiHeader = "solid";
  if (std::string(buf) == asciiHeader) {
    f.seekg(0, std::ios::beg);
    if (f.good()) {
      std::string s0, s1;
      int i = 0;
      while (!f.eof()) {
        f >> s0;
        if (s0 == "facet") {
          STLtriangle<T> tri;
          f >> s1 >> tri.normal[0] >> tri.normal[1] >> tri.normal[2];
          f >> s0 >> s1;
          f >> s0 >> tri.point[0].r[0] >> tri.point[0].r[1]
            >> tri.point[0].r[2];
          f >> s0 >> tri.point[1].r[0] >> tri.point[1].r[1]
            >> tri.point[1].r[2];
          f >> s0 >> tri.point[2].r[0] >> tri.point[2].r[1]
            >> tri.point[2].r[2];
          f >> s0;
          f >> s0;
          for (int k = 0; k < 3; k++) {
            tri.point[0].r[k] *= stlSize;
            tri.point[1].r[k] *= stlSize;
            tri.point[2].r[k] *= stlSize;
          }
          if (i == 0) {
            _min*=T();
            _max*=T();

            _min[0] = tri.point[0].r[0];
            _min[1] = tri.point[0].r[1];
            _min[2] = tri.point[0].r[2];

            _max[0] = tri.point[0].r[0];
            _max[1] = tri.point[0].r[1];
            _max[2] = tri.point[0].r[2];

            _min[0] = std::min(_min[0], (T) tri.point[1].r[0]);
            _min[1] = std::min(_min[1], (T) tri.point[1].r[1]);
            _min[2] = std::min(_min[2], (T) tri.point[1].r[2]);

            _max[0] = std::max(_max[0], (T) tri.point[1].r[0]);
            _max[1] = std::max(_max[1], (T) tri.point[1].r[1]);
            _max[2] = std::max(_max[2], (T) tri.point[1].r[2]);

            _min[0] = std::min(_min[0], (T) tri.point[2].r[0]);
            _min[1] = std::min(_min[1], (T) tri.point[2].r[1]);
            _min[2] = std::min(_min[2], (T) tri.point[2].r[2]);

            _max[0] = std::max(_max[0], (T) tri.point[2].r[0]);
            _max[1] = std::max(_max[1], (T) tri.point[2].r[1]);
            _max[2] = std::max(_max[2], (T) tri.point[2].r[2]);

          } else {
            _min[0] = std::min(_min[0], (T) tri.point[0].r[0]);
            _min[1] = std::min(_min[1], (T) tri.point[0].r[1]);
            _min[2] = std::min(_min[2], (T) tri.point[0].r[2]);

            _max[0] = std::max(_max[0], (T) tri.point[0].r[0]);
            _max[1] = std::max(_max[1], (T) tri.point[0].r[1]);
            _max[2] = std::max(_max[2], (T) tri.point[0].r[2]);

            _min[0] = std::min(_min[0], (T) tri.point[1].r[0]);
            _min[1] = std::min(_min[1], (T) tri.point[1].r[1]);
            _min[2] = std::min(_min[2], (T) tri.point[1].r[2]);

            _max[0] = std::max(_max[0], (T) tri.point[1].r[0]);
            _max[1] = std::max(_max[1], (T) tri.point[1].r[1]);
            _max[2] = std::max(_max[2], (T) tri.point[1].r[2]);

            _min[0] = std::min(_min[0], (T) tri.point[2].r[0]);
            _min[1] = std::min(_min[1], (T) tri.point[2].r[1]);
            _min[2] = std::min(_min[2], (T) tri.point[2].r[2]);

            _max[0] = std::max(_max[0], (T) tri.point[2].r[0]);
            _max[1] = std::max(_max[1], (T) tri.point[2].r[1]);
            _max[2] = std::max(_max[2], (T) tri.point[2].r[2]);
          }

          i++;
          tri.init();
          _triangles.push_back(tri);

          _maxDist2 = std::max(distPoints(tri.point[0], tri.point[1]),
                               _maxDist2);
          _maxDist2 = std::max(distPoints(tri.point[2], tri.point[1]),
                               _maxDist2);
          _maxDist2 = std::max(distPoints(tri.point[0], tri.point[2]),
                               _maxDist2);
        } else if (s0 == "endsolid") {
          break;
        }
      }
    }
  } else {
    f.close();
    f.open(fName.c_str(), std::ios::in | std::ios::binary);
    char comment[80];
    f.read(comment, 80);

    if (!f.good()) {
      throw std::runtime_error("STL File not valid.");
    }

    comment[79] = 0;
    int32_t nFacets;
    f.read(reinterpret_cast<char *>(&nFacets), sizeof(int32_t));

    if (!f.good()) {
      throw std::runtime_error("STL File not valid.");
    }

    float v[12];
    unsigned short uint16;
    for (int32_t i = 0; i < nFacets; ++i) {
      for (unsigned int j = 0; j < 12; ++j) {
        f.read(reinterpret_cast<char *>(&v[j]), sizeof(float));
      }
      f.read(reinterpret_cast<char *>(&uint16), sizeof(unsigned short));
      STLtriangle<T> tri;
      tri.normal[0] = v[0];
      tri.normal[1] = v[1];
      tri.normal[2] = v[2];
      tri.point[0].r[0] = v[3];
      tri.point[0].r[1] = v[4];
      tri.point[0].r[2] = v[5];
      tri.point[1].r[0] = v[6];
      tri.point[1].r[1] = v[7];
      tri.point[1].r[2] = v[8];
      tri.point[2].r[0] = v[9];
      tri.point[2].r[1] = v[10];
      tri.point[2].r[2] = v[11];

      for (int k = 0; k < 3; k++) {
        tri.point[0].r[k] *= stlSize;
        tri.point[1].r[k] *= stlSize;
        tri.point[2].r[k] *= stlSize;
      }
      if (i == 0) {
        _min[0] = tri.point[0].r[0];
        _min[1] = tri.point[0].r[1];
        _min[2] = tri.point[0].r[2];

        _max[0] = tri.point[0].r[0];
        _max[1] = tri.point[0].r[1];
        _max[2] = tri.point[0].r[2];

        _min[0] = std::min(_min[0], (T) tri.point[1].r[0]);
        _min[1] = std::min(_min[1], (T) tri.point[1].r[1]);
        _min[2] = std::min(_min[2], (T) tri.point[1].r[2]);

        _max[0] = std::max(_max[0], (T) tri.point[1].r[0]);
        _max[1] = std::max(_max[1], (T) tri.point[1].r[1]);
        _max[2] = std::max(_max[2], (T) tri.point[1].r[2]);

        _min[0] = std::min(_min[0], (T) tri.point[2].r[0]);
        _min[1] = std::min(_min[1], (T) tri.point[2].r[1]);
        _min[2] = std::min(_min[2], (T) tri.point[2].r[2]);

        _max[0] = std::max(_max[0], (T) tri.point[2].r[0]);
        _max[1] = std::max(_max[1], (T) tri.point[2].r[1]);
        _max[2] = std::max(_max[2], (T) tri.point[2].r[2]);

      } else {
        _min[0] = std::min(_min[0], (T) tri.point[0].r[0]);
        _min[1] = std::min(_min[1], (T) tri.point[0].r[1]);
        _min[2] = std::min(_min[2], (T) tri.point[0].r[2]);

        _max[0] = std::max(_max[0], (T) tri.point[0].r[0]);
        _max[1] = std::max(_max[1], (T) tri.point[0].r[1]);
        _max[2] = std::max(_max[2], (T) tri.point[0].r[2]);

        _min[0] = std::min(_min[0], (T) tri.point[1].r[0]);
        _min[1] = std::min(_min[1], (T) tri.point[1].r[1]);
        _min[2] = std::min(_min[2], (T) tri.point[1].r[2]);

        _max[0] = std::max(_max[0], (T) tri.point[1].r[0]);
        _max[1] = std::max(_max[1], (T) tri.point[1].r[1]);
        _max[2] = std::max(_max[2], (T) tri.point[1].r[2]);

        _min[0] = std::min(_min[0], (T) tri.point[2].r[0]);
        _min[1] = std::min(_min[1], (T) tri.point[2].r[1]);
        _min[2] = std::min(_min[2], (T) tri.point[2].r[2]);

        _max[0] = std::max(_max[0], (T) tri.point[2].r[0]);
        _max[1] = std::max(_max[1], (T) tri.point[2].r[1]);
        _max[2] = std::max(_max[2], (T) tri.point[2].r[2]);
      }
      tri.init();
      _triangles.push_back(tri);

      _maxDist2 = std::max(distPoints(tri.point[0], tri.point[1]), _maxDist2);
      _maxDist2 = std::max(distPoints(tri.point[2], tri.point[1]), _maxDist2);
      _maxDist2 = std::max(distPoints(tri.point[0], tri.point[2]), _maxDist2);
    }
  }
  f.close();
}

template<typename T>
STLmesh<T>::STLmesh(const std::vector<std::vector<T>> meshPoints, T stlSize)
  : _fName("meshPoints.stl"),
    _min(T()),
    _max(T()),
    _maxDist2(0),
    clout(std::cout, "STLmesh")
{
  _triangles.reserve(10000);
  for (size_t i = 0; i < meshPoints.size() / 3; i++) {
    STLtriangle<T> tri;
    tri.point[0].r[0] = meshPoints[i*3 + 0][0];
    tri.point[0].r[1] = meshPoints[i*3 + 0][1];
    tri.point[0].r[2] = meshPoints[i*3 + 0][2];

    tri.point[1].r[0] = meshPoints[i*3 + 1][0];
    tri.point[1].r[1] = meshPoints[i*3 + 1][1];
    tri.point[1].r[2] = meshPoints[i*3 + 1][2];

    tri.point[2].r[0] = meshPoints[i*3 + 2][0];
    tri.point[2].r[1] = meshPoints[i*3 + 2][1];
    tri.point[2].r[2] = meshPoints[i*3 + 2][2];
    for (int k = 0; k < 3; k++) {
      tri.point[0].r[k] *= stlSize;
      tri.point[1].r[k] *= stlSize;
      tri.point[2].r[k] *= stlSize;
    }
    if (i == 0) {
      _min*=T();
      _max*=T();

      _min[0] = tri.point[0].r[0];
      _min[1] = tri.point[0].r[1];
      _min[2] = tri.point[0].r[2];

      _max[0] = tri.point[0].r[0];
      _max[1] = tri.point[0].r[1];
      _max[2] = tri.point[0].r[2];

      _min[0] = std::min(_min[0], (T) tri.point[1].r[0]);
      _min[1] = std::min(_min[1], (T) tri.point[1].r[1]);
      _min[2] = std::min(_min[2], (T) tri.point[1].r[2]);

      _max[0] = std::max(_max[0], (T) tri.point[1].r[0]);
      _max[1] = std::max(_max[1], (T) tri.point[1].r[1]);
      _max[2] = std::max(_max[2], (T) tri.point[1].r[2]);

      _min[0] = std::min(_min[0], (T) tri.point[2].r[0]);
      _min[1] = std::min(_min[1], (T) tri.point[2].r[1]);
      _min[2] = std::min(_min[2], (T) tri.point[2].r[2]);

      _max[0] = std::max(_max[0], (T) tri.point[2].r[0]);
      _max[1] = std::max(_max[1], (T) tri.point[2].r[1]);
      _max[2] = std::max(_max[2], (T) tri.point[2].r[2]);

    } else {
      _min[0] = std::min(_min[0], (T) tri.point[0].r[0]);
      _min[1] = std::min(_min[1], (T) tri.point[0].r[1]);
      _min[2] = std::min(_min[2], (T) tri.point[0].r[2]);

      _max[0] = std::max(_max[0], (T) tri.point[0].r[0]);
      _max[1] = std::max(_max[1], (T) tri.point[0].r[1]);
      _max[2] = std::max(_max[2], (T) tri.point[0].r[2]);

      _min[0] = std::min(_min[0], (T) tri.point[1].r[0]);
      _min[1] = std::min(_min[1], (T) tri.point[1].r[1]);
      _min[2] = std::min(_min[2], (T) tri.point[1].r[2]);

      _max[0] = std::max(_max[0], (T) tri.point[1].r[0]);
      _max[1] = std::max(_max[1], (T) tri.point[1].r[1]);
      _max[2] = std::max(_max[2], (T) tri.point[1].r[2]);

      _min[0] = std::min(_min[0], (T) tri.point[2].r[0]);
      _min[1] = std::min(_min[1], (T) tri.point[2].r[1]);
      _min[2] = std::min(_min[2], (T) tri.point[2].r[2]);

      _max[0] = std::max(_max[0], (T) tri.point[2].r[0]);
      _max[1] = std::max(_max[1], (T) tri.point[2].r[1]);
      _max[2] = std::max(_max[2], (T) tri.point[2].r[2]);
    }

    tri.init();
    _triangles.push_back(tri);

    _maxDist2 = std::max(distPoints(tri.point[0], tri.point[1]),
                         _maxDist2);
    _maxDist2 = std::max(distPoints(tri.point[2], tri.point[1]),
                         _maxDist2);
    _maxDist2 = std::max(distPoints(tri.point[0], tri.point[2]),
                         _maxDist2);
  }
}

template<typename T>
T STLmesh<T>::distPoints(STLpoint<T>& p1, STLpoint<T>& p2)
{
  return std::pow(double(p1.r[0] - p2.r[0]), 2)
         + std::pow(double(p1.r[1] - p2.r[1]), 2)
         + std::pow(double(p1.r[2] - p2.r[2]), 2);
}

template<typename T>
void STLmesh<T>::print(bool full)
{
  if (full) {
    int i = 0;
    clout << "Triangles: " << std::endl;
    typename std::vector<STLtriangle<T> >::iterator it = _triangles.begin();

    for (; it != _triangles.end(); ++it) {
      clout << i++ << ": " << it->point[0].r[0] << " " << it->point[0].r[1]
            << " " << it->point[0].r[2] << " | " << it->point[1].r[0] << " "
            << it->point[1].r[1] << " " << it->point[1].r[2] << " | "
            << it->point[2].r[0] << " " << it->point[2].r[1] << " "
            << it->point[2].r[2] << std::endl;
    }
  }
  clout << "nTriangles=" << _triangles.size() << "; maxDist2=" << _maxDist2
        << std::endl;
  clout << "minPhysR(StlMesh)=(" << getMin()[0] << "," << getMin()[1] << ","
        << getMin()[2] << ")";
  clout << "; maxPhysR(StlMesh)=(" << getMax()[0] << "," << getMax()[1] << ","
        << getMax()[2] << ")" << std::endl;
}

template<typename T>
void STLmesh<T>::write(std::string fName)
{
  int rank = 0;
#ifdef PARALLEL_MODE_MPI
  rank = singleton::mpi().getRank();
#endif
  if (rank == 0) {
    std::string fullName = singleton::directories().getVtkOutDir() + fName
                           + ".stl";
    std::ofstream f(fullName.c_str());
    f << "solid ascii " << fullName << "\n";

    for (unsigned int i = 0; i < _triangles.size(); i++) {
      f << "facet normal " << _triangles[i].normal[0] << " "
        << _triangles[i].normal[1] << " " << _triangles[i].normal[2] << "\n";
      f << "    outer loop\n";
      f << "        vertex " << _triangles[i].point[0].r[0] << " "
        << _triangles[i].point[0].r[1] << " " << _triangles[i].point[0].r[2]
        << "\n";
      f << "        vertex " << _triangles[i].point[1].r[0] << " "
        << _triangles[i].point[1].r[1] << " " << _triangles[i].point[1].r[2]
        << "\n";
      f << "        vertex " << _triangles[i].point[2].r[0] << " "
        << _triangles[i].point[2].r[1] << " " << _triangles[i].point[2].r[2]
        << "\n";
      f << "    endloop\n";
      f << "endfacet\n";
    }
    f.close();
  }
  /*if (_verbose)*/clout << "Write ... OK" << std::endl;
}

template<typename T>
bool STLmesh<T>::testRayIntersect(const std::set<unsigned int>& tris, const Vector<T,3>& pt,const Vector<T,3>& dir, Vector<T,3>& q, T& alpha)
{
  std::set<unsigned int>::iterator it = tris.begin();
  for (; it != tris.end(); ++it) {
    if (_triangles[*it].testRayIntersect(pt, dir, q, alpha) && alpha < 1) {
      return true;
    }
  }
  return false;
}

/*
 * STLReader functions
 */
template<typename T>
STLreader<T>::STLreader(const std::string fName, T voxelSize, T stlSize,
                        unsigned short int method, bool verbose, T overlap, T max)
  : _voxelSize(voxelSize),
    _stlSize(stlSize),
    _overlap(overlap),
    _fName(fName),
    _mesh(fName, stlSize),
    _verbose(verbose),
    clout(std::cout, "STLreader")
{
  this->getName() = "STLreader";

  if (_verbose) {
    clout << "Voxelizing ..." << std::endl;
  }

  Vector<T,3> extension = _mesh.getMax() - _mesh.getMin();
  if ( util::nearZero(max) ) {
    max = std::max(extension[0], std::max(extension[1], extension[2])) + _voxelSize;
  }
  int j = 0;
  for (; _voxelSize * std::pow(2, j) < max; j++)
    ;
  Vector<T,3> center;
  T radius = _voxelSize * std::pow(2, j - 1);

  /// Find center of tree and move by _voxelSize/4.
  for (unsigned i = 0; i < 3; i++) {
    center[i] = (_mesh.getMin()[i] + _mesh.getMax()[i]) / 2. - _voxelSize / 4.;
  }

  /// Create tree
  _tree = new Octree<T>(center, radius, &_mesh, j, _overlap);

  /// Compute _myMin, _myMax such that they are the smallest (greatest) Voxel inside the STL.
  for (int i = 0; i < 3; i++) {
    this->_myMin[i] = center[i] + _voxelSize / 2.;
    this->_myMax[i] = center[i] - _voxelSize / 2.;
  }
  for (int i = 0; i < 3; i++) {
    while (this->_myMin[i] > _mesh.getMin()[i]) {
      this->_myMin[i] -= _voxelSize;
    }
    while (this->_myMax[i] < _mesh.getMax()[i]) {
      this->_myMax[i] += _voxelSize;
    }
    this->_myMax[i] -= _voxelSize;
    this->_myMin[i] += _voxelSize;
  }

  /// Indicate nodes of the tree. (Inside/Outside)
  switch (method) {
  case 1:
    indicate1();
    break;
  case 3:
    indicate3();
    break;
  default:
    indicate2();
    break;
  }

  if (_verbose) {
    print();
  }
  if (_verbose) {
    clout << "Voxelizing ... OK" << std::endl;
  }
}

/*
 * STLReader functions
 */
template<typename T>
STLreader<T>::STLreader(const std::vector<std::vector<T>> meshPoints, T voxelSize, T stlSize,
                        unsigned short int method, bool verbose, T overlap, T max)
  : _voxelSize(voxelSize),
    _stlSize(stlSize),
    _overlap(overlap),
    _fName("meshPoints.stl"),
    _mesh(meshPoints, stlSize),
    _verbose(verbose),
    clout(std::cout, "STLreader")
{
  this->getName() = "STLreader";

  if (_verbose) {
    clout << "Voxelizing ..." << std::endl;
  }

  Vector<T,3> extension = _mesh.getMax() - _mesh.getMin();
  if ( util::nearZero(max) ) {
    max = std::max(extension[0], std::max(extension[1], extension[2])) + _voxelSize;
  }
  int j = 0;
  for (; _voxelSize * std::pow(2, j) < max; j++)
    ;
  Vector<T,3> center;
  T radius = _voxelSize * std::pow(2, j - 1);

  /// Find center of tree and move by _voxelSize/4.
  for (unsigned i = 0; i < 3; i++) {
    center[i] = (_mesh.getMin()[i] + _mesh.getMax()[i]) / 2. - _voxelSize / 4.;
  }

  /// Create tree

  _tree = new Octree<T>(center, radius, &_mesh, j, _overlap);

  /// Compute _myMin, _myMax such that they are the smallest (greatest) Voxel inside the STL.
  for (int i = 0; i < 3; i++) {
    this->_myMin[i] = center[i] + _voxelSize / 2.;
    this->_myMax[i] = center[i] - _voxelSize / 2.;
  }
  for (int i = 0; i < 3; i++) {
    while (this->_myMin[i] > _mesh.getMin()[i]) {
      this->_myMin[i] -= _voxelSize;
    }
    while (this->_myMax[i] < _mesh.getMax()[i]) {
      this->_myMax[i] += _voxelSize;
    }
    this->_myMax[i] -= _voxelSize;
    this->_myMin[i] += _voxelSize;
  }

  // Indicate nodes of the tree. (Inside/Outside)
  switch (method) {
  case 1:
    indicate1();
    break;
  case 3:
    indicate3();
    break;
  default:
    indicate2();
    break;
  }

  if (_verbose) {
    print();
  }
  if (_verbose) {
    clout << "Voxelizing ... OK" << std::endl;
  }
}

template<typename T>
STLreader<T>::~STLreader()
{
  delete _tree;
}

/*
 *  Old indicate function (slower, more stable)
 *  Define three rays (X-, Y-, Z-direction) for each leaf and count intersections
 *  with STL for each ray. Odd number of intersection means inside (Majority vote).
 */

template<typename T>
void STLreader<T>::indicate1()
{
  std::vector<Octree<T>*> leafs;
  _tree->getLeafs(leafs);
  typename std::vector<Octree<T>*>::iterator it = leafs.begin();
  Vector<T,3> dir, pt, s;

  int intersections = 0;
  int inside = 0;
  Octree<T>* node = nullptr;
  T step = 1. / 1000. * _voxelSize;
  for (; it != leafs.end(); ++it) {
    inside = 0;

    pt = (*it)->getCenter();
    intersections = 0;
    s = pt;  // + step;

    /// X+ dir
    dir[0] = 1;
    dir[1] = 0;
    dir[2] = 0;
    while (s[0] < _mesh.getMax()[0] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections += node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
    }
    inside += (intersections % 2);

    /// Y+ Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = 1;
    dir[2] = 0;
    while (s[1] < _mesh.getMax()[1] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections += node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
    }
    inside += (intersections % 2);

    /// Z+ Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = 0;
    dir[2] = 1;
    while (s[2] < _mesh.getMax()[2] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections += node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
    }
    inside += (intersections % 2);
    (*it)->setInside(inside > 1);
  }
}

/*
 *  New indicate function (faster, less stable)
 *  Define ray in Z-direction for each Voxel in XY-layer. Indicate all nodes on the fly.
 */
template<typename T>
void STLreader<T>::indicate2()
{
  T rad = _tree->getRadius();
  Vector<T,3> rayPt = _tree->getCenter() - rad + .5 * _voxelSize;
  Vector<T,3> pt = rayPt;
  Vector<T,3> rayDir;
  rayDir[0] = 0.;
  rayDir[1] = 0.;
  rayDir[2] = 1.;
  //Vector<T,3> maxEdge = _tree->getCenter() + rad;

  T step = 1. / 1000. * _voxelSize;

  Octree<T>* node = nullptr;
  unsigned short rayInside = 0;
  Vector<T,3> nodeInters;
  while (pt[0] < _mesh.getMax()[0] + std::numeric_limits<T>::epsilon()) {
    node = _tree->find(pt);
    nodeInters = pt;
    nodeInters[2] = node->getCenter()[2] - node->getRadius();
    rayInside = 0;
    while (pt[1] < _mesh.getMax()[1] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(pt);
      nodeInters = pt;
      nodeInters[2] = node->getCenter()[2] - node->getRadius();
      rayInside = 0;
      while (pt[2] < _mesh.getMax()[2] + std::numeric_limits<T>::epsilon()) {
        node = _tree->find(pt);
        node->checkRay(nodeInters, rayDir, rayInside);
        node->intersectRayNode(pt, rayDir, nodeInters);
        pt = nodeInters + step * rayDir;
      }
      pt[2] = rayPt[2];
      pt[1] += _voxelSize;
    }
    pt[1] = rayPt[1];
    pt[0] += _voxelSize;
  }
}


/*
 *  Double ray approach: two times (X-, Y-, Z-direction) for each leaf.
 *  Could be use to deal with double layer triangles and face intersections.
 */
template<typename T>
void STLreader<T>::indicate3()
{
  std::vector<Octree<T>*> leafs;
  _tree->getLeafs(leafs);
  typename std::vector<Octree<T>*>::iterator it = leafs.begin();

  Vector<T,3> dir, pt, s;
  Octree<T>* node = nullptr;
  T step = 1. / 1000. * _voxelSize;
  int intersections;
  int sum_intersections;

  for (; it != leafs.end(); ++it) {
    pt = (*it)->getCenter();
    intersections = 0;
    sum_intersections = 0;
    s = pt;  // + step;

    /// X+ dir
    dir[0] = 1;
    dir[1] = 0;
    dir[2] = 0;
    while (s[0] < _mesh.getMax()[0] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }

    /// Y+ Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = 1;
    dir[2] = 0;
    while (s[1] < _mesh.getMax()[1] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }

    /// Z+ Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = 0;
    dir[2] = 1;
    while (s[2] < _mesh.getMax()[2] + std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }

    /// X- dir
    intersections = 0;
    s = pt;  // + step;
    dir[0] = -1;
    dir[1] = 0;
    dir[2] = 0;
    while (s[0] > _mesh.getMin()[0] - std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }

    /// Y- Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = -1;
    dir[2] = 0;
    while (s[1] > _mesh.getMin()[1] - std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }

    /// Z- Test
    intersections = 0;
    s = pt;  // + step;
    dir[0] = 0;
    dir[1] = 0;
    dir[2] = -1;
    while (s[2] > _mesh.getMin()[2] - std::numeric_limits<T>::epsilon()) {
      node = _tree->find(s, (*it)->getMaxdepth());
      intersections = node->testIntersection(pt, dir);
      node->intersectRayNode(pt, dir, s);
      s = s + step * dir;
      if (intersections > 0) {
        sum_intersections++;
        break;
      }
    }
    (*it)->setInside(sum_intersections > 5);
  }
}

template<typename T>
bool STLreader<T>::operator() (bool output[], const T input[])
{

  output[0] = false;
  T r = _tree->getRadius();
  Vector<T,3> c(_tree->getCenter());
  if (c[0] - r < input[0] && input[0] < c[0] + r && c[1] - r < input[1]
      && input[1] < c[1] + r && c[2] - r < input[2] && input[2] < c[2] + r) {
    std::vector<T> tmp(input, input + 3);
    output[0] = _tree->find(tmp)->getInside();
  }
  return true;
}

template<typename T>
bool STLreader<T>::distance(T& distance, const Vector<T,3>& origin,
                            const Vector<T,3>& direction, int iC)
{
  Octree<T>* node = nullptr;
  Vector<T,3> dir(direction);
  dir.normalize();
  Vector<T,3> extends = _mesh.getMax() - _mesh.getMin();
  Vector<T,3> pt(origin);
  Vector<T,3> q;
  Vector<T,3> s;
  Vector<T,3> center = _mesh.getMin() + 1 / 2. * extends;
  T step = _voxelSize / 1000., a = 0;

  for (int i = 0; i < 3; i++) {
    extends[i] /= 2.;
  }

  if (!(_mesh.getMin()[0] < origin[0] && origin[0] < _mesh.getMax()[0]
        && _mesh.getMin()[1] < origin[1] && origin[1] < _mesh.getMax()[1]
        && _mesh.getMin()[2] < origin[2] && origin[2] < _mesh.getMax()[2])) {
    T t = T(), d = T();
    bool foundQ = false;

    if (dir[0] > 0) {
      d = _mesh.getMin()[0];
      t = (d - origin[0]) / dir[0];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];

      if (_mesh.getMin()[1] < pt[1] && pt[1] < _mesh.getMax()[1]
          && _mesh.getMin()[2] < pt[2] && pt[2] < _mesh.getMax()[2]) {
        foundQ = true;
      }
    } else if (dir[0] < 0) {
      d = _mesh.getMax()[0];
      t = (d - origin[0]) / dir[0];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];
      if (_mesh.getMin()[1] < pt[1] && pt[1] < _mesh.getMax()[1]
          && _mesh.getMin()[2] < pt[2] && pt[2] < _mesh.getMax()[2]) {
        foundQ = true;
      }
    }

    if (dir[1] > 0 && !foundQ) {
      d = _mesh.getMin()[1];
      t = (d - origin[1]) / dir[1];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];
      if (_mesh.getMin()[0] < pt[0] && pt[0] < _mesh.getMax()[0]
          && _mesh.getMin()[2] < pt[2] && pt[2] < _mesh.getMax()[2]) {
        foundQ = true;
      }
    } else if (dir[1] < 0 && !foundQ) {
      d = _mesh.getMax()[1];
      t = (d - origin[1]) / dir[1];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];
      if (_mesh.getMin()[0] < pt[0] && pt[0] < _mesh.getMax()[0]
          && _mesh.getMin()[2] < pt[2] && pt[2] < _mesh.getMax()[2]) {
        foundQ = true;
      }
    }

    if (dir[2] > 0 && !foundQ) {
      d = _mesh.getMin()[2];
      t = (d - origin[2]) / dir[2];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];
      if (_mesh.getMin()[0] < pt[0] && pt[0] < _mesh.getMax()[0]
          && _mesh.getMin()[1] < pt[1] && pt[1] < _mesh.getMax()[1]) {
        foundQ = true;
      }
    } else if (dir[2] < 0 && !foundQ) {
      d = _mesh.getMax()[2];
      t = (d - origin[2]) / dir[2];
      pt[0] = origin[0] + (t + step) * dir[0];
      pt[1] = origin[1] + (t + step) * dir[1];
      pt[2] = origin[2] + (t + step) * dir[2];
      if (_mesh.getMin()[0] < pt[0] && pt[0] < _mesh.getMax()[0]
          && _mesh.getMin()[1] < pt[1] && pt[1] < _mesh.getMax()[1]) {
        foundQ = true;
      }
    }

    if (!foundQ) {
      return false;
    }
  }

  while ((std::fabs(pt[0] - center[0]) < extends[0])
         && (std::fabs(pt[1] - center[1]) < extends[1])
         && (std::fabs(pt[2] - center[2]) < extends[2])) {
    node = _tree->find(pt);
    if (node->closestIntersection(Vector<T,3>(origin), dir, q, a)) {
      Vector<T,3> vek(q - Vector<T,3>(origin));
      distance = vek.norm();
      return true;
    } else {
      Octree<T>* tmpNode = _tree->find(pt);
      tmpNode->intersectRayNode(pt, dir, s);
      for (int i = 0; i < 3; i++) {
        pt[i] = s[i] + step * dir[i];
      }
    }
  }

  if (_verbose) {
    clout << "Returning false" << std::endl;
  }
  return false;
}

template<typename T>
void STLreader<T>::print()
{
  _mesh.print();
  _tree->print();
  clout << "voxelSize=" << _voxelSize << "; stlSize=" << _stlSize << std::endl;
  clout << "minPhysR(VoxelMesh)=(" << this->_myMin[0] << "," << this->_myMin[1]
        << "," << this->_myMin[2] << ")";
  clout << "; maxPhysR(VoxelMesh)=(" << this->_myMax[0] << ","
        << this->_myMax[1] << "," << this->_myMax[2] << ")" << std::endl;
}

template<typename T>
void STLreader<T>::writeOctree()
{
  _tree->write(_fName);
}

template<typename T>
void STLreader<T>::writeSTL(std::string stlName)
{
  if (stlName == "") {
    _mesh.write(_fName);
  } else {
    _mesh.write(stlName);
  }
}

template<typename T>
void STLreader<T>::setNormalsOutside()
{
  unsigned int noTris = _mesh.triangleSize();
  Vector<T,3> center;
  //Octree<T>* node = nullptr;
  for (unsigned int i = 0; i < noTris; i++) {
    center[0] = (_mesh.getTri(i).point[0].r[0] + _mesh.getTri(i).point[1].r[0]
                 + _mesh.getTri(i).point[2].r[0]) / 3.;
    center[1] = (_mesh.getTri(i).point[0].r[1] + _mesh.getTri(i).point[1].r[1]
                 + _mesh.getTri(i).point[2].r[1]) / 3.;
    center[2] = (_mesh.getTri(i).point[0].r[2] + _mesh.getTri(i).point[1].r[2]
                 + _mesh.getTri(i).point[2].r[2]) / 3.;
    if (_tree->find(
          center + _mesh.getTri(i).normal * std::sqrt(3.) * _voxelSize)->getInside()) {
      //      cout << "Wrong direction" << std::endl;
      Vector<T,3> pt(_mesh.getTri(i).point[0].r);
      _mesh.getTri(i).point[0].r = _mesh.getTri(i).point[2].r;
      _mesh.getTri(i).point[2].r = pt;
      _mesh.getTri(i).init();
      //      _mesh.getTri(i).getNormal()[0] *= -1.;
      //      _mesh.getTri(i).getNormal()[1] *= -1.;
      //      _mesh.getTri(i).getNormal()[2] *= -1.;
    }
  }
}

}  // namespace olb

#endif
