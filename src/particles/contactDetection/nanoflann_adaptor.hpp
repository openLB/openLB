/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Thomas Henn, Mathias J. Krause
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

#ifndef NanoflannParticleAdaptor_H_
#define NanoflannParticleAdaptor_H_

#include "nanoflann.hpp"

using namespace std;
using namespace nanoflann;

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
class ContactDetection;

template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D;

template<typename coord_t, typename Derived>
class NanoflannParticleAdaptor {
 public:
//		typedef typename T coord_t;
  const Derived &obj;  //!< A const ref to the data set origin

  /// The constructor that sets the data set source
  NanoflannParticleAdaptor(const Derived &obj_)
      : obj(obj_) {
  }

  /// CRTP helper method
  inline const Derived& derived() const {
    return obj;
  }

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const {
    return derived().sizeInclShadow();
  }

  // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
  inline coord_t kdtree_distance(const coord_t *p1, const size_t idx_p2,
                                 size_t size) const {
    return std::pow(p1[0] - derived()[idx_p2].getPos()[0], 2) +
           std::pow(p1[1] - derived()[idx_p2].getPos()[1], 2) +
           std::pow(p1[2] - derived()[idx_p2].getPos()[2], 2);
//    return (p1[0] - derived()[idx_p2].getPos()[0]) * (p1[0] - derived()[idx_p2].getPos()[0]) +
//           (p1[1] - derived()[idx_p2].getPos()[1]) * (p1[1] - derived()[idx_p2].getPos()[1]) +
//           (p1[2] - derived()[idx_p2].getPos()[2]) * (p1[2] - derived()[idx_p2].getPos()[2]);
  }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate value, the
  //  "if/else's" are actually solved at compile time.
  inline coord_t kdtree_get_pt(const size_t idx, int dim) const {
    return derived()[idx].getPos()[dim];
  }

  // Optional bounding-box computation: return false to default to a standard bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
  //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
  template<class BBOX>
  bool kdtree_get_bbox(BBOX &bb) const {
    return false;
  }
};

template<typename T, template<typename U> class PARTICLETYPE>
class NanoflannContact : public ContactDetection<T, PARTICLETYPE> {
 public:

  NanoflannContact(ParticleSystem3D<T, PARTICLETYPE>& pSys, T sRad) : ContactDetection<T, PARTICLETYPE> (pSys, "Nanoflann"),
  _pc2kd(pSys),
  _index(3,_pc2kd, KDTreeSingleIndexAdaptorParams(10) ),
  _sRad2(sRad*sRad),
  _sRad(sRad)
 {
    _index.init();
    _params.sorted = false;
  }
  virtual ContactDetection<T, PARTICLETYPE>* generate(ParticleSystem3D<T, PARTICLETYPE>& pSys) {
    //std::cout << "calling NanoflannContact.generate()" << std::endl;
    return new NanoflannContact(pSys,_sRad);
  }

  void sort() {
    _index.buildIndex();
  }

  int getMatches(int pInt, std::vector<std::pair<size_t, T> >& matches) {
   _index.radiusSearch(&this->_pSys[pInt].getPos()[0], _sRad2, matches, _params);
   return matches.size();
  }

 private:
    typedef NanoflannParticleAdaptor<T, ParticleSystem3D<T, PARTICLETYPE> > PC2KD;
    const PC2KD _pc2kd;
    typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<T, PC2KD > ,PC2KD, 3> kd_tree;
    kd_tree _index;
    nanoflann::SearchParams _params;
    T _sRad2;
    T _sRad;
};

}  // namespace olb

#endif
