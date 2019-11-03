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

#ifndef PLATTICE_H_
#define PLATTICE_H_

#include "sortAlgorithms.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
class PLattice : public ContactDetection<T, PARTICLETYPE> {
 public:
  PLattice(ParticleSystem3D<T, PARTICLETYPE>& pSys, T overlap, T spacing);

  void sort();
  int getMatches(int pInt,
                 std::vector<std::pair<size_t, T> >& matches);

 private:
  PLattice<T, PARTICLETYPE>* generate(ParticleSystem3D<T, PARTICLETYPE>& pSys);


  std::vector<T> _physPos, _physExtend;
  std::vector<int> _intExtend;
  T _overlap;
  T _spacing;

  std::vector<std::vector<std::vector<std::list<int>> > > _pLattice;

};

template<typename T, template<typename U> class PARTICLETYPE>
PLattice<T, PARTICLETYPE>* PLattice<T, PARTICLETYPE>::generate(ParticleSystem3D<T, PARTICLETYPE>& pSys) {
  //std::cout << "calling NanoflannContact.generate()" << std::endl;
  return new PLattice(pSys,_overlap, _spacing);
}

template<typename T, template<typename U> class PARTICLETYPE>
PLattice<T, PARTICLETYPE>::PLattice(ParticleSystem3D<T, PARTICLETYPE>& pSys,
                                    T overlap, T spacing)
    : ContactDetection<T, PARTICLETYPE>(pSys, "PLattice"),
      _physPos(3, T()),
      _physExtend(3, T()),
      _intExtend(3, 0),
      _overlap(overlap),
      _spacing(spacing) {

  _physPos = this->_pSys.getPhysPos();
  _physExtend = this->_pSys.getPhysExtend();

  int intOverlap = std::ceil(overlap / _spacing);
  cout << "intOverlap: " << intOverlap << std::endl;
  for (int i = 0; i < 3; i++) {
    _physPos[i] -= overlap;
    _intExtend[i] = std::ceil(_physExtend[i] / _spacing + 2 * intOverlap + 1);
    cout << "intExtend[" << i << "]: " << _physExtend[i] << std::endl;
  }

  cout << "intOverlap: " << intOverlap << std::endl;

  _pLattice.resize(_intExtend[0]);
  for (int iX = 0; iX < _intExtend[0]; ++iX) {
    _pLattice[iX].resize(_intExtend[1]);
    for (int iY = 0; iY < _intExtend[1]; ++iY) {
      _pLattice[iX][iY].resize(_intExtend[2]);
      for (int iZ = 0; iZ < _intExtend[2]; ++iZ) {
      }
    }
  }
//  cout << "intExtend " << _intExtend[0] << " "<< _intExtend[1] << " "<< _intExtend[2] << " " << std::endl;
}

template<typename T, template<typename U> class PARTICLETYPE>
void PLattice<T, PARTICLETYPE>::sort() {
//  _pLattice.clear();
  for (int i = 0; i < _intExtend[0]; i++) {
    for (int j = 0; j < _intExtend[1]; j++) {
      for (int k = 0; k < _intExtend[2]; k++) {
        _pLattice[i][j][k].clear();
      }
    }
  }

  std::vector<T> pos(3, T());
  for (unsigned int i = 0; i < this->_pSys.sizeInclShadow(); i++) {
    pos = this->_pSys[i].getPos();
#if OLB_DEBUG
    int aa = (pos[0] - _physPos[0]) / _spacing;
    int bb = (pos[1] - _physPos[1]) / _spacing;
    int cc = (pos[2] - _physPos[2]) / _spacing;

    OLB_PRECONDITION( aa >= 0);
    OLB_PRECONDITION( bb >= 0);
    OLB_PRECONDITION( cc >= 0);
    OLB_PRECONDITION( aa < _intExtend[0]);
    OLB_PRECONDITION( bb < _intExtend[1]);
    OLB_PRECONDITION( cc < _intExtend[2]);
#endif

    _pLattice[(int) floor((pos[0] - _physPos[0]) / _spacing)][(int) floor(
        (pos[1] - _physPos[1]) / _spacing)][(int) floor(
        (pos[2] - _physPos[2]) / _spacing)].push_back(i);
  }

  int iX = 0, iY = 0, iZ = 0;
  int x = 0, y = 0, z = 0;
  for (unsigned int s = 0; s < this->_pSys.size(); s++) {
    pos = this->_pSys[s].getPos();
    iX = floor((pos[0] - _physPos[0]) / _spacing);
    iY = floor((pos[1] - _physPos[1]) / _spacing);
    iZ = floor((pos[2] - _physPos[2]) / _spacing);
    //
    if (_pLattice[iX][iY][iZ].back() != -1) {

      int a = 1, b = 1, c = 1;
      int d = -1, e = -1, f = -1;
      if (iX <= 0)
        d = 0;
      if (iY <= 0)
        e = 0;
      if (iZ <= 0)
        f = 0;
      if (iX >= _intExtend[0])
        a = 0;
      if (iY >= _intExtend[1])
        b = 0;
      if (iZ >= _intExtend[2])
        c = 0;

      typename std::list<int>::iterator it;
      for (int i = d; i <= a; i++) {
        for (int j = e; j <= b; j++) {
          for (int k = f; k <= c; k++) {
            if (i != 0 || j != 0 || k != 0) {
              it = _pLattice[iX + i][iY + j][iZ + k].begin();
              int size = _pLattice[iX + i][iY + j][iZ + k].size();
              for (int lauf = 0; lauf < size; lauf++, it++) {
                if (*it != -1) {
                  pos = this->_pSys[*it].getPos();
                  T rad = this->_pSys[*it].getRad();
                  x = floor((pos[0] - _physPos[0] - i * rad) / _spacing);
                  y = floor((pos[1] - _physPos[1] - j * rad) / _spacing);
                  z = floor((pos[2] - _physPos[2] - k * rad) / _spacing);
                  if (x != iX || y != iY || z != iZ) {
                    _pLattice[iX][iY][iZ].push_back(*it);
                  }
                }
              }
            }
          }
        }
      }
      _pLattice[iX][iY][iZ].push_back(-1);
    }
  }
}

template<typename T, template<typename U> class PARTICLETYPE>
int PLattice<T, PARTICLETYPE>::getMatches(int pInt,
    std::vector<std::pair<size_t, T> >& matches) {
  matches.clear();
  PARTICLETYPE<T>& p = this->_pSys[pInt];
  int iX = floor((p.getPos()[0] - _physPos[0]) / _spacing);
  int iY = floor((p.getPos()[1] - _physPos[1]) / _spacing);
  int iZ = floor((p.getPos()[2] - _physPos[2]) / _spacing);

  typename std::list<int>::iterator it = _pLattice[iX][iY][iZ].begin();
  std::vector<T>& pos = p.getPos();
  for (; it != _pLattice[iX][iY][iZ].end(); ++it) {
    if (*it != -1) {
      std::vector<T>& pos2 = this->_pSys[*it].getPos();
      matches.push_back(
          make_pair<size_t, T>(
              *it,
              std::pow(pos[0] - pos2[0], 2) + std::pow(pos[1] - pos2[1], 2)
                  + std::pow(pos[2] - pos2[2], 2)));
    }
  }
  return matches.size();
}

}  // namespace olb

#endif
