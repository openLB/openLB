/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013-2018 Mathias Krause, Albert Mink, Adrian Kummerlaender
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

#ifndef BLOCK_GEOMETRY_FACES_3D_HH
#define BLOCK_GEOMETRY_FACES_3D_HH

#include "blockGeometryFaces3D.h"

namespace olb {


template <typename T>
BlockGeometryFaces3D<T>::BlockGeometryFaces3D(BlockIndicatorF3D<T>& indicatorF, T latticeL)
  : BlockF3D<T>(indicatorF.getBlockStructure(), 7),
    _indicatorF(indicatorF),
    _latticeL(latticeL)
{
  this->getName() = "blockGeometryFaces";
}

template <typename T>
bool BlockGeometryFaces3D<T>::operator() (T output[], const int input[])
{
  for (int i=0; i<7; ++i) {
    output[i] = T();
  }

  std::size_t counter[7] = {0};

  if (!_indicatorF.isEmpty()) {
    auto& blockGeometry = _indicatorF.getBlockGeometryStructure();
    const Vector<int,3> min = _indicatorF.getMin();
    const Vector<int,3> max = _indicatorF.getMax();

    // Iterate over all cells and count the cells of the face
    for (int iX = min[0]; iX <= max[0]; ++iX) {
      for (int iY = min[1]; iY <= max[1]; ++iY) {
        for (int iZ = min[2]; iZ <= max[2]; ++iZ) {
          // Lock at solid nodes only
          if (_indicatorF(iX, iY, iZ)) {
            if (blockGeometry.getMaterial(iX-1, iY, iZ) == 1) {
              counter[0]++;
            }
            if (blockGeometry.getMaterial(iX, iY-1, iZ) == 1) {
              counter[1]++;
            }
            if (blockGeometry.getMaterial(iX, iY, iZ-1) == 1) {
              counter[2]++;
            }
            if (blockGeometry.getMaterial(iX+1, iY, iZ) == 1) {
              counter[3]++;
            }
            if (blockGeometry.getMaterial(iX, iY+1, iZ) == 1) {
              counter[4]++;
            }
            if (blockGeometry.getMaterial(iX, iY, iZ+1) == 1) {
              counter[5]++;
            }
          }
        }
      }
    }

    const T dx2 = _latticeL*_latticeL;
    for (int i=0; i<6; ++i) {
      output[i]  = (T) counter[i] * dx2;
      output[6] += (T) counter[i] * dx2;
    }
  }

  return true;
}


}

#endif
