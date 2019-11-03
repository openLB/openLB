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
 * Representation of a 2d block geometry structure -- generic implementation.
 */

#ifndef BLOCK_GEOMETRY_STRUCTURE_2D_HH
#define BLOCK_GEOMETRY_STRUCTURE_2D_HH

#include <math.h>
#include "geometry/blockGeometryStructure2D.h"
#include "functors/lattice/indicator/blockIndicatorBaseF2D.h"


namespace olb {

template<typename T>
BlockGeometryStructure2D<T>::BlockGeometryStructure2D(int iCglob)
  : _iCglob(iCglob), _statistics(this), clout(std::cout,"BlockGeometryStructure2D") { }

template<typename T>
int const& BlockGeometryStructure2D<T>::getIcGlob() const
{
  return _iCglob;
}

template<typename T>
Vector<int,2> const BlockGeometryStructure2D<T>::getExtend() const
{
  return  Vector<int,2> (getNx(), getNy());
}

template<typename T>
void BlockGeometryStructure2D<T>::getPhysR(T physR[2], const int latticeR[2]) const
{
  getPhysR(physR, latticeR[0], latticeR[1]);
  return;
}

template<typename T>
int&  BlockGeometryStructure2D<T>::get(std::vector<int> latticeR)
{
  return get(latticeR[0], latticeR[1]);
}

template<typename T>
int const&  BlockGeometryStructure2D<T>::get(std::vector<int> latticeR) const
{
  return get(latticeR[0], latticeR[1]);
}

template<typename T>
int BlockGeometryStructure2D<T>::clean(bool verbose)
{
  int counter=0;
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (get(iX, iY) != 1 && get(iX, iY)!= 0) {
        if (   getMaterial(iX, iY) != 1
               && getMaterial(iX + 1, iY) != 1
               && getMaterial(iX - 1, iY) != 1
               && getMaterial(iX, iY + 1) != 1
               && getMaterial(iX + 1, iY + 1) != 1
               && getMaterial(iX - 1, iY + 1) != 1
               && getMaterial(iX, iY - 1) != 1
               && getMaterial(iX + 1, iY - 1) != 1
               && getMaterial(iX - 1, iY - 1) != 1 ) {
          get(iX, iY) = 0;
          counter++;
        }
      }
    }
  }
  if (verbose) {
    clout << "cleaned "<< counter << " outer boundary voxel(s)" << std::endl;
  }
  return counter;
}

template<typename T>
int BlockGeometryStructure2D<T>::outerClean(bool verbose)
{
  int counter=0;
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (get(iX, iY) == 1) {
        if (   getMaterial(iX + 1, iY) == 0
               || getMaterial(iX - 1, iY) == 0
               || getMaterial(iX, iY + 1) == 0
               || getMaterial(iX + 1, iY + 1) == 0
               || getMaterial(iX - 1, iY + 1) == 0
               || getMaterial(iX, iY - 1) == 0
               || getMaterial(iX + 1, iY - 1) == 0
               || getMaterial(iX - 1, iY - 1) == 0 ) {
          get(iX, iY) = 0;
          counter++;
        }
      }
    }
  }
  if (verbose) {
    clout << "cleaned "<< counter << " outer fluid voxel(s)" << std::endl;
  }
  return counter;
}

template<typename T>
int BlockGeometryStructure2D<T>::innerClean(bool verbose)
{
  int count = 0;
  int count2 = 0;

  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (get(iX, iY) != 1 && get(iX, iY) != 0) {
        count++;

        if (getMaterial(iX - 1, iY) == 1) {
          if (getMaterial(iX, iY + 1) == 1) {
            if (getMaterial(iX, iY - 1) == 1) {
              get(iX, iY) = 1;
              count2++;
            }
          }
        }
        if (getMaterial(iX + 1, iY) == 1) {
          if (getMaterial(iX, iY + 1) == 1) {
            if (getMaterial(iX, iY - 1) == 1) {
              get(iX, iY) = 1;
              count2++;
            }
          }
        }
        if (getMaterial(iX, iY + 1) == 1) {
          if (getMaterial(iX + 1, iY) == 1) {
            if (getMaterial(iX - 1, iY) == 1) {
              get(iX, iY) = 1;
              count2++;
            }
          }
        }
        if (getMaterial(iX, iY - 1) == 1) {
          if (getMaterial(iX + 1, iY) == 1) {
            if (getMaterial(iX - 1, iY) == 1) {
              get(iX, iY) = 1;
              count2++;
            }
          }
        }
      }
    }
  }
  if (verbose) {
    this->clout << "cleaned "<< count2 << " inner boundary voxel(s)" << std::endl;
  }
  return count2;
}

template<typename T>
int BlockGeometryStructure2D<T>::innerClean(int fromM, bool verbose)
{
  int count = 0;
  int count2 = 0;

  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (get(iX, iY) != 1 && get(iX, iY)!= 0 && get(iX, iY) == fromM) {
        count++;
        if (getMaterial(iX - 1, iY) == 1) {
          if (getMaterial(iX, iY + 1) == 1) {
            if (getMaterial(iX, iY - 1) == 1) {
              get(iX, iY) = 1;
              count2++;
            }
          }
        }
        if (getMaterial(iX + 1, iY) == 1) {
          if (getMaterial(iX, iY + 1) == 1) {
            if (getMaterial(iX, iY - 1) == 1) {
              get(iX, iY) = 1;
              count2++;
            }
          }
        }
        if (getMaterial(iX, iY + 1) == 1) {
          if (getMaterial(iX + 1, iY) == 1) {
            if (getMaterial(iX - 1, iY) == 1) {
              get(iX, iY) = 1;
              count2++;
            }
          }
        }
        if (getMaterial(iX, iY - 1) == 1) {
          if (getMaterial(iX + 1, iY) == 1) {
            if (getMaterial(iX - 1, iY) == 1) {
              get(iX, iY) = 1;
              count2++;
            }
          }
        }
      }
    }
  }
  if (verbose)
    this->clout << "cleaned "<< count2
                << " inner boundary voxel(s) of Type " << fromM << std::endl;
  return count2;
}

template<typename T>
void BlockGeometryStructure2D<T>::reset(IndicatorF2D<T>& domain)
{
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      T physR[2] { };
      getPhysR(physR, iX, iY);
      if (domain(physR)) {
        get(iX, iY) = 0;
      }
    }
  }
}

template<typename T>
template<typename V, typename DESCRIPTOR >
bool BlockGeometryStructure2D<T>::findStreamDirections(
  int iX, int iY,
  BlockIndicatorF2D<T>& boundaryIndicator, BlockIndicatorF2D<T>& bulkIndicator,
  bool streamDirections[])
{
  if (boundaryIndicator(iX, iY)) {
    bool found = false;
    streamDirections[0] = false;
    for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
      streamDirections[iPop] = false;
      if (bulkIndicator(iX + descriptors::c<DESCRIPTOR >(iPop,0), iY + descriptors::c<DESCRIPTOR >(iPop,1))) {
        streamDirections[iPop] = true;
        found = true;
      }
    }
    return found;
  }
  else {
    return false;
  }
}

template<typename T>
template<typename V, typename DESCRIPTOR >
bool BlockGeometryStructure2D<T>::findStreamDirections(int iX, int iY, int material, std::list<int> bulkMaterials, bool streamDirections[])
{
  bool found = false;
  if (getMaterial(iX, iY) != material) {
    return false;
  }
  else {
    std::list<int>::iterator mat;
    streamDirections[0] = false;
    for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
      streamDirections[iPop] = false;
      for (mat=bulkMaterials.begin(); !streamDirections[iPop] && mat!=bulkMaterials.end(); ++mat) {
        if (getMaterial(iX + descriptors::c<DESCRIPTOR >(iPop,0), iY + descriptors::c<DESCRIPTOR >(iPop,1)) == *mat ) {
          streamDirections[iPop] = true;
          found = true;
        }
      }
    }
    return found;
  }
}

template<typename T>
bool BlockGeometryStructure2D<T>::find(int material, unsigned offsetX, unsigned offsetY,
                                       int& foundX, int& foundY)
{

  bool found = false;
  for (foundX = 0; foundX < getNx(); foundX++) {
    for (foundY = 0; foundY < getNy(); foundY++) {
      found = check(material, foundX, foundY, offsetX, offsetY);
      if (found) {
        return found;
      }
    }
  }
  return found;
}

template<typename T>
bool BlockGeometryStructure2D<T>::check(int material, int iX, int iY,
                                        unsigned offsetX, unsigned offsetY)
{
  bool found = true;
  for (int iOffsetX = -offsetX; iOffsetX <= (int) offsetX; ++iOffsetX) {
    for (int iOffsetY = -offsetY; iOffsetY <= (int) offsetY; ++iOffsetY) {
      if (getMaterial(iX + iOffsetX, iY + iOffsetY) != material) {
        found = false;
      }
    }
  }
  return found;
}

template<typename T>
bool BlockGeometryStructure2D<T>::checkForErrors(bool verbose) const
{
  bool error = false;

  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (get(iX, iY) == 0) {
        if (   getMaterial(iX + 1, iY) == 1
               || getMaterial(iX - 1, iY) == 1
               || getMaterial(iX, iY + 1) == 1
               || getMaterial(iX + 1, iY + 1) == 1
               || getMaterial(iX - 1, iY + 1) == 1
               || getMaterial(iX, iY - 1) == 1
               || getMaterial(iX + 1, iY - 1) == 1
               || getMaterial(iX - 1, iY - 1) == 1 ) {
          error = true;
        }
      }
    }
  }
  if (verbose) {
    if (error) {
      this->clout << "error!" << std::endl;
    }
    else {
      this->clout << "the model is correct!" << std::endl;
    }
  }
  return error;
}

template<typename T>
void BlockGeometryStructure2D<T>::rename(int fromM, int toM)
{
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (get(iX, iY) == fromM) {
        get(iX, iY) = toM;
      }
    }
  }
}

template<typename T>
void BlockGeometryStructure2D<T>::rename(int fromM, int toM, IndicatorF2D<T>& condition)
{
  T physR[2];
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (get(iX, iY) == fromM) {
        getPhysR(physR, iX,iY);
        bool inside[1];
        condition(inside, physR);
        if (inside[0]) {
          get(iX, iY) = toM;
        }
      }
    }
  }
}

template<typename T>
void BlockGeometryStructure2D<T>::rename(int fromM, int toM, unsigned offsetX,
    unsigned offsetY)
{
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (get(iX, iY) == fromM) {
        bool found = true;
        for (int iOffsetX = -offsetX; iOffsetX <= (int) offsetX; ++iOffsetX) {
          for (int iOffsetY = -offsetY; iOffsetY <= (int) offsetY; ++iOffsetY) {
            if (getMaterial(iX + iOffsetX, iY + iOffsetY) != fromM) {
              if (getMaterial(iX + iOffsetX, iY + iOffsetY) != 1245) {
                found = false;
              }
            }
          }
        }
        if (found) {
          get(iX, iY) = 1245;
        }
      }
    }
  }
  rename(1245,toM);
}

template<typename T>
void BlockGeometryStructure2D<T>::rename(int fromM, int toM, int testM,
    std::vector<int> testDirection)
{
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (get(iX, iY) == fromM) {

        // flag that indicates the renaming of the current voxel, valid voxels are not renamed
        bool isValid = true;
        for (int iOffsetX = std::min(testDirection[0],0); iOffsetX <= std::max(testDirection[0],0); ++iOffsetX) {
          for (int iOffsetY = std::min(testDirection[1],0); iOffsetY <= std::max(testDirection[1],0); ++iOffsetY) {
            if (iOffsetX!=0 || iOffsetY!=0) {
              if (getMaterial(iX + iOffsetX, iY + iOffsetY) != testM) {
                isValid = false;
              }
            }
          }
        }
        if (!isValid) {
          get(iX, iY) = toM;
        }
      }
    }
  }
}


template<typename T>
void BlockGeometryStructure2D<T>::rename(int fromM, int toM, int fluidM,
    IndicatorF2D<T>& condition, std::vector<int> discreteNormal)
{
  rename(fromM, toM, condition);
  std::vector<int> testDirection(discreteNormal);
  T physR[2];
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (get(iX, iY) == toM) {
        getPhysR(physR, iX,iY);
        bool inside[1];
        condition(inside, physR);
        if (inside[0]) {
          if (getMaterial(iX+testDirection[0],iY+testDirection[1])!=fluidM ||
              getMaterial(iX+2*testDirection[0],iY+2*testDirection[1])!=fluidM ||
              getMaterial(iX-testDirection[0],iY-testDirection[1])!=0 ) {
            get(iX, iY) = fromM;
          }
        }
      }
    }
  }
}

template<typename T>
void BlockGeometryStructure2D<T>::rename(int fromM, int toM, int fluidM,
    IndicatorF2D<T>& condition)
{
  rename(fromM, toM, condition);
  std::vector<int> testDirection = getStatistics().computeDiscreteNormal(toM);
  T physR[3];
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      if (get(iX, iY) == toM) {
        getPhysR(physR, iX,iY);
        bool inside[1];
        condition(inside, physR);
        if (inside[0]) {
          if (getMaterial(iX+testDirection[0],iY+testDirection[1])!=fluidM ||
              getMaterial(iX+2*testDirection[0],iY+2*testDirection[1])!=fluidM ||
              getMaterial(iX-testDirection[0],iY-testDirection[1])!=0 ) {
            get(iX, iY) = fromM;
          }
        }
      }
    }
  }
}

template<typename T>
void BlockGeometryStructure2D<T>::regionGrowing(int fromM, int toM, int seedX, int seedY,
    int offsetX, int offsetY, std::map<std::vector<int>, int>* tmp)
{
  std::map<std::vector<int>, int> tmp2;
  bool firstCall = false;
  if (tmp == nullptr) {
    tmp = &tmp2;
    firstCall = true;
  }

  if (getMaterial(seedX, seedY) == fromM) {
    std::vector<int> found;
    found.push_back(seedX);
    found.push_back(seedY);
    if (tmp->count(found) == 0) {
      (*tmp)[found] = 2;
      if (offsetX != 0) {
        regionGrowing(fromM, toM, seedX + 1, seedY, offsetX,
                      offsetY, tmp);
        regionGrowing(fromM, toM, seedX - 1, seedY, offsetX,
                      offsetY, tmp);
      }
      if (offsetY != 0) {
        regionGrowing(fromM, toM, seedX, seedY + 1, offsetX,
                      offsetY, tmp);
        regionGrowing(fromM, toM, seedX, seedY - 1, offsetX,
                      offsetY, tmp);
      }
    }
  }
  if (firstCall) {
    std::map<std::vector<int>, int>::iterator iter;
    for (iter = tmp->begin(); iter != tmp->end(); iter++) {
      get((iter->first)[0],(iter->first)[1]) = toM;
    }
  }
  return;
}


} // namespace olb

#endif
