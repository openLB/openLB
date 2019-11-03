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
 * Representation of a 3d block geometry structure -- generic implementation.
 */

#ifndef BLOCK_GEOMETRY_STRUCTURE_3D_HH
#define BLOCK_GEOMETRY_STRUCTURE_3D_HH

#include <math.h>
#include "geometry/blockGeometryStructure3D.h"
#include "functors/analytical/indicator/indicatorF3D.h"


namespace olb {

template<typename T>
BlockGeometryStructure3D<T>::BlockGeometryStructure3D(int iCglob)
  : _iCglob(iCglob), _statistics(this), clout(std::cout,"BlockGeometryStructure3D") { }

template<typename T>
int const& BlockGeometryStructure3D<T>::getIcGlob() const
{
  return _iCglob;
}

template<typename T>
Vector<int,3> const BlockGeometryStructure3D<T>::getExtend() const
{
  return Vector<int,3> (getNx(), getNy(), getNz());
}

template<typename T>
void BlockGeometryStructure3D<T>::getPhysR(T physR[3], const int latticeR[3]) const
{
  getPhysR(physR, latticeR[0], latticeR[1], latticeR[2]);
  return;
}

template<typename T>
bool BlockGeometryStructure3D<T>::isInside(int iX, int iY, int iZ) const
{
  return 0 <= iX && iX < getNx() &&
         0 <= iY && iY < getNy() &&
         0 <= iZ && iZ < getNz();
}

template<typename T>
int&  BlockGeometryStructure3D<T>::get(std::vector<int> latticeR)
{
  return get(latticeR[0], latticeR[1], latticeR[2]);
}

template<typename T>
int const&  BlockGeometryStructure3D<T>::get(std::vector<int> latticeR) const
{
  return get(latticeR[0], latticeR[1], latticeR[2]);
}

template<typename T>
int BlockGeometryStructure3D<T>::clean(bool verbose)
{
  int counter=0;
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {

        if (get(iX, iY, iZ) != 1 && get(iX, iY, iZ) != 0) {

          if (
            getMaterial(iX - 1, iY, iZ) != 1
            && getMaterial(iX, iY - 1, iZ) != 1
            && getMaterial(iX, iY, iZ - 1) != 1
            && getMaterial(iX - 1, iY - 1, iZ) != 1
            && getMaterial(iX, iY - 1, iZ - 1) != 1
            && getMaterial(iX - 1, iY, iZ - 1) != 1
            && getMaterial(iX - 1, iY - 1, iZ - 1) != 1
            && getMaterial(iX + 1, iY, iZ) != 1
            && getMaterial(iX, iY + 1, iZ) != 1
            && getMaterial(iX, iY, iZ + 1) != 1
            && getMaterial(iX + 1, iY + 1, iZ) != 1
            && getMaterial(iX, iY + 1, iZ + 1) != 1
            && getMaterial(iX + 1, iY, iZ + 1) != 1
            && getMaterial(iX + 1, iY + 1, iZ + 1) != 1
            && getMaterial(iX - 1, iY + 1, iZ) != 1
            && getMaterial(iX + 1, iY - 1, iZ) != 1
            && getMaterial(iX, iY - 1, iZ + 1) != 1
            && getMaterial(iX, iY + 1, iZ - 1) != 1
            && getMaterial(iX - 1, iY, iZ + 1) != 1
            && getMaterial(iX + 1, iY, iZ - 1) != 1
            && getMaterial(iX + 1, iY + 1, iZ - 1) != 1
            && getMaterial(iX + 1, iY - 1, iZ - 1) != 1
            && getMaterial(iX + 1, iY - 1, iZ + 1) != 1
            && getMaterial(iX - 1, iY + 1, iZ + 1) != 1
            && getMaterial(iX - 1, iY - 1, iZ + 1) != 1
            && getMaterial(iX - 1, iY + 1, iZ - 1) != 1) {

            get(iX, iY, iZ) = 0;
            counter++;
          }
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
int BlockGeometryStructure3D<T>::clean(int material, bool verbose)
{
  int counter=0;
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {

        if (get(iX, iY, iZ) != 1 && get(iX, iY, iZ) != 0 && get(iX, iY, iZ) != material) {

          if (
            getMaterial(iX - 1, iY, iZ) != 1
            && getMaterial(iX, iY - 1, iZ) != 1
            && getMaterial(iX, iY, iZ - 1) != 1
            && getMaterial(iX - 1, iY - 1, iZ) != 1
            && getMaterial(iX, iY - 1, iZ - 1) != 1
            && getMaterial(iX - 1, iY, iZ - 1) != 1
            && getMaterial(iX - 1, iY - 1, iZ - 1) != 1
            && getMaterial(iX + 1, iY, iZ) != 1
            && getMaterial(iX, iY + 1, iZ) != 1
            && getMaterial(iX, iY, iZ + 1) != 1
            && getMaterial(iX + 1, iY + 1, iZ) != 1
            && getMaterial(iX, iY + 1, iZ + 1) != 1
            && getMaterial(iX + 1, iY, iZ + 1) != 1
            && getMaterial(iX + 1, iY + 1, iZ + 1) != 1
            && getMaterial(iX - 1, iY + 1, iZ) != 1
            && getMaterial(iX + 1, iY - 1, iZ) != 1
            && getMaterial(iX, iY - 1, iZ + 1) != 1
            && getMaterial(iX, iY + 1, iZ - 1) != 1
            && getMaterial(iX - 1, iY, iZ + 1) != 1
            && getMaterial(iX + 1, iY, iZ - 1) != 1
            && getMaterial(iX + 1, iY + 1, iZ - 1) != 1
            && getMaterial(iX + 1, iY - 1, iZ - 1) != 1
            && getMaterial(iX + 1, iY - 1, iZ + 1) != 1
            && getMaterial(iX - 1, iY + 1, iZ + 1) != 1
            && getMaterial(iX - 1, iY - 1, iZ + 1) != 1
            && getMaterial(iX - 1, iY + 1, iZ - 1) != 1) {

            get(iX, iY, iZ) = 0;
            counter++;
          }
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
int BlockGeometryStructure3D<T>::outerClean(bool verbose)
{
  int counter=0;
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {

        if (get(iX, iY, iZ) == 1) {
          if (   getMaterial(iX - 1, iY,     iZ    ) == 0
                 || getMaterial(iX,     iY - 1, iZ    ) == 0
                 || getMaterial(iX,     iY,     iZ - 1) == 0
                 || getMaterial(iX - 1, iY - 1, iZ    ) == 0
                 || getMaterial(iX,     iY - 1, iZ - 1) == 0
                 || getMaterial(iX - 1, iY,     iZ - 1) == 0
                 || getMaterial(iX - 1, iY - 1, iZ - 1) == 0
                 || getMaterial(iX + 1, iY,     iZ    ) == 0
                 || getMaterial(iX,     iY + 1, iZ    ) == 0
                 || getMaterial(iX,     iY,     iZ + 1) == 0
                 || getMaterial(iX + 1, iY + 1, iZ    ) == 0
                 || getMaterial(iX,     iY + 1, iZ + 1) == 0
                 || getMaterial(iX + 1, iY,     iZ + 1) == 0
                 || getMaterial(iX + 1, iY + 1, iZ + 1) == 0
                 || getMaterial(iX - 1, iY + 1, iZ    ) == 0
                 || getMaterial(iX + 1, iY - 1, iZ    ) == 0
                 || getMaterial(iX,     iY - 1, iZ + 1) == 0
                 || getMaterial(iX,     iY + 1, iZ - 1) == 0
                 || getMaterial(iX - 1, iY,     iZ + 1) == 0
                 || getMaterial(iX + 1, iY,     iZ - 1) == 0
                 || getMaterial(iX + 1, iY + 1, iZ - 1) == 0
                 || getMaterial(iX + 1, iY - 1, iZ - 1) == 0
                 || getMaterial(iX + 1, iY - 1, iZ + 1) == 0
                 || getMaterial(iX - 1, iY + 1, iZ + 1) == 0
                 || getMaterial(iX - 1, iY - 1, iZ + 1) == 0
                 || getMaterial(iX - 1, iY + 1, iZ - 1) == 0) {
            get(iX, iY, iZ) = 0;
            counter++;
          }
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
int BlockGeometryStructure3D<T>::innerClean(bool verbose)
{
  int count = 0;
  int count2 = 0;

  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {

        if (get(iX, iY, iZ) != 1 && get(iX, iY, iZ) != 0) {
          count++;

          if (getMaterial(iX - 1, iY, iZ) == 1) {
            if (getMaterial(iX, iY + 1, iZ) == 1) {
              if (getMaterial(iX, iY - 1, iZ) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX + 1, iY, iZ) == 1) {
            if (getMaterial(iX, iY + 1, iZ) == 1) {
              if (getMaterial(iX, iY - 1, iZ) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX, iY, iZ - 1) == 1) {
            if (getMaterial(iX, iY + 1, iZ) == 1) {
              if (getMaterial(iX, iY - 1, iZ) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX, iY, iZ + 1) == 1) {
            if (getMaterial(iX, iY + 1, iZ) == 1) {
              if (getMaterial(iX, iY - 1, iZ) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX, iY + 1, iZ) == 1) {
            if (getMaterial(iX + 1, iY, iZ) == 1) {
              if (getMaterial(iX - 1, iY, iZ) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX, iY - 1, iZ) == 1) {
            if (getMaterial(iX + 1, iY, iZ) == 1) {
              if (getMaterial(iX - 1, iY, iZ) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX, iY, iZ + 1) == 1) {
            if (getMaterial(iX + 1, iY, iZ) == 1) {
              if (getMaterial(iX - 1, iY, iZ) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX, iY, iZ - 1) == 1) {
            if (getMaterial(iX + 1, iY, iZ) == 1) {
              if (getMaterial(iX - 1, iY, iZ) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX, iY + 1, iZ) == 1) {
            if (getMaterial(iX, iY, iZ + 1) == 1) {
              if (getMaterial(iX, iY, iZ - 1) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX, iY - 1, iZ) == 1) {
            if (getMaterial(iX, iY, iZ + 1) == 1) {
              if (getMaterial(iX, iY, iZ - 1) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX + 1, iY, iZ) == 1) {
            if (getMaterial(iX, iY, iZ + 1) == 1) {
              if (getMaterial(iX, iY, iZ - 1) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX - 1, iY, iZ) == 1) {
            if (getMaterial(iX, iY, iZ + 1) == 1) {
              if (getMaterial(iX, iY, iZ - 1) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
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
int BlockGeometryStructure3D<T>::innerClean(int fromM, bool verbose)
{
  int count = 0;
  int count2 = 0;

  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {

        if (get(iX, iY, iZ) != 1 && get(iX, iY, iZ)
            != 0 && get(iX, iY, iZ) == fromM) {
          count++;

          if (getMaterial(iX - 1, iY, iZ) == 1) {
            if (getMaterial(iX, iY + 1, iZ) == 1) {
              if (getMaterial(iX, iY - 1, iZ) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX + 1, iY, iZ) == 1) {
            if (getMaterial(iX, iY + 1, iZ) == 1) {
              if (getMaterial(iX, iY - 1, iZ) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX, iY, iZ - 1) == 1) {
            if (getMaterial(iX, iY + 1, iZ) == 1) {
              if (getMaterial(iX, iY - 1, iZ) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX, iY, iZ + 1) == 1) {
            if (getMaterial(iX, iY + 1, iZ) == 1) {
              if (getMaterial(iX, iY - 1, iZ) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX, iY + 1, iZ) == 1) {
            if (getMaterial(iX + 1, iY, iZ) == 1) {
              if (getMaterial(iX - 1, iY, iZ) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX, iY - 1, iZ) == 1) {
            if (getMaterial(iX + 1, iY, iZ) == 1) {
              if (getMaterial(iX - 1, iY, iZ) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX, iY, iZ + 1) == 1) {
            if (getMaterial(iX + 1, iY, iZ) == 1) {
              if (getMaterial(iX - 1, iY, iZ) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX, iY, iZ - 1) == 1) {
            if (getMaterial(iX + 1, iY, iZ) == 1) {
              if (getMaterial(iX - 1, iY, iZ) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX, iY + 1, iZ) == 1) {
            if (getMaterial(iX, iY, iZ + 1) == 1) {
              if (getMaterial(iX, iY, iZ - 1) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX, iY - 1, iZ) == 1) {
            if (getMaterial(iX, iY, iZ + 1) == 1) {
              if (getMaterial(iX, iY, iZ - 1) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX + 1, iY, iZ) == 1) {
            if (getMaterial(iX, iY, iZ + 1) == 1) {
              if (getMaterial(iX, iY, iZ - 1) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
            }
          }
          if (getMaterial(iX - 1, iY, iZ) == 1) {
            if (getMaterial(iX, iY, iZ + 1) == 1) {
              if (getMaterial(iX, iY, iZ - 1) == 1) {
                get(iX, iY, iZ) = 1;
                count2++;
              }
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
bool BlockGeometryStructure3D<T>::find(int material, unsigned offsetX, unsigned offsetY,
                                       unsigned offsetZ, int& foundX, int& foundY, int& foundZ)
{

  bool found = false;
  for (foundX = 0; foundX < getNx(); foundX++) {
    for (foundY = 0; foundY < getNy(); foundY++) {
      for (foundZ = 0; foundZ < getNz(); foundZ++) {
        found = check(material, foundX, foundY, foundZ, offsetX, offsetY, offsetZ);
        if (found) {
          return found;
        }
      }
    }
  }
  return found;
}

template<typename T>
bool BlockGeometryStructure3D<T>::check(int material, int iX, int iY, int iZ,
                                        unsigned offsetX, unsigned offsetY, unsigned offsetZ)
{
  bool found = true;
  for (int iOffsetX = -offsetX; iOffsetX <= (int) offsetX; ++iOffsetX) {
    for (int iOffsetY = -offsetY; iOffsetY <= (int) offsetY; ++iOffsetY) {
      for (int iOffsetZ = -offsetZ; iOffsetZ <= (int) offsetZ; ++iOffsetZ) {
        if (getMaterial(iX + iOffsetX, iY + iOffsetY, iZ + iOffsetZ) != material) {
          found = false;
        }
      }
    }
  }
  return found;
}

template<typename T>
bool BlockGeometryStructure3D<T>::checkForErrors(bool verbose) const
{
  bool error = false;

  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        if (get(iX, iY, iZ) == 0) {
          if (   getMaterial(iX - 1, iY,     iZ    ) == 1
                 || getMaterial(iX,     iY - 1, iZ    ) == 1
                 || getMaterial(iX,     iY,     iZ - 1) == 1
                 || getMaterial(iX - 1, iY - 1, iZ    ) == 1
                 || getMaterial(iX,     iY - 1, iZ - 1) == 1
                 || getMaterial(iX - 1, iY,     iZ - 1) == 1
                 || getMaterial(iX - 1, iY - 1, iZ - 1) == 1
                 || getMaterial(iX + 1, iY,     iZ    ) == 1
                 || getMaterial(iX,     iY + 1, iZ    ) == 1
                 || getMaterial(iX,     iY,     iZ + 1) == 1
                 || getMaterial(iX + 1, iY + 1, iZ    ) == 1
                 || getMaterial(iX,     iY + 1, iZ + 1) == 1
                 || getMaterial(iX + 1, iY,     iZ + 1) == 1
                 || getMaterial(iX + 1, iY + 1, iZ + 1) == 1
                 || getMaterial(iX - 1, iY + 1, iZ    ) == 1
                 || getMaterial(iX + 1, iY - 1, iZ    ) == 1
                 || getMaterial(iX,     iY - 1, iZ + 1) == 1
                 || getMaterial(iX,     iY + 1, iZ - 1) == 1
                 || getMaterial(iX - 1, iY,     iZ + 1) == 1
                 || getMaterial(iX + 1, iY,     iZ - 1) == 1
                 || getMaterial(iX + 1, iY + 1, iZ - 1) == 1
                 || getMaterial(iX + 1, iY - 1, iZ - 1) == 1
                 || getMaterial(iX + 1, iY - 1, iZ + 1) == 1
                 || getMaterial(iX - 1, iY + 1, iZ + 1) == 1
                 || getMaterial(iX - 1, iY - 1, iZ + 1) == 1
                 || getMaterial(iX - 1, iY + 1, iZ - 1) == 1) {
            error = true;
          }
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
void BlockGeometryStructure3D<T>::rename(int fromM, int toM)
{

  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        if (get(iX, iY, iZ) == fromM) {
          get(iX, iY, iZ) = toM;
        }
      }
    }
  }
}

template<typename T>
void BlockGeometryStructure3D<T>::rename(int fromM, int toM, IndicatorF3D<T>& condition)
{
  T physR[3];
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        if (get(iX, iY, iZ) == fromM) {
          getPhysR(physR, iX,iY,iZ);
          bool inside[1];
          condition(inside, physR);
          if (inside[0]) {
            get(iX, iY, iZ) = toM;
          }
        }
      }
    }
  }
}

template<typename T>
void BlockGeometryStructure3D<T>::rename(int fromM, int toM, unsigned offsetX,
    unsigned offsetY, unsigned offsetZ)
{
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        if (get(iX, iY, iZ) == fromM) {
          bool found = true;
          for (int iOffsetX = -offsetX; iOffsetX <= (int) offsetX; ++iOffsetX) {
            for (int iOffsetY = -offsetY; iOffsetY <= (int) offsetY; ++iOffsetY) {
              for (int iOffsetZ = -offsetZ; iOffsetZ <= (int) offsetZ; ++iOffsetZ) {
                if (getMaterial(iX + iOffsetX, iY + iOffsetY, iZ + iOffsetZ) != fromM) {
                  if (getMaterial(iX + iOffsetX, iY + iOffsetY, iZ + iOffsetZ) != 1245) {
                    found = false;
                  }
                }
              }
            }
          }
          if (found) {
            get(iX, iY, iZ) = 1245;
          }
        }
      }
    }
  }
  rename(1245,toM);
}

template<typename T>
void BlockGeometryStructure3D<T>::rename(int fromM, int toM, int testM,
    std::vector<int> testDirection)
{
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        if (get(iX, iY, iZ) == fromM) {

          // flag that indicates the renaming of the current voxel, valid voxels are not renamed
          bool isValid = true;
          for (int iOffsetX = std::min(testDirection[0],0); iOffsetX <= std::max(testDirection[0],0); ++iOffsetX) {
            for (int iOffsetY = std::min(testDirection[1],0); iOffsetY <= std::max(testDirection[1],0); ++iOffsetY) {
              for (int iOffsetZ = std::min(testDirection[2],0); iOffsetZ <= std::max(testDirection[2],0); ++iOffsetZ) {
                if (iOffsetX!=0 || iOffsetY!=0 || iOffsetZ!=0) {
                  if (getMaterial(iX + iOffsetX, iY + iOffsetY, iZ + iOffsetZ) != testM) {
                    isValid = false;
                  }
                }
              }
            }
          }
          if (isValid) {
            get(iX, iY, iZ) = toM;
          }
        }
      }
    }
  }
}


template<typename T>
void BlockGeometryStructure3D<T>::rename(int fromM, int toM, int fluidM, IndicatorF3D<T>& condition,
    std::vector<int> discreteNormal)
{
  rename(fromM, toM, condition);
  std::vector<int> testDirection(discreteNormal);
  T physR[3];
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        if (get(iX, iY, iZ) == toM) {
          getPhysR(physR, iX,iY,iZ);
          bool inside[1];
          condition(inside, physR);
          if (inside[0]) {
            if (getMaterial(iX+testDirection[0],iY+testDirection[1],iZ+testDirection[2])!=fluidM ||
                getMaterial(iX+2*testDirection[0],iY+2*testDirection[1],iZ+2*testDirection[2])!=fluidM ||
                getMaterial(iX-testDirection[0],iY-testDirection[1],iZ-testDirection[2])!=0 ) {
              get(iX, iY, iZ) = fromM;
            }
          }
        }
      }
    }
  }
}

template<typename T>
void BlockGeometryStructure3D<T>::rename(int fromM, int toM, int fluidM, IndicatorF3D<T>& condition)
{
  rename(fromM, toM, condition);
  std::vector<int> testDirection = getStatistics().computeDiscreteNormal(toM);
  T physR[3];
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        if (get(iX, iY, iZ) == toM) {
          getPhysR(physR, iX,iY,iZ);
          bool inside[1];
          condition(inside, physR);
          if (inside[0]) {
            if (getMaterial(iX+testDirection[0],iY+testDirection[1],iZ+testDirection[2])!=fluidM ||
                getMaterial(iX+2*testDirection[0],iY+2*testDirection[1],iZ+2*testDirection[2])!=fluidM ||
                getMaterial(iX-testDirection[0],iY-testDirection[1],iZ-testDirection[2])!=0 ) {
              get(iX, iY, iZ) = fromM;
            }
          }
        }
      }
    }
  }
}

template<typename T>
void BlockGeometryStructure3D<T>::copyMaterialLayer(IndicatorF3D<T>& condition, int discreteNormal[3], int numberOfLayers)
{
  T physR[3];
  for (int iX = 0; iX < getNx(); iX++) {
    for (int iY = 0; iY < getNy(); iY++) {
      for (int iZ = 0; iZ < getNz(); iZ++) {
        getPhysR(physR, iX,iY,iZ);
        bool inside[1];
        condition(inside, physR);
        if (inside[0]) {
          for (int i = 0; i < numberOfLayers; i++) {
            if (0 <= iX + i * discreteNormal[0] && iX + i * discreteNormal[0] < getNx() &&
                0 <= iY + i * discreteNormal[1] && iY + i * discreteNormal[1] < getNy() &&
                0 <= iZ + i * discreteNormal[2] && iZ + i * discreteNormal[2] < getNz()) {
              get(iX + i * discreteNormal[0], iY + i * discreteNormal[1], iZ + i * discreteNormal[2]) = get(iX, iY, iZ);
            }
          }
        }
      }
    }
  }
}


template<typename T>
void BlockGeometryStructure3D<T>::regionGrowing(int fromM, int toM, int seedX, int seedY,
    int seedZ, int offsetX, int offsetY, int offsetZ, std::map<std::vector<int>, int>* tmp)
{
  std::map<std::vector<int>, int> tmp2;
  bool firstCall = false;
  if (tmp == nullptr) {
    tmp = &tmp2;
    firstCall = true;
  }

  if (getMaterial(seedX, seedY, seedZ) == fromM) {
    std::vector<int> found;
    found.push_back(seedX);
    found.push_back(seedY);
    found.push_back(seedZ);
    if (tmp->count(found) == 0) {
      (*tmp)[found] = 2;
      if (offsetX != 0) {
        regionGrowing(fromM, toM, seedX + 1, seedY, seedZ, offsetX,
                      offsetY, offsetZ, tmp);
        regionGrowing(fromM, toM, seedX - 1, seedY, seedZ, offsetX,
                      offsetY, offsetZ, tmp);
      }
      if (offsetY != 0) {
        regionGrowing(fromM, toM, seedX, seedY + 1, seedZ, offsetX,
                      offsetY, offsetZ, tmp);
        regionGrowing(fromM, toM, seedX, seedY - 1, seedZ, offsetX,
                      offsetY, offsetZ, tmp);
      }
      if (offsetZ != 0) {
        regionGrowing(fromM, toM, seedX, seedY, seedZ + 1, offsetX,
                      offsetY, offsetZ, tmp);
        regionGrowing(fromM, toM, seedX, seedY, seedZ - 1, offsetX,
                      offsetY, offsetZ, tmp);
      }
    }
  }
  if (firstCall) {
    std::map<std::vector<int>, int>::iterator iter;
    for (iter = tmp->begin(); iter != tmp->end(); iter++) {
      get((iter->first)[0],(iter->first)[1],(iter->first)[2]) = toM;
    }
  }
  return;
}

} // namespace olb

#endif
