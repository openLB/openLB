/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007, 2008 Jonas Latt
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

#ifndef LATTICE_STATISTICS_HH
#define LATTICE_STATISTICS_HH

#include "blockLattice2D.h"
#include "blockLattice3D.h"
#include "latticeStatistics.h"
#include <cmath>
#include <numeric>
#include <limits>
#include "util.h"

namespace olb {

////////////////////// Class LatticeStatistics /////////////////

template<typename T>
LatticeStatistics<T>::LatticeStatistics() : clout(std::cout,"LatticeStatistics")
{
  initialize();
}

template<typename T>
LatticeStatistics<T>::~LatticeStatistics()
{
}

template<typename T>
void LatticeStatistics<T>::reset()
{
  // avoid division by zero
  if (tmpNumCells == 0) {
    for (unsigned iVect=0; iVect<averageVect.size(); ++iVect) {
      averageVect[iVect] = T();
    }
    for (unsigned iVect=0; iVect<sumVect.size(); ++iVect) {
      sumVect[iVect] = T();
    }
    for (unsigned iVect=0; iVect<minVect.size(); ++iVect) {
      minVect[iVect] = T();
    }
    for (unsigned iVect=0; iVect<maxVect.size(); ++iVect) {
      maxVect[iVect] = T();
    }
    numCells = 0;
    firstCall = false;
  } else {
    // The average density is actually used in the "ConstRhoBgk" model.
    // Depending on the simulation setup, it is possible that it has
    // a nonsensical value before the simulation is started. For this
    // and similar cases, averages are initialized to 1.
    if (firstCall) {
      for (unsigned iVect=0; iVect<averageVect.size(); ++iVect) {
        averageVect[iVect] = (T)1;
      }
      firstCall = false;
    } else {
      for (unsigned iVect=0; iVect<averageVect.size(); ++iVect) {
        averageVect[iVect] = tmpAv[iVect] / (T)tmpNumCells;
      }
    }
    for (unsigned iVect=0; iVect<sumVect.size(); ++iVect) {
      sumVect[iVect] = tmpSum[iVect];
    }
    for (unsigned iVect=0; iVect<minVect.size(); ++iVect) {
      minVect[iVect] = tmpMin[iVect];
    }
    for (unsigned iVect=0; iVect<maxVect.size(); ++iVect) {
      maxVect[iVect] = tmpMax[iVect];
    }
    averageVect[avEnergy] *= (T)0.5; // energy is 0.5 *uSqr
    maxVect[maxU]         = sqrt(maxVect[maxU]); // u is sqrt(uSqr)
    numCells              = tmpNumCells;
  }

  for (unsigned iVect=0; iVect<averageVect.size(); ++iVect) {
    tmpAv[iVect]    = T();
  }
  for (unsigned iVect=0; iVect<sumVect.size(); ++iVect) {
    tmpSum[iVect]   = T();
  }
  for (unsigned iVect=0; iVect<minVect.size(); ++iVect) {
    tmpMin[iVect]   = std::numeric_limits<T>::max();
  }
  for (unsigned iVect=0; iVect<maxVect.size(); ++iVect) {
    tmpMax[iVect]   = std::numeric_limits<T>::min();
  }

  tmpNumCells     = 0;
}

template<typename T>
void LatticeStatistics<T>::reset (
  T average_rho_, T average_energy_, T maxU_, size_t numCells_ )
{
  averageVect[avRho]    = average_rho_;
  averageVect[avEnergy] = average_energy_;
  maxVect[maxU]         = maxU_;
  numCells              = numCells_;

  tmpAv[avRho]    = T();
  tmpAv[avEnergy] = T();
  tmpMax[maxU]    = T();
  tmpNumCells     = 0;
}

template<typename T>
void LatticeStatistics<T>::initialize()
{
  tmpAv.resize(2);
  averageVect.resize(2);
  tmpMax.resize(1);
  maxVect.resize(1);

  tmpAv[avRho]    = T();
  tmpAv[avEnergy] = T();
  tmpMax[maxU]    = T();
  tmpNumCells     = 0;

  averageVect[avRho]    = (T)1;
  averageVect[avEnergy] = T();
  maxVect[maxU]         = T();

  firstCall = true;

  resetTime();
}

template<typename T>
int LatticeStatistics<T>::subscribeAverage()
{
  int newSize = tmpAv.size()+1;
  tmpAv.resize(newSize);
  averageVect.resize(newSize);
  return newSize-1;
}

template<typename T>
int LatticeStatistics<T>::subscribeSum()
{
  int newSize = tmpSum.size()+1;
  tmpSum.resize(newSize);
  sumVect.resize(newSize);
  return newSize-1;
}

template<typename T>
int LatticeStatistics<T>::subscribeMin()
{
  int newSize = tmpMin.size()+1;
  tmpMin.resize(newSize);
  minVect.resize(newSize);
  return newSize-1;
}

template<typename T>
int LatticeStatistics<T>::subscribeMax()
{
  int newSize = tmpMax.size()+1;
  tmpMax.resize(newSize);
  maxVect.resize(newSize);
  return newSize-1;
}

template<typename T>
void LatticeStatistics<T>::incrementStats(T rho, T uSqr)
{
  tmpAv[avRho]    += rho;
  tmpAv[avEnergy] += uSqr;
  if (uSqr > tmpMax[maxU]) {
    tmpMax[maxU] = uSqr;
  }
  ++tmpNumCells;
}

template<typename T>
void LatticeStatistics<T>::gatherAverage(int whichAverage, T value)
{
  OLB_PRECONDITION( whichAverage < (int) tmpAv.size() );
  tmpAv[whichAverage] += value;
}

template<typename T>
void LatticeStatistics<T>::gatherSum(int whichSum, T value)
{
  OLB_PRECONDITION( whichSum < (int) tmpSum.size() );
  tmpSum[whichSum] += value;
}

template<typename T>
void LatticeStatistics<T>::gatherMin(int whichMin, T value)
{
  OLB_PRECONDITION( whichMin < (int) tmpMin.size() );
  if (value < tmpMin[whichMin]) {
    tmpMin[whichMin] = value;
  }
}

template<typename T>
void LatticeStatistics<T>::gatherMax(int whichMax, T value)
{
  OLB_PRECONDITION( whichMax < (int) tmpMax.size() );
  if (value > tmpMax[whichMax]) {
    tmpMax[whichMax] = value;
  }
}

template<typename T>
void LatticeStatistics<T>::incrementStats()
{
  ++tmpNumCells;
}

template<typename T>
T LatticeStatistics<T>::getAverageRho() const
{
  return averageVect[avRho];
}

template<typename T>
T LatticeStatistics<T>::getAverageEnergy() const
{
  return averageVect[avEnergy];
}

template<typename T>
T LatticeStatistics<T>::getMaxU() const
{
  return maxVect[maxU];
}

template<typename T>
size_t const& LatticeStatistics<T>::getNumCells() const
{
  return numCells;
}

template<typename T>
T LatticeStatistics<T>::getAverage(int whichAverage) const
{
  OLB_PRECONDITION( whichAverage < (int) tmpAv.size() );
  return averageVect[whichAverage];
}

template<typename T>
T LatticeStatistics<T>::getSum(int whichSum) const
{
  OLB_PRECONDITION( whichSum < (int) tmpSum.size() );
  return sumVect[whichSum];
}

template<typename T>
T LatticeStatistics<T>::getMin(int whichMin) const
{
  OLB_PRECONDITION( whichMin < (int) tmpMin.size() );
  return minVect[whichMin];
}

template<typename T>
T LatticeStatistics<T>::getMax(int whichMax) const
{
  OLB_PRECONDITION( whichMax < (int) tmpMax.size() );
  return maxVect[whichMax];
}

template<typename T>
std::vector<T>& LatticeStatistics<T>::getAverageVect()
{
  return averageVect;
}

template<typename T>
std::vector<T>& LatticeStatistics<T>::getSumVect()
{
  return sumVect;
}

template<typename T>
std::vector<T>& LatticeStatistics<T>::getMinVect()
{
  return minVect;
}

template<typename T>
std::vector<T>& LatticeStatistics<T>::getMaxVect()
{
  return maxVect;
}

template<typename T>
void LatticeStatistics<T>::incrementTime()
{
  ++latticeTime;
}

template<typename T>
void LatticeStatistics<T>::resetTime(size_t value)
{
  latticeTime=value;
}

template<typename T>
size_t LatticeStatistics<T>::getTime() const
{
  return latticeTime;
}

template<typename T>
void LatticeStatistics<T>::print(int iterationStep, T physicalTime) const
{
  clout
      << "step=" << iterationStep << "; "
      << "t=" << physicalTime << "; "
      << "uMax=" << getMaxU() << "; "
      << "avEnergy=" << getAverageEnergy() << "; "
      << "avRho=" << getAverageRho()
      << std::endl;
}

}  // namespace olb

#endif
