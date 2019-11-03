/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
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

#ifndef ADVECTION_DIFFUSION_BOUNDARIES_HH
#define ADVECTION_DIFFUSION_BOUNDARIES_HH

#include "advectionDiffusionBoundaries.h"
#include "dynamics/latticeDescriptors.h"
#include "core/util.h"
#include "dynamics/lbHelpers.h"

namespace olb {

 

//==================================================================================================
//==================== For regularized Advection Diffusion Boundary Condition ======================
//============================================================================================


// For flat Walls

template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
AdvectionDiffusionBoundariesDynamics<T,DESCRIPTOR,Dynamics,direction,orientation>::
AdvectionDiffusionBoundariesDynamics( T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BasicDynamics<T,DESCRIPTOR>(momenta_), boundaryDynamics(omega_, momenta_)
{
}

template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
T AdvectionDiffusionBoundariesDynamics<T,DESCRIPTOR,Dynamics,direction,orientation>::
computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return lbHelpers<T,DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
}


template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
void AdvectionDiffusionBoundariesDynamics<T,DESCRIPTOR,Dynamics,direction,orientation>::
collide(Cell<T,DESCRIPTOR>& cell,LatticeStatistics<T>& statistics)
{
  typedef DESCRIPTOR L;
  typedef lbHelpers<T,DESCRIPTOR> lbH;

  T dirichletTemperature = this->_momenta.computeRho(cell);
  T* u = cell.template getFieldPointer<descriptors::VELOCITY>();

  std::vector<int> unknownIndexes = util::subIndexOutgoing<L, direction,
      orientation>();
  std::vector<int> knownIndexes = util::remainingIndexes<L>(unknownIndexes);

  int missingNormal = 0;

  if ((L::d == 3 && L::q == 7)||(L::d == 2 && L::q == 5)) {
    T sum = T();
    for (unsigned i = 0; i < knownIndexes.size(); ++i) {
      sum += cell[knownIndexes[i]];
    }

    T difference = dirichletTemperature - (T) 1 - sum; // on cell there are non-shiftet values -> temperature has to be changed

    // here I know all missing and non missing f_i
    for (unsigned i = 0; i < unknownIndexes.size(); ++i) {
      int numOfNonNullComp = 0;
      for (int iDim = 0; iDim < L::d; ++iDim) {
        numOfNonNullComp += abs(descriptors::c<L>(unknownIndexes[i],iDim));
      }
      if (numOfNonNullComp == 1) {
        missingNormal = unknownIndexes[i];
        // here missing diagonal directions are erased
        // just the normal direction stays (D3Q7)
        unknownIndexes.erase(unknownIndexes.begin() + i);
        break;

      }
    }
    cell[missingNormal] = difference; // on cell there are non-shiftet values -> temperature has to be changed
    boundaryDynamics.collide(cell, statistics); // only for D3Q7
  } else {
    // part for q=19 copied from AdvectionDiffusionEdgesDynamics.collide()
    // but here just all missing directions, even at border of inlet area
    // has to be checked!
    for (unsigned iteratePop = 0; iteratePop < unknownIndexes.size();
        ++iteratePop) {
      cell[unknownIndexes[iteratePop]] =
          lbH::equilibriumFirstOrder(unknownIndexes[iteratePop], dirichletTemperature, u)
              - (cell[util::opposite<L>(unknownIndexes[iteratePop])]
                  - lbH::equilibriumFirstOrder(
                      util::opposite<L>(unknownIndexes[iteratePop]),
                      dirichletTemperature, u));
    }
  }
}

template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
T AdvectionDiffusionBoundariesDynamics<T,DESCRIPTOR,Dynamics,direction,orientation>::
getOmega() const
{
  return boundaryDynamics.getOmega();
}

template<typename T, typename DESCRIPTOR, typename Dynamics, int direction, int orientation>
void AdvectionDiffusionBoundariesDynamics<T,DESCRIPTOR,Dynamics,direction,orientation>::
setOmega(T omega_)
{
  boundaryDynamics.setOmega(omega_);
}

//=================================================================
// For 2D Corners with regularized Dynamic ==============================================
//=================================================================
template<typename T, typename DESCRIPTOR, typename Dynamics, int xNormal, int yNormal>
AdvectionDiffusionCornerDynamics2D<T,DESCRIPTOR,Dynamics,xNormal,yNormal>::AdvectionDiffusionCornerDynamics2D(
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BasicDynamics<T,DESCRIPTOR>(momenta_),
    boundaryDynamics(omega_, momenta_)
{
}

template<typename T, typename DESCRIPTOR, typename Dynamics,  int xNormal, int yNormal>
T AdvectionDiffusionCornerDynamics2D<T,DESCRIPTOR,Dynamics,xNormal,yNormal>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return lbHelpers<T,DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
}


template<typename T, typename DESCRIPTOR, typename Dynamics,  int xNormal, int yNormal>
void AdvectionDiffusionCornerDynamics2D<T,DESCRIPTOR,Dynamics,xNormal,yNormal>::collide(Cell<T,DESCRIPTOR>& cell,LatticeStatistics<T>& statistics)
{
  typedef DESCRIPTOR L;
  typedef lbHelpers<T,DESCRIPTOR> lbH;

  T temperature = this->_momenta.computeRho(cell);
  T* u = cell.template getFieldPointer<descriptors::VELOCITY>();
  // I need to get Missing information on the corners !!!!
  std::vector<int> unknownIndexes = util::subIndexOutgoing2DonCorners<L,xNormal,yNormal>();
  // here I know all missing and non missing f_i


  // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
  // Given the rule f_i_neq = -f_opposite(i)_neq
  // I have the right number of equations for the number of unknowns using these lattices

  for (unsigned iPop = 0; iPop < unknownIndexes.size(); ++iPop) {
    cell[unknownIndexes[iPop]] = lbH::equilibriumFirstOrder(unknownIndexes[iPop], temperature, u)
                                 -(cell[util::opposite<L>(unknownIndexes[iPop])]
                                 - lbH::equilibriumFirstOrder(util::opposite<L>(unknownIndexes[iPop]), temperature, u) ) ;
  }

  // Once all the f_i are known, I can call the collision for the Regularized Model.
  boundaryDynamics.collide(cell, statistics);

}

template<typename T, typename DESCRIPTOR, typename Dynamics, int xNormal, int yNormal>
T AdvectionDiffusionCornerDynamics2D<T,DESCRIPTOR,Dynamics,xNormal,yNormal>::getOmega() const
{
  return boundaryDynamics.getOmega();
}

template<typename T, typename DESCRIPTOR, typename Dynamics, int xNormal, int yNormal>
void AdvectionDiffusionCornerDynamics2D<T,DESCRIPTOR,Dynamics,xNormal,yNormal>::setOmega(T omega_)
{
  boundaryDynamics.setOmega(omega_);
}



//=================================================================
// For 3D Corners with regularized Dynamic ==============================================
//=================================================================
template<typename T, typename DESCRIPTOR, typename Dynamics, int xNormal, int yNormal, int zNormal>
AdvectionDiffusionCornerDynamics3D<T,DESCRIPTOR,Dynamics,xNormal,yNormal,zNormal>::AdvectionDiffusionCornerDynamics3D(
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BasicDynamics<T,DESCRIPTOR>(momenta_),
    boundaryDynamics(omega_, momenta_)
{
}

template<typename T, typename DESCRIPTOR, typename Dynamics,  int xNormal, int yNormal, int zNormal>
T AdvectionDiffusionCornerDynamics3D<T,DESCRIPTOR,Dynamics,xNormal,yNormal,zNormal>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return lbHelpers<T, DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
}


template<typename T, typename DESCRIPTOR, typename Dynamics,  int xNormal, int yNormal, int zNormal>
void AdvectionDiffusionCornerDynamics3D<T,DESCRIPTOR,Dynamics,xNormal,yNormal,zNormal>::collide(Cell<T,DESCRIPTOR>& cell,LatticeStatistics<T>& statistics)
{
  typedef DESCRIPTOR L;
  typedef lbHelpers<T,DESCRIPTOR> lbH;

  T temperature = this->_momenta.computeRho(cell);
  T* u = cell.template getFieldPointer<descriptors::VELOCITY>();
  // I need to get Missing information on the corners !!!!
  std::vector<int> unknownIndexes = util::subIndexOutgoing3DonCorners<L,xNormal,yNormal,zNormal>();
  // here I know all missing and non missing f_i


  // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
  // Given the rule f_i_neq = -f_opposite(i)_neq
  // I have the right number of equations for the number of unknowns using these lattices

  for (unsigned iPop = 0; iPop < unknownIndexes.size(); ++iPop) {
    cell[unknownIndexes[iPop]] = lbH::equilibriumFirstOrder(unknownIndexes[iPop], temperature, u)
                                 -(cell[util::opposite<L>(unknownIndexes[iPop])]
                                   - lbH::equilibriumFirstOrder(util::opposite<L>(unknownIndexes[iPop]), temperature, u) ) ;
  }

  // Once all the f_i are known, I can call the collision for the Regularized Model.
  boundaryDynamics.collide(cell, statistics);

}

template<typename T, typename DESCRIPTOR, typename Dynamics, int xNormal, int yNormal, int zNormal>
T AdvectionDiffusionCornerDynamics3D<T,DESCRIPTOR,Dynamics,xNormal,yNormal,zNormal>::getOmega() const
{
  return boundaryDynamics.getOmega();
}

template<typename T, typename DESCRIPTOR, typename Dynamics, int xNormal, int yNormal, int zNormal>
void AdvectionDiffusionCornerDynamics3D<T,DESCRIPTOR,Dynamics,xNormal,yNormal,zNormal>::setOmega(T omega_)
{
  boundaryDynamics.setOmega(omega_);
}

//=================================================================
// For 3D Edges with regularized Dynamic ==============================================
//=================================================================
template<typename T, typename DESCRIPTOR, typename Dynamics, int plane, int normal1, int normal2>
AdvectionDiffusionEdgesDynamics<T,DESCRIPTOR,Dynamics,plane,normal1, normal2>::AdvectionDiffusionEdgesDynamics(
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BasicDynamics<T,DESCRIPTOR>(momenta_),
    boundaryDynamics(omega_, momenta_)
{
}

template<typename T, typename DESCRIPTOR, typename Dynamics, int plane, int normal1, int normal2>
T AdvectionDiffusionEdgesDynamics<T,DESCRIPTOR,Dynamics,plane,normal1, normal2>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return lbHelpers<T,DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
}


template<typename T, typename DESCRIPTOR, typename Dynamics,  int plane, int normal1, int normal2>
void AdvectionDiffusionEdgesDynamics<T,DESCRIPTOR,Dynamics,plane,normal1, normal2>::collide(Cell<T,DESCRIPTOR>& cell,LatticeStatistics<T>& statistics)
{
  typedef DESCRIPTOR L;
  typedef lbHelpers<T,DESCRIPTOR> lbH;

  T temperature = this->_momenta.computeRho(cell);
  T* u = cell.template getFieldPointer<descriptors::VELOCITY>();
  // I need to get Missing information on the corners !!!!
  std::vector<int> unknownIndexes = util::subIndexOutgoing3DonEdges<L,plane,normal1, normal2>();
  // here I know all missing and non missing f_i


  // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
  // Given the rule f_i_neq = -f_opposite(i)_neq
  // I have the right number of equations for the number of unknowns using these lattices

  for (unsigned iPop = 0; iPop < unknownIndexes.size(); ++iPop) {
    cell[unknownIndexes[iPop]] = lbH::equilibriumFirstOrder(unknownIndexes[iPop], temperature, u)
                                 -(cell[util::opposite<L>(unknownIndexes[iPop])]
                                 - lbH::equilibriumFirstOrder(util::opposite<L>(unknownIndexes[iPop]), temperature, u) ) ;
  }

  // Once all the f_i are known, I can call the collision for the Regularized Model.
  boundaryDynamics.collide(cell, statistics);

}

template<typename T, typename DESCRIPTOR, typename Dynamics, int plane, int normal1, int normal2>
T AdvectionDiffusionEdgesDynamics<T,DESCRIPTOR,Dynamics,plane,normal1, normal2>::getOmega() const
{
  return boundaryDynamics.getOmega();
}

template<typename T, typename DESCRIPTOR, typename Dynamics, int plane, int normal1, int normal2>
void AdvectionDiffusionEdgesDynamics<T,DESCRIPTOR,Dynamics,plane,normal1, normal2>::setOmega(T omega_)
{
  boundaryDynamics.setOmega(omega_);
}




}  // namespace olb




#endif
