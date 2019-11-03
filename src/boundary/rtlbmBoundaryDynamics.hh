/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Albert Mink
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

#ifndef RTLBM_BOUNDARY_DYNAMICS_HH
#define RTLBM_BOUNDARY_DYNAMICS_HH

#include "rtlbmBoundaryDynamics.h"
#include "dynamics/lbHelpers.h"

namespace olb {



// flat diffuse
template<typename T, typename DESCRIPTOR, int direction, int orientation>
RtlbmDiffuseBoundaryDynamics<T,DESCRIPTOR,direction,orientation>::RtlbmDiffuseBoundaryDynamics( T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BasicDynamics<T,DESCRIPTOR>(momenta_)
{
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
T RtlbmDiffuseBoundaryDynamics<T,DESCRIPTOR,direction,orientation>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return descriptors::t<T,DESCRIPTOR>(iPop)*rho - descriptors::t<T,DESCRIPTOR>(iPop);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void RtlbmDiffuseBoundaryDynamics<T,DESCRIPTOR,direction,orientation>::collide(Cell<T,DESCRIPTOR>& cell,LatticeStatistics<T>& statistics)
{
  typedef DESCRIPTOR L;
  T dirichletTemperature = this->_momenta.computeRho(cell);
  std::vector<int> const missing_iPop = util::subIndexOutgoing<L,direction,orientation>();
  // compute summ of weights for all missing directions
  double sumWeights = 0;
  for ( int i : missing_iPop ) {
    sumWeights += descriptors::t<T,L>(i);
  }
  // construct missing directions such that 0th moment equals emposed dirichletTemperature
  for ( int i : missing_iPop ) {
    cell[i] = descriptors::t<T,L>(i)*dirichletTemperature/sumWeights - descriptors::t<T,L>(i);
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
T RtlbmDiffuseBoundaryDynamics<T,DESCRIPTOR,direction,orientation>::getOmega() const
{
  return T(-1);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void RtlbmDiffuseBoundaryDynamics<T,DESCRIPTOR,direction,orientation>::setOmega(T omega_)
{
}

// edge diffuse
template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
RtlbmDiffuseEdgeBoundaryDynamics<T,DESCRIPTOR,plane,normal1,normal2>::RtlbmDiffuseEdgeBoundaryDynamics( T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BasicDynamics<T,DESCRIPTOR>(momenta_)
{
}

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
T RtlbmDiffuseEdgeBoundaryDynamics<T,DESCRIPTOR,plane,normal1,normal2>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return descriptors::t<T,DESCRIPTOR>(iPop)*rho - descriptors::t<T,DESCRIPTOR>(iPop);
}

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
void RtlbmDiffuseEdgeBoundaryDynamics<T,DESCRIPTOR,plane,normal1,normal2>::collide(Cell<T,DESCRIPTOR>& cell,LatticeStatistics<T>& statistics)
{
  typedef DESCRIPTOR L;
  T dirichletTemperature = this->_momenta.computeRho(cell);
  std::vector<int> missing_iPop = util::subIndexOutgoing3DonEdges<L,plane,normal1,normal2>();
  // compute summ of weights for all missing directions
  double sumWeights = 0;
  for ( int i : missing_iPop ) {
    sumWeights += descriptors::t<T,L>(i);
  }
  // construct missing directions such that 0th moment equals emposed dirichletTemperature
  for ( int i : missing_iPop ) {
    cell[i] = descriptors::t<T,L>(i)*dirichletTemperature/sumWeights - descriptors::t<T,L>(i);
  }
}

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
T RtlbmDiffuseEdgeBoundaryDynamics<T,DESCRIPTOR,plane,normal1,normal2>::getOmega() const
{
  return T(-1);
}

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
void RtlbmDiffuseEdgeBoundaryDynamics<T,DESCRIPTOR,plane,normal1,normal2>::setOmega(T omega_)
{
}

// corner diffuse
template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
RtlbmDiffuseCornerBoundaryDynamics<T,DESCRIPTOR,xNormal,yNormal,zNormal>::RtlbmDiffuseCornerBoundaryDynamics( T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BasicDynamics<T,DESCRIPTOR>(momenta_)
{
}

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
T RtlbmDiffuseCornerBoundaryDynamics<T,DESCRIPTOR,xNormal,yNormal,zNormal>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return descriptors::t<T,DESCRIPTOR>(iPop)*rho - descriptors::t<T,DESCRIPTOR>(iPop);
}

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
void RtlbmDiffuseCornerBoundaryDynamics<T,DESCRIPTOR,xNormal,yNormal,zNormal>::collide(Cell<T,DESCRIPTOR>& cell,LatticeStatistics<T>& statistics)
{
  typedef DESCRIPTOR L;
  T dirichletTemperature = this->_momenta.computeRho(cell);
  std::vector<int> const missing_iPop = util::subIndexOutgoing3DonCorners<L,xNormal,yNormal,zNormal>();
  // compute summ of weights for all missing directions
  double sumWeights = 0;
  for ( int i : missing_iPop ) {
    sumWeights += descriptors::t<T,L>(i);
  }
  // construct missing directions such that 0th moment equals emposed dirichletTemperature
  for ( int i : missing_iPop ) {
    cell[i] = descriptors::t<T,L>(i)*dirichletTemperature/sumWeights - descriptors::t<T,L>(i);
  }
}

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
T RtlbmDiffuseCornerBoundaryDynamics<T,DESCRIPTOR,xNormal,yNormal,zNormal>::getOmega() const
{
  return T(-1);
}

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
void RtlbmDiffuseCornerBoundaryDynamics<T,DESCRIPTOR,xNormal,yNormal,zNormal>::setOmega(T omega_)
{
}




// flat diffuse constant density
template<typename T, typename DESCRIPTOR, int direction, int orientation>
RtlbmDiffuseConstBoundaryDynamics<T,DESCRIPTOR,direction,orientation>::RtlbmDiffuseConstBoundaryDynamics( T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BasicDynamics<T,DESCRIPTOR>(momenta_)
{
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
T RtlbmDiffuseConstBoundaryDynamics<T,DESCRIPTOR,direction,orientation>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return descriptors::t<T,DESCRIPTOR>(iPop)*rho - descriptors::t<T,DESCRIPTOR>(iPop);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void RtlbmDiffuseConstBoundaryDynamics<T,DESCRIPTOR,direction,orientation>::collide(Cell<T,DESCRIPTOR>& cell,LatticeStatistics<T>& statistics)
{
  // For direction i \in I_in define
  // cell_i = w_i * dirichlet/sumWeights - w_i
  // For direction i \in I_out defube
  // cell_i = - w_i
  // This construction yields
  // sum_{i=0}^{q-1} cell_i == dirichlet - 1

  typedef DESCRIPTOR L;
  // shift all: cell_i = f_i - weight_i
  for ( int iPop = 0; iPop < L::q; ++iPop ) {
    cell[iPop] = - descriptors::t<T,L>(iPop);
  }

  std::vector<int> const missing_iPop = util::subIndexOutgoing<L,direction,orientation>();
  // compute summ of weights for all missing directions
  double sumWeights = 0;
  for ( int i : missing_iPop ) {
    sumWeights += descriptors::t<T,L>(i);
  }
  // construct missing directions such that 0th moment equals emposed dirichletTemperature
  T dirichletTemperature = this->_momenta.computeRho(cell);
  for ( int i : missing_iPop ) {
    cell[i] = descriptors::t<T,L>(i)*dirichletTemperature/sumWeights - descriptors::t<T,L>(i);
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
T RtlbmDiffuseConstBoundaryDynamics<T,DESCRIPTOR,direction,orientation>::getOmega() const
{
  return T(-1);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void RtlbmDiffuseConstBoundaryDynamics<T,DESCRIPTOR,direction,orientation>::setOmega(T omega_)
{
}



// edge diffuse with constant density
template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
RtlbmDiffuseConstEdgeBoundaryDynamics<T,DESCRIPTOR,plane,normal1,normal2>::RtlbmDiffuseConstEdgeBoundaryDynamics( T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BasicDynamics<T,DESCRIPTOR>(momenta_)
{
}

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
T RtlbmDiffuseConstEdgeBoundaryDynamics<T,DESCRIPTOR,plane,normal1,normal2>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return descriptors::t<T,DESCRIPTOR>(iPop)*rho - descriptors::t<T,DESCRIPTOR>(iPop);
}

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
void RtlbmDiffuseConstEdgeBoundaryDynamics<T,DESCRIPTOR,plane,normal1,normal2>::collide(Cell<T,DESCRIPTOR>& cell,LatticeStatistics<T>& statistics)
{
    // For direction i \in I_in define
    // cell_i = w_i * dirichlet/sumWeights - w_i
    // For direction i \in I_out defube
    // cell_i = - w_i
    // This construction yields
    // sum_{i=0}^{q-1} cell_i == dirichlet - 1

  typedef DESCRIPTOR L;

  // shift all: cell_i = f_i - weight_i
  for ( int iPop = 0; iPop < L::q; ++iPop ) {
    cell[iPop] = - descriptors::t<T,L>(iPop);
  }

  std::vector<int> missing_iPop = util::subIndexOutgoing3DonEdges<L,plane,normal1,normal2>();
  double sumWeights = 0;
  for ( int i : missing_iPop ) {
    sumWeights += descriptors::t<T,L>(i);
  }

  T dirichletTemperature = this->_momenta.computeRho(cell);
  for ( int i : missing_iPop ) {
    cell[i] = descriptors::t<T,L>(i)*dirichletTemperature/sumWeights - descriptors::t<T,L>(i);
  }
}

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
T RtlbmDiffuseConstEdgeBoundaryDynamics<T,DESCRIPTOR,plane,normal1,normal2>::getOmega() const
{
  return T(-1);
}

template<typename T, typename DESCRIPTOR, int plane, int normal1, int normal2>
void RtlbmDiffuseConstEdgeBoundaryDynamics<T,DESCRIPTOR,plane,normal1,normal2>::setOmega(T omega_)
{
}



// corner diffuse with constant density
template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
RtlbmDiffuseConstCornerBoundaryDynamics<T,DESCRIPTOR,xNormal,yNormal,zNormal>::RtlbmDiffuseConstCornerBoundaryDynamics( T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BasicDynamics<T,DESCRIPTOR>(momenta_)
{
}

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
T RtlbmDiffuseConstCornerBoundaryDynamics<T,DESCRIPTOR,xNormal,yNormal,zNormal>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return descriptors::t<T,DESCRIPTOR>(iPop)*rho - descriptors::t<T,DESCRIPTOR>(iPop);
}

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
void RtlbmDiffuseConstCornerBoundaryDynamics<T,DESCRIPTOR,xNormal,yNormal,zNormal>::collide(Cell<T,DESCRIPTOR>& cell,LatticeStatistics<T>& statistics)
{
    // For direction i \in I_in define
    // cell_i = w_i * dirichlet/sumWeights - w_i
    // For direction i \in I_out defube
    // cell_i = - w_i
    // This construction yields
    // sum_{i=0}^{q-1} cell_i == dirichlet - 1

  typedef DESCRIPTOR L;

  // shift all: cell_i = f_i - weight_i
  for ( int iPop = 0; iPop < L::q; ++iPop ) {
    cell[iPop] = - descriptors::t<T,L>(iPop);
  }

  std::vector<int> const missing_iPop = util::subIndexOutgoing3DonCorners<L,xNormal,yNormal,zNormal>();
  double sumWeights = 0;
  for ( int i : missing_iPop ) {
    sumWeights += descriptors::t<T,L>(i);
  }

  T dirichletTemperature = this->_momenta.computeRho(cell);
  for ( int i : missing_iPop ) {
    cell[i] = descriptors::t<T,L>(i)*dirichletTemperature/sumWeights - descriptors::t<T,L>(i);
  }
}

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
T RtlbmDiffuseConstCornerBoundaryDynamics<T,DESCRIPTOR,xNormal,yNormal,zNormal>::getOmega() const
{
  return T(-1);
}

template<typename T, typename DESCRIPTOR, int xNormal, int yNormal, int zNormal>
void RtlbmDiffuseConstCornerBoundaryDynamics<T,DESCRIPTOR,xNormal,yNormal,zNormal>::setOmega(T omega_)
{
}


// directed wall
template<typename T, typename DESCRIPTOR, int direction, int orientation>
RtlbmDirectedBoundaryDynamics<T,DESCRIPTOR,direction,orientation>::RtlbmDirectedBoundaryDynamics( T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BasicDynamics<T,DESCRIPTOR>(momenta_)
{
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
T RtlbmDirectedBoundaryDynamics<T,DESCRIPTOR,direction,orientation>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return lbHelpers<T,DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void RtlbmDirectedBoundaryDynamics<T,DESCRIPTOR,direction,orientation>::collide(Cell<T,DESCRIPTOR>& cell,LatticeStatistics<T>& statistics)
{
  typedef DESCRIPTOR L;
  T dirichletTemperature = this->_momenta.computeRho(cell);

  for( int iPop = 0; iPop < L::q; ++iPop ) {
    cell[iPop] = - descriptors::t<T,L>(iPop);
  }

  std::vector<int> const missingDiagonal = util::subIndexOutgoing<L,direction,orientation>();
  for ( int i : missingDiagonal ) {
    // compute norm of c_iPopMissing
    // is direction axis parallel
    if ( util::normSqr<int,L::d>(descriptors::c<L>(i)) == 1 ) {
      cell[i] = dirichletTemperature - descriptors::t<T,L>(i);
    }
  }
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
T RtlbmDirectedBoundaryDynamics<T,DESCRIPTOR,direction,orientation>::getOmega() const
{
  return T(-1);
}

template<typename T, typename DESCRIPTOR, int direction, int orientation>
void RtlbmDirectedBoundaryDynamics<T,DESCRIPTOR,direction,orientation>::setOmega(T omega_)
{
}







}  // namespace olb


#endif
