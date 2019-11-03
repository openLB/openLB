/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Marc Haussmann
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

#ifndef WALLFUNCTION_BOUNDARY_POST_PROCESSORS_3D_HH
#define WALLFUNCTION_BOUNDARY_POST_PROCESSORS_3D_HH

#include "wallFunctionBoundaryPostProcessors3D.h"
#include "core/finiteDifference3D.h"
#include "core/blockLattice3D.h"
#include "dynamics/firstOrderLbHelpers.h"
#include "core/util.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template <typename T, typename S>
Musker<T,S>::Musker(T nu, T y, T rho) : AnalyticalF1D<T,S>(1), _nu(nu), _y(y),_rho(rho)
{
  this->getName() = "Musker";
}

template <typename T, typename S>
bool Musker<T,S>::operator()(T output[], const S tau_w[])
{
  T y_plus = _y*sqrt(tau_w[0]/_rho)/_nu;

  T a = 5.424;
  T b = 0.119760479041916168;
  T c = 0.488023952095808383;
  T d = 0.434;
  T e = 3.50727901936264842;

  output[0] = sqrt(tau_w[0]/_rho)*(a*atan(b*y_plus - c) +
                                   d*log(pow(y_plus+10.6, 9.6)/pow(pow(y_plus, 2) - 8.15*y_plus + 86, 2)) - e);

  // Condition for the sub-viscous layer : TODO MARC H
  if (output[0] < 0) {
    output[0] = y_plus * sqrt(tau_w[0]/_rho);
  }

  return true;
}

template <typename T, typename S>
PowerLawProfile<T,S>::PowerLawProfile(T nu, T u2, T y2, T y1, T rho) : AnalyticalF1D<T,S>(2), _nu(nu), _u2(u2), _y2(y2), _y1(y1), _rho(rho)
{
  this->getName() = "PowerLawProfile";
}

template <typename T, typename S>
bool PowerLawProfile<T,S>::operator()(T output[], const S input[])
{
  T tau_w = 0.0246384 * _rho * pow(_nu, 0.25) * pow(_u2, 1.75) / pow(_y2, 0.25);
  T u_tau = sqrt(tau_w/_rho);
  T y_plus = _y1 * u_tau / _nu;

  if (y_plus > 30.0) {
    output[0] = tau_w;
    output[1] = u_tau * 8.3 * pow(y_plus, 1./7.);
  } else if (y_plus < 30.0  && y_plus > 5.) {
    output[0] = tau_w;
    output[1] = u_tau * (-2.81742 + 4.85723 * log(y_plus));
  } else {
    output[0] = 2.*_u2 * _rho * _nu / _y2;
    if (output[0] < 0.) {
      output[0]*=-1.;
    }
    u_tau = sqrt(output[0]/_rho);
    y_plus = _y1 * u_tau / _nu;
    output[1] = u_tau * y_plus;
  }
  return true;
}


////////  WallFunctionBoundaryProcessor3D ////////////////////////////////
template<typename T, typename DESCRIPTOR>
WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::WallFunctionBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                                                                            BlockGeometryStructure3D<T>& blockGeometryStructure,
                                                                            std::vector<int> discreteNormal, std::vector<int> missingIndices,
                                                                            UnitConverter<T, DESCRIPTOR> const& converter, wallFunctionParam<T> const& wallFunctionParam,
                                                                            IndicatorF3D<T>* geoIndicator)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_),
    _blockGeometryStructure(blockGeometryStructure),
    _discreteNormal(discreteNormal), _missingIndices(missingIndices),
    _converter(converter), _wallFunctionParam(wallFunctionParam)
{
  // Define Direction and orientatation
  discreteNormalX = _discreteNormal[0];
  discreteNormalY = _discreteNormal[1];
  discreteNormalZ = _discreteNormal[2];

  int Type_BC = discreteNormalX*discreteNormalX + discreteNormalY*discreteNormalY + discreteNormalZ*discreteNormalZ;
  T normal_norm = sqrt(Type_BC); // l2 norm : magnitude of the vector
  if (Type_BC == 1) { // Straight plane
    if (discreteNormalX != 0) {
      orientation = discreteNormalX;
      direction = 0;
    } else if (discreteNormalY != 0) {
      orientation = discreteNormalY;
      direction = 1;
    } else if (discreteNormalZ != 0) {
      orientation = discreteNormalZ;
      direction = 2;
    }
  } else if (Type_BC == 2) { // Edge
    orientation = 0;
    direction = 0;
  } else if (Type_BC == 3) { // Corner
    orientation = 0;
    direction = 0;
  }

  if (geoIndicator == NULL) {
    /// === STEP 1 :  Define distance for boundary and neighbor nodes
    y_1 = 0.5; // Default half-way Bounce-Back distance - [m] DESCRIPTOR units
    if (_wallFunctionParam.latticeWalldistance > 0.) {
      y_1 = _wallFunctionParam.latticeWalldistance; // [m] DESCRIPTOR units
    }
    y_1 *= normal_norm; // DESCRIPTOR units

    // Define the unit normal vector
    unit_normal[0] = discreteNormalX / normal_norm;
    unit_normal[1] = discreteNormalY / normal_norm;
    unit_normal[2] = discreteNormalZ / normal_norm;
  } else {
    calculateWallDistances(geoIndicator);
    if (y_1 > 2. || std::isnan(unit_normal[0]) || std::isnan(unit_normal[1]) || std::isnan(unit_normal[2])) {
      y_1 = 0.5; // Default half-way Bounce-Back distance - [m] DESCRIPTOR units
      if (_wallFunctionParam.latticeWalldistance > 0.) {
        y_1 = _wallFunctionParam.latticeWalldistance; // [m] DESCRIPTOR units
      }
      y_1 *= normal_norm; // DESCRIPTOR units

      // Define the unit normal vector
      unit_normal[0] = discreteNormalX / normal_norm;
      unit_normal[1] = discreteNormalY / normal_norm;
      unit_normal[2] = discreteNormalZ / normal_norm;
    }
  }
  y_2 = y_1 + normal_norm; // DESCRIPTOR units

  // Computation of Rho - Zou an He - Speed up
  getIndices(direction, 0, onWallIndices);
  getIndices(direction, orientation, normalIndices);

  // Vector needed for Pi Computation - FNeq Bounce-Back - Speep up
  // Malaspinas condition c*n < zero
  // Definition of directions pointing towards the fluid
  getIndices(direction, -orientation, normalInwardsIndices);
}

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    int iX;
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          ComputeWallFunction(blockLattice, iX,iY,iZ);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::getIndices(int index, int value , std::vector<int>& indices)
{
  for (int iVel=0; iVel<DESCRIPTOR::q; ++iVel) {
    if (descriptors::c<DESCRIPTOR>(iVel,index)==value) {
      indices.push_back(iVel);
    }
  }
}

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::calculateWallDistances(IndicatorF3D<T>* indicator)
{
  T physR[3];
  int iX = x0;
  int iY = y1;
  int iZ = z1;
  T scaling = _converter.getConversionFactorLength() * 0.1;
  _blockGeometryStructure.getPhysR(physR,iX, iY, iZ);
  Vector<T,3> origin(physR[0],physR[1],physR[2]);
  Vector<T,3> normal(0.,0.,0.);
  T distance = 0.;
  Vector<T,3> direction(0.,0.,0.);
  int smallestDistance_i = 0;
  T smallestDistance = 0.;
  bool firstHit = true;
  origin[0] = physR[0];
  origin[1] = physR[1];
  origin[2] = physR[2];
  int discreteDirection[6][3];
  for (int i=0; i < 6; i++) {
    for (int j=0; j < 3; j++) {
      discreteDirection[i][j] = 0;
    }
  }
  discreteDirection[0][0] = 1;
  discreteDirection[1][0] = -1;
  discreteDirection[2][1] = 1;
  discreteDirection[3][1] = -1;
  discreteDirection[4][2] = 1;
  discreteDirection[5][2] = -1;
  for (int i=0; i < 6; i++) {
    direction[0] = discreteDirection[i][0] * scaling;
    direction[1] = discreteDirection[i][1] * scaling;
    direction[2] = discreteDirection[i][2] * scaling;
    if (indicator->distance(distance, origin, direction)) {
      if (firstHit) {
        smallestDistance = distance;
        smallestDistance_i = i;
        firstHit = false;
      } else {
        if (distance < smallestDistance) {
          smallestDistance = distance;
          smallestDistance_i = i;
        }
      }
    }
  }

  direction[0] = discreteDirection[smallestDistance_i][0] * scaling;
  direction[1] = discreteDirection[smallestDistance_i][1] * scaling;
  direction[2] = discreteDirection[smallestDistance_i][2] * scaling;

  Vector<T,3> direction2(direction);
  Vector<T,3> direction3(direction);
  if (smallestDistance_i == 0 || smallestDistance_i == 1 ) {
    direction2[1] = direction2[0] * scaling;
    direction3[2] = direction3[0] * scaling;
  } else if (smallestDistance_i == 2 || smallestDistance_i == 3 ) {
    direction2[0] = direction2[1] * scaling;
    direction3[2] = direction3[1] * scaling;
  } else {
    direction2[0] = direction2[2] * scaling;
    direction3[1] = direction3[2] * scaling;
  }

  Vector<S,3> directionN = direction*(1/const_cast<Vector<S,3>&> (direction).norm());
  Vector<S,3> POS(origin + smallestDistance*directionN); //Point on Surface

  indicator->distance(distance, origin, direction2);
  Vector<S,3> direction2N = direction2*(1/const_cast<Vector<S,3>&> (direction2).norm());
  Vector<S,3> POS2(origin + distance*direction2N); //Point on Surface

  indicator->distance(distance, origin, direction3);

  Vector<S,3> direction3N = direction3*(1/const_cast<Vector<S,3>&> (direction3).norm());
  Vector<S,3> POS3(origin + distance*direction3N); //Point on Surface

  Vector<S,3> vec1 (POS - POS2);
  Vector<S,3> vec2 (POS - POS3);

  normal[0] = -(vec1[1]*vec2[2] - vec1[2]*vec2[1]);
  normal[1] = -(vec1[2]*vec2[0] - vec1[0]*vec2[2]);
  normal[2] = -(vec1[0]*vec2[1] - vec1[1]*vec2[0]);

  T normalMagnitude = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  normal[0] /= normalMagnitude;
  normal[1] /= normalMagnitude;
  normal[2] /= normalMagnitude;

  unit_normal[0] = normal[0];
  unit_normal[1] = normal[0];
  unit_normal[2] = normal[0];

  direction[0] = normal[0] * scaling;
  direction[1] = normal[1] * scaling;
  direction[2] = normal[2] * scaling;

  indicator->distance(distance, origin, direction);
  y_1 = distance / _converter.getConversionFactorLength();
}

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::VelGradFromSecondOrderFD(bool NormalGradient, T Vel_BC[DESCRIPTOR::d], T Vel_1[DESCRIPTOR::d], T Vel_2[DESCRIPTOR::d], T VelGrad[DESCRIPTOR::d])
{
  if (NormalGradient) {
    if (orientation == 1) {
      for (int Dim=0; Dim<DESCRIPTOR::d; Dim++) {
        VelGrad[Dim] = fd::BSGradient(Vel_BC[Dim], Vel_1[Dim], Vel_2[Dim]); // Backward second-order velocity gradient
      }
    } else { // orientation == -1
      for (int Dim=0; Dim<DESCRIPTOR::d; Dim++) {
        VelGrad[Dim] = fd::FSGradient(Vel_BC[Dim], Vel_1[Dim], Vel_2[Dim]); // Forward second-order velocity gradient
      }
    }
  } else {
    if (orientation == 1) {
      for (int Dim=0; Dim<DESCRIPTOR::d; Dim++) {
        VelGrad[Dim] = fd::centralGradient(Vel_1[Dim], Vel_2[Dim]); // central second-order velocity gradient: centralGradient( y(x+1), y(x-1) )
      }
    } else { // orientation == -1
      for (int Dim=0; Dim<DESCRIPTOR::d; Dim++) {
        VelGrad[Dim] = fd::centralGradient(Vel_2[Dim], Vel_1[Dim]); // central second-order velocity gradient: centralGradient( y(x+1), y(x-1) )
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::computeVelocityGradientTensor(T u_bc[DESCRIPTOR::d], T u_x1[DESCRIPTOR::d], T u_x2[DESCRIPTOR::d], T u_y1[DESCRIPTOR::d], T u_y2[DESCRIPTOR::d], T u_z1[DESCRIPTOR::d], T u_z2[DESCRIPTOR::d],  T VelGrad[DESCRIPTOR::d][DESCRIPTOR::d])
{
  T dx_u[DESCRIPTOR::d], dy_u[DESCRIPTOR::d], dz_u[DESCRIPTOR::d];
  VelGradFromSecondOrderFD(direction == 0, u_bc, u_x1, u_x2, dx_u);
  VelGradFromSecondOrderFD(direction == 1, u_bc, u_y1, u_y2, dy_u);
  VelGradFromSecondOrderFD(direction == 2, u_bc, u_z1, u_z2, dz_u);

  // du/dx  du/dy  du/dz
  // dv/dx  dv/dy  dv/dz
  // dw/dx  dw/dy  dw/dz
  for (int Dim = 0; Dim < DESCRIPTOR::d; Dim++) {
    VelGrad[Dim][0] = dx_u[Dim];
    VelGrad[Dim][1] = dy_u[Dim];
    VelGrad[Dim][2] = dz_u[Dim];
  }
}

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::computeNeighborsU(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x, int y, int z,
                                                                             T u_x1[DESCRIPTOR::d], T u_x2[DESCRIPTOR::d], T u_y1[DESCRIPTOR::d], T u_y2[DESCRIPTOR::d], T u_z1[DESCRIPTOR::d], T u_z2[DESCRIPTOR::d])
{
  using namespace olb::util::tensorIndices3D;
  // Computation of local external velocity around lattice node (x,y,z)
  if (direction == 0) {
    blockLattice.get(x -   discreteNormalX, y, z).computeU(u_x1);
    blockLattice.get(x - 2*discreteNormalX, y, z).computeU(u_x2);
  } else {
    ComputeUWall(blockLattice, x + discreteNormalY + discreteNormalZ, y, z, u_x1);
    ComputeUWall(blockLattice, x - discreteNormalY - discreteNormalZ, y, z, u_x2);
  }
  if (direction == 1) {
    blockLattice.get(x, y -   discreteNormalY, z).computeU(u_y1);
    blockLattice.get(x, y - 2*discreteNormalY, z).computeU(u_y2);
  } else {
    ComputeUWall(blockLattice, x, y + discreteNormalX + discreteNormalZ, z, u_y1);
    ComputeUWall(blockLattice, x, y - discreteNormalX - discreteNormalZ, z, u_y2);
  }
  if (direction == 2) {
    blockLattice.get(x, y, z -   discreteNormalZ).computeU(u_z1);
    blockLattice.get(x, y, z - 2*discreteNormalZ).computeU(u_z2);
  } else {
    ComputeUWall(blockLattice, x, y, z + discreteNormalX + discreteNormalY, u_z1);
    ComputeUWall(blockLattice, x, y, z - discreteNormalX - discreteNormalY, u_z2);
  }
}

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::computeVelocityGradient(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x, int y, int z, T u_bc[DESCRIPTOR::d],  T VelGrad[DESCRIPTOR::d][DESCRIPTOR::d])
{
  // Computation of neighbor velocity around lattice node (x,y,z)
  T u_x1[3], u_x2[3], u_y1[3], u_y2[3], u_z1[3], u_z2[3];
  computeNeighborsU(blockLattice, x, y, z, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2);
  // Computation of Velocity Gradient Tensor
  computeVelocityGradientTensor(u_bc, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2, VelGrad);
}

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::computeNeighborsRho(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x, int y, int z,
                                                                               T u_x1[DESCRIPTOR::d], T u_x2[DESCRIPTOR::d], T u_y1[DESCRIPTOR::d], T u_y2[DESCRIPTOR::d], T u_z1[DESCRIPTOR::d], T u_z2[DESCRIPTOR::d],
                                                                               T& rho_x1, T& rho_x2, T& rho_y1, T& rho_y2, T& rho_z1, T& rho_z2)
{
  using namespace olb::util::tensorIndices3D;
  // Computation of local external velocity around lattice node (x,y,z)
  if (direction == 0) {
    rho_x1 = blockLattice.get(x -   discreteNormalX, y, z).computeRho();
    rho_x2 = blockLattice.get(x - 2*discreteNormalX, y, z).computeRho();
  } else {
    ComputeRhoWall(blockLattice, blockLattice.get(x + discreteNormalY + discreteNormalZ, y, z), x + discreteNormalY + discreteNormalZ, y, z, u_x1, rho_x1);
    ComputeRhoWall(blockLattice, blockLattice.get(x - discreteNormalY - discreteNormalZ, y, z), x - discreteNormalY - discreteNormalZ, y, z, u_x2, rho_x2);
  }
  if (direction == 1) {
    rho_y1 = blockLattice.get(x, y -   discreteNormalY, z).computeRho();
    rho_y2 = blockLattice.get(x, y - 2*discreteNormalY, z).computeRho();
  } else {
    ComputeRhoWall(blockLattice, blockLattice.get(x, y - discreteNormalX - discreteNormalZ, z), x, y - discreteNormalX - discreteNormalZ, z, u_y1, rho_y1);
    ComputeRhoWall(blockLattice, blockLattice.get(x, y - discreteNormalX - discreteNormalZ, z), x, y - discreteNormalX - discreteNormalZ, z, u_y2, rho_y2);
  }
  if (direction == 2) {
    rho_z1 = blockLattice.get(x, y, z -   discreteNormalZ).computeRho();
    rho_z2 = blockLattice.get(x, y, z - 2*discreteNormalZ).computeRho();
  } else {
    ComputeRhoWall(blockLattice, blockLattice.get(x, y, z + discreteNormalX + discreteNormalY), x, y, z + discreteNormalX + discreteNormalY, u_z1, rho_z1);
    ComputeRhoWall(blockLattice, blockLattice.get(x, y, z - discreteNormalX - discreteNormalY), x, y, z - discreteNormalX - discreteNormalY, u_z2, rho_z2);
  }
}

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::computeNeighborsRhoU(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x, int y, int z,
                                                                                T u_x1[DESCRIPTOR::d], T u_x2[DESCRIPTOR::d], T u_y1[DESCRIPTOR::d], T u_y2[DESCRIPTOR::d], T u_z1[DESCRIPTOR::d], T u_z2[DESCRIPTOR::d],
                                                                                T& rho_x1, T& rho_x2, T& rho_y1, T& rho_y2, T& rho_z1, T& rho_z2)
{
  // Computation of neighbor velocity around lattice node (x,y,z)
  computeNeighborsU(blockLattice, x, y, z, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2);
  // Computation of neighbor density around lattice node (x,y,z)
  computeNeighborsRho(blockLattice, x, y, z, u_x1, u_x2, u_y1, u_y2, u_z1, u_z2, rho_x1, rho_x2, rho_y1, rho_y2, rho_z1, rho_z2);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::ComputeUWall(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x, int y, int z, T u[DESCRIPTOR::d])
{
  /// === Computation of velocity with Musker Profile - Malaspinas and Sagaut (2014) === ///
  using namespace olb::util::tensorIndices3D;
  Cell<T,DESCRIPTOR>& cell = blockLattice.get(x,y,z);

  /// === STEP 2 : compute u_2 and rho_2 (neighbor node) -> [m/s] - [m] DESCRIPTOR units
  T u_2[DESCRIPTOR::d]; // Velocity to the neighbord lattice in the normal inward direction
  T rho_2; // Density in the neighbor lattice in the inwards direction
  blockLattice.get(x - discreteNormalX, y - discreteNormalY, z - discreteNormalZ).computeRhoU(rho_2,u_2);

  /// === STEP 3 : get local basis vector and neighbor wall parallel velocity
  T u_2_dot_unit_normal = u_2[0] * unit_normal[0] +
                          u_2[1] * unit_normal[1] +
                          u_2[2] * unit_normal[2];  //scalar Product of the u2 and uniy normal
  T u_2_parallel[3];
  u_2_parallel[0] = u_2[0] - (u_2_dot_unit_normal * unit_normal[0]);
  u_2_parallel[1] = u_2[1] - (u_2_dot_unit_normal * unit_normal[1]);
  u_2_parallel[2] = u_2[2] - (u_2_dot_unit_normal * unit_normal[2]);

  T u_2_parallel_norm = sqrt(u_2_parallel[0] * u_2_parallel[0] +
                             u_2_parallel[1] * u_2_parallel[1] +
                             u_2_parallel[2] * u_2_parallel[2]); // l2 norm : magnitude of the vector

  T e_x_loc[3] = {0., 0., 0.}; // Streamwise direction
  if (u_2_parallel_norm != 0) {
    T inv_u_2_parallel_norm = (1. /u_2_parallel_norm);
    e_x_loc[0] = u_2_parallel[0] * inv_u_2_parallel_norm;
    e_x_loc[1] = u_2_parallel[1] * inv_u_2_parallel_norm;
    e_x_loc[2] = u_2_parallel[2] * inv_u_2_parallel_norm;
  }
  T u2 = e_x_loc[0] * u_2[0] +
         e_x_loc[1] * u_2[1] +
         e_x_loc[2] * u_2[2];  //scalar Product of the basis vector e_x_loc and u_2

  // STEP 5 :  u1 can only be 0 if u2==0 for musker profile
  //T* tau_w = new T [1];
  T tau_w[1] = {0.};
  T u1[1];
  if (u2 != 0) {

    T rho_phy = _converter.getPhysDensity(); // [kg/mÂ³] SI units
    T rho_lat = _converter.getLatticeDensity(rho_phy); // [kg/mÂ³] DESCRIPTOR units
    T Mol_Lat_Nu = _converter.getLatticeViscosity(); // [mÂ²/s] DESCRIPTOR units

    if (_wallFunctionParam.wallProfile == 0) {
      // Musker profile for boundary node using the molecular viscosity and the characteristic density in lattice units
      Musker<T,T> musker_u2_Lat_MolVis(Mol_Lat_Nu, y_2, rho_lat);
      util::Newton1D<T> newton(musker_u2_Lat_MolVis,u2);
      T tau_w_guess = *(cell.template getFieldPointer<descriptors::TAU_W>());
      tau_w[0] = newton.solve(tau_w_guess,false); // Computation of local wall shear stress

      if (std::isnan(tau_w[0])) {
        //OstreamManager clout(std::cout, "First NAN");
        //clout << "Position = [" << x << ", " << y <<  "," << z << "]" << std::endl;
        //clout << "Tau guess  = " << tau_w_guess << std::endl;
        tau_w[0] = newton.solve(0.,false);
      }

      if (std::isnan(tau_w[0])) {
        OstreamManager clout(std::cout, "Second NAN");
        clout << "Position = [" << x << ", " << y <<  "," << z << "]" << std::endl;
        clout << "Tau guess  = " << tau_w_guess << std::endl;
        tau_w[0] = 0.;
      }

      // compute u1 - boundary node
      Musker<T,T> musker_u1_Lat_MolVis(Mol_Lat_Nu, y_1, rho_lat);
      musker_u1_Lat_MolVis(u1,tau_w);

    } else {
      // Werner and Wengle profile for boundary node using the molecular viscosity and the characteristic density in lattice units
      PowerLawProfile<T,T> PLP (Mol_Lat_Nu, u2, y_2, y_1, rho_lat);
      T input [1];
      T output[2];
      PLP(output,input);
      tau_w[0] = output[0];
      u1[0] = output[1];
    }
  } else {
    u1[0] = 0.;
  }
  // save tau_w for next step
  cell.template defineField<descriptors::TAU_W>(&(tau_w[0]));
  // STEP 6 : compute velocity vector at the boundary
  u[0] = e_x_loc[0] * u1[0];
  u[1] = e_x_loc[1] * u1[0];
  u[2] = e_x_loc[2] * u1[0];

  if (_wallFunctionParam.bodyForce) {
    for (int iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
      u[iDim] -= cell.template getFieldPointer<descriptors::FORCE>()[iDim]/2.;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::computeVanDriestTauEff(T y_bc, T tau_w, T u_bc, T u_1, T u_2, T& tau_eff)
{

  T rho_phy = _converter.getPhysDensity(); // [kg/mÂ³] SI units
  T rho = _converter.getLatticeDensity(rho_phy); // [kg/mÂ³] DESCRIPTOR units
  T nu_mol = _converter.getLatticeViscosity(); // [mÂ²/s] DESCRIPTOR units

  T y_plus = y_bc*sqrt(tau_w/rho)/nu_mol;
  T uxdz_abs = std::abs(fd::boundaryGradient(u_bc, u_1, u_2));

  T vanDriest = 1 - std::exp(-y_plus / 26.);

  T nu_turb = pow(_wallFunctionParam.vonKarman * y_bc * vanDriest, 2.) * uxdz_abs;

  T nu_eff = nu_turb + nu_mol;

  tau_eff = (3. * nu_eff) + 0.5; // Rectangular integration rule - Collision
}
////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::ComputeTauEff(BlockLattice3D<T,DESCRIPTOR>& blockLattice, Cell<T,DESCRIPTOR>& cell, int x, int y, int z, T u_bc[DESCRIPTOR::d])
{

    T u_z[3];
    T u_z2[3];
    blockLattice.get(x -   discreteNormalX, y -   discreteNormalY, z -   discreteNormalZ).computeU(u_z);
    blockLattice.get(x - 2*discreteNormalX, y - 2*discreteNormalY, z - 2*discreteNormalZ).computeU(u_z2);

    T y_bc = 0.5; // Default half-way Bounce-Back distance - [m] DESCRIPTOR units
    if (_wallFunctionParam.latticeWalldistance > 0.) {
      y_bc = _wallFunctionParam.latticeWalldistance; // [m] DESCRIPTOR units
    }
    T normal_norm = sqrt(discreteNormalX * discreteNormalX +
                         discreteNormalY * discreteNormalY +
                         discreteNormalZ * discreteNormalZ); // l2 norm : magnitude of the vector
    y_bc *= normal_norm;
    T tau_w = cell.template getFieldPointer<descriptors::TAU_W>()[0];

    T tau_eff;

  if (_wallFunctionParam.curved == true || orientation == 0) {
    T inv_normal_norm = 1. / normal_norm;
    T normal[3];
    normal[0] = discreteNormalX * (inv_normal_norm);
    normal[1] = discreteNormalY * (inv_normal_norm);
    normal[2] = discreteNormalZ * (inv_normal_norm);

    T u_2_parallel[3];
    T sP_normal_u2 = normal[0] * u_z[0] +
                     normal[1] * u_z[1] +
                     normal[2] * u_z[2];  //scalar Product of the normal and u2

    u_2_parallel[0] = u_z[0] - sP_normal_u2 * normal[0];
    u_2_parallel[1] = u_z[1] - sP_normal_u2 * normal[1];
    u_2_parallel[2] = u_z[2] - sP_normal_u2 * normal[2];

    T u_2_parallel_norm = sqrt(u_2_parallel[0] * u_2_parallel[0] +
                               u_2_parallel[1] * u_2_parallel[1] +
                               u_2_parallel[2] * u_2_parallel[2]); // l2 norm : magnitude of the vector

    T e_x_loc[3];
    T inv_u_2_parallel_norm = (1. /u_2_parallel_norm );

    if (u_2_parallel_norm != 0) {
      e_x_loc[0] = u_2_parallel[0] * inv_u_2_parallel_norm;
      e_x_loc[1] = u_2_parallel[1] * inv_u_2_parallel_norm;
      e_x_loc[2] = u_2_parallel[2] * inv_u_2_parallel_norm;
    } else {
      e_x_loc[0] = T(-discreteNormalX);
      e_x_loc[1] = T(-discreteNormalY);
      e_x_loc[2] = T(-discreteNormalZ);
    }

    T u1;
    T u2;
    T u3;
    u1 = e_x_loc[0] * u_bc[0] +
         e_x_loc[1] * u_bc[1] +
         e_x_loc[2] * u_bc[2];
    u2 = e_x_loc[0] * u_z[0] +
         e_x_loc[1] * u_z[1] +
         e_x_loc[2] * u_z[2];
    u3 = e_x_loc[0] * u_z2[0] +
         e_x_loc[1] * u_z2[1] +
         e_x_loc[2] * u_z2[2];

    computeVanDriestTauEff(y_bc, tau_w, u1, u2, u3, tau_eff);
  } else {
    computeVanDriestTauEff(y_bc, tau_w, u_bc[direction], u_z[direction], u_z2[direction], tau_eff);
  }
    cell.template defineField<descriptors::TAU_EFF>(&(tau_eff));

}

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::ComputeRhoWall(BlockLattice3D<T,DESCRIPTOR>& blockLattice, Cell<T,DESCRIPTOR>& cell, int x, int y, int z, T u_bc[DESCRIPTOR::d], T& rho_bc)
{

/// === Computation of density - Finite Difference Scheme === ///
  if (_wallFunctionParam.rhoMethod == 0) {
    /// === Computation of density - Zou and He (1997) === ///
    // The code have been copied from VelocityBM -> compute rho method. This implementation it is only valid for straight boundaries
    T u_bc_Perpendicular = 0.;
    if (_wallFunctionParam.bodyForce) {
      T u_bc_tmp[DESCRIPTOR::d];
      for (int iDim = 0; iDim<DESCRIPTOR::d; ++iDim){
        u_bc_tmp[iDim] = u_bc[iDim] - cell.template getFieldPointer<descriptors::FORCE>()[iDim] / 2.;
      }
      for (int iDim = 0; iDim<DESCRIPTOR::d; ++iDim){
        u_bc_Perpendicular += u_bc_tmp[iDim]*unit_normal[iDim];
      }
    } else {
      for (int iDim = 0; iDim<DESCRIPTOR::d; ++iDim){
        u_bc_Perpendicular += u_bc[iDim]*unit_normal[iDim];
      }
    }

    T rhoOnWall = T();
    for (unsigned fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
      rhoOnWall += cell[onWallIndices[fIndex]];
    }

    T rhoNormal = T();
    for (unsigned fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
      rhoNormal += cell[normalIndices[fIndex]];
    }

    rho_bc = ((T)2*rhoNormal+rhoOnWall+(T)1) /
          ((T)1+u_bc_Perpendicular);

  } else if (_wallFunctionParam.rhoMethod == 1) { // Extrapolation Fluid

    rho_bc = blockLattice.get(x - discreteNormalX, y - discreteNormalY, z - discreteNormalZ).computeRho();

  } else if (_wallFunctionParam.rhoMethod == 2) { // constant rho
    rho_bc = 1.;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::computeRFneqfromFneq(T fneq_bc[DESCRIPTOR::q])
{
  T sum_cell_fneq = 0.;
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    sum_cell_fneq += fneq_bc[iPop];
  }
  T pi_bc[util::TensorVal< DESCRIPTOR >::n];
  int iPi = 0;
  for (int iAlpha=0; iAlpha < DESCRIPTOR::d; ++iAlpha) {
    for (int iBeta=iAlpha; iBeta < DESCRIPTOR::d; ++iBeta) {
      pi_bc[iPi] = T();
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        pi_bc[iPi] += descriptors::c<DESCRIPTOR>(iPop)[iAlpha]*descriptors::c<DESCRIPTOR>(iPop)[iBeta]*fneq_bc[iPop];
      }
      if (iAlpha==iBeta) {
        pi_bc[iPi] -= (1./descriptors::invCs2<T,DESCRIPTOR>())*sum_cell_fneq;
      }
      ++iPi;
    }
  }

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    fneq_bc[iPop] = firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi_bc);
  }
}

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::computeFneqRNEBB(Cell<T,DESCRIPTOR>& cell, T u_bc[DESCRIPTOR::d], T rho_bc, T fneq_bc[DESCRIPTOR::q])
{
  Dynamics<T,DESCRIPTOR>* dynamics = cell.getDynamics();
  T uSqr_bc = util::normSqr<T,DESCRIPTOR::d>(u_bc);
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    fneq_bc[iPop] = cell[iPop] - dynamics -> computeEquilibrium(iPop,rho_bc,u_bc,uSqr_bc);
  }
  for (unsigned fIndex=0; fIndex<normalInwardsIndices.size(); ++fIndex) {
    fneq_bc[normalInwardsIndices[fIndex]] = fneq_bc[DESCRIPTOR::opposite[normalInwardsIndices[fIndex]]];
  }

  computeRFneqfromFneq(fneq_bc);
}

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::computeFneqENeq(BlockLattice3D<T,DESCRIPTOR>& blockLattice, Cell<T,DESCRIPTOR>& cell, int x, int y, int z, T u_bc[DESCRIPTOR::d], T rho_bc, T fneq_bc[DESCRIPTOR::q])
{
  Cell<T,DESCRIPTOR>& cell_fluid =  blockLattice.get(x - discreteNormalX, y - discreteNormalY, z - discreteNormalZ);
  Dynamics<T,DESCRIPTOR>* dynamics_fluid = cell_fluid.getDynamics();
  T rho_fluid, u_fluid[DESCRIPTOR::d];
  cell_fluid.computeRhoU(rho_fluid,u_fluid);
  T uSqr_fluid = util::normSqr<T,DESCRIPTOR::d>(u_fluid);
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    fneq_bc[iPop] = cell_fluid[iPop] - dynamics_fluid -> computeEquilibrium(iPop,rho_fluid,u_fluid,uSqr_fluid);
  }
}

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::computeFneqRSOFD(BlockLattice3D<T,DESCRIPTOR>& blockLattice, Cell<T,DESCRIPTOR>& cell, int x, int y, int z, T u_bc[DESCRIPTOR::d], T rho_bc, T fneq_bc[DESCRIPTOR::q])
{
  T pi_bc[util::TensorVal< DESCRIPTOR >::n];
  T Velocity_Grad[DESCRIPTOR::d][DESCRIPTOR::d];
  computeVelocityGradient(blockLattice, x, y, z, u_bc, Velocity_Grad);
  // Creation of strain Rate
  T tau_eff = cell.template getFieldPointer<descriptors::TAU_EFF>()[0];
  T Factor = -1.*tau_eff*rho_bc/descriptors::invCs2<T,DESCRIPTOR>();
  int iPi = 0;
  for (int Alpha=0; Alpha<DESCRIPTOR::d; ++Alpha) {
    for (int Beta=Alpha; Beta<DESCRIPTOR::d; ++Beta) {
      pi_bc[iPi] = Factor*(Velocity_Grad[Alpha][Beta] + Velocity_Grad[Beta][Alpha]);
      ++iPi;
    }
  }

  if (_wallFunctionParam.bodyForce) {
    // Force tensor
    //Creation of body force tensor (rho/2.)*(G_alpha*U_beta + U_alpha*G_Beta)
    T F[DESCRIPTOR::d];
    for (int iDim = 0; iDim < DESCRIPTOR::d; ++iDim) {
      F[iDim] = cell.template getFieldPointer<descriptors::FORCE>()[iDim];
    }

    iPi = 0;
    for (int Alpha=0; Alpha<DESCRIPTOR::d; ++Alpha) {
      for (int Beta=Alpha; Beta<DESCRIPTOR::d; ++Beta) {
        pi_bc[iPi] += -1.*(rho_bc/2.)*(F[Alpha]*u_bc[Beta] + u_bc[Alpha]*F[Beta]);
        ++iPi;
      }
    }
  }

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    fneq_bc[iPop] = firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi_bc);
  }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::ComputeFneqWall(BlockLattice3D<T,DESCRIPTOR>& blockLattice, Cell<T,DESCRIPTOR>& cell, int x, int y, int z, T u_bc[DESCRIPTOR::d], T rho_bc, T fneq_bc[DESCRIPTOR::q])
{

  //regularized NEBB (Latt)
  if (_wallFunctionParam.fneqMethod == 0) {
    computeFneqRNEBB(cell, u_bc, rho_bc, fneq_bc);

  //extrapolation NEQ (Guo Zhaoli)
  } else if (_wallFunctionParam.fneqMethod == 1) {
    computeFneqENeq(blockLattice, cell, x, y, z, u_bc, rho_bc, fneq_bc);

  //Regularized second order finite Differnce
  } else if (_wallFunctionParam.fneqMethod == 2) {
    computeFneqRSOFD(blockLattice, cell, x, y, z, u_bc, rho_bc, fneq_bc);

  //equilibrium scheme
  } else if (_wallFunctionParam.fneqMethod == 3) {
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      fneq_bc[iPop] = 0.;
    }
  }

}

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::ComputeWallFunction(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x, int y, int z)
{
  Cell<T,DESCRIPTOR>& cell_bc = blockLattice.get(x,y,z);
  T rho_bc = 0.;
  T u_bc[DESCRIPTOR::d];
  T fneq_bc[DESCRIPTOR::q];

  // Computation of boundary velocity from the wall function
  ComputeUWall(blockLattice, x, y, z, u_bc);

  if (_wallFunctionParam.useVanDriest) {
    // Computation of effective relaxation time from the vanDriest damping function
    ComputeTauEff(blockLattice, cell_bc, x, y, z, u_bc);
  }

  // Computation of density at the boundary node
  ComputeRhoWall(blockLattice, cell_bc, x, y, z, u_bc, rho_bc);

  // Computation of the second-order moment of non-equilibrium distribution function
  ComputeFneqWall(blockLattice, cell_bc, x, y, z, u_bc, rho_bc, fneq_bc);

  // Computation of the particle distribution functions according to the regularized formula
  Dynamics<T,DESCRIPTOR>* dynamics_bc = cell_bc.getDynamics();
  T uSqr_bc = util::normSqr<T,DESCRIPTOR::d>(u_bc);

  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    cell_bc[iPop] = dynamics_bc -> computeEquilibrium(iPop,rho_bc,u_bc,uSqr_bc) + fneq_bc[iPop];
    if (std::isnan(cell_bc[iPop])) {
      OstreamManager clout(std::cout, "Slip Musker Profile");
      clout << "Musker Profile Computation" << std::endl;
      clout << "Position = [" << x << ", " << y <<  "," << z << "]" << std::endl;
      clout << "Normal outwards = [" << discreteNormalX << ", " << discreteNormalY <<  "," << discreteNormalZ << "]" << std::endl;
      clout << "Velocity at boundary u_bc = [" << u_bc[0] << ", " << u_bc[1] <<  "," << u_bc[2] << "]" << std::endl;
      clout << "Density at boundary rho_bc = " << rho_bc << std::endl;
      exit(1);
    }
  }

}

template<typename T, typename DESCRIPTOR>
void WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>::process(BlockLattice3D<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

////////  WallFunctionBoundaryProcessorGenerator3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
WallFunctionBoundaryProcessorGenerator3D<T,DESCRIPTOR>::WallFunctionBoundaryProcessorGenerator3D(int x0, int x1, int y0, int y1, int z0, int z1, BlockGeometryStructure3D<T>& blockGeometryStructure,
                                                                                              std::vector<int> discreteNormal, std::vector<int> missingIndices,
                                                                                              UnitConverter<T, DESCRIPTOR> const& converter, wallFunctionParam<T> const& wallFunctionParam,
                                                                                              IndicatorF3D<T>* geoIndicator)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0, x1, y0, y1, z0, z1),
    _blockGeometryStructure(blockGeometryStructure),
    _discreteNormal(discreteNormal), _missingIndices(missingIndices),
    _converter(converter), _wallFunctionParam(wallFunctionParam), _geoIndicator(geoIndicator)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
WallFunctionBoundaryProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new WallFunctionBoundaryProcessor3D<T,DESCRIPTOR>( this->x0, this->x1, this->y0, this->y1, this->z0, this->z1,
                                                         _blockGeometryStructure, _discreteNormal, _missingIndices,
                                                         _converter, _wallFunctionParam, _geoIndicator);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
WallFunctionBoundaryProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new WallFunctionBoundaryProcessorGenerator3D<T,DESCRIPTOR> (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1,
                                                                  _blockGeometryStructure, _discreteNormal, _missingIndices,
                                                                  _converter, _wallFunctionParam, _geoIndicator);
}

}  // namespace olb

#endif
