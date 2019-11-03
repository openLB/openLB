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

#ifndef WALLFUNCTION_BOUNDARY_POST_PROCESSORS_3D_H
#define WALLFUNCTION_BOUNDARY_POST_PROCESSORS_3D_H

#include "core/postProcessing.h"
#include "momentaOnBoundaries.h"
#include "core/blockLattice3D.h"

namespace olb {

template <typename T>
struct wallFunctionParam
{
  /*  Used method for density reconstruction
   *  0: Zou-He
   *  1: extrapolation
   *  2: constant
   */
  int rhoMethod = 1;

  /*  Used method for non-equilibrium particle distribution reconstruction
   *  0: regularized NEBB (Latt)
   *  1: extrapolation NEQ (Guo Zhaoli)
   *  2: regularized second order finite Differnce
   *  3: equilibrium scheme
   */
  int fneqMethod = 1;

  /*  Used wall profile
   *  0: Musker profile
   *  1: power law profile
   */
  int wallProfile = 0;

  /// check if descriptor with body force is used
  bool bodyForce;

  /// special formulation for straight boundaries
  bool curved = true;

  /// use van Driest damping function in boundary cell
  bool useVanDriest = true;

  /// von Karman constant for van Driest model
  T vonKarman = 0.375;

  /// wall distance in lattice units
  T latticeWalldistance = 0.5;
};

/// Musker profile
template <typename T, typename S>
class Musker : public AnalyticalF1D<T,S> {
private:
  T _nu;
  T _y;
  T _rho;
public:
  Musker(T nu, T y, T rho);
  bool operator() (T output[], const S tau_w[]);
};

/// PowerLaw profile
template <typename T, typename S>
class PowerLawProfile : public AnalyticalF1D<T,S> {
private:
  T _nu;
  T _u2;
  T _y2;
  T _y1;
  T _rho;
public:
  PowerLawProfile(T nu, T u2, T y2, T y1, T rho);
  bool operator() (T output[], const S tau_w[]);
};

template<typename T, typename DESCRIPTOR>
class WallFunctionBoundaryProcessor3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  WallFunctionBoundaryProcessor3D(int x0, int x1, int y0, int y1, int z0, int z1, BlockGeometryStructure3D<T>& blockGeometryStructure,
                                            std::vector<int> discreteNormal, std::vector<int> missingIndices,
                                            UnitConverter<T, DESCRIPTOR> const& converter, wallFunctionParam<T> const& wallFunctionParam,
                                            IndicatorF3D<T>* geoIndicator);
  virtual int extent() const {
    return 2;
  }
  virtual int extent(int whichDirection) const {
    return 2;
  }
  virtual void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice);
  virtual void processSubDomain ( BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ );
  virtual void ComputeWallFunction(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x, int y, int z);
private:
  void getIndices(int index, int value, std::vector<int>& indices);
  void calculateWallDistances(IndicatorF3D<T>* indicator);
  // FD Difference Methods
  void VelGradFromSecondOrderFD(bool NormalGradient, T Vel_BC[DESCRIPTOR::d], T Vel_1[DESCRIPTOR::d], T Vel_2[DESCRIPTOR::d], T VelGrad[DESCRIPTOR::d]);
  void computeNeighborsU(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x, int y, int z,
                         T u_x1[DESCRIPTOR::d], T u_x2[DESCRIPTOR::d], T u_y1[DESCRIPTOR::d], T u_y2[DESCRIPTOR::d], T u_z1[DESCRIPTOR::d], T u_z2[DESCRIPTOR::d]);
  void computeVelocityGradientTensor(T u_bc[DESCRIPTOR::d], T u_x1[DESCRIPTOR::d], T u_x2[DESCRIPTOR::d], T u_y1[DESCRIPTOR::d],
                                     T u_y2[DESCRIPTOR::d], T u_z1[DESCRIPTOR::d], T u_z2[DESCRIPTOR::d],  T VelGrad[DESCRIPTOR::d][DESCRIPTOR::d]);
  void computeVelocityGradient(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x, int y, int z, T u_BC[DESCRIPTOR::d],  T VelGrad[DESCRIPTOR::d][DESCRIPTOR::d]);
  void ComputeUWallNeighbor(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x, int y, int z, T (&u)[DESCRIPTOR::d]);
  void computeNeighborsRho(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x, int y, int z,
                           T u_x1[DESCRIPTOR::d], T u_x2[DESCRIPTOR::d], T u_y1[DESCRIPTOR::d], T u_y2[DESCRIPTOR::d], T u_z1[DESCRIPTOR::d], T u_z2[DESCRIPTOR::d],
                           T& rho_x1, T& rho_x2, T& rho_y1, T& rho_y2, T& rho_z1, T& rho_z2);
  void computeNeighborsRhoU(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x, int y, int z,
                           T u_x1[DESCRIPTOR::d], T u_x2[DESCRIPTOR::d], T u_y1[DESCRIPTOR::d], T u_y2[DESCRIPTOR::d], T u_z1[DESCRIPTOR::d], T u_z2[DESCRIPTOR::d],
                           T& rho_x1, T& rho_x2, T& rho_y1, T& rho_y2, T& rho_z1, T& rho_z2);

  // Van Driest Method
  void computeVanDriestTauEff(T y_bc, T tau_w, T u_bc, T u_1, T u_2, T& tau_eff);
  //
  void ComputeUWall(BlockLattice3D<T,DESCRIPTOR>& blockLattice, int x, int y, int z, T u[DESCRIPTOR::d]);
  void ComputeTauEff(BlockLattice3D<T,DESCRIPTOR>& blockLattice, Cell<T,DESCRIPTOR>& cell, int x, int y, int z, T u_bc[DESCRIPTOR::d]);
  void ComputeRhoWall(BlockLattice3D<T,DESCRIPTOR>& blockLattice, Cell<T,DESCRIPTOR>& cell, int x, int y, int z, T u_bc[DESCRIPTOR::d], T& rho_bc);

  // Methods FneqWall
  void computeRFneqfromFneq(T fneq_bc[DESCRIPTOR::q]);
  void computeFneqRNEBB(Cell<T,DESCRIPTOR>& cell, T u_bc[DESCRIPTOR::d], T rho_bc, T fneq_bc[DESCRIPTOR::q]);
  void computeFneqENeq(BlockLattice3D<T,DESCRIPTOR>& blockLattice, Cell<T,DESCRIPTOR>& cell, int x, int y, int z, T u_bc[DESCRIPTOR::d], T rho_bc, T fneq_bc[DESCRIPTOR::q]);
  void computeFneqRSOFD(BlockLattice3D<T,DESCRIPTOR>& blockLattice, Cell<T,DESCRIPTOR>& cell, int x, int y, int z, T u_bc[DESCRIPTOR::d], T rho_bc, T fneq_bc[DESCRIPTOR::q]);

  void ComputeFneqWall(BlockLattice3D<T,DESCRIPTOR>& blockLattice, Cell<T,DESCRIPTOR>& cell, int x, int y, int z, T u_bc[DESCRIPTOR::d], T rho_bc, T fneq_bc[DESCRIPTOR::q]);

  int x0, x1, y0, y1, z0, z1;
  BlockGeometryStructure3D<T>& _blockGeometryStructure;
  std::vector<int> _discreteNormal;
  std::vector<int> _missingIndices;
  UnitConverter<T, DESCRIPTOR> const& _converter;
  wallFunctionParam<T> const& _wallFunctionParam;

  T y_1;
  T y_2;

  int discreteNormalX;
  int discreteNormalY;
  int discreteNormalZ;

  int direction;
  int orientation;
  T unit_normal[3];
  std::vector<int> onWallIndices;
  std::vector<int> normalIndices;
  std::vector<int> normalInwardsIndices;
};


template<typename T, typename DESCRIPTOR>
class WallFunctionBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  WallFunctionBoundaryProcessorGenerator3D(int x0, int x1, int y0, int y1, int z0, int z1, BlockGeometryStructure3D<T>& blockGeometryStructure,
                                           std::vector<int> discreteNormal, std::vector<int> missingIndices,
                                           UnitConverter<T, DESCRIPTOR> const& converter, wallFunctionParam<T> const& wallFunctionParam, IndicatorF3D<T>* geoIndicator);
  virtual PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  virtual PostProcessorGenerator3D<T,DESCRIPTOR>*  clone() const override;
private:
  BlockGeometryStructure3D<T>& _blockGeometryStructure;
  std::vector<int> _discreteNormal;
  std::vector<int> _missingIndices;
  UnitConverter<T, DESCRIPTOR> const& _converter;
  wallFunctionParam<T> const& _wallFunctionParam;
  IndicatorF3D<T>* _geoIndicator;
};

}



#endif
