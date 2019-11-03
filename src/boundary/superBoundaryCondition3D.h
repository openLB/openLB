/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Mathias J. Krause
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
 * A helper for initialising 3D boundaries -- header file.
 */

#ifndef SUPER_BOUNDARY_CONDITION_3D_H
#define SUPER_BOUNDARY_CONDITION_3D_H

#include <vector>
#include "io/ostreamManager.h"
#include "utilities/functorPtr.h"
#include "extendedFiniteDifferenceBoundary3D.h"

/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T, typename DESCRIPTOR> class OnLatticeAdvectionDiffusionBoundaryCondition3D;
template<typename T, typename DESCRIPTOR> class OnLatticeBoundaryCondition3D;
template<typename T, typename DESCRIPTOR> class SuperLattice3D;
template<typename T> class SuperGeometry3D;
template<typename T> class SuperIndicatorF3D;

/// A helper for initialising 3D boundaries for super lattices.
/** Here we have methods that initializes the local postprocessors and the
 * communicator (_commBC in SuperLattice) for boundary conditions
 * for a given global point or global range.
 *
 * This class is not intended to be derived from.
 */
template<typename T, typename DESCRIPTOR>
class sOnLatticeBoundaryCondition3D {
public:
  sOnLatticeBoundaryCondition3D(SuperLattice3D<T, DESCRIPTOR>& sLattice);
  sOnLatticeBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> const& rhs);
  sOnLatticeBoundaryCondition3D operator=(sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> rhs);
  ~sOnLatticeBoundaryCondition3D();

  void addVelocityBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T omega);
  void addVelocityBoundary(SuperGeometry3D<T>& superGeometry, int material, T omega);

  void addSlipBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator);
  void addSlipBoundary(SuperGeometry3D<T>& superGeometry, int material);

  void addPartialSlipBoundary(T tuner, FunctorPtr<SuperIndicatorF3D<T>>&& indicator);
  void addPartialSlipBoundary(T tuner, SuperGeometry3D<T>& superGeometry, int material);

  void addPressureBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T omega);
  void addPressureBoundary(SuperGeometry3D<T>& superGeometry, int material, T omega);

  void addConvectionBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                             T omega, T* uAv=NULL);
  void addConvectionBoundary(SuperGeometry3D<T>& superGeometry, int material,
                             T omega, T* uAv=NULL);
  void addConvectionBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator);
  void addConvectionBoundary(SuperGeometry3D<T>& superGeometry, int material);

  void addWallFunctionBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
                               UnitConverter<T, DESCRIPTOR> const& converter,
                               wallFunctionParam<T> const& wallFunctionParam,
                               IndicatorF3D<T>* geoIndicator=NULL);
  void addWallFunctionBoundary(SuperGeometry3D<T>& superGeometry, int material,
                               UnitConverter<T, DESCRIPTOR> const& converter,
                               wallFunctionParam<T> const& wallFunctionParam,
                               IndicatorF3D<T>* geoIndicator=NULL);

  void addTemperatureBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T omega);
  void addTemperatureBoundary(SuperGeometry3D<T>& superGeometry, int material, T omega);

  void addExtFieldBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator, int offset);
  void addExtFieldBoundary(SuperGeometry3D<T>& superGeometry, int material, int offset);

  void addZeroDistributionBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator);
  void addZeroDistributionBoundary(SuperGeometry3D<T>& superGeometry, int material);

  /// Implementation of a wetting boundary condition for the binary free energy model, consisting of a BounceBack
  /// dynamics and an FreeEnergyWall PostProcessor.
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to surface tension. [lattice units]
  /// \param[in] kappa2_ - Parameter related to surface tension. [lattice units]
  /// \param[in] h1_ - Parameter related to resulting contact angle of the boundary. [lattice units]
  /// \param[in] h2_ - Parameter related to resulting contact angle of the boundary. [lattice units]
  /// \param[in] latticeNumber - determines the number of the free energy lattice to set the boundary accordingly
  void addFreeEnergyWallBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
		  T alpha, T kappa1, T kappa2, T h1, T h2, int latticeNumber);
  /// Implementation of a wetting boundary condition for the binary free energy model, consisting of a BounceBack
  /// dynamics and an FreeEnergyWall PostProcessor.
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to surface tension. [lattice units]
  /// \param[in] kappa2_ - Parameter related to surface tension. [lattice units]
  /// \param[in] h1_ - Parameter related to resulting contact angle of the boundary. [lattice units]
  /// \param[in] h2_ - Parameter related to resulting contact angle of the boundary. [lattice units]
  /// \param[in] latticeNumber - determines the number of the free energy lattice to set the boundary accordingly
  void addFreeEnergyWallBoundary(SuperGeometry3D<T>& superGeometry, int material,
		  T alpha, T kappa1, T kappa2, T h1, T h2, int latticeNumber);

  /// Implementation of a wetting boundary condition for the ternary free energy model, consisting of a BounceBack
  /// dynamics and an FreeEnergyWall PostProcessor.
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to surface tension. [lattice units]
  /// \param[in] kappa2_ - Parameter related to surface tension. [lattice units]
  /// \param[in] kappa3_ - Parameter related to surface tension. [lattice units]
  /// \param[in] h1_ - Parameter related to resulting contact angle of the boundary. [lattice units]
  /// \param[in] h2_ - Parameter related to resulting contact angle of the boundary. [lattice units]
  /// \param[in] h3_ - Parameter related to resulting contact angle of the boundary. [lattice units]
  /// \param[in] latticeNumber - determines the number of the free energy lattice to set the boundary accordingly
  void addFreeEnergyWallBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
		  T alpha, T kappa1, T kappa2, T kappa3, T h1, T h2, T h3, int latticeNumber);
  /// Implementation of a wetting boundary condition for the ternary free energy model, consisting of a BounceBack
  /// dynamics and an FreeEnergyWall PostProcessor.
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to surface tension. [lattice units]
  /// \param[in] kappa2_ - Parameter related to surface tension. [lattice units]
  /// \param[in] kappa3_ - Parameter related to surface tension. [lattice units]
  /// \param[in] h1_ - Parameter related to resulting contact angle of the boundary. [lattice units]
  /// \param[in] h2_ - Parameter related to resulting contact angle of the boundary. [lattice units]
  /// \param[in] h3_ - Parameter related to resulting contact angle of the boundary. [lattice units]
  /// \param[in] latticeNumber - determines the number of the free energy lattice to set the boundary accordingly
  void addFreeEnergyWallBoundary(SuperGeometry3D<T>& superGeometry, int material,
		  T alpha, T kappa1, T kappa2, T kappa3, T h1, T h2, T h3, int latticeNumber);
 
  /// Implementation of a inlet boundary condition for the partner lattices of the binary or the ternary free energy model.
  void addFreeEnergyInletBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T omega, std::string type, int latticeNumber);
  /// Implementation of a inlet boundary condition for the partner lattices of the binary or the ternary free energy model.
  void addFreeEnergyInletBoundary(SuperGeometry3D<T>& superGeometry, int material, T omega, std::string type, int latticeNumber);
 
  /// Implementation of a outlet boundary condition for the partner lattices of the binary or the ternary free energy model.
  void addFreeEnergyOutletBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T omega, std::string type, int latticeNumber);
  /// Implementation of a outlet boundary condition for the partner lattices of the binary or the ternary free energy model.
  void addFreeEnergyOutletBoundary(SuperGeometry3D<T>& superGeometry, int material, T omega, std::string type, int latticeNumber);

  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  void addPoints2CommBC(FunctorPtr<SuperIndicatorF3D<T>>&& indicator);
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  void addPoints2CommBC(SuperGeometry3D<T>& superGeometry, int material);

  SuperLattice3D<T, DESCRIPTOR>& getSuperLattice();
  std::vector<OnLatticeBoundaryCondition3D<T, DESCRIPTOR>*>& getBlockBCs();
  std::vector<OnLatticeAdvectionDiffusionBoundaryCondition3D<T, DESCRIPTOR>*>& getADblockBCs();
  int getOverlap();
  void setOverlap(int overlap);

  void outputOn();
  void outputOff();

private:
  mutable OstreamManager clout;
  SuperLattice3D<T, DESCRIPTOR>& _sLattice;
  std::vector<OnLatticeBoundaryCondition3D<T, DESCRIPTOR>*> _blockBCs;
  std::vector<OnLatticeAdvectionDiffusionBoundaryCondition3D<T, DESCRIPTOR>*> _ADblockBCs;
  int _overlap;
  bool _output;
};


////////////////// Factory functions //////////////////////////////////

template<typename T, typename DESCRIPTOR>
void createLocalBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& sBC);

template<typename T, typename DESCRIPTOR, typename MixinDynamics=BGKdynamics<T,DESCRIPTOR> >
void createInterpBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& sBC);

template<typename T, typename DESCRIPTOR, typename MixinDynamics=BGKdynamics<T,DESCRIPTOR> >
void createExtFdBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T, DESCRIPTOR>& sBC);


} // namespace olb

#endif
