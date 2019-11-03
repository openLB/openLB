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
 * A helper for initialising 2D boundaries -- header file.
 */


#ifndef SUPER_BOUNDARY_CONDITION_2D_H
#define SUPER_BOUNDARY_CONDITION_2D_H

#include <vector>
#include "boundaryCondition2D.h"
#include "advectionDiffusionBoundaryCondition2D.h"     // -> for AdvectionDiffusion
#include "geometry/blockGeometryStatistics2D.h"
#include "core/superLattice2D.h"
#include "io/ostreamManager.h"
#include "geometry/superGeometry2D.h"
#include "utilities/functorPtr.h"


/// All OpenLB code is contained in this namespace.
namespace olb {

/// A helper for initialising 2D boundaries for super lattices.
/** Here we have methods that initializes for a given global
 * point or global range the local postprocessors and the
 * communicator (_commBC in SuperLattice) for boundary conditions.
 *
 * This class is not intended to be derived from.
 */
template<typename T, typename DESCRIPTOR>
class sOnLatticeBoundaryCondition2D {
public:
  /// Constructor
  sOnLatticeBoundaryCondition2D(SuperLattice2D<T,DESCRIPTOR>& sLattice);
  /// Copy construction
  sOnLatticeBoundaryCondition2D(sOnLatticeBoundaryCondition2D<T,DESCRIPTOR> const& rhs);
  /// Copy assignment
  sOnLatticeBoundaryCondition2D operator=(sOnLatticeBoundaryCondition2D<T,DESCRIPTOR> rhs);
  /// Destructor
  ~sOnLatticeBoundaryCondition2D();

  void addVelocityBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T omega);
  void addVelocityBoundary(SuperGeometry2D<T>& superGeometry, int material, T omega);

  void addSlipBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& indicator);
  void addSlipBoundary(SuperGeometry2D<T>& superGeometry, int material);

  void addPartialSlipBoundary(T tuner, FunctorPtr<SuperIndicatorF2D<T>>&& indicator);
  void addPartialSlipBoundary(T tuner, SuperGeometry2D<T>& superGeometry, int material);

  void addTemperatureBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T omega);
  void addTemperatureBoundary(SuperGeometry2D<T>& superGeometry, int material, T omega);

  void addRegularizedTemperatureBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T omega);
  void addRegularizedTemperatureBoundary(SuperGeometry2D<T>& superGeometry, int material, T omega);

  void addRegularizedHeatFluxBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T omega, T *heatFlux=nullptr);
  void addRegularizedHeatFluxBoundary(SuperGeometry2D<T>& superGeometry, int material, T omega, T *heatFlux=nullptr);

  void addPressureBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T omega);
  void addPressureBoundary(SuperGeometry2D<T>& superGeometry, int material, T omega);

  void addConvectionBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T omega, T* uAv=NULL);
  void addConvectionBoundary(SuperGeometry2D<T>& superGeometry, int material, T omega, T* uAv=NULL);
 
  /// Implementation of a wetting boundary condition for the binary free energy model, consisting of a BounceBack
  /// dynamics and an FreeEnergyWall PostProcessor.
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to surface tension. [lattice units]
  /// \param[in] kappa2_ - Parameter related to surface tension. [lattice units]
  /// \param[in] h1_ - Parameter related to resulting contact angle of the boundary. [lattice units]
  /// \param[in] h2_ - Parameter related to resulting contact angle of the boundary. [lattice units]
  /// \param[in] latticeNumber - determines the number of the free energy lattice to set the boundary accordingly
  void addFreeEnergyWallBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T alpha, T kappa1, T kappa2, T h1, T h2, int latticeNumber);
  /// Implementation of a wetting boundary condition for the binary free energy model, consisting of a BounceBack
  /// dynamics and an FreeEnergyWall PostProcessor.
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to surface tension. [lattice units]
  /// \param[in] kappa2_ - Parameter related to surface tension. [lattice units]
  /// \param[in] h1_ - Parameter related to resulting contact angle of the boundary. [lattice units]
  /// \param[in] h2_ - Parameter related to resulting contact angle of the boundary. [lattice units]
  /// \param[in] latticeNumber - determines the number of the free energy lattice to set the boundary accordingly
  void addFreeEnergyWallBoundary(SuperGeometry2D<T>& superGeometry, int material, T alpha, T kappa1, T kappa2, T h1, T h2, int latticeNumber);
 
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
  void addFreeEnergyWallBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T alpha, T kappa1, T kappa2, T kappa3, T h1, T h2, T h3, int latticeNumber);
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
  void addFreeEnergyWallBoundary(SuperGeometry2D<T>& superGeometry, int material, T alpha, T kappa1, T kappa2, T kappa3, T h1, T h2, T h3, int latticeNumber);

  /// Implementation of a inlet boundary condition for the partner lattices of the binary or ternary free energy model.
  void addFreeEnergyInletBoundary(SuperGeometry2D<T>& superGeometry, int material, T omega, std::string type, int latticeNumber);
  /// Implementation of a inlet boundary condition for the partner lattices of the binary or ternary free energy model.
  void addFreeEnergyInletBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T omega, std::string type, int latticeNumber);

  /// Implementation of a outlet boundary condition for the partner lattices of the binary or ternary free energy model.
  void addFreeEnergyOutletBoundary(SuperGeometry2D<T>& superGeometry, int material, T omega, std::string type, int latticeNumber);
  /// Implementation of a outlet boundary condition for the partner lattices of the binary or ternary free energy model.
  void addFreeEnergyOutletBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& indicator, T omega, std::string type, int latticeNumber);

  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  void addPoints2CommBC(FunctorPtr<SuperIndicatorF2D<T>>&& indicator);
  void addPoints2CommBC(SuperGeometry2D<T>& superGeometry, int material);

  SuperLattice2D<T,DESCRIPTOR>& getSuperLattice()
  {
    return _sLattice;
  };
  std::vector<OnLatticeBoundaryCondition2D<T,DESCRIPTOR>* >& getBlockBCs()
  {
    return _blockBCs;
  };
  std::vector<OnLatticeAdvectionDiffusionBoundaryCondition2D<T, DESCRIPTOR>*>& getADblockBCs()
  {
    return _ADblockBCs;
  };
  int getOverlap()
  {
    return _overlap;
  };
  void setOverlap(int overlap)
  {
    _overlap = overlap;
  };

  void outputOn();
  void outputOff();

private:
  mutable OstreamManager clout;
  SuperLattice2D<T,DESCRIPTOR>& _sLattice;
  std::vector<OnLatticeBoundaryCondition2D<T,DESCRIPTOR>* > _blockBCs;
  std::vector<OnLatticeAdvectionDiffusionBoundaryCondition2D<T, DESCRIPTOR>*> _ADblockBCs;       // -> for AdvectionDiffusion
  int _overlap;
  bool _output;
};


template<typename T, typename DESCRIPTOR, typename MixinDynamics=RLBdynamics<T,DESCRIPTOR> >
void createLocalBoundaryCondition2D(sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& sBC);

template<typename T, typename DESCRIPTOR, typename MixinDynamics=BGKdynamics<T,DESCRIPTOR> >
void createInterpBoundaryCondition2D(sOnLatticeBoundaryCondition2D<T,DESCRIPTOR>& sBC);

template<typename T, typename DESCRIPTOR, typename MixinDynamics=BGKdynamics<T,DESCRIPTOR> >
void createExtFdBoundaryCondition2D(sOnLatticeBoundaryCondition2D<T, DESCRIPTOR>& sBC);

}  // namespace olb

#endif
