/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2016 Jonas Kratzke, Mathias J. Krause
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


#ifndef SUPER_OFF_BOUNDARY_CONDITION_2D_H
#define SUPER_OFF_BOUNDARY_CONDITION_2D_H

#include <vector>
#include <list>

#include "offBoundaryCondition2D.h"
#include "geometry/superGeometry2D.h"
#include "core/superLattice2D.h"
#include "io/ostreamManager.h"
#include "functors/analytical/analyticalF.h"
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
class sOffLatticeBoundaryCondition2D {
public:
  /// Constructor
  sOffLatticeBoundaryCondition2D(SuperLattice2D<T,DESCRIPTOR>& sLattice, T epsFraction_ = 0.0001);
  /// Copy construction
  sOffLatticeBoundaryCondition2D(sOffLatticeBoundaryCondition2D<T,DESCRIPTOR> const& rhs);
  /// Copy assignment
  sOffLatticeBoundaryCondition2D operator=(sOffLatticeBoundaryCondition2D<T,DESCRIPTOR> rhs);
  /// Destructor
  ~sOffLatticeBoundaryCondition2D();

  /// Set offDynamics with boundary links and post processors using indicators
  /**
   * Add offDynamics with initialisation of boundary links and the corresponding
   * post processors
   * Note: Uses information of the second neighbours of the cell (x,y)
   * Add post processors. Ensure that offDynamics are defined!
   *
   * \param boundaryIndicator Indicator describing boundary cells
   * \param bulkIndicator     Indicator describing bulk cells
   * \param geometryIndicator Indicator describing the geometry to be bounded
   **/
  void addVelocityBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& boundaryIndicator,
                           FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
                           IndicatorF2D<T>&                   geometryIndicator);
  void addVelocityBoundary(SuperGeometry2D<T>& superGeometry, int material,
                           IndicatorF2D<T>& geometryIndicator,
                           std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void addVelocityBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& boundaryIndicator,
                           FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator);
  void addVelocityBoundary(SuperGeometry2D<T>& superGeometry, int material,
                           std::vector<int> bulkMaterials = std::vector<int>(1,1));

  void addZeroVelocityBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& boundaryIndicator,
                               FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
                               IndicatorF2D<T>&                   geometryIndicator);
  void addZeroVelocityBoundary(SuperGeometry2D<T>& superGeometry, int material,
                               IndicatorF2D<T>& geometryIndicator,
                               std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void addZeroVelocityBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& boundaryIndicator,
                               FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator);
  void addZeroVelocityBoundary(SuperGeometry2D<T>& superGeometry, int material,
                               std::vector<int> bulkMaterials = std::vector<int>(1,1));

  void addPressureBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& boundaryIndicator,
                           FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
                           IndicatorF2D<T>& geometryIndicator);
  void addPressureBoundary(SuperGeometry2D<T>& superGeometry, int material,
                           IndicatorF2D<T>& geometryIndicator,
                           std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void addPressureBoundary(FunctorPtr<SuperIndicatorF2D<T>>&& boundaryIndicator,
                           FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator);
  void addPressureBoundary(SuperGeometry2D<T>& superGeometry, int material,
                           std::vector<int> bulkMaterials = std::vector<int>(1,1));

  void defineU(FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
               FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
               AnalyticalF2D<T,T>& u);
  void defineU(SuperGeometry2D<T>& superGeometry, int material,
               AnalyticalF2D<T,T>& u,
               std::vector<int> bulkMaterials = std::vector<int>(1,1));

  void defineRho(FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
                 FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
                 AnalyticalF2D<T,T>&                rho);
  void defineRho(SuperGeometry2D<T>& superGeometry, int material,
                 AnalyticalF2D<T,T>& rho,
                 std::vector<int> bulkMaterials = std::vector<int>(1,1));

  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  void addPoints2CommBC(FunctorPtr<SuperIndicatorF2D<T>>&& indicator);
  void addPoints2CommBC(SuperGeometry2D<T>& superGeometry, int material);

  SuperLattice2D<T,DESCRIPTOR>& getSuperLattice()
  {
    return _sLattice;
  };
  std::vector<OffLatticeBoundaryCondition2D<T,DESCRIPTOR>* >& getBlockBCs()
  {
    return _blockBCs;
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
  std::vector<OffLatticeBoundaryCondition2D<T,DESCRIPTOR>* > _blockBCs;
  T _epsFraction;
  int _overlap;
  bool _output;
};

template<typename T, typename DESCRIPTOR, typename MixinDynamics=BGKdynamics<T,DESCRIPTOR> >
void createBouzidiBoundaryCondition2D(sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>& sBC);

template<typename T, typename DESCRIPTOR>
void createBounceBackBoundaryCondition2D(sOffLatticeBoundaryCondition2D<T,DESCRIPTOR>& sBC);

}  // namespace olb

#endif
