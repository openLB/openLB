/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Jonas Kratzke, Mathias J. Krause
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


#ifndef SUPER_OFF_BOUNDARY_CONDITION_3D_H
#define SUPER_OFF_BOUNDARY_CONDITION_3D_H

#include <vector>
#include <list>
#include "offBoundaryCondition3D.h"
#include "geometry/superGeometry3D.h"
#include "core/superLattice3D.h"
#include "io/ostreamManager.h"
#include "functors/analytical/analyticalF.h"
#include "utilities/functorPtr.h"


/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T> class SuperIndicatorF3D;

/// A helper for initialising 3D boundaries for super lattices.
/**
 * Here we have methods that initializes for a given global
 * point or global range the local postprocessors and the
 * communicator (_commBC in SuperLattice) for boundary conditions.
 *
 * This class is not intended to be derived from.
 */
template<typename T, typename DESCRIPTOR>
class sOffLatticeBoundaryCondition3D {
public:
  /// Constructor
  sOffLatticeBoundaryCondition3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, T epsFraction_ = 0.0001);
  /// Copy construction
  sOffLatticeBoundaryCondition3D(sOffLatticeBoundaryCondition3D<T,DESCRIPTOR> const& rhs);
  /// Copy assignment
  sOffLatticeBoundaryCondition3D operator=(sOffLatticeBoundaryCondition3D<T,DESCRIPTOR> rhs);
  /// Destructor
  ~sOffLatticeBoundaryCondition3D();

  /// Set offDynamics with boundary links and post processors using indicators
  /**
   * Add offDynamics with initialisation of boundary links and the corresponding
   * post processors
   * Note: Uses information of the second neighbours of the cell (x,y,z)
   * Add post processors. Ensure that offDynamics are defined!
   *
   * \param boundaryIndicator Indicator describing boundary cells
   * \param bulkIndicator     Indicator describing bulk cells
   * \param geometryIndicator Indicator describing the geometry to be bounded
   **/
  void addVelocityBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator,
                           FunctorPtr<SuperIndicatorF3D<T>>&& bulkIndicator,
                           IndicatorF3D<T>&                   geometryIndicator);
  void addVelocityBoundary(SuperGeometry3D<T>& superGeometry, int material,
                           IndicatorF3D<T>& indicator,
                           std::vector<int> bulkMaterials = std::vector<int>(1,1));

  void addZeroVelocityBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator,
                               FunctorPtr<SuperIndicatorF3D<T>>&& bulkIndicator,
                               IndicatorF3D<T>&                   geometryIndicator);
  void addZeroVelocityBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator,
                               IndicatorF3D<T>&                   geometryIndicator,
                               std::vector<int> bulkMaterials = std::vector<int>(1,1));
  void addZeroVelocityBoundary(SuperGeometry3D<T>& superGeometry, int material,
                               IndicatorF3D<T>& geometryIndicator,
                               std::vector<int> bulkMaterials = std::vector<int>(1,1));

  void defineU(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
               FunctorPtr<SuperIndicatorF3D<T>>&& bulkIndicator,
               AnalyticalF3D<T,T>& u);
  void defineU(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
               AnalyticalF3D<T,T>& u,
               std::vector<int> bulkMaterials = std::vector<int>(1,1) );
  void defineU(SuperGeometry3D<T>& superGeometry, int material,
               AnalyticalF3D<T,T>& u,
               std::vector<int> bulkMaterials = std::vector<int>(1,1) );

  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  void addPoints2CommBC(FunctorPtr<SuperIndicatorF3D<T>>&& indicator);
  /// Adds needed Cells to the Communicator _commBC in SuperLattice
  void addPoints2CommBC(SuperGeometry3D<T>& superGeometry, int material);

  SuperLattice3D<T,DESCRIPTOR>& getSuperLattice();
  std::vector<OffLatticeBoundaryCondition3D<T,DESCRIPTOR>* >& getBlockBCs();
  int getOverlap();
  void setOverlap(int overlap);

  void outputOn();
  void outputOff();

private:
  mutable OstreamManager clout;
  SuperLattice3D<T,DESCRIPTOR>& _sLattice;
  std::vector<OffLatticeBoundaryCondition3D<T,DESCRIPTOR>* > _blockBCs;
  T _epsFraction;
  int _overlap;
  bool _output;
};

template<typename T, typename DESCRIPTOR, typename MixinDynamics=BGKdynamics<T,DESCRIPTOR> >
void createBouzidiBoundaryCondition3D(sOffLatticeBoundaryCondition3D<T,DESCRIPTOR>& sBC);


}  // namespace olb

#endif
