/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
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

#ifndef BOUNDARY_CONDITION_2D_H
#define BOUNDARY_CONDITION_2D_H

#include "core/blockLatticeStructure2D.h"
#include "momentaOnBoundaries2D.h"
#include "boundaryPostProcessors2D.h"
#include "dynamics/dynamics.h"
#include "geometry/blockGeometry2D.h"
#include "functors/lattice/indicator/blockIndicatorF2D.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
class OnLatticeBoundaryCondition2D {
public:
  virtual ~OnLatticeBoundaryCondition2D() { }

  virtual void addVelocityBoundary0N(int x0, int x1, int y0, int y1, T omega) =0;
  virtual void addVelocityBoundary0P(int x0, int x1, int y0, int y1, T omega) =0;
  virtual void addVelocityBoundary1N(int x0, int x1, int y0, int y1, T omega) =0;
  virtual void addVelocityBoundary1P(int x0, int x1, int y0, int y1, T omega) =0;

  virtual void addConvectionBoundary0N(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) =0;
  virtual void addConvectionBoundary0P(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) =0;
  virtual void addConvectionBoundary1N(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) =0;
  virtual void addConvectionBoundary1P(int x0, int x1, int y0, int y1, T omega, T* uAv=NULL) =0;

  virtual void addExternalVelocityCornerNN(int x, int y, T omega) =0;
  virtual void addExternalVelocityCornerNP(int x, int y, T omega) =0;
  virtual void addExternalVelocityCornerPN(int x, int y, T omega) =0;
  virtual void addExternalVelocityCornerPP(int x, int y, T omega) =0;

  virtual void addInternalVelocityCornerNN(int x, int y, T omega) =0;
  virtual void addInternalVelocityCornerNP(int x, int y, T omega) =0;
  virtual void addInternalVelocityCornerPN(int x, int y, T omega) =0;
  virtual void addInternalVelocityCornerPP(int x, int y, T omega) =0;

  /**
   * \anchor BCimplI2D
   * \name Indicator-accepting boundary condition interfaces
   *
   * \param indicator Block indicator defining boundary cells
   * \param x0,x1,y0,y1 Range of cells to be traversed
   * \{
   **/

  /// Add velocity boundary for indicated cells
  /**
   * \param omega Omega value of velocity BC
   **/
  virtual void addVelocityBoundary(BlockIndicatorF2D<T>& indicator,
                                   int x0, int x1, int y0, int y1,
                                   T omega) =0;
  /// Add slip boundary for indicated cells
  virtual void addSlipBoundary(BlockIndicatorF2D<T>& indicator,
                               int x0, int x1, int y0, int y1) =0;
  /// Add partial slip boundary for indicated cells
  /**
   * \param tuner Value between 0 (=no slip) and 1(=free slip)
   **/
  virtual void addPartialSlipBoundary(T tuner, BlockIndicatorF2D<T>& indicator, int x0,
                                      int x1, int y0, int y1) =0;
  /// Add pressure boundary for indicated cells
  /**
   * \param omega Omega value of pressure BC
   **/
  virtual void addPressureBoundary(BlockIndicatorF2D<T>& indicator,
                                   int x0, int x1, int y0, int y1,
                                   T omega) =0;
  /// Add convection boundary for indicated cells
  /**
   * \param omega Omega value of convection BC
   * \param uAv   Optional param for post processor
   **/
  virtual void addConvectionBoundary(BlockIndicatorF2D<T>& indicator,
                                     int x0, int x1, int y0, int y1,
                                     T omega, T* uAv=NULL) =0;

  virtual void addFreeEnergyWallBoundary(BlockIndicatorF2D<T>& indicator,
                                   int x0, int x1, int y0, int y1,
                                   T addend, int latticeNumber) =0;

  virtual void addFreeEnergyInletBoundary(BlockIndicatorF2D<T>& indicator,
                                   int x0, int x1, int y0, int y1, T omega,
                                   std::string type, int latticeNumber) =0;

  virtual void addFreeEnergyOutletBoundary(BlockIndicatorF2D<T>& indicator,
                                   int x0, int x1, int y0, int y1, T omega,
                                   std::string type, int latticeNumber) =0;
  ///\}

  /**
   * \name Convenience wrappers for boundary functions
   * In practice it is often preferable to define a boundary on a single material
   * number instead of instantiating an appropriate indicator by hand.
   *
   * These convenience functions are simple wrappers around the actual
   * \ref BCimplI2D "boundary condition interfaces".
   * \{
   **/

  /// Add velocity boundary for a single material number
  void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
                           int x0, int x1, int y0, int y1,
                           T omega);
  /// Add velocity boundary for any indicated cells inside the block domain
  void addVelocityBoundary(BlockIndicatorF2D<T>& indicator, T omega, bool includeOuterCells=false);
  /// Add velocity boundary for all cells of a material number inside the block domain
  void addVelocityBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, T omega, bool includeOuterCells=false);

  /// Add slip boundary for a single material number
  void addSlipBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
                       int x0, int x1, int y0, int y1);
  /// Add slip boundary for any indicated cells inside the block domain
  void addSlipBoundary(BlockIndicatorF2D<T>& indicator, bool includeOuterCells=false);
  /// Add slip boundary for all cells of a material number inside the block domain
  void addSlipBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, bool includeOuterCells=false);

  /// Add partial slip boundary for a single material number
  void addPartialSlipBoundary(T tuner, BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
                              int x0, int x1, int y0, int y1);
  /// Add partial slip boundary for any indicated cells inside the block domain
  void addPartialSlipBoundary(T tuner, BlockIndicatorF2D<T>& indicator, bool includeOuterCells=false);
  /// Add partial slip boundary for all cells of a material number inside the block domain
  void addPartialSlipBoundary(T tuner, BlockGeometryStructure2D<T>& blockGeometryStructure, int material, bool includeOuterCells=false);

  /// Add pressure boundary for a single material number
  void addPressureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
                           int x0, int x1, int y0, int y1,
                           T omega);
  /// Add pressure boundary for any indicated cells inside the block domain
  void addPressureBoundary(BlockIndicatorF2D<T>& indicator, T omega, bool includeOuterCells=false);
  /// Add pressure boundary for all cells of a material number inside the block domain
  void addPressureBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
                           T omega, bool includeOuterCells=false);

  /// Add convection boundary for a single material number
  void addConvectionBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
                             int x0, int x1, int y0, int y1,
                             T omega, T* uAv=NULL);
  /// Add convection boundary for any indicated cells inside the block domain
  void addConvectionBoundary(BlockIndicatorF2D<T>& indicator,
                             T omega, T* uAv=NULL, bool includeOuterCells=false);
  /// Add convection boundary for all cells of a material number inside the block domain
  void addConvectionBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
                             T omega, T* uAv=NULL, bool includeOuterCells=false);

  void addFreeEnergyWallBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
                                 int x0, int x1, int y0, int y1, T addend, int latticeNumber);
  void addFreeEnergyWallBoundary(BlockIndicatorF2D<T>& indicator, 
                                 T addend, int latticeNumber, bool includeOuterCells=false);
  void addFreeEnergyWallBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
                                 T addend, int latticeNumber, bool includeOuterCells=false);

  void addFreeEnergyInletBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0,
                                  int x1, int y0, int y1, T omega, std::string type, int latticeNumber);
  void addFreeEnergyInletBoundary(BlockIndicatorF2D<T>& indicator, T omega, std::string type,
                                  int latticeNumber, bool includeOuterCells=false);
  void addFreeEnergyInletBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
                                  T omega, std::string type, int latticeNumber, bool includeOuterCells=false);

  void addFreeEnergyOutletBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material, int x0,
                                   int x1, int y0, int y1, T omega, std::string type, int latticeNumber);
  void addFreeEnergyOutletBoundary(BlockIndicatorF2D<T>& indicator, T omega, std::string type,
                                   int latticeNumber, bool includeOuterCells=false);
  void addFreeEnergyOutletBoundary(BlockGeometryStructure2D<T>& blockGeometryStructure, int material,
                                   T omega, std::string type, int latticeNumber, bool includeOuterCells=false);
  ///\}

  virtual BlockLatticeStructure2D<T,DESCRIPTOR>& getBlock() =0;
  virtual BlockLatticeStructure2D<T,DESCRIPTOR> const& getBlock() const =0;

  virtual void outputOn() =0;
  virtual void outputOff() =0;
};

////////// Factory functions //////////////////////////////////////////////////

template<typename T, typename DESCRIPTOR, typename MixinDynamics=RLBdynamics<T,DESCRIPTOR> >
OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
createLocalBoundaryCondition2D(BlockLatticeStructure2D<T,DESCRIPTOR>& block);

template<typename T, typename DESCRIPTOR, typename MixinDynamics=BGKdynamics<T,DESCRIPTOR> >
OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
createInterpBoundaryCondition2D(BlockLatticeStructure2D<T,DESCRIPTOR>& block);


}

#endif
