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
 * A helper for initialising 3D boundaries -- header file.
 */

#ifndef BOUNDARY_CONDITION_3D_H
#define BOUNDARY_CONDITION_3D_H

#include "dynamics/dynamics.h"
#include "core/unitConverter.h"
#include "wallFunctionBoundaryPostProcessors3D.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
class OnLatticeBoundaryCondition3D {
public:
  virtual ~OnLatticeBoundaryCondition3D() { }

  virtual void addVelocityBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addVelocityBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addVelocityBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addVelocityBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addVelocityBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addVelocityBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;

  virtual void addPressureBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addPressureBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addPressureBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addPressureBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addPressureBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addPressureBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;

  virtual void addConvectionBoundary0N(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) =0;
  virtual void addConvectionBoundary0P(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) =0;
  virtual void addConvectionBoundary1N(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) =0;
  virtual void addConvectionBoundary1P(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) =0;
  virtual void addConvectionBoundary2N(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) =0;
  virtual void addConvectionBoundary2P(int x0, int x1, int y0, int y1, int z0, int z1, T omega, T* uAv=NULL) =0;

  virtual void addExternalVelocityEdge0NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addExternalVelocityEdge0NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addExternalVelocityEdge0PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addExternalVelocityEdge0PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addExternalVelocityEdge1NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addExternalVelocityEdge1NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addExternalVelocityEdge1PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addExternalVelocityEdge1PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addExternalVelocityEdge2NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addExternalVelocityEdge2NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addExternalVelocityEdge2PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addExternalVelocityEdge2PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;

  virtual void addInternalVelocityEdge0NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addInternalVelocityEdge0NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addInternalVelocityEdge0PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addInternalVelocityEdge0PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addInternalVelocityEdge1NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addInternalVelocityEdge1NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addInternalVelocityEdge1PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addInternalVelocityEdge1PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addInternalVelocityEdge2NN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addInternalVelocityEdge2NP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addInternalVelocityEdge2PN(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;
  virtual void addInternalVelocityEdge2PP(int x0, int x1, int y0, int y1, int z0, int z1, T omega) =0;

  virtual void addExternalVelocityCornerNNN(int x, int y, int z, T omega) =0;
  virtual void addExternalVelocityCornerNNP(int x, int y, int z, T omega) =0;
  virtual void addExternalVelocityCornerNPN(int x, int y, int z, T omega) =0;
  virtual void addExternalVelocityCornerNPP(int x, int y, int z, T omega) =0;
  virtual void addExternalVelocityCornerPNN(int x, int y, int z, T omega) =0;
  virtual void addExternalVelocityCornerPNP(int x, int y, int z, T omega) =0;
  virtual void addExternalVelocityCornerPPN(int x, int y, int z, T omega) =0;
  virtual void addExternalVelocityCornerPPP(int x, int y, int z, T omega) =0;

  virtual void addInternalVelocityCornerNNN(int x, int y, int z, T omega) =0;
  virtual void addInternalVelocityCornerNNP(int x, int y, int z, T omega) =0;
  virtual void addInternalVelocityCornerNPN(int x, int y, int z, T omega) =0;
  virtual void addInternalVelocityCornerNPP(int x, int y, int z, T omega) =0;
  virtual void addInternalVelocityCornerPNN(int x, int y, int z, T omega) =0;
  virtual void addInternalVelocityCornerPNP(int x, int y, int z, T omega) =0;
  virtual void addInternalVelocityCornerPPN(int x, int y, int z, T omega) =0;
  virtual void addInternalVelocityCornerPPP(int x, int y, int z, T omega) =0;

  /**
   * \anchor BCimplI3D
   * \name Indicator-accepting boundary condition interfaces
   *
   * \param indicator Block indicator defining boundary cells
   * \param x0,x1,y0,y1,z0,z1 Range of cells to be traversed
   * \{
   **/

  /// Add velocity boundary for indicated cells
  /**
   * \param omega Omega value of velocity BC
   **/
  virtual void addVelocityBoundary(BlockIndicatorF3D<T>& indicator,
                                   int x0, int x1, int y0, int y1, int z0, int z1,
                                   T omega) =0;
  /// Add slip boundary for indicated cells
  virtual void addSlipBoundary(BlockIndicatorF3D<T>& indicator,
                               int x0, int x1, int y0, int y1, int z0, int z1) =0;
  /// Add partial slip boundary for indicated cells
  /**
   * \param tuner Value between 0 (=no slip) and 1(=free slip)
   **/
  virtual void addPartialSlipBoundary(T tuner, BlockIndicatorF3D<T>& indicator,
                               int x0, int x1, int y0, int y1, int z0, int z1) =0;
  /// Add pressure boundary for indicated cells
  /**
   * \param omega Omega value of pressure BC
   **/
  virtual void addPressureBoundary(BlockIndicatorF3D<T>& indicator,
                                   int x0, int x1, int y0, int y1, int z0, int z1,
                                   T omega) =0;
  /// Add convection boundary for indicated cells
  /**
   * \param omega Omega value of convection BC
   * \param uAv   Optional param for post processor
   **/
  virtual void addConvectionBoundary(BlockIndicatorF3D<T>& indicator,
                                     int x0, int x1, int y0, int y1, int z0, int z1,
                                     T omega, T* uAv=NULL) =0;
  /// Add wall function boundary for indicated cells
  /**
   * converter, wallFunctionParam and geoIndicator are post processor parameters.
   **/
  virtual void addWallFunctionBoundary(BlockIndicatorF3D<T>& indicator,
                                       int x0, int x1, int y0, int y1, int z0, int z1,
                                       UnitConverter<T, DESCRIPTOR> const& converter,
                                       wallFunctionParam<T> const& wallFunctionParam,
                                       IndicatorF3D<T>* geoIndicator=NULL) =0;

  virtual void addFreeEnergyWallBoundary(BlockIndicatorF3D<T>& indicator,
                                         int x0, int x1, int y0, int y1, int z0, int z1,
                                         T addend, int latticeNumber) =0;

  virtual void addFreeEnergyInletBoundary(BlockIndicatorF3D<T>& indicator,
                                         int x0, int x1, int y0, int y1, int z0, int z1,
                                         T omega, std::string type, int latticeNumber) =0;

  virtual void addFreeEnergyOutletBoundary(BlockIndicatorF3D<T>& indicator,
                                         int x0, int x1, int y0, int y1, int z0, int z1,
                                         T omega, std::string type, int latticeNumber) =0;

  ///\}

  /**
   * \name Convenience wrappers for boundary functions
   * In practice it is often preferable to define a boundary on a single material
   * number instead of instantiating an appropriate indicator by hand.
   *
   * These convenience functions are simple wrappers around the actual
   * \ref BCimplI3D "boundary condition interfaces".
   * \{
   **/

  /// Add velocity boundary for a single material number
  void addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                           int x0, int x1, int y0, int y1, int z0, int z1,
                           T omega);
  /// Add velocity boundary for any indicated cells inside the block domain
  void addVelocityBoundary(BlockIndicatorF3D<T>& indicator, T omega, bool includeOuterCells=false);
  /// Add velocity boundary for any cells of a material number inside the block domain
  void addVelocityBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, T omega, bool includeOuterCells=false);

  /// Add slip boundary for a single material number
  void addSlipBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                       int x0, int x1, int y0, int y1, int z0, int z1) ;
  /// Add slip boundary for any indicated cells inside the block domain
  void addSlipBoundary(BlockIndicatorF3D<T>& indicator, bool includeOuterCells=false);
  /// Add slip boundary for all cells of a material number inside the block domain
  void addSlipBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material, bool includeOuterCells=false);

  /// Add partial slip boundary for a single material number
  void addPartialSlipBoundary(T tuner, BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                       int x0, int x1, int y0, int y1, int z0, int z1);
  /// Add partial slip boundary for any indicated cells inside the block domain
  void addPartialSlipBoundary(T tuner, BlockIndicatorF3D<T>& indicator, bool includeOuterCells=false);
  /// Add partial slip boundary for all cells of a material number inside the block domain
  void addPartialSlipBoundary(T tuner, BlockGeometryStructure3D<T>& blockGeometryStructure, int material, bool includeOuterCells=false);

  /// Add pressure boundary for a single material number
  void addPressureBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                           int x0, int x1, int y0, int y1, int z0, int z1,
                           T omega);
  /// Add pressure boundary for any indicated cells inside the block domain
  void addPressureBoundary(BlockIndicatorF3D<T>& indicator, T omega, bool includeOuterCells=false);
  /// Add pressure boundary for all cells of a material number inside the block domain
  void addPressureBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                           T omega, bool includeOuterCells=false);

  /// Add convection boundary for a single material number
  void addConvectionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                             int x0, int x1, int y0, int y1, int z0, int z1,
                             T omega, T* uAv=NULL);
  /// Add convection boundary for any indicated cells inside the block domain
  void addConvectionBoundary(BlockIndicatorF3D<T>& indicator, T omega, T* uAv=NULL, bool includeOuterCells=false);
  /// Add convection boundary for all cells of a material number inside the block domain
  void addConvectionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                             T omega, T* uAv=NULL, bool includeOuterCells=false);

  /// Add wall function boundary for a single material number
  void addWallFunctionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                               int x0, int x1, int y0, int y1, int z0, int z1,
                               UnitConverter<T, DESCRIPTOR> const& converter,
                               wallFunctionParam<T> const&      wallFunctionParam,
                               IndicatorF3D<T>*                 geoIndicator=NULL);
  /// Add wall function boundary for any indicated cells inside the block domain
  void addWallFunctionBoundary(BlockIndicatorF3D<T>& indicator,
                               UnitConverter<T, DESCRIPTOR> const& converter,
                               wallFunctionParam<T> const& wallFunctionParam,
                               IndicatorF3D<T>*            geoIndicator=NULL,
                               bool includeOuterCells=false);
  /// Add wall function boundary for all cells of a material number inside the block domain
  void addWallFunctionBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                               UnitConverter<T, DESCRIPTOR> const& converter,
                               wallFunctionParam<T> const& wallFunctionParam,
                               IndicatorF3D<T>*            geoIndicator=NULL,
                               bool includeOuterCells=false);

  void addFreeEnergyWallBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                                 int x0, int x1, int y0, int y1, int z0, int z1, T addend, int latticeNumber);
  void addFreeEnergyWallBoundary(BlockIndicatorF3D<T>& indicator, 
                                 T addend, int latticeNumber, bool includeOuterCells=false);
  void addFreeEnergyWallBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                                 T addend, int latticeNumber, bool includeOuterCells=false);

  void addFreeEnergyInletBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                                  int x0, int x1, int y0, int y1, int z0, int z1,
                                  T omega, std::string type, int latticeNumber);
  void addFreeEnergyInletBoundary(BlockIndicatorF3D<T>& indicator, T omega, std::string type,
                                  int latticeNumber, bool includeOuterCells=false);
  void addFreeEnergyInletBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                                  T omega, std::string type, int latticeNumber, bool includeOuterCells=false);

  void addFreeEnergyOutletBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                                   int x0, int x1, int y0, int y1, int z0, int z1,
                                   T omega, std::string type, int latticeNumber);
  void addFreeEnergyOutletBoundary(BlockIndicatorF3D<T>& indicator, T omega, std::string type,
                                   int latticeNumber, bool includeOuterCells=false);
  void addFreeEnergyOutletBoundary(BlockGeometryStructure3D<T>& blockGeometryStructure, int material,
                                   T omega, std::string type, int latticeNumber, bool includeOuterCells=false);

  ///\}

  virtual void outputOn() =0;
  virtual void outputOff() =0;

};


////////// Factory functions //////////////////////////////////////////////////

template<typename T, typename DESCRIPTOR, typename MixinDynamics=RLBdynamics<T,DESCRIPTOR>>
OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* createLocalBoundaryCondition3D(BlockLatticeStructure3D<T,DESCRIPTOR>& block);

template<typename T, typename DESCRIPTOR, typename MixinDynamics=BGKdynamics<T,DESCRIPTOR>>
OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* createInterpBoundaryCondition3D(BlockLatticeStructure3D<T,DESCRIPTOR>& block);

}

#endif
