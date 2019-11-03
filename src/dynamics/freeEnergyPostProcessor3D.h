/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Sam Avis, Robin Trunk
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

#ifndef FREE_ENERGY_POST_PROCESSOR_3D_H
#define FREE_ENERGY_POST_PROCESSOR_3D_H

#include "core/spatiallyExtendedObject3D.h"
#include "core/postProcessing.h"
#include "core/blockLattice3D.h"

/* \file
 * PostProcessor classes organising the coupling between the lattices for the free energy 
 * model. The PostProcessor for the calculation of the chemical potential needs to be 
 * applied first, as the force relies on its results.
 */

namespace olb {

/// This class calculates the chemical potential and stores it in the external field of
/// the respective lattice.
template<typename T, typename DESCRIPTOR>
class FreeEnergyChemicalPotentialCoupling3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to the surface tension (needs to be >0). [lattice units]
  /// \param[in] kappa2_ - Parameter related to the surface tension (needs to be >0). [lattice units]
  /// \param[in] kappa3_ - Parameter related to the surface tension (needs to be >0). [lattice units]
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  FreeEnergyChemicalPotentialCoupling3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                                             T alpha_, T kappa1_, T kappa2_, T kappa3_,
                                             std::vector<SpatiallyExtendedObject3D*> partners_);
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to the surface tension (needs to be >0). [lattice units]
  /// \param[in] kappa2_ - Parameter related to the surface tension (needs to be >0). [lattice units]
  /// \param[in] kappa3_ - Parameter related to the surface tension (needs to be >0). [lattice units]
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  FreeEnergyChemicalPotentialCoupling3D(T alpha_, T kappa1_, T kappa2_, T kappa3_,
                                             std::vector<SpatiallyExtendedObject3D*> partners_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
private:
  int x0, x1, y0, y1, z0, z1;
  T alpha, kappa1, kappa2, kappa3;
  std::vector<SpatiallyExtendedObject3D*> partners;
};

/// PostProcessor calculating the interfacial force in the free energy model. On the fist
/// lattice the force is stored for the Guo forcing scheme. On the other lattices a velocity,
/// calculated from that force, is stored which is used in the equilibrium distribution function.
template<typename T, typename DESCRIPTOR>
class FreeEnergyForceCoupling3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  FreeEnergyForceCoupling3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                                 std::vector<SpatiallyExtendedObject3D*> partners_);
  FreeEnergyForceCoupling3D(std::vector<SpatiallyExtendedObject3D*> partners_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
private:
  int x0, x1, y0, y1, z0, z1;
  std::vector<SpatiallyExtendedObject3D*> partners;
};

/// PostProcessor for assigning the velocity at inlet and outlets to lattice two and three.
/// This should be assigned to the second lattice after FreeEnergyForcePostProcessor.
/// The first lattice should be the first partner lattice.
template<typename T, typename DESCRIPTOR>
class FreeEnergyInletOutletCoupling3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  FreeEnergyInletOutletCoupling3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                                       std::vector<SpatiallyExtendedObject3D*> partners_);
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  FreeEnergyInletOutletCoupling3D(std::vector<SpatiallyExtendedObject3D*> partners_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
private:
  int x0, x1, y0, y1, z0, z1;
  std::vector<SpatiallyExtendedObject3D*> partners;
};

/// PostProcessor for setting a constant density outlet.
/// This should be used before the bulk chemical potential post-
/// processor because it depends upon the result of this.
template<typename T, typename DESCRIPTOR>
class FreeEnergyDensityOutletCoupling3D : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  /// \param[in] rho_ - Gives the value of the constraint.
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  FreeEnergyDensityOutletCoupling3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                                  T rho_, std::vector<SpatiallyExtendedObject3D*> partners_);
  /// \param[in] rho_ - Gives the value of the constraint.
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  FreeEnergyDensityOutletCoupling3D(T rho_, std::vector<SpatiallyExtendedObject3D*> partners_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice3D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
private:
  int x0, x1, y0, y1, z0, z1;
  T rho;
  std::vector<SpatiallyExtendedObject3D*> partners;
};


/// Generator class for the PostProcessors calculating the chemical potential.
template<typename T, typename DESCRIPTOR>
class FreeEnergyChemicalPotentialGenerator3D : public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  /// Two component free energy model
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to the surface tension (need to be >0). [lattice units]
  /// \param[in] kappa2_ - Parameter related to the surface tension (need to be >0). [lattice units]
  FreeEnergyChemicalPotentialGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                                         T alpha_, T kappa1_, T kappa2_);
  /// Two component free energy model
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to the surface tension (need to be >0). [lattice units]
  /// \param[in] kappa2_ - Parameter related to the surface tension (need to be >0). [lattice units]
  FreeEnergyChemicalPotentialGenerator3D(T alpha_, T kappa1_, T kappa2_);
  /// Three component free energy model
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to the surface tension (need to be >0). [lattice units]
  /// \param[in] kappa2_ - Parameter related to the surface tension (need to be >0). [lattice units]
  /// \param[in] kappa3_ - Parameter related to the surface tension (need to be >0). [lattice units]
  FreeEnergyChemicalPotentialGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                                         T alpha_, T kappa1_, T kappa2_, T kappa3_);
  /// Three component free energy model
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to the surface tension (need to be >0). [lattice units]
  /// \param[in] kappa2_ - Parameter related to the surface tension (need to be >0). [lattice units]
  /// \param[in] kappa3_ - Parameter related to the surface tension (need to be >0). [lattice units]
  FreeEnergyChemicalPotentialGenerator3D(T alpha_, T kappa1_, T kappa2_, T kappa3_);
  PostProcessor3D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject3D*> partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;
private:
  T alpha, kappa1, kappa2, kappa3;
};

/// Generator class for the PostProcessors calculating the interfacial force.
template<typename T, typename DESCRIPTOR>
class FreeEnergyForceGenerator3D : public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  FreeEnergyForceGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);
  FreeEnergyForceGenerator3D( );
  PostProcessor3D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject3D*> partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;
};

/// Generator class for the PostProcessors assigning the velocity at the outlet to lattice two and three.
template<typename T, typename DESCRIPTOR>
class FreeEnergyInletOutletGenerator3D : public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  FreeEnergyInletOutletGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);
  FreeEnergyInletOutletGenerator3D( );
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  PostProcessor3D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject3D*> partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;
};

/// Generator class for the PostProcessors assigning the density boundary condition at the outlet.
template<typename T, typename DESCRIPTOR>
class FreeEnergyDensityOutletGenerator3D : public LatticeCouplingGenerator3D<T,DESCRIPTOR> {
public:
  /// \param[in] rho_ - Gives the value of the constraint.
  FreeEnergyDensityOutletGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T rho_);
  /// \param[in] rho_ - Gives the value of the constraint.
  FreeEnergyDensityOutletGenerator3D(T rho_);
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  PostProcessor3D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject3D*> partners) const override;
  LatticeCouplingGenerator3D<T,DESCRIPTOR>* clone() const override;
private:
  T rho;
};

}

#endif
