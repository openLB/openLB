/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Robin Trunk, Sam Avis
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

#ifndef FREE_ENERGY_POST_PROCESSOR_2D_H
#define FREE_ENERGY_POST_PROCESSOR_2D_H

#include "core/spatiallyExtendedObject2D.h"
#include "core/postProcessing.h"
#include "core/blockLattice2D.h"

/* \file
 * PostProcessor classes organising the coupling between the lattices for the free energy 
 * model.
 *
 * The PostProcessor for the calculation of the chemical potential needs to be applied first,
 * as the force relies on its results. This post processor should be assigned to the first
 * lattice with the second (and third) lattices given as partner lattices.
 * Then the force post processor should be assigned to the second lattice with the first (and
 * third) given as partner lattices. Between the execution of these post processors the
 * chemical potential should be communicated.
 */

namespace olb {

/// This class calculates the chemical potential and stores it in the external field of
/// the respective lattice.
template<typename T, typename DESCRIPTOR>
class FreeEnergyChemicalPotentialCoupling2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to the surface tension (needs to be >0). [lattice units]
  /// \param[in] kappa2_ - Parameter related to the surface tension (needs to be >0). [lattice units]
  /// \param[in] kappa3_ - Parameter related to the surface tension (needs to be >0). [lattice units]
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  FreeEnergyChemicalPotentialCoupling2D(int x0_, int x1_, int y0_, int y1_, T alpha_, 
                                T kappa1_, T kappa2_, T kappa3_,
                                std::vector<SpatiallyExtendedObject2D*> partners_);
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to the surface tension (needs to be >0). [lattice units]
  /// \param[in] kappa2_ - Parameter related to the surface tension (needs to be >0). [lattice units]
  /// \param[in] kappa3_ - Parameter related to the surface tension (needs to be >0). [lattice units]
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  FreeEnergyChemicalPotentialCoupling2D(T alpha_, T kappa1_, T kappa2_, T kappa3_,
                                std::vector<SpatiallyExtendedObject2D*> partners_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_) override;
private:
  int x0, x1, y0, y1;
  T alpha, kappa1, kappa2, kappa3;
  std::vector<SpatiallyExtendedObject2D*> partners;
};

/// PostProcessor calculating the interfacial force in the free energy model. On the fist
/// lattice the force is stored for the Guo forcing scheme. On the other lattices a velocity,
/// calculated from that force, is stored which is used in the equilibrium distribution function.
/// This should be assigned to the second lattice, with the first lattice being the first partner lattice.
template<typename T, typename DESCRIPTOR>
class FreeEnergyForceCoupling2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  FreeEnergyForceCoupling2D(int x0_, int x1_, int y0_, int y1_,
                                std::vector<SpatiallyExtendedObject2D*> partners_);
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  FreeEnergyForceCoupling2D(std::vector<SpatiallyExtendedObject2D*> partners_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_) override;
private:
  int x0, x1, y0, y1;
  std::vector<SpatiallyExtendedObject2D*> partners;
};

/// PostProcessor for assigning the velocity at inlet and outlets to lattice two and three.
/// This should be assigned to the second lattice after FreeEnergyForcePostProcessor.
/// The first lattice should be the first partner lattice.
template<typename T, typename DESCRIPTOR>
class FreeEnergyInletOutletCoupling2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  FreeEnergyInletOutletCoupling2D(int x0_, int x1_, int y0_, int y1_,
                                       std::vector<SpatiallyExtendedObject2D*> partners_);
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  FreeEnergyInletOutletCoupling2D(std::vector<SpatiallyExtendedObject2D*> partners_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_) override;
private:
  int x0, x1, y0, y1;
  std::vector<SpatiallyExtendedObject2D*> partners;
};

/// PostProcessor for setting a constant density outlet.
/// This should be used before the bulk chemical potential post-
/// processor because it depends upon the result of this.
template<typename T, typename DESCRIPTOR>
class FreeEnergyDensityOutletCoupling2D : public LocalPostProcessor2D<T,DESCRIPTOR> {
public:
  /// \param[in] rho_ - Gives the value of the density constraint.
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  FreeEnergyDensityOutletCoupling2D(int x0_, int x1_, int y0_, int y1_, T rho_,
                                  std::vector<SpatiallyExtendedObject2D*> partners_);
  /// \param[in] rho_ - Gives the value of the density constraint.
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  FreeEnergyDensityOutletCoupling2D(T rho_, std::vector<SpatiallyExtendedObject2D*> partners_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice2D<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_) override;
private:
  int x0, x1, y0, y1;
  T rho;
  std::vector<SpatiallyExtendedObject2D*> partners;
};


/// Generator class for the PostProcessors calculating the chemical potential.
template<typename T, typename DESCRIPTOR>
class FreeEnergyChemicalPotentialGenerator2D : public LatticeCouplingGenerator2D<T,DESCRIPTOR> {
public:
  /// Two component free energy model
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to the surface tension (need to be >0). [lattice units]
  /// \param[in] kappa2_ - Parameter related to the surface tension (need to be >0). [lattice units]
  FreeEnergyChemicalPotentialGenerator2D(int x0_, int x1_, int y0_, int y1_, T alpha_,
                            T kappa1_, T kappa2_);
  /// Two component free energy model
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to the surface tension (need to be >0). [lattice units]
  /// \param[in] kappa2_ - Parameter related to the surface tension (need to be >0). [lattice units]
  FreeEnergyChemicalPotentialGenerator2D(T alpha_, T kappa1_, T kappa2_);
  /// Three component free energy model
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to the surface tension (need to be >0). [lattice units]
  /// \param[in] kappa2_ - Parameter related to the surface tension (need to be >0). [lattice units]
  /// \param[in] kappa3_ - Parameter related to the surface tension (need to be >0). [lattice units]
  FreeEnergyChemicalPotentialGenerator2D(int x0_, int x1_, int y0_, int y1_, T alpha_,
                                         T kappa1_, T kappa2_, T kappa3_);
  /// Three component free energy model
  /// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
  /// \param[in] kappa1_ - Parameter related to the surface tension (need to be >0). [lattice units]
  /// \param[in] kappa2_ - Parameter related to the surface tension (need to be >0). [lattice units]
  /// \param[in] kappa3_ - Parameter related to the surface tension (need to be >0). [lattice units]
  FreeEnergyChemicalPotentialGenerator2D(T alpha_, T kappa1_, T kappa2_, T kappa3_);
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  PostProcessor2D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject2D*> partners) const override;
  LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const override;
private:
  T alpha, kappa1, kappa2, kappa3;
};

/// Generator class for the PostProcessors calculating the interfacial force.
template<typename T, typename DESCRIPTOR>
class FreeEnergyForceGenerator2D : public LatticeCouplingGenerator2D<T,DESCRIPTOR> {
public:
  FreeEnergyForceGenerator2D(int x0_, int x1_, int y0_, int y1_ );
  FreeEnergyForceGenerator2D( );
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  PostProcessor2D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject2D*> partners) const override;
  LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const override;
};

/// Generator class for the PostProcessors assigning the velocity at the outlet to lattice two and three.
template<typename T, typename DESCRIPTOR>
class FreeEnergyInletOutletGenerator2D : public LatticeCouplingGenerator2D<T,DESCRIPTOR> {
public:
  FreeEnergyInletOutletGenerator2D(int x0_, int x1_, int y0_, int y1_);
  FreeEnergyInletOutletGenerator2D( );
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  PostProcessor2D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject2D*> partners) const override;
  LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const override;
};

/// Generator class for the PostProcessors assigning the density boundary condition at the outlet.
template<typename T, typename DESCRIPTOR>
class FreeEnergyDensityOutletGenerator2D : public LatticeCouplingGenerator2D<T,DESCRIPTOR> {
public:
  /// \param[in] rho_ - Gives the value of the density constraint.
  FreeEnergyDensityOutletGenerator2D(int x0_, int x1_, int y0_, int y1_, T rho_);
  /// \param[in] rho_ - Gives the value of the density constraint.
  FreeEnergyDensityOutletGenerator2D(T rho_);
  /// \param[in] partners_ - Contains one partner lattice for two fluid components, or two lattices for three components.
  PostProcessor2D<T,DESCRIPTOR>* generate(std::vector<SpatiallyExtendedObject2D*> partners) const override;
  LatticeCouplingGenerator2D<T,DESCRIPTOR>* clone() const override;
private:
  T rho;
};

}

#endif
