/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#ifndef SUPER_LATTICE_REFINEMENT_METRIC_F_2D_H
#define SUPER_LATTICE_REFINEMENT_METRIC_F_2D_H

#include "superBaseF2D.h"
#include "core/superLattice2D.h"


namespace olb {


/// SuperLatticeKnudsen2D measures cell-local ratio between non-equilibrium and equilibrium distribution
/**
 * This value approximates the dimensionless Knudsen number in sufficiently well-resolved grids.
 *
 * As such it may be used to infer information on the quality of the simulation.
 *
 * See SuperLatticeRefinementMetricKnudsen2D.
 **/
template <typename T, typename DESCRIPTOR>
class SuperLatticeKnudsen2D final : public SuperLatticeF2D<T, DESCRIPTOR> {
public:
  SuperLatticeKnudsen2D(SuperLattice2D<T, DESCRIPTOR>& lattice);

  bool operator() (T output[], const int input[]) override;
};

/// SuperLatticeRefinementMetricKnudsen2D suggests a per-block grid refinement factor
/**
 *  This functor is a direct implementation of the refinement criterion described in
 *  "Automatic grid refinement criterion for lattice Boltzmann method" [arXiv:1507.06767]
 *  by D. Lagrava et al.
 **/
template <typename T, typename DESCRIPTOR>
class SuperLatticeRefinementMetricKnudsen2D final : public SuperLatticeF2D<T, DESCRIPTOR> {
public:
  /// Constructor
  /**
   * \param converter Required to get the theoretical Knudsen number as reference
   **/
  SuperLatticeRefinementMetricKnudsen2D(
    SuperLattice2D<T, DESCRIPTOR>&      lattice,
    const UnitConverter<T, DESCRIPTOR>& converter);

  /**
   * \param output
   *        output[0] == 0 suggests sufficient grid resolution,
   *        output[0]  > 0 suggests the grid to be refined output[0] times
   * \param glob
   *        cuboid to be considered
   **/
  bool operator() (T output[], int glob);
  /**
   * \param output
   *        output[0] == 0 suggests sufficient grid resolution,
   *        output[0]  > 0 suggests the grid to be refined output[0] times
   * \input input
   *        grid cell to be considered
   **/
  bool operator() (T output[], const int input[]) override;

  /// Convenience method for printing per-block refinement factors
  void print();
};


}

#endif
