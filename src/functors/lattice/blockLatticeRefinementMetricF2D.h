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

#ifndef BLOCK_LATTICE_REFINEMENT_METRIC_F_2D_H
#define BLOCK_LATTICE_REFINEMENT_METRIC_F_2D_H

#include "blockBaseF2D.h"
#include "core/blockLattice2D.h"
#include "core/blockLatticeStructure2D.h"


namespace olb {


template <typename T, typename DESCRIPTOR>
class BlockLatticeKnudsen2D : public BlockLatticeF2D<T, DESCRIPTOR> {
public:
  BlockLatticeKnudsen2D(BlockLatticeStructure2D<T, DESCRIPTOR>& blockLattice);

  bool operator() (T output[], const int input[]) override;
};

template <typename T, typename DESCRIPTOR>
class BlockLatticeRefinementMetricKnudsen2D final : public BlockLatticeKnudsen2D<T, DESCRIPTOR> {
private:
  const T _knudsen;
public:
  BlockLatticeRefinementMetricKnudsen2D(
    BlockLatticeStructure2D<T, DESCRIPTOR>& blockLattice,
    const UnitConverter<T, DESCRIPTOR>&     converter);

  bool operator() (T output[]);
  bool operator() (T output[], const int input[]) override;
};


}

#endif
