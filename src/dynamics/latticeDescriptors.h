/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
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
 * Descriptor for all types of 2D and 3D lattices. In principle, thanks
 * to the fact that the OpenLB code is generic, it is sufficient to
 * write a new descriptor when a new type of lattice is to be used.
 *  -- header file
 */
#ifndef LATTICE_DESCRIPTORS_H
#define LATTICE_DESCRIPTORS_H

#include "descriptorBase.h"
#include "descriptorField.h"
#include "descriptorFunction.h"

namespace olb {

/// Descriptors for the 2D and 3D lattices.
/** \warning Attention: The lattice directions must always be ordered in
 * such a way that c[i] = -c[i+(q-1)/2] for i=1..(q-1)/2, and c[0] = 0 must
 * be the rest velocity. Furthermore, the velocities c[i] for i=1..(q-1)/2
 * must verify
 *  - in 2D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *  - in 3D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *                       || (c[i][0]==0 && c[i][1]==0 && c[i][2]<0)
 * Otherwise some of the code will work erroneously, because the
 * aformentioned relations are taken as given to enable a few
 * optimizations.
*/
namespace descriptors {


/// D2Q9 lattice
template <typename... FIELDS>
struct D2Q9 : public DESCRIPTOR_BASE<2,9,POPULATION,FIELDS...> {
  typedef D2Q9<FIELDS...> BaseDescriptor;
  D2Q9() = delete;
};

namespace data {

template <>
constexpr int vicinity<2,9> = 1;

template <>
constexpr int c<2,9>[9][2] = {
  { 0, 0},
  {-1, 1}, {-1, 0}, {-1,-1}, { 0,-1},
  { 1,-1}, { 1, 0}, { 1, 1}, { 0, 1}
};

template <>
constexpr int opposite<2,9>[9] = {
  0, 5, 6, 7, 8, 1, 2, 3, 4
};

template <>
constexpr Fraction t<2,9>[9] = {
  {4, 9}, {1, 36}, {1, 9}, {1, 36}, {1, 9},
  {1, 36}, {1, 9}, {1, 36}, {1, 9}
};

template <>
constexpr Fraction cs2<2,9> = {1, 3};

}

using D2Q9Descriptor                   = D2Q9<>;
using AdvectionDiffusionD2Q9Descriptor = D2Q9<VELOCITY>;
using ForcedD2Q9Descriptor             = D2Q9<FORCE>;
using SmagorinskyForcedD2Q9Descriptor  = D2Q9<FORCE,TAU_EFF>;
using MixedScaleForcedD2Q9Descriptor   = D2Q9<FORCE,TAU_EFF,CUTOFF_KIN_ENERGY>;
using FreeEnergyD2Q9Descriptor         = D2Q9<CHEM_POTENTIAL,FORCE>;
using V6ForcedD2Q9Descriptor           = D2Q9<FORCE,V6>;


/// D2Q5 lattice
template <typename... FIELDS>
struct D2Q5 : public DESCRIPTOR_BASE<2,5,POPULATION,FIELDS...> {
  typedef D2Q5<FIELDS...> BaseDescriptor;
  D2Q5() = delete;
};

namespace data {

template <>
constexpr int vicinity<2,5> = 1;

template <>
constexpr int c<2,5>[5][2] = {
  { 0, 0},
  {-1, 0}, {0, -1}, {1,0}, { 0,1}
};

template <>
constexpr int opposite<2,5>[5] = {
  0, 3, 4, 1, 2
};

template <>
constexpr Fraction t<2,5>[5] = {
  {1, 3},
  {1, 6}, {1, 6},
  {1, 6}, {1, 6}
};

template <>
constexpr Fraction cs2<2,5> = {1, 3};

}

using D2Q5Descriptor                   = D2Q5<>;
using AdvectionDiffusionD2Q5Descriptor = D2Q5<VELOCITY>;
using SourcedAdvectionDiffusionD2Q5Descriptor = D2Q5<SOURCE,VELOCITY>;
using SmagorinskyAdvectionDiffusionD2Q5Descriptor = D2Q5<VELOCITY,TAU_EFF>;
using MixedScaleAdvectionDiffusionD2Q5Descriptor  = D2Q5<VELOCITY,TAU_EFF,CUTOFF_HEAT_FLUX>;


/// D3Q19 lattice
template <typename... FIELDS>
struct D3Q19 : public DESCRIPTOR_BASE<3,19,POPULATION,FIELDS...> {
  typedef D3Q19<FIELDS...> BaseDescriptor;
  D3Q19() = delete;
};

namespace data {

template <>
constexpr int vicinity<3,19> = 1;

template <>
constexpr int c<3,19>[19][3] = {
  { 0, 0, 0},

  {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
  {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
  {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},

  { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
  { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
  { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1}
};

template <>
constexpr int opposite<3,19>[19] = {
  0, 10, 11, 12, 13, 14, 15, 16, 17, 18, 1, 2, 3, 4, 5, 6, 7, 8, 9
};

template <>
constexpr Fraction t<3,19>[19] = {
  {1, 3},

  {1, 18}, {1, 18}, {1, 18},
  {1, 36}, {1, 36}, {1, 36},
  {1, 36}, {1, 36}, {1, 36},

  {1, 18}, {1, 18}, {1, 18},
  {1, 36}, {1, 36}, {1, 36},
  {1, 36}, {1, 36}, {1, 36}
};

template <>
constexpr Fraction cs2<3,19> = {1, 3};

}

using D3Q19Descriptor                   = D3Q19<>;
using ForcedD3Q19Descriptor             = D3Q19<FORCE>;
using SmagorinskyForcedD3Q19Descriptor  = D3Q19<FORCE,TAU_EFF>;
using V12ForcedD3Q19Descriptor          = D3Q19<FORCE,V12>;
using FreeEnergyD3Q19Descriptor         = D3Q19<CHEM_POTENTIAL,FORCE>;
using ParticleAdvectionDiffusionD3Q19Descriptor = D3Q19<VELOCITY>;

/// D3Q7 lattice
template <typename... FIELDS>
struct D3Q7 : public DESCRIPTOR_BASE<3,7,POPULATION,FIELDS...> {
  typedef D3Q7<FIELDS...> BaseDescriptor;
  D3Q7() = delete;
};

namespace data {

template <>
constexpr int vicinity<3,7> = 1;

template <>
constexpr int c<3,7>[7][3] = {
  { 0, 0, 0},

  {-1, 0, 0}, {0,-1, 0},
  { 0, 0,-1}, {1, 0, 0},
  { 0, 1, 0}, {0, 0, 1},
};

template <>
constexpr int opposite<3,7>[7] = {
  0, 4, 5, 6, 1, 2, 3
};

template <>
constexpr Fraction cs2<3,7> = {1, 4};

template <>
constexpr Fraction t<3,7>[7] = {
  {1, 4},

  {1, 8}, {1, 8}, {1, 8},
  {1, 8}, {1, 8}, {1, 8}
};

}

using D3Q7Descriptor                              = D3Q7<>;
using AdvectionDiffusionD3Q7Descriptor            = D3Q7<VELOCITY>;
using SourcedAdvectionDiffusionD3Q7Descriptor     = D3Q7<SOURCE,VELOCITY>;
using SmagorinskyAdvectionDiffusionD3Q7Descriptor = D3Q7<VELOCITY,TAU_EFF>;
using ParticleAdvectionDiffusionD3Q7Descriptor    = D3Q7<VELOCITY,VELOCITY2>;
using ParticleAdvectionDiffusionMRTD3Q7Descriptor = D3Q7<VELOCITY,VELOCITY2>;

/// D3Q13 lattice
template <typename... FIELDS>
struct D3Q13 : public DESCRIPTOR_BASE<3,13,POPULATION,FIELDS...> {
  typedef D3Q13<FIELDS...> BaseDescriptor;
  D3Q13() = delete;
};

namespace data {

template <>
constexpr int vicinity<3,13> = 1;

template <>
constexpr int c<3,13>[13][3] = {
  { 0, 0, 0},

  {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
  {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},

  { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
  { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1}
};

template <>
constexpr int opposite<3,13>[13] = {
  0, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6
};

template <>
constexpr Fraction cs2<3,13> = {1, 3};

template <>
constexpr Fraction t<3,13>[13] = {
  {1, 2},

  {1, 24}, {1, 24}, {1, 24},
  {1, 24}, {1, 24}, {1, 24},

  {1, 24}, {1, 24}, {1, 24},
  {1, 24}, {1, 24}, {1, 24}
};

template <>
constexpr Fraction lambda_e<3,13> = {3, 2};

template <>
constexpr Fraction lambda_h<3,13> = {9, 5};

}

using D3Q13Descriptor       = D3Q13<>;
using ForcedD3Q13Descriptor = D3Q13<FORCE>;


/// D3Q15 lattice
template <typename... FIELDS>
struct D3Q15 : public DESCRIPTOR_BASE<3,15,POPULATION,FIELDS...> {
  typedef D3Q15<FIELDS...> BaseDescriptor;
  D3Q15() = delete;
};

namespace data {

template <>
constexpr int vicinity<3,15> = 1;

template <>
constexpr int c<3,15>[15][3] = {
  { 0, 0, 0},

  {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
  {-1,-1,-1}, {-1,-1, 1}, {-1, 1,-1}, {-1, 1, 1},

  { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
  { 1, 1, 1}, { 1, 1,-1}, { 1,-1, 1}, { 1,-1,-1}
};

template <>
constexpr int opposite<3,15>[15] = {
  0, 8, 9, 10, 11, 12, 13, 14, 1, 2, 3, 4, 5, 6, 7
};

template <>
constexpr Fraction cs2<3,15> = {1, 3};

template <>
constexpr Fraction t<3,15>[15] = {
  {2, 9},

  {1, 9}, {1, 9}, {1, 9},
  {1, 72}, {1, 72}, {1, 72}, {1, 72},

  {1, 9}, {1, 9}, {1, 9},
  {1, 72}, {1, 72}, {1, 72}, {1, 72}
};

}

using D3Q15Descriptor = D3Q15<>;
using ForcedD3Q15Descriptor = D3Q15<FORCE>;


/// D3Q27 lattice
template <typename... FIELDS>
struct D3Q27 : public DESCRIPTOR_BASE<3,27,POPULATION,FIELDS...> {
  typedef D3Q27<FIELDS...> BaseDescriptor;
  D3Q27() = delete;
};

namespace data {

template <>
constexpr int vicinity<3,27> = 1;

template <>
constexpr int c<3,27>[27][3] = {
  { 0, 0, 0},

  {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
  {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
  {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},
  {-1,-1,-1}, {-1,-1, 1}, {-1, 1,-1}, {-1, 1, 1},

  { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
  { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
  { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1},
  { 1, 1, 1}, { 1, 1,-1}, { 1,-1, 1}, { 1,-1,-1}
};

template <>
constexpr int opposite<3,27>[27] = {
  0, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
};

template <>
constexpr Fraction cs2<3,27> = {1, 3};

template <>
constexpr Fraction t<3,27>[27] = {
  {8, 27},

  {2, 27},  {2, 27},  {2, 27},
  {1, 54},  {1, 54},  {1, 54},
  {1, 54},  {1, 54},  {1, 54},
  {1, 216}, {1, 216}, {1, 216}, {1, 216},

  {2, 27},  {2, 27},  {2, 27},
  {1, 54},  {1, 54},  {1, 54},
  {1, 54},  {1, 54},  {1, 54},
  {1, 216}, {1, 216}, {1, 216}, {1, 216}
};

}

using D3Q27Descriptor       = D3Q27<>;
using ForcedD3Q27Descriptor = D3Q27<FORCE>;


}  // namespace descriptors

}  // namespace olb

#endif
