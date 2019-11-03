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
 * DESCRIPTORBASE for all types of 2D and 3D lattices. In principle, thanks
 * to the fact that the OpenLB code is generic, it is sufficient to
 * write a new descriptor when a new type of lattice is to be used.
 *  -- header file
 */
#ifndef MRT_LATTICE_DESCRIPTORS_H
#define MRT_LATTICE_DESCRIPTORS_H

#include "latticeDescriptors.h"

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

namespace tag {

struct MRT : public CATEGORY, public DESCRIPTOR_TAG { };

}

using MRTD2Q5Descriptor = D2Q5<tag::MRT>;

/// Advection Diffusion MRT D2Q5
/**
 * Based on: Liu, Q., & He, Y. L. (2015). Double multiple-relaxation-time lattice Boltzmann model
 *           for solid–liquid phase change with natural convection in porous media.
 *           Physica A: Statistical Mechanics and its Applications, 438, 94-106.
 **/
using AdvectionDiffusionMRTD2Q5Descriptor = D2Q5<tag::MRT,VELOCITY>;

/// MRT D2Q9 lattice. The numbering follows the one in "Viscous flow computations
/// with the method of lattice Boltzmann equation", D. Yu, L.-S. Luo, W. Shi,
/// Progress in Aerospace Sciences 39, (2003), p. 329-367
using MRTD2Q9Descriptor = D2Q9<tag::MRT>;

using ForcedMRTD2Q9Descriptor = D2Q9<tag::MRT,FORCE>;

using MRTD3Q7Descriptor = D3Q7<tag::MRT>;

/// Advection Diffusion MRT D3Q7
/**
 * Based on: Wu, H., Wang, J., & Tao, Z. (2011). Passive heat transfer in a turbulent
 *           channel flow simulation using large eddy simulation based on the lattice
 *           Boltzmann method framework.
 *           International Journal of Heat and Fluid Flow, 32(6), 1111-1119.
 *
 * There are some differences in respect to the order of the columns based on the lattice directions
 **/
using AdvectionDiffusionMRTD3Q7Descriptor = D3Q7<tag::MRT,VELOCITY>;

/// MRT D3Q19 lattice. The numbering follows the one in "Multiple-relaxation-
/// time lattice Boltzmann models in three dimensions", D. D'Humières,
/// I. Ginzburg, M. Krafzcyk, P. Lallemand, L.-S. Luo,
/// Phil. Trans. R. Soc. Lond. A (2002) 660, p. 437-451
using MRTD3Q19Descriptor = D3Q19<tag::MRT>;

using ForcedMRTD3Q19Descriptor = D3Q19<tag::MRT,FORCE>;


namespace mrt_data {

using utilities::Fraction;

// Matrix of base change between f and moments : moments=M.f
template <unsigned D, unsigned Q>
constexpr Fraction M[Q][Q] = {};

// inverse of base change matrix : f=invM.moments
template <unsigned D, unsigned Q>
constexpr Fraction invM[Q][Q] = {};

// relaxation times
template <unsigned D, unsigned Q>
constexpr Fraction s[Q] = {};

// relaxation times
template <unsigned D, unsigned Q>
constexpr Fraction s_2[Q] = {};

template <unsigned D, unsigned Q>
constexpr int shearIndexes = {};

// relevant indexes of r. t. for shear viscosity
template <unsigned D, unsigned Q>
constexpr int shearViscIndexes[shearIndexes<D,Q>] = {};

// relevant index of r. t. for bulk viscosity
template <unsigned D, unsigned Q>
constexpr int bulkViscIndex = {};

template <>
constexpr Fraction M<2,5>[5][5] = {
  { 1, 1, 1, 1, 1},
  { 0,-1, 0, 1, 0},
  { 0, 0,-1, 0, 1},
  {-4, 1, 1, 1, 1},
  { 0, 1,-1, 1,-1}
};

template <>
constexpr Fraction M<2,9>[9][9] = {
  { 1,  1,  1,  1,  1,  1,  1,  1,  1},
  {-4,  2, -1,  2, -1,  2, -1,  2, -1},
  { 4,  1, -2,  1, -2,  1, -2,  1, -2},
  { 0, -1, -1, -1,  0,  1,  1,  1,  0},
  { 0, -1,  2, -1,  0,  1, -2,  1,  0},
  { 0,  1,  0, -1, -1, -1,  0,  1,  1},
  { 0,  1,  0, -1,  2, -1,  0,  1, -2},
  { 0,  0,  1,  0, -1,  0,  1,  0, -1},
  { 0, -1,  0,  1,  0, -1,  0,  1,  0}
};

template <>
constexpr Fraction M<3,7>[7][7] = {
  //  Li, Yang et al 2016: The directions are modified for the OpenLB definition
  {1,  1,  1,  1,  1,  1,  1},
  {0,  1,  0, -1,  0,  0,  0},
  {0,  0, -1,  0,  0,  1,  0},
  {0,  0,  0,  0,  1,  0, -1},
  {6, -1, -1, -1, -1, -1, -1},
  {0,  2, -1,  2, -1, -1, -1},
  {0,  0,  1,  0, -1,  1, -1}
};

template <>
constexpr Fraction M<3,19>[19][19] = {
  {  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
  {-30,-11,-11,-11,  8,  8,  8,  8,  8,  8,-11,-11,-11,  8,  8,  8,  8,  8,  8},
  { 12, -4, -4, -4,  1,  1,  1,  1,  1,  1, -4, -4, -4,  1,  1,  1,  1,  1,  1},
  {  0, -1,  0,  0, -1, -1, -1, -1,  0,  0,  1,  0,  0,  1,  1,  1,  1,  0,  0},
  {  0,  4,  0,  0, -1, -1, -1, -1,  0,  0, -4,  0,  0,  1,  1,  1,  1,  0,  0},
  {  0,  0, -1,  0, -1,  1,  0,  0, -1, -1,  0,  1,  0,  1, -1,  0,  0,  1,  1},
  {  0,  0,  4,  0, -1,  1,  0,  0, -1, -1,  0, -4,  0,  1, -1,  0,  0,  1,  1},
  {  0,  0,  0, -1,  0,  0, -1,  1, -1,  1,  0,  0,  1,  0,  0,  1, -1,  1, -1},
  {  0,  0,  0,  4,  0,  0, -1,  1, -1,  1,  0,  0, -4,  0,  0,  1, -1,  1, -1},
  {  0,  2, -1, -1,  1,  1,  1,  1, -2, -2,  2, -1, -1,  1,  1,  1,  1, -2, -2},
  {  0, -4,  2,  2,  1,  1,  1,  1, -2, -2, -4,  2,  2,  1,  1,  1,  1, -2, -2},
  {  0,  0,  1, -1,  1,  1, -1, -1,  0,  0,  0,  1, -1,  1,  1, -1, -1,  0,  0},
  {  0,  0, -2,  2,  1,  1, -1, -1,  0,  0,  0, -2,  2,  1,  1, -1, -1,  0,  0},
  {  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0},
  {  0,  0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  1, -1},
  {  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0},
  {  0,  0,  0,  0, -1, -1,  1,  1,  0,  0,  0,  0,  0,  1,  1, -1, -1,  0,  0},
  {  0,  0,  0,  0,  1, -1,  0,  0, -1, -1,  0,  0,  0, -1,  1,  0,  0,  1,  1},
  {  0,  0,  0,  0,  0,  0, -1,  1,  1, -1,  0,  0,  0,  0,  0,  1, -1, -1,  1}
};

template <>
constexpr Fraction invM<2,5>[5][5] = {
  {{1, 5},       0,       0, {-1,  5},       0},
  {{1, 5}, {-1, 2},       0, { 1, 20}, { 1, 4}},
  {{1, 5},       0, {-1, 2}, { 1, 20}, {-1, 4}},
  {{1, 5}, { 1, 2},       0, { 1, 20}, { 1, 4}},
  {{1, 5},       0, { 1, 2}, { 1, 20}, {-1, 4}}
};

template <>
constexpr Fraction invM<2,9>[9][9] = {
  {{1, 9}, {-1,  9}, { 1,  9},       0,        0,       0,        0,       0,       0},
  {{1, 9}, { 1, 18}, { 1, 36}, {-1, 6}, {-1, 12}, { 1, 6}, { 1, 12},       0, {-1, 4}},
  {{1, 9}, {-1, 36}, {-1, 18}, {-1, 6}, { 1,  6},       0,        0, { 1, 4},       0},
  {{1, 9}, { 1, 18}, { 1, 36}, {-1, 6}, {-1, 12}, {-1, 6}, {-1, 12},       0, { 1, 4}},
  {{1, 9}, {-1, 36}, {-1, 18},       0,        0, {-1, 6}, { 1,  6}, {-1, 4},       0},
  {{1, 9}, { 1, 18}, { 1, 36}, { 1, 6}, { 1, 12}, {-1, 6}, {-1, 12},       0, {-1, 4}},
  {{1, 9}, {-1, 36}, {-1, 18}, { 1, 6}, {-1,  6},       0,        0, { 1, 4},       0},
  {{1, 9}, { 1, 18}, { 1, 36}, { 1, 6}, { 1, 12}, { 1, 6}, { 1, 12},       0, { 1, 4}},
  {{1, 9}, {-1, 36}, {-1, 18},       0,        0, { 1, 6}, {-1,  6}, {-1, 4},       0}
};

template <>
constexpr Fraction invM<3,7>[7][7] = {
  //  Li, Yang et al 2016: The directions are modified for the OpenLB definition
  {{1, 7},       0,       0,       0, { 1,  7},        0,       0},
  {{1, 7}, { 1, 2},       0,       0, {-1, 42}, { 1,  6},       0},
  {{1, 7},       0, {-1, 2},       0, {-1, 42}, {-1, 12}, { 1, 4}},
  {{1, 7}, {-1, 2},       0,       0, {-1, 42}, { 1,  6},       0},
  {{1, 7},       0,       0, { 1, 2}, {-1, 42}, {-1, 12}, {-1, 4}},
  {{1, 7},       0, { 1, 2},       0, {-1, 42}, {-1, 12}, { 1, 4}},
  {{1, 7},       0,       0, {-1, 2}, {-1, 42}, {-1, 12}, {-1, 4}}
};

template <>
constexpr Fraction invM<3,19>[19][19] = {
  {{1,19}, { -5,  399}, { 1,  21},        0,        0,        0,        0,
                     0,         0,        0,        0,        0,        0,
                     0,         0,        0,        0,        0,        0},/*0*/
  {{1,19}, {-11, 2394}, {-1,  63}, {-1, 10}, { 1, 10},        0,        0,
                     0,         0, { 1, 18}, {-1, 18},        0,        0,
                     0,         0,        0,        0,        0,        0},/*1*/
  {{1,19}, {-11, 2394}, {-1,  63},        0,        0, {-1, 10}, { 1, 10},
                     0,         0, {-1, 36}, { 1, 36}, { 1, 12}, {-1, 12},
                     0,         0,        0,        0,        0,        0},/*2*/
  {{1,19}, {-11, 2394}, {-1,  63},        0,        0,        0,        0,
           { -1,   10}, { 1,  10}, {-1, 36}, { 1, 36}, {-1, 12}, { 1, 12},
                     0,         0,        0,        0,        0,        0},/*3*/
  {{1,19}, {  4, 1197}, { 1, 252}, {-1, 10}, {-1, 40}, {-1, 10}, {-1, 40},
                     0,         0, { 1, 36}, { 1, 72}, { 1, 12}, { 1, 24},
           {  1,    4},         0,        0, {-1,  8}, { 1,  8},        0},/*4*/
  {{1,19}, {  4, 1197}, { 1, 252}, {-1, 10}, {-1, 40}, { 1, 10}, { 1, 40},
                     0,         0, { 1, 36}, { 1, 72}, { 1, 12}, { 1, 24},
           { -1,    4},         0,        0, {-1,  8}, {-1,  8},        0},/*5*/
  {{1,19}, {  4, 1197}, { 1, 252}, {-1, 10}, {-1, 40},        0,        0,
           { -1,   10}, {-1,  40}, { 1, 36}, { 1, 72}, {-1, 12}, {-1, 24},
                     0,         0, { 1,  4}, { 1,  8},        0, {-1,  8}},/*6*/
  {{1,19}, {  4, 1197}, { 1, 252}, {-1, 10}, {-1, 40},        0,        0,
           {  1,   10}, { 1,  40}, { 1, 36}, { 1, 72}, {-1, 12}, {-1, 24},
                     0,         0, {-1,  4}, { 1,  8},        0, { 1,  8}},/*7*/
  {{1,19}, {  4, 1197}, { 1, 252},        0,        0, {-1, 10}, {-1, 40},
           { -1,   10}, {-1,  40}, {-1, 18}, {-1, 36},        0,        0,
                     0, { 1,   4},        0,        0, {-1,  8}, { 1,  8}},/*8*/
  {{1,19}, {  4, 1197}, { 1, 252},        0,        0, {-1, 10}, {-1, 40},
           {  1,   10}, { 1,  40}, {-1, 18}, {-1, 36},        0,        0,
                     0, {-1,   4},        0,        0, {-1,  8}, {-1,  8}},/*9*/
  {{1,19}, {-11, 2394}, {-1,  63}, { 1, 10}, {-1, 10},        0,        0,
                     0,         0, { 1, 18}, {-1, 18},        0,        0,
                     0,         0,        0,        0,        0,        0},/*10*/
  {{1,19}, {-11, 2394}, {-1,  63},        0,        0, { 1, 10}, {-1, 10},
                     0,         0, {-1, 36}, { 1, 36}, { 1, 12}, {-1, 12},
                     0,         0,        0,        0,        0,        0},/*11*/
  {{1,19}, {-11, 2394}, {-1,  63},        0,        0,        0,        0,
           {  1,   10}, {-1,  10}, {-1, 36}, { 1, 36}, {-1, 12}, { 1, 12},
                     0,         0,        0,        0,        0,        0},/*12*/
  {{1,19}, {  4, 1197}, { 1, 252}, { 1, 10}, { 1, 40}, { 1, 10}, { 1, 40},
                     0,         0, { 1, 36}, { 1, 72}, { 1, 12}, { 1, 24},
           {  1,    4},         0,        0, { 1,  8}, {-1,  8},        0},/*13*/
  {{1,19}, {  4, 1197}, { 1, 252}, { 1, 10}, { 1, 40}, {-1, 10}, {-1, 40},
                     0,         0, { 1, 36}, { 1, 72}, { 1, 12}, { 1, 24},
           { -1,    4},         0,        0, { 1,  8}, { 1,  8},        0},/*14*/
  {{1,19}, {  4, 1197}, { 1, 252}, { 1, 10}, { 1, 40},        0,        0,
           {  1,   10}, { 1,  40}, { 1, 36}, { 1, 72}, {-1, 12}, {-1, 24},
                     0,         0, { 1,  4}, {-1,  8},        0, { 1,  8}},/*15*/
  {{1,19}, {  4, 1197}, { 1, 252}, { 1, 10}, { 1, 40},        0,        0,
           { -1,   10}, {-1,  40}, { 1, 36}, { 1, 72}, {-1, 12}, {-1, 24},
                     0,         0, {-1,  4}, {-1,  8},        0, {-1,  8}},/*16*/
  {{1,19}, {  4, 1197}, { 1, 252},        0,        0, { 1, 10}, { 1, 40},
           {  1,   10}, { 1,  40}, {-1, 18}, {-1, 36},        0,        0,
                     0, { 1,   4},        0,        0, { 1,  8}, {-1,  8}},/*17*/
  {{1,19}, {  4, 1197}, { 1, 252},        0,        0, { 1, 10}, { 1, 40},
           { -1,   10}, {-1,  40}, {-1, 18}, {-1, 36},        0,        0,
                     0, {-1,   4},        0,        0, { 1,  8}, { 1,  8}}/*18*/
};

template <>
constexpr Fraction s<2,5>[5] = {
  0, 0, 0, {3, 2}, {3, 2}
};

// s7=s8 to have a shear viscosity nu
// and the bulk viscosity depends on s2.
template <>
constexpr Fraction s<2,9>[9] = {
  0, {11, 10}, {11, 10}, 0, {11, 10}, 0, {11, 10}, 0, 0
};

template <>
constexpr Fraction s<3,7>[7] = {
  // Original MRT Relaxation times
  /*s0*/  0,  // rho (conserved)
  /*s1*/  0,  // Function of the thermal diffusivity: S_a = 1/t_a = 1/(4*a + 1/2)
  /*s2*/  0,  // Function of the thermal diffusivity: S_a = 1/t_a = 1/(4*a + 1/2)
  /*s3*/  0,  // Function of the thermal diffusivity: S_a = 1/t_a = 1/(4*a + 1/2)
  /*s4*/  {19, 10},
  /*s5*/  {19, 10},
  /*s6*/  {19, 10}
};

// Original MRT Relaxation times
template <>
constexpr Fraction s<3,19>[19] = {
  /*s0*/            0, // rho (conserved)
  /*s1*/  { 119, 100},
  /*s2*/  {   7,   5},
  /*s3*/            0, // rho*ux (conserved)
  /*s4*/  {   6,   5},
  /*s5*/            0, // rho*uy (conserved)
  /*s6*/  {   6,   5}, // = s4
  /*s7*/            0, // rho*uz (conserved)
  /*s8*/  {   6,   5}, // = s4
  /*s9*/            0, //should be equal to s13, used to define nu
  /*s10*/ {   7,   5},
  /*s11*/           0, // = s9,
  /*s12*/ {   7,   5},
  /*s13*/           0, //should be equal to s9, used to define nu
  /*s14*/           0, // = s13,
  /*s15*/           0, // = s13,
  /*s16*/ {  99,  50},
  /*s17*/ {  99,  50}, // = s16,
  /*s18*/ {  99,  50} // = s16,
};

//  Use these relaxation time for higher stability
template <>
constexpr Fraction s_2<3,19>[19] = {
  /*s0*/  0, // rho (conserved)
  /*s1*/  1,
  /*s2*/  1,
  /*s3*/  0, // rho*ux (conserved)
  /*s4*/  1,
  /*s5*/  0, // rho*uy (conserved)
  /*s6*/  1, // = s4
  /*s7*/  0, // rho*uz (conserved)
  /*s8*/  1, // = s4
  /*s9*/  0, //should be equal to s13, used to define nu
  /*s10*/ 1,
  /*s11*/ 0, // = s9,
  /*s12*/ 1,
  /*s13*/ 0, //should be equal to s9, used to define nu
  /*s14*/ 0, // = s13,
  /*s15*/ 0, // = s13,
  /*s16*/ 1,
  /*s17*/ 1, // = s16,
  /*s18*/ 1  // = s16,
};

template <>
constexpr int shearIndexes<2,5> = 2;

template <>
constexpr int shearIndexes<2,9> = 2;

template <>
constexpr int shearIndexes<3,7> = 3;

template <>
constexpr int shearIndexes<3,19> = 5;

template <>
constexpr int shearViscIndexes<2,5>[shearIndexes<2,5>] = { 1, 2 };

template <>
constexpr int shearViscIndexes<2,9>[shearIndexes<2,9>] = { 7, 8 };

template <>
constexpr int shearViscIndexes<3,7>[shearIndexes<3,7>] = { 1, 2, 3 };

template <>
constexpr int shearViscIndexes<3,19>[shearIndexes<3,19>] = { 9, 11, 13, 14, 15 };

template <>
constexpr int bulkViscIndex<2,9> = 2;

template <>
constexpr int bulkViscIndex<3,19> = 1;

} // mrt_data

template <typename T, unsigned D, unsigned Q>
constexpr T t(unsigned iPop, tag::MRT)
{
  return data::t<D,Q>[iPop].template as<T>();
}

template <typename T, unsigned D, unsigned Q>
constexpr T m(unsigned iPop, unsigned jPop, tag::MRT)
{
  return mrt_data::M<D,Q>[iPop][jPop].template as<T>();
}

template <typename T, typename DESCRIPTOR>
constexpr T m(unsigned iPop, unsigned jPop)
{
  return m<T, DESCRIPTOR::d, DESCRIPTOR::q>(iPop, jPop, typename DESCRIPTOR::category_tag());
}

template <typename T, unsigned D, unsigned Q>
constexpr T invM(unsigned iPop, unsigned jPop, tag::MRT)
{
  return mrt_data::invM<D,Q>[iPop][jPop].template as<T>();
}

template <typename T, typename DESCRIPTOR>
constexpr T invM(unsigned iPop, unsigned jPop)
{
  return invM<T, DESCRIPTOR::d, DESCRIPTOR::q>(iPop, jPop, typename DESCRIPTOR::category_tag());
}

template <typename T, unsigned D, unsigned Q>
constexpr T s(unsigned iPop, tag::MRT)
{
  return mrt_data::s<D,Q>[iPop].template as<T>();
}

template <typename T, typename DESCRIPTOR>
constexpr T s(unsigned iPop)
{
  return s<T, DESCRIPTOR::d, DESCRIPTOR::q>(iPop, typename DESCRIPTOR::category_tag());
}

template <typename T, unsigned D, unsigned Q>
constexpr T s_2(unsigned iPop, tag::MRT)
{
  return mrt_data::s_2<D,Q>[iPop].template as<T>();
}

template <typename T, typename DESCRIPTOR>
constexpr T s_2(unsigned iPop)
{
  return s_2<T, DESCRIPTOR::d, DESCRIPTOR::q>(iPop, typename DESCRIPTOR::category_tag());
}

template <unsigned D, unsigned Q>
constexpr int shearIndexes(tag::MRT)
{
  return mrt_data::shearIndexes<D,Q>;
}

template <typename DESCRIPTOR>
constexpr int shearIndexes()
{
  return shearIndexes<DESCRIPTOR::d, DESCRIPTOR::q>(typename DESCRIPTOR::category_tag());
}

template <unsigned D, unsigned Q>
constexpr int shearViscIndexes(unsigned iPop, tag::MRT)
{
  return mrt_data::shearViscIndexes<D,Q>[iPop];
}

template <typename DESCRIPTOR>
constexpr int shearViscIndexes(unsigned iPop)
{
  return shearViscIndexes<DESCRIPTOR::d, DESCRIPTOR::q>(iPop, typename DESCRIPTOR::category_tag());
}

template <unsigned D, unsigned Q>
constexpr int bulkViscIndex(tag::MRT)
{
  return mrt_data::bulkViscIndex<D,Q>;
}

template <typename DESCRIPTOR>
constexpr int bulkViscIndex()
{
  return bulkViscIndex<DESCRIPTOR::d, DESCRIPTOR::q>(typename DESCRIPTOR::category_tag());
}

}  // namespace descriptors

}  // namespace olb

#endif
