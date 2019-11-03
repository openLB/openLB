/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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
 * Efficient 3D specializations for lbHelpers3D.h
 */
#ifndef FIRST_ORDER_LB_HELPERS_3D_H
#define FIRST_ORDER_LB_HELPERS_3D_H

namespace olb {
/*

/// Compute Pi tensor efficiently on D3Q15 lattice
namespace neqPiD3Q15 {

using namespace olb::util::tensorIndices3D;

template<typename T>
T fromPiToFneq0(const T pi[6]) {
    return -(T)1/(T)3*pi[xx] - (T)1/(T)3*pi[yy] - (T)1/(T)3*pi[zz];
}

template<typename T>
T fromPiToFneq1(const T pi[6]) {
    return (T)1/(T)2 * (
              (T)2/(T)3*pi[xx] - (T)1/(T)3*pi[yy] - (T)1/(T)3*pi[zz]
           );
}

template<typename T>
T fromPiToFneq2(const T pi[6]) {
    return (T)1/(T)2 * (
             -(T)1/(T)3*pi[xx] + (T)2/(T)3*pi[yy] - (T)1/(T)3*pi[zz]
           );
}

template<typename T>
T fromPiToFneq3(const T pi[6]) {
    return (T)1/(T)2 * (
             -(T)1/(T)3*pi[xx] - (T)1/(T)3*pi[yy] + (T)2/(T)3*pi[zz]
           );
}

template<typename T>
T fromPiToFneq4(const T pi[6]) {
    return (T)1/(T)16 * (
              (T)2/(T)3*pi[xx] + (T)2/(T)3*pi[yy] + (T)2/(T)3*pi[zz]
              + pi[xy] + pi[xz] + pi[yz]
           );
}

template<typename T>
T fromPiToFneq5(const T pi[6]) {
    return (T)1/(T)16 * (
              (T)2/(T)3*pi[xx] + (T)2/(T)3*pi[yy] + (T)2/(T)3*pi[zz]
              + pi[xy] - pi[xz] - pi[yz]
           );
}

template<typename T>
T fromPiToFneq6(const T pi[6]) {
    return (T)1/(T)16 * (
              (T)2/(T)3*pi[xx] + (T)2/(T)3*pi[yy] + (T)2/(T)3*pi[zz]
              - pi[xy] + pi[xz] - pi[yz]
           );
}

template<typename T>
T fromPiToFneq7(const T pi[6]) {
    return (T)1/(T)16 * (
              (T)2/(T)3*pi[xx] + (T)2/(T)3*pi[yy] + (T)2/(T)3*pi[zz]
              - pi[xy] - pi[xz] + pi[yz]
           );
}

}  // namespace neqPiD3Q15


/// Compute RLB collision efficiently on D3Q15 lattice
template<typename T>
struct rlbHelpers<T, descriptors::D3Q15<>> {

static T rlbCollision (
    Cell<T,descriptors::D3Q15<>>& cell,
    T rho, const T u[3], const T pi[6], T omega )
{
    typedef lbHelpers<T, descriptors::D3Q15<>> LH;
    const T uSqr = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];

    cell[0]  = LH::equilibrium(0, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q15::fromPiToFneq0<T>(pi);
    cell[1]  = LH::equilibrium(1, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q15::fromPiToFneq1<T>(pi);
    cell[2]  = LH::equilibrium(2, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q15::fromPiToFneq2<T>(pi);
    cell[3]  = LH::equilibrium(3, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q15::fromPiToFneq3<T>(pi);
    cell[4]  = LH::equilibrium(4, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q15::fromPiToFneq4<T>(pi);
    cell[5]  = LH::equilibrium(5, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q15::fromPiToFneq5<T>(pi);
    cell[6]  = LH::equilibrium(6, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q15::fromPiToFneq6<T>(pi);
    cell[7]  = LH::equilibrium(7, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q15::fromPiToFneq7<T>(pi);
    cell[8]  = LH::equilibrium(8, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q15::fromPiToFneq1<T>(pi);
    cell[9]  = LH::equilibrium(9, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q15::fromPiToFneq2<T>(pi);
    cell[10]  = LH::equilibrium(10, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q15::fromPiToFneq3<T>(pi);
    cell[11]  = LH::equilibrium(11, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q15::fromPiToFneq4<T>(pi);
    cell[12]  = LH::equilibrium(12, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q15::fromPiToFneq5<T>(pi);
    cell[13]  = LH::equilibrium(13, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q15::fromPiToFneq6<T>(pi);
    cell[14]  = LH::equilibrium(14, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q15::fromPiToFneq7<T>(pi);
    return uSqr;
}

};  // struct rlbHelpers<D3Q15<>>


namespace neqPiD3Q19 {

using namespace olb::util::tensorIndices3D;

template<typename T>
T fromPiToFneq0(const T pi[6]) {
    return (T)3/(T)2 * (
             -(T)1/(T)3*pi[xx] - (T)1/(T)3*pi[yy] - (T)1/(T)3*pi[zz]
           );
}

template<typename T>
T fromPiToFneq1(const T pi[6]) {
    return (T)1/(T)4 * (
              (T)2/(T)3*pi[xx] - (T)1/(T)3*pi[yy] - (T)1/(T)3*pi[zz]
           );
}

template<typename T>
T fromPiToFneq2(const T pi[6]) {
    return (T)1/(T)4 * (
             -(T)1/(T)3*pi[xx] + (T)2/(T)3*pi[yy] - (T)1/(T)3*pi[zz]
           );
}

template<typename T>
T fromPiToFneq3(const T pi[6]) {
    return (T)1/(T)4 * (
             -(T)1/(T)3*pi[xx] - (T)1/(T)3*pi[yy] + (T)2/(T)3*pi[zz]
           );
}

template<typename T>
T fromPiToFneq4(const T pi[6]) {
    return (T)1/(T)8 * (
              (T)2/(T)3*pi[xx] + (T)2/(T)3*pi[yy] - (T)1/(T)3*pi[zz]
              + (T)2*pi[xy]
           );
}

template<typename T>
T fromPiToFneq5(const T pi[6]) {
    return (T)1/(T)8 * (
              (T)2/(T)3*pi[xx] + (T)2/(T)3*pi[yy] - (T)1/(T)3*pi[zz]
              - (T)2*pi[xy]
           );
}

template<typename T>
T fromPiToFneq6(const T pi[6]) {
    return (T)1/(T)8 * (
              (T)2/(T)3*pi[xx] - (T)1/(T)3*pi[yy] + (T)2/(T)3*pi[zz]
              + (T)2*pi[xz]
           );
}

template<typename T>
T fromPiToFneq7(const T pi[6]) {
    return (T)1/(T)8 * (
              (T)2/(T)3*pi[xx] - (T)1/(T)3*pi[yy] + (T)2/(T)3*pi[zz]
              - (T)2*pi[xz]
           );
}

template<typename T>
T fromPiToFneq8(const T pi[6]) {
    return (T)1/(T)8 * (
             -(T)1/(T)3*pi[xx] + (T)2/(T)3*pi[yy] + (T)2/(T)3*pi[zz]
              + (T)2*pi[yz]
           );
}

template<typename T>
T fromPiToFneq9(const T pi[6]) {
    return (T)1/(T)8 * (
             -(T)1/(T)3*pi[xx] + (T)2/(T)3*pi[yy] + (T)2/(T)3*pi[zz]
              - (T)2*pi[yz]
           );
}


}  // namespace neqPiD3Q19

template<typename T>
struct rlbHelpers<T, descriptors::D3Q19<>> {
static T rlbCollision (
    Cell<T,descriptors::D3Q19<>>& cell,
    T rho, const T u[3], const T pi[6], T omega )
{
    typedef lbHelpers<T, descriptors::D3Q19<>> LH;
    const T uSqr = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];

    cell[0]  = LH::equilibrium(0, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq0<T>(pi);
    cell[1]  = LH::equilibrium(1, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq1<T>(pi);
    cell[2]  = LH::equilibrium(2, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq2<T>(pi);
    cell[3]  = LH::equilibrium(3, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq3<T>(pi);
    cell[4]  = LH::equilibrium(4, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq4<T>(pi);
    cell[5]  = LH::equilibrium(5, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq5<T>(pi);
    cell[6]  = LH::equilibrium(6, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq6<T>(pi);
    cell[7]  = LH::equilibrium(7, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq7<T>(pi);
    cell[8]  = LH::equilibrium(8, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq8<T>(pi);
    cell[9]  = LH::equilibrium(9, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq9<T>(pi);
    cell[10]  = LH::equilibrium(10, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq1<T>(pi);
    cell[11]  = LH::equilibrium(11, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq2<T>(pi);
    cell[12]  = LH::equilibrium(12, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq3<T>(pi);
    cell[13]  = LH::equilibrium(13, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq4<T>(pi);
    cell[14]  = LH::equilibrium(14, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq5<T>(pi);
    cell[15]  = LH::equilibrium(15, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq6<T>(pi);
    cell[16]  = LH::equilibrium(16, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq7<T>(pi);
    cell[17]  = LH::equilibrium(17, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq8<T>(pi);
    cell[18]  = LH::equilibrium(18, rho, u, uSqr) +
                   ((T)1-omega) * neqPiD3Q19::fromPiToFneq9<T>(pi);
    return uSqr;
}

};  // struct rlbHelpers<D3Q19<>>
*/

}  // namespace olb


#endif
