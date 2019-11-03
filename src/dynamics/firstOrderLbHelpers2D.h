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
 * Efficient 2D specializations for lbHelpersD2Q9.h
 */
#ifndef FIRST_ORDER_LB_HELPERS_2D_H
#define FIRST_ORDER_LB_HELPERS_2D_H

namespace olb {

/// Compute Pi tensor efficiently on D2Q9 lattice
namespace neqPiD2Q9 {

using namespace olb::util::tensorIndices2D;

template<typename T>
T fromPiToFneq0(const T pi[3])
{
  return (T)2 * (-(T)1/(T)3*pi[xx] - (T)1/(T)3*pi[yy]);
}

template<typename T>
T fromPiToFneq1(const T pi[3])
{
  return (T)1/(T)4 * ((T)1/(T)3*pi[xx] + (T)1/(T)3*pi[yy] - pi[xy]);
}

template<typename T>
T fromPiToFneq2(const T pi[3])
{
  return (T)1/(T)2 * ((T)2/(T)3*pi[xx] - (T)1/(T)3*pi[yy]);
}

template<typename T>
T fromPiToFneq3(const T pi[3])
{
  return (T)1/(T)4 * ((T)1/(T)3*pi[xx] + (T)1/(T)3*pi[yy] + pi[xy]);
}

template<typename T>
T fromPiToFneq4(const T pi[3])
{
  return (T)1/(T)2 * (-(T)1/(T)3*pi[xx] + (T)2/(T)3*pi[yy]);
}

}  // namespace neqPiD2Q9

template<typename T>
struct rlbHelpers<T, descriptors::D2Q9<>> {

  static T rlbCollision (
    Cell<T,descriptors::D2Q9<>>& cell,
    T rho, const T u[2], const T pi[3], T omega )
  {
    typedef lbHelpers<T, descriptors::D2Q9<>> LH;
    const T uSqr = u[0]*u[0] + u[1]*u[1];

    cell[0]  = LH::equilibrium(0, rho, u, uSqr) +
               ((T)1-omega) * neqPiD2Q9::fromPiToFneq0(pi);
    cell[1]  = LH::equilibrium(1, rho, u, uSqr) +
               ((T)1-omega) * neqPiD2Q9::fromPiToFneq1(pi);
    cell[2]  = LH::equilibrium(2, rho, u, uSqr) +
               ((T)1-omega) * neqPiD2Q9::fromPiToFneq2(pi);
    cell[3]  = LH::equilibrium(3, rho, u, uSqr) +
               ((T)1-omega) * neqPiD2Q9::fromPiToFneq3(pi);
    cell[4]  = LH::equilibrium(4, rho, u, uSqr) +
               ((T)1-omega) * neqPiD2Q9::fromPiToFneq4(pi);
    cell[5]  = LH::equilibrium(5, rho, u, uSqr) +
               ((T)1-omega) * neqPiD2Q9::fromPiToFneq1(pi);
    cell[6]  = LH::equilibrium(6, rho, u, uSqr) +
               ((T)1-omega) * neqPiD2Q9::fromPiToFneq2(pi);
    cell[7]  = LH::equilibrium(7, rho, u, uSqr) +
               ((T)1-omega) * neqPiD2Q9::fromPiToFneq3(pi);
    cell[8]  = LH::equilibrium(8, rho, u, uSqr) +
               ((T)1-omega) * neqPiD2Q9::fromPiToFneq4(pi);
    return uSqr;
  }

};  // struct rlbHelpers<D2Q9<>>

}  // namespace olb


#endif
