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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */
#ifndef ENTROPIC_LB_HELPERS_H
#define ENTROPIC_LB_HELPERS_H

#include <cmath>

namespace olb {

template<typename T, typename DESCRIPTOR>
struct entropicLbHelpers {
  /// Computation of equilibrium distribution
  static T equilibrium( int iPop, T rho, const T u[DESCRIPTOR::d])
  {
    typedef DESCRIPTOR L;
    const T invCs = sqrt(descriptors::invCs2<T,L>());
    const T sqt3 = sqrt(3.0);
    T prod = (T)1;
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      T uc = u[iD] * invCs; // u[iD] / c_s

      prod *= ((T)2 - sqrt(1.0+uc*uc)) *
              pow((2.0 / sqt3 * uc +
                   sqrt(1.0+uc*uc))/(1.0-uc/sqt3),
                  descriptors::c<L>(iPop,iD)/sqt3*invCs);
    }
    return rho*descriptors::t<T,L>(iPop)*prod-descriptors::t<T,L>(iPop);
  }

  /// Computation of equilibrium distribution
  static T equilibriumApprox( int iPop, T rho, const T u[DESCRIPTOR::d])
  {
    typedef DESCRIPTOR L;

    T uSqr = util::normSqr<T,L::d>(u);
    T cu = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      cu += descriptors::c<L>(iPop,iD)*u[iD];
    }

    return rho * descriptors::t<T,L>(iPop) * (1.0 +
                               cu*descriptors::invCs2<T,L>() - 0.5 * uSqr*descriptors::invCs2<T,L>() + 0.5*pow(descriptors::invCs2<T,L>(),2)*cu*cu
                               - 0.5*pow(descriptors::invCs2<T,L>(),2)*cu*uSqr + pow(cu,3)*pow(descriptors::invCs2<T,L>(),3)/6.0
                               + 0.125*uSqr*uSqr*pow(descriptors::invCs2<T,L>(),2) - 0.25*cu*cu*uSqr*pow(descriptors::invCs2<T,L>(),3)
                               + pow(cu,4)*pow(descriptors::invCs2<T,L>(),4)/24.0)-descriptors::t<T,L>(iPop);
  }
};

}

#include "entropicLbHelpers2D.h"
#include "entropicLbHelpers3D.h"

#endif
