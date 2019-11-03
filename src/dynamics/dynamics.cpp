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
 * can be instantiated -- template instantiation.
 */
#include "dynamics.h"
#include "dynamics.hh"
#include "latticeDescriptors.h"
 

namespace olb {

template struct Dynamics<double, descriptors::D2Q9<>>;
template struct Momenta<double, descriptors::D2Q9<>>;
template class BasicDynamics<double, descriptors::D2Q9<>>;
template class BGKdynamics<double, descriptors::D2Q9<>>;
template class ConstRhoBGKdynamics<double, descriptors::D2Q9<>>;
template class IncBGKdynamics<double, descriptors::D2Q9<>>;
template class RLBdynamics<double, descriptors::D2Q9<>>;
template class CombinedRLBdynamics<double, descriptors::D2Q9<>,
                                   RLBdynamics<double, descriptors::D2Q9<>> >;
template class CombinedRLBdynamics<double, descriptors::D2Q9<>,
                                   BGKdynamics<double, descriptors::D2Q9<>> >;
template class CombinedRLBdynamics<double, descriptors::D2Q9<>,
                                   ConstRhoBGKdynamics<double, descriptors::D2Q9<>> >;
template struct BulkMomenta<double, descriptors::D2Q9<>>;
template class BounceBack<double, descriptors::D2Q9<>>;
template class BounceBackVelocity<double, descriptors::D2Q9<>>;
template class BounceBackAnti<double, descriptors::D2Q9<>>;
template class NoDynamics<double, descriptors::D2Q9<>>;
template class OffDynamics<double, descriptors::D2Q9<>>;
template class ZeroDistributionDynamics<double, descriptors::D2Q9<>>;

template struct Dynamics<double, descriptors::D3Q19<>>;
template struct Momenta<double, descriptors::D3Q19<>>;
template class BasicDynamics<double, descriptors::D3Q19<>>;
template class BGKdynamics<double, descriptors::D3Q19<>>;
template class ConstRhoBGKdynamics<double, descriptors::D3Q19<>>;
template class IncBGKdynamics<double, descriptors::D3Q19<>>;
template class RLBdynamics<double, descriptors::D3Q19<>>;
template class CombinedRLBdynamics<double, descriptors::D3Q19<>,
                                   RLBdynamics<double, descriptors::D3Q19<>> >;
template class CombinedRLBdynamics<double, descriptors::D3Q19<>,
                                   BGKdynamics<double, descriptors::D3Q19<>> >;
template class CombinedRLBdynamics<double, descriptors::D3Q19<>,
                                   ConstRhoBGKdynamics<double, descriptors::D3Q19<>> >;
template struct BulkMomenta<double, descriptors::D3Q19<>>;
template class BounceBack<double, descriptors::D3Q19<>>;
template class BounceBackVelocity<double, descriptors::D3Q19<>>;
template class BounceBackAnti<double, descriptors::D3Q19<>>;
template class NoDynamics<double, descriptors::D3Q19<>>;
template class OffDynamics<double, descriptors::D3Q19<>>;
template class ZeroDistributionDynamics<double, descriptors::D3Q19<>>;


namespace instances {

template BulkMomenta<double, descriptors::D2Q9<>>& getBulkMomenta();

template BounceBack<double, descriptors::D2Q9<>>& getBounceBack();

template BounceBackAnti<double, descriptors::D2Q9<>>& getBounceBackAnti(const double rho);

template BounceBackVelocity<double, descriptors::D2Q9<>>& getBounceBackVelocity(const double rho, const double u[2]);

template NoDynamics<double, descriptors::D2Q9<>>& getNoDynamics(double rho);

template ZeroDistributionDynamics<double, descriptors::D2Q9<>>& getZeroDistributionDynamics();

template BulkMomenta<double, descriptors::D3Q19<>>& getBulkMomenta();

template BounceBack<double, descriptors::D3Q19<>>& getBounceBack();

template BounceBackVelocity<double, descriptors::D3Q19<>>& getBounceBackVelocity(const double rho, const double u[3]);

template BounceBackAnti<double, descriptors::D3Q19<>>& getBounceBackAnti(const double rho);

template NoDynamics<double, descriptors::D3Q19<>>& getNoDynamics(double rho);

template ZeroDistributionDynamics<double, descriptors::D3Q19<>>& getZeroDistributionDynamics();
}

}
