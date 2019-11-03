/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause
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

#ifndef FRINGE_2D_H
#define FRINGE_2D_H

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm

#include "analyticalF.h"
#include "io/ostreamManager.h"

/** A fringe zone is a spatial domain, which is included into the
 *  computational domain to aim a transition from turbulent into
 * laminar flow. It is often applied to improve the use of
 * artificial boundary conditions.
 *
 *                 lambda_max
 *               _____________
 *              /             \
 *             /               \
 *            /                 \
 *           /                   \
 *          /                     \
 * ________/                       \________
 *         |x_start                | x_ende
 *         |-----|           |-----|
 *             b_rise         b_fall
 *
 * fringe function:
 *   lambda(x) = lambda_max*[S((x-x_start)/(b_rise)) - S((x-x_end)/(b_fall)+1)]
 * The fringe function sets the strength of the damping. The
 * maximum strength is lambda_max and the shape of the function is defined by
 * the stepfunction S(x) and the parameters b_rise and b_fall.
 *
 * lambda_max: maximal damping force
 * x_start: begin of the fringe zone
 * x_end: end of the fringe zone
 * b_rise: rise distance
 * b_fall: fall distance
 *
 * S is a smooth step function:
 * S(x)=0,     for x<=0
 * S(x)=1/( 1 + exp( (1/(x-1)) + (1/x) ) ),  for 0<x<1,
 * S(x)=1,     for x>=1.
 *
 * --> G = lambda*(U - u) is the volume force, which is added to
 * the flow equations in order to transfer the actual velocity into the wanted
 * velocity field
 *
 * u: actual velocity
 * U: wanted velocity - either average profil or poiseuille profile
 *
 * BibteX listing of the main paper:
 *
 * @TECHREPORT{lundbladh:99,
 *   author = {Anders Lundbladh and Stellan Berlin and Martin
 *             Skote and Casper Hildings and Jaisig Choi and
 *             John Kim and Dan Henningson},
 *   title = {An Efficient Spectral Method for Simulation of
 *            Incompressible Flow Over a Flat Plate},
 *   institution = {not set},
 *   url = {http://www.fluidosol.se/thesismod/paper9.pdf},
 *   year = {1999} }
 */


namespace olb {

template <typename T, typename S>
class Fringe2D : public AnalyticalF2D<T,S> {

private:
  /// Wanted lattice velocity (in lattice units)
  AnalyticalF2D<T,S>& _wantedVelocity;
  /// Start of the fringe zone in _direction in physR (Si units)
  T _start;
  /// End of the fringe zone in _direction in physR (Si units)
  T _end;
  /// Direction of the fringe zone (parallel to x-axis means _direction = 0, parallel to y-axis means _direction = 1)
  int _direction;
  /// Maximum damping (TODO physical? depent on dt?)
  T _lambdaMax;
  /// Percentage of the length of the fringezone at the beginning, aim: softer increase of lambda
  T _rise;
  /// Percentage of the length of the fringezone at the end
  T _fall;

public:
  Fringe2D(AnalyticalF2D<T,S>& wantedVelocity, T start, T end, int direction, T lambdaMax=.005, T rise=.4, T fall=.2);

  /// Returns coefficients (a1, a2, B11, B12, B21, B22)
  /// for a linear velocity force modell where F = a + B*lattticeVelocity
  /// a = -scaleLambda*wanted, B = (scaleLamda 0 0 scaleLamda), so that
  /// F = scaleLambda*(latticeVelocity - wantedVelocity)
  bool operator()(T output[6], const S input[2]);

private:
  /// Step function
  T s(const T& x) const;
  /// Damping function
  T lambda(const T& x) const;
};

} // end namespace olb

#endif


