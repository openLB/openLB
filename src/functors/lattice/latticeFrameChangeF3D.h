/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2015 Gilles Zahnd, Mathias J. Krause
 *  Marie-Luise Maier
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

#ifndef LATTICE_FRAME_CHANGE_F_3D_H
#define LATTICE_FRAME_CHANGE_F_3D_H

#include <vector>
#include <string>

#include "functors/analytical/analyticalF.h"
#include "functors/analytical/frameChangeF2D.h"
#include "superBaseF3D.h"
#include "superLatticeLocalF3D.h"

/** To enable simulations in a rotating frame, the axis is set in the
  * constructor with axisPoint and axisDirection. The axisPoint can be the
  * coordinate of any point on the axis. The axisDirection has to be a normed to
  * 1. The pulse w is in rad/s. It determines the pulse speed by its norm while
  * the trigonometric or clockwise direction is determined by its sign: When the
  * axisDirection is pointing "towards you", a positive pulse makes it turn in
  * the trigonometric way. It has to be noticed that putting both axisDirection
  * into -axisDirection and w into -w yields an exactly identical situation.
  */


namespace olb {


template<typename T, typename DESCRIPTOR> class SuperLatticeF3D;

/**
  * This functor gives a parabolic profile for a given point x as it computes
  * the distance between x and the axis.

    The forces set are the fake forces caused by a non-Galilean frame, here rotating around an axis with a pulse w.
    There are:
          - The centrifuge force F = rho * w * w * (x-hx, y-hy, z-hz)         | where (hx, hy, hz) is the projection of the point on the axis.
          - The coriolis force F = -2 * rho * w * (vx, vy, vz)°(Dx, Dy, Dz)   | where ° is the vector product and (Dx, Dy, Dz) is the direction vector of the axis normed to 1 (named axisDirection in the code)

     The boolean terms enable to choose having centrifugue or coriolis forces too or not.
*/
/// Functor for the rotation of forces.
template <typename T, typename DESCRIPTOR>
class RotatingForceField3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
protected:
  SuperGeometry3D<T>& sg;
  const UnitConverter<T,DESCRIPTOR>& converter;
  std::vector<T> axisPoint;
  std::vector<T> axisDirection;
  T w;
  bool centrifugeForceOn;
  bool coriolisForceOn;
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity;
  SuperLatticeDensity3D<T, DESCRIPTOR> rho;

public:
  RotatingForceField3D(SuperLattice3D<T,DESCRIPTOR>& sLattice_,
                       SuperGeometry3D<T>& superGeometry_,
                       const UnitConverter<T,DESCRIPTOR>& converter_,
                       std::vector<T> axisPoint_,
                       std::vector<T> axisDirection_,
                       T w_,
                       bool centrifugeForceOn_ = true,
                       bool coriolisForceOn_ = true);

  bool operator() (T output[], const int x[]) override;
};

template<typename T, typename DESCRIPTOR> class SuperLatticeF3D;

/**
  * This functor gives a parabolic profile for a given point x as it computes
  * the distance between x and the axis.

    The forces set are the fake forces caused by a non-Galilean frame, here harmonic oscillation around an axis.
    Oscillation is determined by an amplitude and a frequency

*/
/// Functor for the rotation of forces.
template <typename T, typename DESCRIPTOR>
class HarmonicOscillatingRotatingForceField3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
protected:
  SuperGeometry3D<T>& sg;
  const UnitConverter<T,DESCRIPTOR>& converter;
  std::vector<T> axisPoint;
  std::vector<T> axisDirection;
  T amplitude;
  T resonanceFrequency;
  T w;
  T dwdt;
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity;

public:
  HarmonicOscillatingRotatingForceField3D(SuperLattice3D<T,DESCRIPTOR>& sLattice_,
                                           SuperGeometry3D<T>& superGeometry_,
                                           const UnitConverter<T,DESCRIPTOR>& converter_,
                                           std::vector<T> axisPoint_,
                                           std::vector<T> axisDirection_,
                                           T amplitude_,
                                           T frequency_);
  void updateTimeStep(int iT);

  bool operator() (T output[], const int x[]) override;
};


} // end namespace olb

#endif
