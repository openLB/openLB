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

#ifndef LATTICE_FRAME_CHANGE_F_3D_HH
#define LATTICE_FRAME_CHANGE_F_3D_HH

#include<cmath>

#include "latticeFrameChangeF3D.h"
#include "functors/analytical/frameChangeF2D.h"
#include "core/superLattice3D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "utilities/vectorHelpers.h"  // for normalize
#include "geometry/superGeometry3D.h"

namespace olb {


template <typename T, typename DESCRIPTOR>
RotatingForceField3D<T,DESCRIPTOR>::RotatingForceField3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice_, SuperGeometry3D<T>& superGeometry_,
 const UnitConverter<T,DESCRIPTOR>& converter_, std::vector<T> axisPoint_,
 std::vector<T> axisDirection_, T w_, bool centrifugeForceOn_,
 bool coriolisForceOn_)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice_,3), sg(superGeometry_),
    converter(converter_), axisPoint(axisPoint_), axisDirection(axisDirection_),
    w(w_), centrifugeForceOn(centrifugeForceOn_), coriolisForceOn(coriolisForceOn_),
    velocity(sLattice_,converter_), rho(sLattice_)
{
  this->getName() = "rotatingForce";
}


template <typename T, typename DESCRIPTOR>
bool RotatingForceField3D<T,DESCRIPTOR>::operator()(T output[], const int x[])
{
  std::vector<T> F_centri(3,0);
  std::vector<T> F_coriolis(3,0);

  if ( this->_sLattice.getLoadBalancer().rank(x[0]) == singleton::mpi().getRank() ) {
    // local coords are given, fetch local cell and compute value(s)
    std::vector<T> physR(3,T());
    this->sg.getCuboidGeometry().getPhysR(&(physR[0]),&(x[0]));

    T scalar =  (physR[0]-axisPoint[0])*axisDirection[0]
                +(physR[1]-axisPoint[1])*axisDirection[1]
                +(physR[2]-axisPoint[2])*axisDirection[2];

    if (centrifugeForceOn) {
      F_centri[0] = w*w*(physR[0]-axisPoint[0]-scalar*axisDirection[0]);
      F_centri[1] = w*w*(physR[1]-axisPoint[1]-scalar*axisDirection[1]);
      F_centri[2] = w*w*(physR[2]-axisPoint[2]-scalar*axisDirection[2]);
    }
    if (coriolisForceOn) {
      T _vel[3];
      (velocity)(_vel,x);
      F_coriolis[0] = -2*w*(axisDirection[1]*_vel[2]-axisDirection[2]*_vel[1]);
      F_coriolis[1] = -2*w*(axisDirection[2]*_vel[0]-axisDirection[0]*_vel[2]);
      F_coriolis[2] = -2*w*(axisDirection[0]*_vel[1]-axisDirection[1]*_vel[0]);
    }
    // return latticeForce
    output[0] = (F_coriolis[0]+F_centri[0]) * converter.getConversionFactorTime() / converter.getConversionFactorVelocity();
    output[1] = (F_coriolis[1]+F_centri[1]) * converter.getConversionFactorTime() / converter.getConversionFactorVelocity();
    output[2] = (F_coriolis[2]+F_centri[2]) * converter.getConversionFactorTime() / converter.getConversionFactorVelocity();
  }
  return true;
}


template <typename T, typename DESCRIPTOR>
HarmonicOscillatingRotatingForceField3D<T,DESCRIPTOR>::HarmonicOscillatingRotatingForceField3D
(SuperLattice3D<T,DESCRIPTOR>& sLattice_, SuperGeometry3D<T>& superGeometry_,
 const UnitConverter<T,DESCRIPTOR>& converter_, std::vector<T> axisPoint_,
 std::vector<T> axisDirection_, T amplitude_, T frequency_)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice_,3), sg(superGeometry_),
    converter(converter_), axisPoint(axisPoint_), axisDirection(axisDirection_),
    amplitude(amplitude_), resonanceFrequency(2.*4.*std::atan(1.0)*frequency_), w(0.0), dwdt(0.0),
    velocity(sLattice_,converter_)
{
  this->getName() = "harmonicOscillatingrotatingForce";
}

template <typename T, typename DESCRIPTOR>
void HarmonicOscillatingRotatingForceField3D<T,DESCRIPTOR>::updateTimeStep(int iT) {
  w = resonanceFrequency * amplitude * cos(resonanceFrequency*converter.getPhysTime(iT));
  dwdt = -resonanceFrequency * resonanceFrequency * amplitude * sin(resonanceFrequency*converter.getPhysTime(iT));
}


template <typename T, typename DESCRIPTOR>
bool HarmonicOscillatingRotatingForceField3D<T,DESCRIPTOR>::operator()(T output[], const int x[])
{


  std::vector<T> F_centri(3,0);
  std::vector<T> F_coriolis(3,0);
  std::vector<T> F_euler(3,0);

  if ( this->_sLattice.getLoadBalancer().rank(x[0]) == singleton::mpi().getRank() ) {
    // local coords are given, fetch local cell and compute value(s)
    std::vector<T> physR(3,T());
    this->sg.getCuboidGeometry().getPhysR(&(physR[0]),&(x[0]));

    T scalar =  (physR[0]-axisPoint[0])*axisDirection[0]
                +(physR[1]-axisPoint[1])*axisDirection[1]
                +(physR[2]-axisPoint[2])*axisDirection[2];

    T r[3];
    r[0] = physR[0]-axisPoint[0]-scalar*axisDirection[0];
    r[1] = physR[1]-axisPoint[1]-scalar*axisDirection[1];
    r[2] = physR[2]-axisPoint[2]-scalar*axisDirection[2];

    F_centri[0] = w*w*(r[0]);
    F_centri[1] = w*w*(r[1]);
    F_centri[2] = w*w*(r[2]);

    T _vel[3];
    (velocity)(_vel,x);
    F_coriolis[0] = -2*w*(axisDirection[1]*_vel[2]-axisDirection[2]*_vel[1]);
    F_coriolis[1] = -2*w*(axisDirection[2]*_vel[0]-axisDirection[0]*_vel[2]);
    F_coriolis[2] = -2*w*(axisDirection[0]*_vel[1]-axisDirection[1]*_vel[0]);

    F_euler[0] = -dwdt*(axisDirection[1]*r[2]-axisDirection[2]*r[1]);
    F_euler[1] = -dwdt*(axisDirection[2]*r[0]-axisDirection[0]*r[2]);
    F_euler[2] = -dwdt*(axisDirection[0]*r[1]-axisDirection[1]*r[0]);


    // return latticeForce
    output[0] = (F_coriolis[0]+F_centri[0]+F_euler[0]) * converter.getConversionFactorTime() / converter.getConversionFactorVelocity();
    output[1] = (F_coriolis[1]+F_centri[1]+F_euler[1]) * converter.getConversionFactorTime() / converter.getConversionFactorVelocity();
    output[2] = (F_coriolis[2]+F_centri[2]+F_euler[2]) * converter.getConversionFactorTime() / converter.getConversionFactorVelocity();
  }
  return true;
}


} // end namespace olb

#endif
