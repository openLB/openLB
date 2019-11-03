/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Robin Trunk
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

#ifndef ADVECTION_DIFFUSION_FORCES_HH
#define ADVECTION_DIFFUSION_FORCES_HH

#include "advectionDiffusionForces.h"

namespace olb {

template<typename T, typename DESCRIPTOR,
typename ADLattice>
AdvDiffDragForce3D<T,DESCRIPTOR,ADLattice>::AdvDiffDragForce3D(UnitConverter<T,DESCRIPTOR> const& converter_, T St_)
{
  initArg = 8;
  dragCoeff = (converter_.getCharPhysVelocity()*converter_.getConversionFactorTime()) / (St_ * converter_.getCharPhysLength());
}

template<typename T, typename DESCRIPTOR,
typename ADLattice>
AdvDiffDragForce3D<T,DESCRIPTOR,ADLattice>::AdvDiffDragForce3D(UnitConverter<T,DESCRIPTOR> const& converter_, T pRadius_, T pRho_)
{
  initArg = 8;
  dragCoeff = (9.*converter_.getPhysViscosity()*converter_.getPhysDensity()*converter_.getConversionFactorTime()) / (2.*pRho_*pRadius_*pRadius_);
}

template<typename T, typename DESCRIPTOR,
typename ADLattice>
void AdvDiffDragForce3D<T,DESCRIPTOR,ADLattice>::applyForce(T force[], Cell<T,DESCRIPTOR> *nsCell, Cell<T,ADLattice> *adCell, T vel[], int latticeR[])
{
  T velF[3] = {0.,0.,0.};
  nsCell->computeU(velF);
  for (int i=0; i < DESCRIPTOR::d; i++) {
    force[i] += dragCoeff*(velF[i]-vel[i]);
  }
}

template<typename T, typename DESCRIPTOR,
typename ADLattice>
AdvDiffRotatingForce3D<T,DESCRIPTOR,ADLattice>::AdvDiffRotatingForce3D(SuperGeometry3D<T>& superGeometry_,
    const UnitConverter<T,DESCRIPTOR>& converter_, std::vector<T> axisPoint_, std::vector<T> axisDirection_,
    T w_, T* frac_, bool centrifugeForceOn_, bool coriolisForceOn_) :
  sg(superGeometry_), axisPoint(axisPoint_), axisDirection(axisDirection_),
  w(w_), frac(frac_), centrifugeForceOn(centrifugeForceOn_), coriolisForceOn(coriolisForceOn_)
{
  invMassLessForce = converter_.getConversionFactorTime() * converter_.getConversionFactorTime() / converter_.getConversionFactorLength();
}

template<typename T, typename DESCRIPTOR,
typename ADLattice>
void AdvDiffRotatingForce3D<T,DESCRIPTOR,ADLattice>::applyForce(T force[], Cell<T,DESCRIPTOR> *nsCell, Cell<T,ADLattice> *adCell, T vel[], int latticeR[])
{
  std::vector<T> F_centri(3,0);
  std::vector<T> F_coriolis(3,0);
  T wf = w*(*frac);
//  if ( this->_sLattice.getLoadBalancer().rank(latticeR[0]) == singleton::mpi().getRank() ) {
  // local coords are given, fetch local cell and compute value(s)
  std::vector<T> physR(3,T());
  this->sg.getCuboidGeometry().getPhysR(&(physR[0]),&(latticeR[0]));

  T scalar =  (physR[0]-axisPoint[0])*axisDirection[0]
              +(physR[1]-axisPoint[1])*axisDirection[1]
              +(physR[2]-axisPoint[2])*axisDirection[2];

  if (centrifugeForceOn) {
    F_centri[0] = wf*wf*(physR[0]-axisPoint[0]-scalar*axisDirection[0]);
    F_centri[1] = wf*wf*(physR[1]-axisPoint[1]-scalar*axisDirection[1]);
    F_centri[2] = wf*wf*(physR[2]-axisPoint[2]-scalar*axisDirection[2]);
  }
  if (coriolisForceOn) {
    F_coriolis[0] = -2*wf*(axisDirection[1]*vel[2]-axisDirection[2]*vel[1]);
    F_coriolis[1] = -2*wf*(axisDirection[2]*vel[0]-axisDirection[0]*vel[2]);
    F_coriolis[2] = -2*wf*(axisDirection[0]*vel[1]-axisDirection[1]*vel[0]);
  }
  force[0] += (F_coriolis[0]+F_centri[0])*invMassLessForce;
  force[1] += (F_coriolis[1]+F_centri[1])*invMassLessForce;
  force[2] += (F_coriolis[2]+F_centri[2])*invMassLessForce;
//  }
}

template<typename T, typename DESCRIPTOR,
typename ADLattice>
AdvDiffMagneticWireForce3D<T,DESCRIPTOR,ADLattice>::AdvDiffMagneticWireForce3D(SuperGeometry3D<T>& superGeometry_, UnitConverter<T,DESCRIPTOR> const& converter_, T pMass, AnalyticalF3D<T, T>& getMagForce) : sg(superGeometry_), _getMagForce(getMagForce)
{
  initArg = 8;
  _pMass = converter_.getConversionFactorTime() / pMass;
  _conversionVelocity = converter_.getConversionFactorVelocity();
}

template<typename T, typename DESCRIPTOR,
typename ADLattice>
void AdvDiffMagneticWireForce3D<T,DESCRIPTOR,ADLattice>::applyForce(T force[], Cell<T,DESCRIPTOR> *nsCell, Cell<T,ADLattice> *adCell, T vel[], int latticeR[])
{
  std::vector<T> physR(3,T());
  this->sg.getCuboidGeometry().getPhysR(&(physR[0]),&(latticeR[0]));
  T pos[3] = { T(), T(), T() };
  pos[0] = physR[0];
  pos[1] = physR[1];
  pos[2] = physR[2];
  T forceHelp[3] = { T(), T(), T() };
  _getMagForce(forceHelp, pos);
  for (int i=0; i < DESCRIPTOR::d; i++) {
    force[i] += forceHelp[i] * _pMass / _conversionVelocity;
//    std::cout << "----->>>>> force " << forceHelp[i] << std::endl;
  }
}


}
#endif
