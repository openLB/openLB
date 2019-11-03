/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2015 Mathias J. Krause, Vojtech Cvrcekt, Davide Dapelo
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
 * Porous-particle BGK Dynamics with adjusted omega
 * and Smagorinsky turbulence model -- generic implementation.
 * Strain rate similar to "J.Boyd, J. Buick and S.Green: A second-order accurate lattice Boltzmann non-Newtonian flow model"
 * Power Law similar to "Huidan Yu, Sharath S. Girimaji, Li-Shi Luo - DNS and LES of decaying isotropic turbulence with and without frame rotation using lattice Boltzmann method"
 */
#ifndef SMAGORINSKY_POWER_LAW_POROUS_BGK_DYNAMICS_HH
#define SMAGORINSKY_POWER_LAW_POROUS_BGK_DYNAMICS_HH

#include "SmagorinskyPowerLawPorousBGKdynamics.h"
#include "SmagorinskyPorousParticleBGKdynamics.hh"
#include "math.h"

namespace olb {

////////////////////// Class SmagorinskyPowerLawPorousParticleBGKdynamics //////////////////////////

/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
SmagorinskyPowerLawPorousParticleBGKdynamics<T,DESCRIPTOR>::SmagorinskyPowerLawPorousParticleBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_, T m_, T n_ , T dtPL_, T nuMin, T nuMax, T smagoConst_)
  : SmagorinskyPorousParticleBGKdynamics<T,DESCRIPTOR>(omega_,momenta_,smagoConst_),
    m(m_),
    n(n_),
    dtPL(dtPL_)
    //preFactor(computePreFactor(omega_,smagoConstPL_) )
{
  omegaMin = 2./(nuMax*2.*descriptors::invCs2<T,DESCRIPTOR>() + 1.);
  omegaMax = 2./(nuMin*2.*descriptors::invCs2<T,DESCRIPTOR>() + 1.);
}

template<typename T, typename DESCRIPTOR>
void SmagorinskyPowerLawPorousParticleBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  // load old omega from dyn. omega descriptor
//  T oldOmega = this->getOmega(); //compute with constant omega
  T oldOmega = cell.template getFieldPointer<descriptors::OMEGA>()[0]; //compute with dynamic omega
  T OmegaPL = computeOmegaPL(oldOmega, rho, pi);
  T* velDenominator = cell.template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>();
  T* velNumerator = cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>();
  T* porosity = cell.template getFieldPointer<descriptors::POROSITY>();
  if (*velDenominator > std::numeric_limits<T>::epsilon()) {
    *porosity = 1.-*porosity; // 1-prod(1-smoothInd)
    for (int i=0; i < DESCRIPTOR::d; i++)  {
      u[i] += *porosity * (*(velNumerator+i) / *velDenominator - u[i]);
    }
  }

  T newOmega = this->computeOmega(OmegaPL, this->preFactor, rho, pi); 

  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  // save new omega to dyn. omega descriptor
  cell.template getFieldPointer<descriptors::OMEGA>()[0] = newOmega; //compute with dynamic omega
  statistics.incrementStats(rho, uSqr);

  cell.template defineField<descriptors::POROSITY>(1.0);
  cell.template defineField<descriptors::VELOCITY_DENOMINATOR>(0.0);
  cell.template defineField<descriptors::VELOCITY_NUMERATOR>(0.0);
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyPowerLawPorousParticleBGKdynamics<T,DESCRIPTOR>::computeOmegaPL(T omega0, T rho, T pi[util::TensorVal<DESCRIPTOR >::n] )
{

  // strain rate tensor without prefactor
  T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];
  }

  T pre2 = pow(descriptors::invCs2<T,DESCRIPTOR>()/2./dtPL* omega0/rho,2.); // prefactor to the strain rate tensor
  T D = pre2*PiNeqNormSqr; // Strain rate tensor
  T gamma = sqrt(2.*D); // shear rate

  T nuNew = m*pow(gamma,n-1.); //nu for non-Newtonian fluid
  //T newOmega = 2./(nuNew*6.+1.);
  T newOmega = 2./(nuNew*2.*descriptors::invCs2<T,DESCRIPTOR>() + 1.);

  /*
     * problem if newOmega too small or too big is see e.g. "Daniel Conrad , Andreas Schneider, Martin Böhle:
     * A viscosity adaption method for Lattice Boltzmann simulations"
    */
  //if (newOmega>1.965) {
  //  newOmega = 1.965;  //std::cout << newOmega << std::endl;
  //}
  //if (newOmega<0.1) {
  //  newOmega = 0.1;  //std::cout << newOmega << std::endl;
  //}
  if (newOmega>omegaMax) {
    newOmega = omegaMax;  //std::cout << newOmega << std::endl;
  }
  if (newOmega<omegaMin) {
    newOmega = omegaMin;  //std::cout << newOmega << std::endl;
  }
//  std::cout << newOmega << std::endl;
  return newOmega;
  //return omega0;
}


}

#endif
