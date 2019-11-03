/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2015 Mathias J. Krause, Jonas Latt, Patrick Nathen
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
 * BGK Dynamics with adjusted omega -- generic implementation.
 */
#ifndef SMAGORINSKY_BGK_DYNAMICS_HH
#define SMAGORINSKY_BGK_DYNAMICS_HH

#include "smagorinskyBGKdynamics.h"
#include "core/cell.h"
#include "core/util.h"
#include "lbHelpers.h"
#include "math.h"

#include <complex> // For shear kalman Smagorinsky - Populations


namespace olb {

/// Smagorinsky Dynamics
template<typename T, typename DESCRIPTOR>
SmagorinskyDynamics<T,DESCRIPTOR>::SmagorinskyDynamics(T smagoConst_)
  : smagoConst(smagoConst_), preFactor(computePreFactor())
{ }

template<typename T, typename DESCRIPTOR>
T SmagorinskyDynamics<T,DESCRIPTOR>::computePreFactor()
{
  return (T)smagoConst*smagoConst*descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>()*2*sqrt(2);
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyDynamics<T,DESCRIPTOR>::getPreFactor()
{
  return preFactor;
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyDynamics<T,DESCRIPTOR>::getSmagoConst()
{
  return smagoConst;
}

///////////////////////// ADM BGK /////////////////////////////

/*template<typename T, typename DESCRIPTOR>
ADMBGKdynamics<T,DESCRIPTOR>::ADMBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_ )
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_), omega(omega_)
{ }

template<typename T, typename DESCRIPTOR>
void ADMBGKdynamics<T,DESCRIPTOR>::collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics )
{
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, *cell., cell[velocityBeginsAt], omega);
  statistics.incrementStats(*cell[rhoIsAt], uSqr);
}*/

///////////////////////// ForcedADM BGK /////////////////////////////

template<typename T, typename DESCRIPTOR>
ForcedADMBGKdynamics<T,DESCRIPTOR>::ForcedADMBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_ )
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_),
    omega(omega_)
{ }

template<typename T, typename DESCRIPTOR>
void ForcedADMBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  OstreamManager clout(std::cout,"Forced ADM collide:");
  T rho, u[DESCRIPTOR::d], utst[DESCRIPTOR::d];

// this->momenta.computeAllMomenta(cell, rho, utst, pi);

  T* rho_fil = cell.template getFieldPointer<descriptors::FIL_RHO>();
  T* u_filX = cell.template getFieldPointer<descriptors::LOCAL_FIL_VEL_X>();
  T* u_filY = cell.template getFieldPointer<descriptors::LOCAL_FIL_VEL_Y>();
  T* u_filZ = cell.template getFieldPointer<descriptors::LOCAL_FIL_VEL_Z>();

  u[0] = *u_filX;/// *rho_fil;
  u[1] = *u_filY;/// *rho_fil;
  u[2] = *u_filZ;/// *rho_fil;

  T* force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }

  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, *rho_fil, u, omega);

  lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, omega);

  statistics.incrementStats(rho, uSqr);
}

///////////////////// Class ConSmagorinskyBGKdynamics //////////////////////////

/////////////////// Consistent Smagorinsky BGK --> Malaspinas/Sagaut //////////////

template<typename T, typename DESCRIPTOR>
ConSmagorinskyBGKdynamics<T,DESCRIPTOR>::ConSmagorinskyBGKdynamics(T omega_,
    Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_)
  : SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega_, momenta_, smagoConst_)
{ }

template<typename T, typename DESCRIPTOR>
T ConSmagorinskyBGKdynamics<T,DESCRIPTOR>::computeEffectiveOmega (Cell<T,DESCRIPTOR>& cell)
{

  T H[util::TensorVal<DESCRIPTOR >::n];
  T conSmagoR[DESCRIPTOR::q];
  T S[util::TensorVal<DESCRIPTOR >::n];
  T tau_mol = 1./this->getOmega();
  T cs2 = 1./descriptors::invCs2<T,DESCRIPTOR>();
  T smagoConst = this->getSmagoConst();

  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(2*PiNeqNormSqr);

  //Strain rate Tensor
  if (PiNeqNorm != 0) {
    T Phi = (-0.5*(-rho*tau_mol*cs2+sqrt(rho*rho*tau_mol*tau_mol*cs2*cs2+2.0*(smagoConst*smagoConst)*rho*PiNeqNorm))/(smagoConst*smagoConst*rho*PiNeqNorm));
    for (int n = 0; n < util::TensorVal<DESCRIPTOR >::n; ++n) {
      S[n] = Phi*pi[n];
    }
  } else {
    for (int n = 0; n < util::TensorVal<DESCRIPTOR >::n; ++n) {
      S[n] = 0;
    }
  }

  //Strain rate Tensor Norm
  T SNormSqr = S[0]*S[0] + 2.0*S[1]*S[1] + S[2]*S[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    SNormSqr += S[2]*S[2] + S[3]*S[3] + 2.0*S[4]*S[4] + S[5]*S[5];
  }
  T SNorm    = sqrt(2*SNormSqr);

  //consistent Samagorinsky additional R term
  for (int q = 0; q < DESCRIPTOR::q; ++q) {
    T t = descriptors::t<T,DESCRIPTOR>(q); //lattice weights

    //Hermite-Polynom H = c*c-cs^2*kronDelta
    H[0] = descriptors::c<DESCRIPTOR>(q,0)*descriptors::c<DESCRIPTOR>(q,0)-cs2;
    H[1] = descriptors::c<DESCRIPTOR>(q,0)*descriptors::c<DESCRIPTOR>(q,1);
    H[2] = descriptors::c<DESCRIPTOR>(q,1)*descriptors::c<DESCRIPTOR>(q,1)-cs2;//2D
    if (util::TensorVal<DESCRIPTOR >::n == 6) {
      H[2] = descriptors::c<DESCRIPTOR>(q,0)*descriptors::c<DESCRIPTOR>(q,2);//3D
      H[3] = descriptors::c<DESCRIPTOR>(q,1)*descriptors::c<DESCRIPTOR>(q,1)-cs2;
      H[4] = descriptors::c<DESCRIPTOR>(q,1)*descriptors::c<DESCRIPTOR>(q,2);
      H[5] = descriptors::c<DESCRIPTOR>(q,2)*descriptors::c<DESCRIPTOR>(q,2)-cs2;
    }

    //contraction or scalar product H*S
    T contractHS = H[0]*S[0] + 2.0*H[1]*S[1] + H[2]*S[2];
    if (util::TensorVal<DESCRIPTOR >::n == 6) {
      contractHS += H[2]*S[2] + H[3]*S[3] + 2.0*H[4]*S[4] + H[5]*S[5];
    }

    //additional term
    conSmagoR[q] = t*this->getPreFactor()*SNorm*contractHS;
  }
  return conSmagoR[0];
}

/////////////////// Consistent Strain Smagorinsky BGK --> Malaspinas/Sagaut //////////////

template<typename T, typename DESCRIPTOR>
ConStrainSmagorinskyBGKdynamics<T,DESCRIPTOR>::ConStrainSmagorinskyBGKdynamics(T omega_,
    Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_)
  : SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega_, momenta_, smagoConst_)
{ }

template<typename T, typename DESCRIPTOR>
T ConStrainSmagorinskyBGKdynamics<T,DESCRIPTOR>::computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell)
{
  T S[util::TensorVal<DESCRIPTOR >::n];
  T cs2 = 1./descriptors::invCs2<T,DESCRIPTOR>();
  T tau_mol = 1./this->getOmega();
  T smagoConst_ = this->getSmagoConst();

  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.0*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(2*PiNeqNormSqr);


  //Strain Tensor
  if ( !util::nearZero(PiNeqNorm) ) {
    T Phi = (-0.5*(-rho*tau_mol*cs2+sqrt(rho*rho*tau_mol*tau_mol*cs2*cs2+2.0*(smagoConst_*smagoConst_)*rho*PiNeqNorm))/(smagoConst_*smagoConst_*rho*PiNeqNorm));
    for (int n = 0; n < util::TensorVal<DESCRIPTOR >::n; ++n) {
      S[n] = Phi*pi[n];
    }
  } else {
    for (int n = 0; n < util::TensorVal<DESCRIPTOR >::n; ++n) {
      S[n] = 0;
    }
  }

  //Strain Tensor Norm
  T SNormSqr = S[0]*S[0] + 2.0*S[1]*S[1] + S[2]*S[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    SNormSqr += S[2]*S[2] + S[3]*S[3] + 2.0*S[4]*S[4] + S[5]*S[5];
  }
  T SNorm    = sqrt(2*SNormSqr);

  /// Turbulent realaxation time
  T tau_turb = pow(smagoConst_,2)*SNorm/cs2;
  /// Effective realaxation time
  T tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;
}

///////////////////////// DYNAMIC SMAGO BGK /////////////////////////////
template<typename T, typename DESCRIPTOR>
DynSmagorinskyBGKdynamics<T,DESCRIPTOR>::DynSmagorinskyBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega_, momenta_, T(0.))
{ }

template<typename T, typename DESCRIPTOR>
T DynSmagorinskyBGKdynamics<T,DESCRIPTOR>::computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell)
{
  // computation of the relaxation time
  T v_t = 0;
  T* dynSmago = cell.template getFieldPointer<descriptors::SMAGO_CONST>();

  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  //v_t = *dynSmago*dx*dx*PiNeqNorm;
  v_t = *dynSmago*PiNeqNorm;
  T tau_t = 3*v_t;
  T tau_0 = 1/this->getOmega();
  T omega_new = 1/(tau_t+tau_0);
  return omega_new;
}

///////////////////SHEAR IMPROVED SMAGORINSKY//////////////////////////
template<typename T, typename DESCRIPTOR>
ShearSmagorinskyBGKdynamics<T,DESCRIPTOR>::ShearSmagorinskyBGKdynamics(T omega_,
    Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_)
  : SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega_, momenta_, smagoConst_)
{ }

template<typename T, typename DESCRIPTOR>
void ShearSmagorinskyBGKdynamics<T,DESCRIPTOR>::collide(Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  T newOmega = computeEffectiveOmega(cell,statistics.getTime());
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T ShearSmagorinskyBGKdynamics<T,DESCRIPTOR>::getEffectiveOmega(Cell<T,DESCRIPTOR>& cell)
{
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  // computation of the relaxation time
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);

  T* avShear = cell.template getFieldPointer<descriptors::AV_SHEAR>();
  T tau_0 = 1./this->getOmega();
  T PiNeqNorm_SISM = PiNeqNorm - *avShear;
  T tau_t = 0.5*(sqrt(tau_0*tau_0+(this->getPreFactor()*PiNeqNorm_SISM/rho))-tau_0);

  T omega_new = 1./(tau_t+tau_0);

  return omega_new;
}

template<typename T, typename DESCRIPTOR>
T ShearSmagorinskyBGKdynamics<T,DESCRIPTOR>::computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell, int iT)
{
  OstreamManager clout(std::cout,"shearImprovedCollide");

  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  // computation of the relaxation time
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);

  T* avShear = cell.template getFieldPointer<descriptors::AV_SHEAR>();
  *avShear = (*avShear*iT + PiNeqNorm)/(iT+1);
  T tau_0 = 1./this->getOmega();
  T PiNeqNorm_SISM = PiNeqNorm - *avShear;
  T tau_t = 0.5*(sqrt(tau_0*tau_0+(this->getPreFactor()*PiNeqNorm_SISM/rho))-tau_0);

  T omega_new = 1./(tau_t+tau_0);

  return omega_new;
}

///////////////////////// FORCED SHEAR SMAGO BGK /////////////////////////////
template<typename T, typename DESCRIPTOR>
ShearSmagorinskyForcedBGKdynamics<T,DESCRIPTOR>::ShearSmagorinskyForcedBGKdynamics(T omega_,
    Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_)
  : SmagorinskyForcedBGKdynamics<T,DESCRIPTOR>(omega_, momenta_, smagoConst_)
{ }

template<typename T, typename DESCRIPTOR>
void ShearSmagorinskyForcedBGKdynamics<T,DESCRIPTOR>::collide(Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  T newOmega = computeEffectiveOmega(cell, statistics.getTime());
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, newOmega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T ShearSmagorinskyForcedBGKdynamics<T,DESCRIPTOR>::getEffectiveOmega(Cell<T,DESCRIPTOR>& cell)
{
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  // computation of the relaxation time
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);

  T* avShear = cell.template getFieldPointer<descriptors::AV_SHEAR>();
  T tau_0 = 1./this->getOmega();
  T PiNeqNorm_SISM = PiNeqNorm - *avShear;
  T tau_t = 0.5*(sqrt(tau_0*tau_0+(this->getPreFactor()*PiNeqNorm_SISM/rho))-tau_0);

  T omega_new = 1./(tau_t+tau_0);

  return omega_new;
}

template<typename T, typename DESCRIPTOR>
T ShearSmagorinskyForcedBGKdynamics<T,DESCRIPTOR>::computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell, int iT)
{
  OstreamManager clout(std::cout,"shearImprovedCollide");

  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  //Creation of body force tensor (rho/2.)*(G_alpha*U_beta + U_alpha*G_Beta)
  T ForceTensor[util::TensorVal<DESCRIPTOR >::n];
  int iPi = 0;
  for (int Alpha=0; Alpha<DESCRIPTOR::d; ++Alpha) {
    for (int Beta=Alpha; Beta<DESCRIPTOR::d; ++Beta) {
      ForceTensor[iPi] = rho/2.*(cell.template getFieldPointer<descriptors::FORCE>()[Alpha]*u[Beta] + u[Alpha]*cell.template getFieldPointer<descriptors::FORCE>()[Beta]);
      ++iPi;
    }
  }
  // Creation of second-order moment off-equilibrium tensor
  for (int iPi=0; iPi < util::TensorVal<DESCRIPTOR >::n; ++iPi){
    pi[iPi] += ForceTensor[iPi];
  }
  // computation of the relaxation time
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);

  T* avShear = cell.template getFieldPointer<descriptors::AV_SHEAR>();
  *avShear = (*avShear*iT+PiNeqNorm)/(iT+1);

  T tau_0 = 1./this->getOmega();
  T PiNeqNorm_SISM = PiNeqNorm - *avShear;
  T tau_t = 0.5*(sqrt(tau_0*tau_0+(this->getPreFactor()*PiNeqNorm_SISM/rho))-tau_0);

  T omega_new = 1./(tau_t+tau_0);

  return omega_new;
}

////////////////////// Class SmagorinskyBGKdynamics //////////////////////////
/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
SmagorinskyBGKdynamics<T,DESCRIPTOR>::SmagorinskyBGKdynamics(T omega_,
    Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_)
  : SmagorinskyDynamics<T,DESCRIPTOR>(smagoConst_),
    BGKdynamics<T,DESCRIPTOR>(omega_,momenta_)
{ }

template<typename T, typename DESCRIPTOR>
void SmagorinskyBGKdynamics<T,DESCRIPTOR>::collide(Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  T newOmega = computeEffectiveOmega(cell);
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyBGKdynamics<T,DESCRIPTOR>::getEffectiveOmega(Cell<T,DESCRIPTOR>& cell)
{
  T newOmega = computeEffectiveOmega(cell);
  return newOmega;
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyBGKdynamics<T,DESCRIPTOR>::computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell)
{
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  /// Molecular realaxation time
  T tau_mol = 1. /this->getOmega();
  /// Turbulent realaxation time
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + this->getPreFactor()/rho*PiNeqNorm) - tau_mol);
  /// Effective realaxation time
  T tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;

}

///////////////////////// FORCED SMAGO BGK /////////////////////////////
template<typename T, typename DESCRIPTOR>
SmagorinskyForcedBGKdynamics<T,DESCRIPTOR>::SmagorinskyForcedBGKdynamics(T omega_,
    Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_)
  : SmagorinskyDynamics<T,DESCRIPTOR>(smagoConst_),
    ForcedBGKdynamics<T,DESCRIPTOR>(omega_,momenta_)
{ }

template<typename T, typename DESCRIPTOR>
void SmagorinskyForcedBGKdynamics<T,DESCRIPTOR>::collide(Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  T newOmega = computeEffectiveOmega(cell);
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, newOmega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyForcedBGKdynamics<T,DESCRIPTOR>::getEffectiveOmega(Cell<T,DESCRIPTOR>& cell)
{
  T newOmega = computeEffectiveOmega(cell);
  return newOmega;
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyForcedBGKdynamics<T,DESCRIPTOR>::computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell)
{
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  //Creation of body force tensor (rho/2.)*(G_alpha*U_beta + U_alpha*G_Beta)
  T ForceTensor[util::TensorVal<DESCRIPTOR >::n];
  int iPi = 0;
  for (int Alpha=0; Alpha<DESCRIPTOR::d; ++Alpha) {
    for (int Beta=Alpha; Beta<DESCRIPTOR::d; ++Beta) {
      ForceTensor[iPi] = rho/2.*(cell.template getFieldPointer<descriptors::FORCE>()[Alpha]*u[Beta] + u[Alpha]*cell.template getFieldPointer<descriptors::FORCE>()[Beta]);
      ++iPi;
    }
  }
  // Creation of second-order moment off-equilibrium tensor
  for (int iPi=0; iPi < util::TensorVal<DESCRIPTOR >::n; ++iPi){
    pi[iPi] += ForceTensor[iPi];
  }
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  /// Molecular realaxation time
  T tau_mol = 1. /this->getOmega();
  /// Turbulent realaxation time
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + this->getPreFactor()/rho*PiNeqNorm) - tau_mol);
  /// Effective realaxation time
  T tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;
}

///////////////////////// External TAU EFF LES BGK /////////////////////////////
template<typename T, typename DESCRIPTOR>
ExternalTauEffLESBGKdynamics<T,DESCRIPTOR>::ExternalTauEffLESBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_)
  : SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega_, momenta_, smagoConst_)
{ }

template<typename T, typename DESCRIPTOR>
void ExternalTauEffLESBGKdynamics<T,DESCRIPTOR>::collide(Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  T newTau = *(cell.template getFieldPointer<descriptors::TAU_EFF>());
  T newOmega = 1./newTau;
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  statistics.incrementStats(rho, uSqr);
}

///////////////////////// External TAU EFF FORCED LES BGK /////////////////////////////
template<typename T, typename DESCRIPTOR>
ExternalTauEffLESForcedBGKdynamics<T,DESCRIPTOR>::ExternalTauEffLESForcedBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_)
  : SmagorinskyForcedBGKdynamics<T,DESCRIPTOR>(omega_, momenta_, smagoConst_)
{ }

template<typename T, typename DESCRIPTOR>
void ExternalTauEffLESForcedBGKdynamics<T,DESCRIPTOR>::collide(Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  T newTau = *(cell.template getFieldPointer<descriptors::TAU_EFF>());
  T newOmega = 1./newTau;
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, newOmega, rho);
  statistics.incrementStats(rho, uSqr);
}

///////////////////////// FORCED Linear Velocity SMAGO BGK /////////////////////////////
template<typename T, typename DESCRIPTOR>
SmagorinskyLinearVelocityForcedBGKdynamics<T,DESCRIPTOR>::SmagorinskyLinearVelocityForcedBGKdynamics(T omega_,
    Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_)
  : SmagorinskyForcedBGKdynamics<T,DESCRIPTOR>(omega_, momenta_, smagoConst_)
{ }

template<typename T, typename DESCRIPTOR>
void SmagorinskyLinearVelocityForcedBGKdynamics<T,DESCRIPTOR>::collide(Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T newOmega = computeEffectiveOmega(cell);
  T* force = cell.template getFieldPointer<descriptors::FORCE>();
  int nDim = DESCRIPTOR::d;
  T forceSave[nDim];
  // adds a+Bu to force, where
  //   d=2: a1=v[0], a2=v[1], B11=v[2], B12=v[3], B21=v[4], B22=v[5]
  //   d=2: a1=v[0], a2=v[1], a3=v[2], B11=v[3], B12=v[4], B13=v[5], B21=v[6], B22=v[7], B23=v[8], B31=v[9], B32=v[10], B33=v[11]
  T* v = cell.template getFieldPointer<descriptors::V12>();
  for (int iDim=0; iDim<nDim; ++iDim) {
    forceSave[iDim] = force[iDim];
    force[iDim] += v[iDim];
    for (int jDim=0; jDim<nDim; ++jDim) {
      force[iDim] += v[jDim + iDim*nDim + nDim]*u[jDim];
    }
  }
  for (int iVel=0; iVel<nDim; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }

  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, newOmega, rho);
  statistics.incrementStats(rho, uSqr);
  // Writing back to froce fector
  for (int iVel=0; iVel<nDim; ++iVel) {
    force[iVel] = forceSave[iVel];
  }
}

////////////////////// Class KrauseBGKdynamics //////////////////////////
/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
KrauseBGKdynamics<T,DESCRIPTOR>::KrauseBGKdynamics(T omega_,
    Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_)
  : SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega_, momenta_, smagoConst_),
    preFactor(computePreFactor() )
{ }

template<typename T, typename DESCRIPTOR>
void KrauseBGKdynamics<T,DESCRIPTOR>::collide(Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  T newOmega[DESCRIPTOR::q];
  this->_momenta.computeRhoU(cell, rho, u);
  computeEffectiveOmega(this->getOmega(), cell, preFactor, rho, u, newOmega);
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T KrauseBGKdynamics<T,DESCRIPTOR>::getEffectiveOmega(Cell<T,DESCRIPTOR>& cell)
{
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  T newOmega[DESCRIPTOR::q];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  computeEffectiveOmega(this->getOmega(), cell, preFactor, rho, uTemp, newOmega);
  T newOmega_average = 0.;
  for (int iPop=0; iPop<DESCRIPTOR::q; iPop++) {
    newOmega_average += newOmega[iPop];
  }
  newOmega_average /= DESCRIPTOR::q;
  return newOmega_average;
}

template<typename T, typename DESCRIPTOR>
T KrauseBGKdynamics<T,DESCRIPTOR>::computePreFactor()
{
  return (T)this->getSmagoConst()*this->getSmagoConst()*3*descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>()*2*sqrt(2);
}

template<typename T, typename DESCRIPTOR>
void KrauseBGKdynamics<T,DESCRIPTOR>::computeEffectiveOmega(T omega0, Cell<T,DESCRIPTOR>& cell, T preFactor_, T rho,
    T u[DESCRIPTOR::d], T newOmega[DESCRIPTOR::q])
{
  T uSqr = u[0]*u[0];
  for (int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
    uSqr += u[iDim]*u[iDim];
  }
  /// Molecular realaxation time
  T tau_mol = 1./omega0;

  for (int iPop=0; iPop<DESCRIPTOR::q; iPop++) {
    T fNeq = std::fabs(cell[iPop] - lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr));
    /// Turbulent realaxation time
    T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + preFactor_/rho*fNeq) - tau_mol);
    /// Effective realaxation time
    T tau_eff = tau_mol + tau_turb;
    newOmega[iPop] = 1./tau_eff;
  }
}

////////////////////// Class WALEBGKdynamics //////////////////////////
/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
WALEBGKdynamics<T,DESCRIPTOR>::WALEBGKdynamics(T omega_,
    Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_)
  : SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega_, momenta_, smagoConst_)
{
  this->preFactor =  this->getSmagoConst()*this->getSmagoConst();
}

template<typename T, typename DESCRIPTOR>
T WALEBGKdynamics<T,DESCRIPTOR>::computePreFactor()
{
  return (T)this->getSmagoConst()*this->getSmagoConst();
}

template<typename T, typename DESCRIPTOR>
T WALEBGKdynamics<T,DESCRIPTOR>::computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell_)
{
  // velocity gradient tensor
  T g[3][3];
  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      g[i][j] = *(cell_.template getFieldPointer<descriptors::VELO_GRAD>()+(i*3 + j));
    }
  }
  // strain rate tensor
  T s[3][3];
  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      s[i][j] = (g[i][j] + g[j][i]) / 2.;
    }
  }
  // traceless symmetric part of the square of the velocity gradient tensor
  T G[3][3];
  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      G[i][j] = 0.;
    }
  }

  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      for ( int k = 0; k < 3; k++) {
        G[i][j] += (g[i][k]*g[k][j] + g[j][k]*g[k][i]) / 2.;  // The change
      }
    }
  }

  T trace = 0.;
  for ( int i = 0; i < 3; i++) {
    trace += (1./3.) * g[i][i] * g[i][i];
  }

  for ( int i = 0; i < 3; i++) {
    G[i][i] -= trace;
  }


  // inner product of the traceless symmetric part of the square of the velocity gradient tensor
  T G_ip = 0;
  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      G_ip = G[i][j] * G[i][j];
    }
  }

  // inner product of the strain rate
  T s_ip = 0;
  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      s_ip = s[i][j] * s[i][j];
    }
  }

  // Turbulent relaxation time
  T tau_turb = 3. * this->getPreFactor() * (pow(G_ip,1.5) / (pow(s_ip,2.5) + pow(G_ip,1.25)));
  if ((pow(s_ip,2.5) + pow(G_ip,1.25)) == 0) {
    tau_turb = 0.;
  }

  // Physical turbulent viscosity must be equal or higher that zero
  if (tau_turb < 0.) {
    tau_turb = 0.;
  }

  /// Molecular relaxation time
  T tau_mol = 1. /this->getOmega();

  /// Effective relaxation time
  T tau_eff = tau_mol + tau_turb;
  T omega_new = 1. / tau_eff;

  return omega_new;

}

////////////////////// Class WALEForcedBGKdynamics //////////////////////////
/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
WALEForcedBGKdynamics<T,DESCRIPTOR>::WALEForcedBGKdynamics(T omega_,
    Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_)
  : SmagorinskyForcedBGKdynamics<T,DESCRIPTOR>(omega_, momenta_, smagoConst_)
{
  this->preFactor =  this->getSmagoConst()*this->getSmagoConst();
}

template<typename T, typename DESCRIPTOR>
T WALEForcedBGKdynamics<T,DESCRIPTOR>::computePreFactor()
{
  return (T)this->getSmagoConst()*this->getSmagoConst();
}

template<typename T, typename DESCRIPTOR>
T WALEForcedBGKdynamics<T,DESCRIPTOR>::computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell_)
{
  // velocity gradient tensor
  T g[3][3];
  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      g[i][j] = *(cell_.template getFieldPointer<descriptors::VELO_GRAD>()+(i*3 + j));
    }
  }
  // strain rate tensor
  T s[3][3];
  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      s[i][j] = (g[i][j] + g[j][i]) / 2.;
    }
  }
  // traceless symmetric part of the square of the velocity gradient tensor
  T G[3][3];
  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      G[i][j] = 0.;
    }
  }

  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      for ( int k = 0; k < 3; k++) {
        G[i][j] += (g[i][k]*g[k][j] + g[j][k]*g[k][i]) / 2.;  // The change
      }
    }
  }

  T trace = 0.;
  for ( int i = 0; i < 3; i++) {
    trace += (1./3.) * g[i][i] * g[i][i];
  }

  for ( int i = 0; i < 3; i++) {
    G[i][i] -= trace;
  }


  // inner product of the traceless symmetric part of the square of the velocity gradient tensor
  T G_ip = 0;
  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      G_ip = G[i][j] * G[i][j];
    }
  }

  // inner product of the strain rate
  T s_ip = 0;
  for ( int i = 0; i < 3; i++) {
    for ( int j = 0; j < 3; j++) {
      s_ip = s[i][j] * s[i][j];
    }
  }

  // Turbulent relaxation time
  T tau_turb = 3. * this->getPreFactor() * (pow(G_ip,1.5) / (pow(s_ip,2.5) + pow(G_ip,1.25)));
  if ((pow(s_ip,2.5) + pow(G_ip,1.25)) == 0) {
    tau_turb = 0.;
  }

  // Physical turbulent viscosity must be equal or higher that zero
  if (tau_turb < 0.) {
    tau_turb = 0.;
  }

  /// Molecular relaxation time
  T tau_mol = 1. /this->getOmega();

  /// Effective relaxation time
  T tau_eff = tau_mol + tau_turb;
  T omega_new = 1. / tau_eff;

  return omega_new;

}

//////////////// Class ShearKalmanFDSmagorinskyBGKdynamics ///////////////////
/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
FDKalmanShearSmagorinskyBGKdynamics<T,DESCRIPTOR>::FDKalmanShearSmagorinskyBGKdynamics(T omega_,
    Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_,  T u_char_lat, T f_char_lat)
  : SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega_, momenta_, smagoConst_),
    VarInVelKal(pow(2.0*(4.0 * std::atan(1.0))*(u_char_lat*f_char_lat)/sqrt(3),2)),
    UCharLat(u_char_lat)
{

  this->preFactor = this->getSmagoConst() * this->getSmagoConst() * descriptors::invCs2<T,DESCRIPTOR>();

}

template<typename T, typename DESCRIPTOR>
T FDKalmanShearSmagorinskyBGKdynamics<T,DESCRIPTOR>::getEffectiveOmega(Cell<T,DESCRIPTOR>& cell)
{
  T FNSR;
  computeNormStrainRate(cell, FNSR);

  T INSR;
  computeNormStrainRate(cell, INSR);

  // Turbulent relaxation time
  T tau_turb = 0.;
  if (INSR > FNSR) {
    tau_turb = this->getPreFactor() * (INSR - FNSR);
  }

  /// Molecular relaxation time
  T tau_mol = 1. /this->getOmega();

  /// Effective Omega
  T omega_new = 1. / tau_mol + tau_turb;

  return omega_new;
}

template<typename T, typename DESCRIPTOR>
T FDKalmanShearSmagorinskyBGKdynamics<T,DESCRIPTOR>::computePreFactor()
{
  return this->getSmagoConst()*this->getSmagoConst()*descriptors::invCs2<T,DESCRIPTOR>();
}

template<typename T, typename DESCRIPTOR>
T FDKalmanShearSmagorinskyBGKdynamics<T,DESCRIPTOR>::computeOmega(Cell<T,DESCRIPTOR>& cell)
{
  OstreamManager clout(std::cout,"shearImprovedKalmanFDCollide");

  // Kalman procedure to update the filtered velocity
  KalmanStep(cell);

  // Norm of filtered Strain Rate
  T FNSR;
  computeNormStrainRate(cell, FNSR);

  // Norm of Instantaneous Strain Rate
  T INSR;
  computeNormStrainRate(cell, INSR);

  T tau_turb = 0.;
  if (INSR > FNSR) {
    tau_turb = this->getPreFactor() * (INSR - FNSR);
  }

  /// Molecular relaxation time
  T tau_mol = 1. /this->getOmega();

  /// Effective Omega
  T omega_new = 1. / tau_mol + tau_turb;

  return omega_new;

}

template<typename T, typename DESCRIPTOR>
void FDKalmanShearSmagorinskyBGKdynamics<T,DESCRIPTOR>::computeNormStrainRate(Cell<T,DESCRIPTOR>& cell, T& NormStrainRate)
{
  int Dim = DESCRIPTOR::d;
  // Velocity gradient in 2D-3D
  T VG[Dim][Dim];
  for ( int i = 0; i < Dim; i++) {
    for ( int j = 0; j < Dim; j++) {
      VG[i][j] = *(cell.template getFieldPointer<descriptors::FILTERED_VEL_GRAD>()+(i*Dim + j));
    }
  }
  // Strain rate tensor
  T S[Dim][Dim];
  for ( int i = 0; i < Dim; i++) {
    for ( int j = 0; j < Dim; j++) {
      S[i][j] = (VG[i][j] + VG[j][i]) / 2.;
    }
  }
  // inner product of the strain tensor
  T SIP = 0;
  for ( int i = 0; i < Dim; i++) {
    for ( int j = 0; j < Dim; j++) {
      SIP = S[i][j] * S[i][j];
    }
  }

  // Norm of the strain rate tensor
  NormStrainRate = sqrt(2. * SIP);

}

template<typename T, typename DESCRIPTOR>
void FDKalmanShearSmagorinskyBGKdynamics<T,DESCRIPTOR>::KalmanStep(Cell<T,DESCRIPTOR>& cell)
{
  // 1. Prediction Step
  T* ErrorCovariance = cell.template getFieldPointer<descriptors::ERROR_COVARIANCE>();
  ErrorCovariance[0] += VarInVelKal;

  // 2. Update Step
  // 2.1. Smooothing Factor : K
  T* Variance = cell.template getFieldPointer<descriptors::VARIANCE>();
  T K = ErrorCovariance[0]/(ErrorCovariance[0] + Variance[0]);

  // 2.2. Kalman filtered Velocity -> Kalman filtered Populations
  //T* KalmanPopulation = cell.template getFieldPointer<descriptors::FILTERED_POPULATION>();
  T u[DESCRIPTOR::d] = {0., 0., 0.};
  cell.computeU(u);

  T* KalmanVel = cell.template getFieldPointer<descriptors::VELOCITY>();
  for (int iVel=0; iVel<DESCRIPTOR::d; iVel++) {
    KalmanVel[iVel] = (KalmanVel[iVel] * (1-K)) + (K * u[iVel]);
  }

  // 2.3. Error covariance : P
  ErrorCovariance[0] *= (1-K);

  // 3. Adapt Step
  T epsilon = 0.1;
  T KalU_InstU[DESCRIPTOR::d] = {0., 0., 0.};
  for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
    KalU_InstU[iVel] = KalmanVel[iVel]-u[iVel];
  }

  Variance[0] = std::max(UCharLat*util::normSqr<T,DESCRIPTOR::d>(KalU_InstU),epsilon*pow(UCharLat,2));
}



//////////////////////////////////////////////////////////////////////////////
//////////// Shear Improved - Kalman Filter - Smagorinsky BGK ////////////////
template<typename T, typename DESCRIPTOR>
ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR>::ShearKalmanSmagorinskyBGKdynamics(T omega_,
    Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_, T u_char_lat, T f_char_lat)
  : SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega_, momenta_, smagoConst_),
    VarInVelKal(pow(2.0*(4.0 * std::atan(1.0))*(u_char_lat*f_char_lat)/sqrt(3),2)),
    UCharLat(u_char_lat)
{ }

template<typename T, typename DESCRIPTOR>
T ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR>::getEffectiveOmega(Cell<T,DESCRIPTOR>& cell)
{
  OstreamManager clout(std::cout,"shearImprovedKalmanCollide");

  // Compute the norm of second moment of non-equilibrium Instantaneous distribution [n+1][n+1]
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  cell.computeAllMomenta(rho, u, pi);
  T PiNeqNormSqrn1;
  computeNormSOM(pi, rho, PiNeqNormSqrn1);

  // Compute the norm of second moment of non-equilibrium filtered distribution <n+1><n+1>
  // Filtered Stress at time n+1
  T KalmanPiNeqNormSqrN1, KalmanUN1[DESCRIPTOR::d], KalmanPiNeqN1[util::TensorVal<DESCRIPTOR >::n];
  computeKalmanUStress(cell,KalmanUN1,KalmanPiNeqN1);
  computeNormSOM(KalmanPiNeqN1, KalmanPiNeqNormSqrN1);

  T tau_mol = 1./this->getOmega();
  T tau_sgs = T(0.);
  if (PiNeqNormSqrn1 > KalmanPiNeqNormSqrN1) {
    tau_sgs = 0.5*( ( pow(tau_mol,2.0) + ( this->getPreFactor()*( sqrt( PiNeqNormSqrn1)-sqrt(KalmanPiNeqNormSqrN1) ) ) ) - tau_mol );
  }

  T EffectiveOmega = 1.0/(tau_mol + tau_sgs);

  return EffectiveOmega;
}

template<typename T, typename DESCRIPTOR>
T ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR>::computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell)
{
  OstreamManager clout(std::cout,"shearImprovedKalmanCollide");

  // Update the filtered velocity wit a Kalman procedure
  KalmanStep(cell);

  // Compute the norm of second moment of non-equilibrium Instantaneous distribution [n+1][n+1]
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  cell.computeAllMomenta(rho, u, pi);
  T PiNeqNormSqrn1;
  computeNormSOM(pi, rho, PiNeqNormSqrn1);

  // Compute the norm of second moment of non-equilibrium filtered distribution <n+1><n+1>
  // Filtered Stress at time n+1
  T KalmanPiNeqNormSqrN1, KalmanUN1[DESCRIPTOR::d], KalmanPiNeqN1[util::TensorVal<DESCRIPTOR >::n];
  computeKalmanUStress(cell,KalmanUN1,KalmanPiNeqN1);
  computeNormSOM(KalmanPiNeqN1, KalmanPiNeqNormSqrN1);

  T tau_mol = 1./this->getOmega();
  T tau_sgs = T(0.);
  if (PiNeqNormSqrn1 > KalmanPiNeqNormSqrN1) {
    tau_sgs = 0.5*( ( pow(tau_mol,2.0) + ( this->getPreFactor()*( sqrt( PiNeqNormSqrn1)-sqrt(KalmanPiNeqNormSqrN1) ) ) ) - tau_mol );
  }

  T EffectiveOmega = 1.0/(tau_mol + tau_sgs);

  return EffectiveOmega;
}

template<typename T, typename DESCRIPTOR>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR>::KalmanStep(Cell<T,DESCRIPTOR>& cell)
{
  // The Kalman filter procedure //
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  cell.computeAllMomenta(rho, u, pi);

  T* KalmanPopulation = cell.template getFieldPointer<descriptors::FILTERED_POPULATION>();
  if (KalmanPopulation[0] == (T)-1.0){
    for (int iPop=0; iPop<DESCRIPTOR::q; iPop++){
      KalmanPopulation[iPop] =  cell[iPop]/rho;
    }
  }

  // 1. Prediction Step
  T* ErrorCovariance = cell.template getFieldPointer<descriptors::ERROR_COVARIANCE>();
  *ErrorCovariance += VarInVelKal;

  // 2. Update Step
  // 2.1. Smoothing Factor : K
  T* Variance = cell.template getFieldPointer<descriptors::VARIANCE>();
  T K = *ErrorCovariance/(*ErrorCovariance + *Variance);

  // 2.2. Kalman filtered Velocity -> Kalman filtered Populations
  for (int iPop=0; iPop<DESCRIPTOR::q; iPop++) {
    KalmanPopulation[iPop] = (KalmanPopulation[iPop] * (1-K)) + (K * cell[iPop]/rho);
  }

  // 2.3. Error covariance : P
  *ErrorCovariance *= (1-K);

  // 3. Adapt Step
  T epsilon = T(0.1);
  T KalU_InstU[DESCRIPTOR::d];
    // Filtered Stress at time n+1
    T KalmanUN1[DESCRIPTOR::d];
    computeKalmanU(cell, KalmanUN1);
  for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
    KalU_InstU[iD] = KalmanUN1[iD]-u[iD];
  }

  *Variance = std::max(UCharLat*util::normSqr<T,DESCRIPTOR::d>(KalU_InstU),epsilon*pow(UCharLat,2));

  // Fin filtering procedure //
}

template<typename T, typename DESCRIPTOR>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR>::computeKalmanUStress(Cell<T,DESCRIPTOR>& cell,
      T (&KalmanU)[DESCRIPTOR::d],T (&KalmanPi)[util::TensorVal<DESCRIPTOR >::n] )
{
  computeKalmanU(cell,KalmanU);
  computeKalmanStress(cell,KalmanU,KalmanPi);
}

template<typename T, typename DESCRIPTOR>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR>::computeKalmanU(Cell<T,DESCRIPTOR>& cell, T (&KalmanU)[DESCRIPTOR::d])
{
  T* KalmanPopulation = cell.template getFieldPointer<descriptors::FILTERED_POPULATION>();
  for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
    KalmanU[iD] = T();
  }
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      KalmanU[iD] += KalmanPopulation[iPop]*descriptors::c<DESCRIPTOR>(iPop,iD);
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR>::computeKalmanStress(Cell<T,DESCRIPTOR>& cell,
      T (&KalmanU)[DESCRIPTOR::d],T (&KalmanPi)[util::TensorVal<DESCRIPTOR >::n] )
{
  T* KalmanPopulation = cell.template getFieldPointer<descriptors::FILTERED_POPULATION>();

  T rhoRelative = T(0.);
  for (int iPop = 0; iPop < DESCRIPTOR::q; iPop++) {
    rhoRelative += KalmanPopulation[iPop];
  }

  int iPi = 0;
  for (int iAlpha=0; iAlpha < DESCRIPTOR::d; ++iAlpha) {
    for (int iBeta=iAlpha; iBeta < DESCRIPTOR::d; ++iBeta) {
      KalmanPi[iPi] = T();
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        KalmanPi[iPi] += descriptors::c<DESCRIPTOR>(iPop,iAlpha)*descriptors::c<DESCRIPTOR>(iPop,iBeta) * KalmanPopulation[iPop];
      }
      // stripe off equilibrium contribution
      KalmanPi[iPi] -= KalmanU[iAlpha]*KalmanU[iBeta];
      if (iAlpha==iBeta) {
        KalmanPi[iPi] -= rhoRelative/descriptors::invCs2<T,DESCRIPTOR>();
      }
      ++iPi;
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR>::computeAndupdateTauSgs(Cell<T,DESCRIPTOR>& cell,
      T rho, T pi[util::TensorVal<DESCRIPTOR >::n], T KalmanPiNeqN[util::TensorVal<DESCRIPTOR >::n],
      T KalmanPiNeqN1[util::TensorVal<DESCRIPTOR >::n], T K, T &tau_sgs)
{
  // Compute the norm of second moment of non-equilibrium filtered distribution <n><n>
  T KalmanPiNeqNormSqrN;
  computeNormSOM(KalmanPiNeqN, KalmanPiNeqNormSqrN);

  // Compute the norm of cross second moment of non-equilibrium filtered-Instantaneous distribution <n>[n+1]
  T KalmanInstPiNeqNormSqrN;
  computeNormSOM(KalmanPiNeqN, pi, rho, KalmanInstPiNeqNormSqrN);

  // Compute the norm of second moment of non-equilibrium Instantaneous distribution [n+1][n+1]
  T PiNeqNormSqrn;
  computeNormSOM(pi, rho, PiNeqNormSqrn);

  computeTauSgs(cell, rho, KalmanPiNeqNormSqrN, KalmanInstPiNeqNormSqrN, PiNeqNormSqrn, K, tau_sgs);

  // Compute the norm of second moment of non-equilibrium filtered distribution <n+1><n+1>
  T KalmanPiNeqNormSqrN1;
  computeNormSOM(KalmanPiNeqN1, KalmanPiNeqNormSqrN1);

  updateTauSgsKalman(cell, KalmanPiNeqNormSqrN, KalmanInstPiNeqNormSqrN, PiNeqNormSqrn, KalmanPiNeqNormSqrN1, K, tau_sgs);
}

template<typename T, typename DESCRIPTOR>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR>::computeNormSOM(T pi[util::TensorVal<DESCRIPTOR >::n], T &piNorm)
{
  // Compute the norm of second moment of non-equilibrium filtered distribution <-><->
  piNorm = pow(pi[0],2) + 2.0*pow(pi[1],2) + pow(pi[2],2);
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    piNorm += pow(pi[2],2) + pow(pi[3],2) + 2*pow(pi[4],2) + pow(pi[5],2);
  }
}

template<typename T, typename DESCRIPTOR>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR>::computeNormSOM(T pi1[util::TensorVal<DESCRIPTOR >::n],
       T pi2[util::TensorVal<DESCRIPTOR >::n], T rho, T &piNorm)
{
  // Compute the norm of cross second moment of non-equilibrium filtered-Instantaneous distribution <->[-]
  piNorm = pi1[0]*pi2[0] + 2.0*pi1[1]*pi2[1] + pi1[2]*pi2[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    piNorm += pi1[2]*pi2[2] + pi1[3]*pi2[3] + 2*pi1[4]*pi2[4] + pi1[5]*pi2[5];
  }
  piNorm /= rho;
}

template<typename T, typename DESCRIPTOR>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR>::computeNormSOM(T pi[util::TensorVal<DESCRIPTOR >::n], T rho, T &piNorm)
{
  // Compute the norm of second moment of non-equilibrium Instantaneous distribution [-][-]
  computeNormSOM(pi, piNorm);
  piNorm /= pow(rho,2.0);
}

template<typename T, typename DESCRIPTOR>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR>::computeTauSgs(Cell<T,DESCRIPTOR>& cell,
      T rho, T KalmanPiNeqNormSqr, T KalmanInstPiNeqNormSqr, T PiNeqNormSqr, T K, T &tau_sgs)
{
  T tau_mol = this->getOmega();
  T tau_eff_n = tau_mol + *(cell.template getFieldPointer<descriptors::TAU_SGS>());
  T F = 1.0/(pow(this->getSmagoConst()/descriptors::invCs2<T,DESCRIPTOR>(),2.0)/(2.0*tau_eff_n));

  // Coefficients of 4th polynomial : Ax^4 + Bx^3 + Cx^2 + Dx + E = 0
  T A = pow(F,2.0);
  T B = 2.0*A*tau_mol;

  T C = F*pow(tau_mol,2.0)
        - (2.0*F*sqrt(2.0*PiNeqNormSqr)*tau_eff_n)
        - (pow((1.0-K),2.0)*2.0*KalmanPiNeqNormSqr);

  T D = -((2.0*F*sqrt(2.0*PiNeqNormSqr)*tau_mol*tau_eff_n)
           + (2.0*pow((1.0-K),2.0)*tau_mol*2.0*KalmanPiNeqNormSqr)
           + (2.0*K*(1.0-K)*2.0*KalmanInstPiNeqNormSqr*tau_eff_n)
         );

  T E = ((1-pow(K,2.0))*2.0*PiNeqNormSqr*pow(tau_eff_n,2.0))
        + (pow((1.0-K)*tau_mol,2.0)*2.0*KalmanPiNeqNormSqr)
        - (2.0*K*(1.0-K)*2.0*KalmanInstPiNeqNormSqr*tau_mol*tau_eff_n);

  std::complex<T> Roots[4];
  computeRoots4thPoly(A, B, C, D, E, Roots);

  tau_sgs = 0.;
  for ( int i = 0; i < 4; i++){
    if (std::imag(Roots[i]) == T(0.)) {
      if (std::real(Roots[i]) > tau_sgs){
        tau_sgs = std::real(Roots[i]);
      }
    }
  }

  /// Update the value of instantaneous effective omega
  //T* EffectiveOmega = cell.template getFieldPointer<descriptors::EFFECTIVE_OMEGA>();
  //*EffectiveOmega = 1.0/(tau_mol + tau_sgs);

}

template<typename T, typename DESCRIPTOR>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR>::computeRoots4thPoly(T A, T B, T C, T D, T E, std::complex<T> (&Roots)[4])
{
  T p = T((8.*A*C - 3.*pow(B,2.0))/(8.0*pow(A,2.0)));
  T q = T((pow(B,3.0) - 4.0*A*B*C + 8.0*pow(A,2.0)*D)/(8.0*pow(A,3.0)));

  T Delta0 = T(pow(C,2.0) - 3.0*B*D + 12.0*A*E);
  T Delta1 = 2.0*pow(C,3.0) - 9.0*B*C*D + 27*pow(B,2.0)*E + 27.0*A*pow(D,2.0) - 72.0*A*C*E;

  // Discriminant
  std::complex<T> Dis = (pow(Delta1,2.0) - 4.0*pow(Delta0,3.0))/(-27.0);

  std::complex<T> Q = pow((Delta1+sqrt(Dis*(-27.0)))/2.0,1.0/3.0);
  std::complex<T> S = 0.5*sqrt((-2.0*p/3.0)+((1.0/(3.0*A))*(Q+(Delta0/Q))));

  std::complex<T> cas1, cas2;
  for ( int i = 0; i < 2; i++) {
    for ( int j = 0; j < 2; j++) {
      cas1 = T(2*i-1);
      cas2 = T(2*j-1);
      Roots[2*i+j] = (-B/4*A) + cas1*S + (cas2*0.5*sqrt((-4.0*S*S)+(2.0*p)+(q/S)));
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ShearKalmanSmagorinskyBGKdynamics<T,DESCRIPTOR>::updateTauSgsKalman(Cell<T,DESCRIPTOR>& cell, T NN, T Nn1, T n1n1, T N1N1, T K, T tau_sgs_n1)
{
  T tau_mol = this->getOmega();
  T* tau_sgs_N = cell.template getFieldPointer<descriptors::TAU_SGS>();

  //T tau_eff_N = tau_mol + *tau_sgs_N;
  //T tau_eff_n1 = *(cell[EffectiveOmegaIsAt]);
  T tau_eff_n1 = *tau_sgs_N;
//  T A;
//  if (!util::nearZero(tau_sgs_n1)) {
//    A = sqrt( pow((1-K) * (*tau_sgs_N) / tau_eff_N,2.0) * sqrt(NN)
//              + 2.0 * K * (1.0-K) * (*tau_sgs_N * tau_sgs_n1) * sqrt(Nn1) / (tau_eff_N * tau_eff_n1)
//              + pow(K * tau_sgs_n1 / tau_eff_n1,2.0) * sqrt(n1n1)
//           );
//  } else {
//    A = sqrt( pow((1-K) * (*tau_sgs_N) / tau_eff_N,2.0) * sqrt(NN));
//  }
//
//  A /= sqrt(N1N1);
//
//  *tau_sgs_N = tau_mol/((1.0/A)-1.0);
  *tau_sgs_N = (sqrt(2.0*N1N1)
                /(sqrt(2.0*n1n1)
                  -(2.0*pow(1./(this->getSmagoConst()*descriptors::invCs2<T,DESCRIPTOR>()),2.0)*tau_sgs_n1*tau_eff_n1)
                 )
               )*tau_eff_n1;

  if ((*tau_sgs_N - tau_mol) > T(0.)){
    *tau_sgs_N -= tau_mol;
  } else {
    *tau_sgs_N = T(0.);
  }

}


}

#endif
