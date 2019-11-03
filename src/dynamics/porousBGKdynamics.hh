/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Mathias J. Krause, Jonas Latt
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
 * BGK Dynamics for porous -- generic implementation.
 */
#ifndef POROUS_BGK_DYNAMICS_HH
#define POROUS_BGK_DYNAMICS_HH

#include "porousBGKdynamics.h"
#include "core/cell.h"
#include "dynamics.h"
#include "core/util.h"
#include "lbHelpers.h"
#include "math.h"

namespace olb {

////////////////////// Class PorousBGKdynamics //////////////////////////

template<typename T, typename DESCRIPTOR>
PorousBGKdynamics<T,DESCRIPTOR>::PorousBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_),
    omega(omega_)
{ }

template<typename T, typename DESCRIPTOR>
void PorousBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* porosity = cell.template getFieldPointer<descriptors::POROSITY>();
  for (int i=0; i<DESCRIPTOR::d; i++)  {
    u[i] *= porosity[0];
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T PorousBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void PorousBGKdynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  omega = omega_;
}


//////////////////// Class ExtendedPorousBGKdynamics ////////////////////

template<typename T, typename DESCRIPTOR>
ExtendedPorousBGKdynamics<T,DESCRIPTOR>::ExtendedPorousBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_),
    omega(omega_)
{
}

template<typename T, typename DESCRIPTOR>
void ExtendedPorousBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* porosity = cell.template getFieldPointer<descriptors::POROSITY>();
  T* localVelocity = cell.template getFieldPointer<descriptors::LOCAL_DRAG>();

  cell.template setField<descriptors::LOCAL_DRAG>(u);

  for (int i=0; i<DESCRIPTOR::d; i++)  {
    u[i] *= porosity[0];
    u[i] += (1.-porosity[0]) * localVelocity[i];
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T ExtendedPorousBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void ExtendedPorousBGKdynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  omega = omega_;
}

//////////////////// Class SubgridParticleBGKdynamics ////////////////////

template<typename T, typename DESCRIPTOR>
SubgridParticleBGKdynamics<T,DESCRIPTOR>::SubgridParticleBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_),
    omega(omega_)
{
  _fieldTmp[0] = T();
  _fieldTmp[1] = T();
  _fieldTmp[2] = T();
  _fieldTmp[3] = T();
}

template<typename T, typename DESCRIPTOR>
void SubgridParticleBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* porosity = cell.template getFieldPointer<descriptors::POROSITY>();
  T* extVelocity = cell.template getFieldPointer<descriptors::LOCAL_DRAG>();
//  if (porosity[0] != 0) {
//    cout << "extVelocity: " << extVelocity[0] << " " <<  extVelocity[1] << " " <<  extVelocity[2] << " " << std::endl;
//    cout << "porosity: " << porosity[0] << std::endl;
//  }
  for (int i=0; i<DESCRIPTOR::d; i++)  {
    u[i] *= (1.-porosity[0]);
    u[i] += extVelocity[i];
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);

  statistics.incrementStats(rho, uSqr);
  cell.template setField<descriptors::POROSITY>(0);
  cell.template setField<descriptors::VELOCITY_NUMERATOR>(0);
  cell.template setField<descriptors::VELOCITY_DENOMINATOR>(0);
}

template<typename T, typename DESCRIPTOR>
T SubgridParticleBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void SubgridParticleBGKdynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  omega = omega_;
}

//////////////////// Class PorousParticleBGKdynamics ////////////////////

template<typename T, typename DESCRIPTOR, bool isStatic>
PorousParticleBGKdynamics<T,DESCRIPTOR,isStatic>::PorousParticleBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_),
    omega(omega_)
{}

template<typename T, typename DESCRIPTOR, bool isStatic>
void PorousParticleBGKdynamics<T,DESCRIPTOR,isStatic>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T tmpMomentumLoss[DESCRIPTOR::d] = { };

  T* const velNumerator   = cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>();
  T* const external = cell.template getFieldPointer<descriptors::POROSITY>(); // TODO: Remove implicit layout requirements in favor of descriptor fields
  T* const velDenominator = cell.template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>();

// use HLBM with Shan-Chen forcing as in previous versions
#ifdef HLBM_FORCING_SHANCHEN
  if (*velDenominator > std::numeric_limits<T>::epsilon()) {
    effectiveVelocity<isStatic>::calculate(velNumerator, u);
    const T tmp_uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    for(int tmp_i=0; tmp_i<DESCRIPTOR::d; tmp_i++) {
      for(int tmp_iPop=0; tmp_iPop<DESCRIPTOR::q; tmp_iPop++) {
        tmpMomentumLoss[tmp_i] -= DESCRIPTOR::c[tmp_iPop][tmp_i] * omega
                * (lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(tmp_iPop, rho, u, tmp_uSqr)-cell[tmp_iPop]);
      }
    }
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);

// use Kuperstokh forcing by default
#else
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  T  uPlus[DESCRIPTOR::d] = { };
  T diff[DESCRIPTOR::q] = {};
  if (*velDenominator > std::numeric_limits<T>::epsilon()) {
    for(int iDim=0; iDim<DESCRIPTOR::d; iDim++)
      uPlus[iDim] = u[iDim];
    effectiveVelocity<isStatic>::calculate(external, uPlus);
    const T uPlusSqr = util::normSqr<T,DESCRIPTOR::d>(uPlus);
    for(int tmp_iPop=0; tmp_iPop<DESCRIPTOR::q; tmp_iPop++) {
      diff[tmp_iPop] += lbHelpers<T,DESCRIPTOR>::equilibrium(tmp_iPop, rho, uPlus, uPlusSqr)
                      - lbHelpers<T,DESCRIPTOR>::equilibrium(tmp_iPop, rho, u, uSqr);
      cell[tmp_iPop] += diff[tmp_iPop];
      for(int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
        tmpMomentumLoss[iDim] -= descriptors::c<DESCRIPTOR>(tmp_iPop,iDim) * diff[tmp_iPop];
      }
    }
    // alternate option to calculate the force
    //for(int iDim=0; iDim<DESCRIPTOR::d; iDim++)
    //  tmpMomentumLoss[iDim] = -rho*(uPlus[iDim]-u[iDim]);
  }
#endif

  statistics.incrementStats(rho, uSqr);
  for(int i_dim=0; i_dim<DESCRIPTOR::d; i_dim++) {
    *(velNumerator+i_dim) = tmpMomentumLoss[i_dim];
  }
}

// TODO: Update for meta descriptor
template<typename T, typename DESCRIPTOR, bool isStatic>
template<bool dummy>
struct PorousParticleBGKdynamics<T,DESCRIPTOR,isStatic>::effectiveVelocity<false, dummy> {
  static void calculate(T* pExternal, T* pVelocity) {
    for (int i=0; i<DESCRIPTOR::d; i++)  {
      pVelocity[i] += (1.-pExternal[DESCRIPTOR::template index<descriptors::POROSITY>() - DESCRIPTOR::q])
          * (pExternal[DESCRIPTOR::template index<descriptors::VELOCITY_NUMERATOR>() - DESCRIPTOR::q +i] / pExternal[DESCRIPTOR::template index<descriptors::VELOCITY_DENOMINATOR>() - DESCRIPTOR::q] - pVelocity[i]);
    }
    pExternal[DESCRIPTOR::template index<descriptors::POROSITY>() - DESCRIPTOR::q] = 1.;
    pExternal[DESCRIPTOR::template index<descriptors::VELOCITY_DENOMINATOR>() - DESCRIPTOR::q] = 0.;
  }
};

// TODO: Update for meta descriptor
template<typename T, typename Lattice, bool isStatic>
template<bool dummy>
struct PorousParticleBGKdynamics<T,Lattice,isStatic>::effectiveVelocity<true, dummy> {
  static void calculate(T* pExternal, T* pVelocity) {
    for (int i=0; i<Lattice::d; i++)  {
      pVelocity[i] -= (1.-pExternal[Lattice::template index<descriptors::POROSITY>() - Lattice::q]) * pVelocity[i];
    }
  }
};

template<typename T, typename DESCRIPTOR, bool isStatic>
T PorousParticleBGKdynamics<T,DESCRIPTOR,isStatic>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR, bool isStatic>
void PorousParticleBGKdynamics<T,DESCRIPTOR,isStatic>::setOmega(T omega_)
{
  omega = omega_;
}

//////////////////// Class KrauseHBGKdynamics ////////////////////

template<typename T, typename DESCRIPTOR>
KrauseHBGKdynamics<T,DESCRIPTOR>::KrauseHBGKdynamics (T omega_,
    Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_, T dx_, T dt_ )
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_), omega(omega_), smagoConst(smagoConst_),
    preFactor(computePreFactor(omega_,smagoConst_) )
{
  _fieldTmp[0] = T(1);
  _fieldTmp[1] = T();
  _fieldTmp[2] = T();
  _fieldTmp[3] = T();
}

template<typename T, typename DESCRIPTOR>
void KrauseHBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  T newOmega[DESCRIPTOR::q];
  this->_momenta.computeRhoU(cell, rho, u);
  computeOmega(this->getOmega(), cell, preFactor, rho, u, newOmega);

  T vel_denom = *cell.template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>();
  if (vel_denom > std::numeric_limits<T>::epsilon()) {
    T porosity = *cell.template getFieldPointer<descriptors::POROSITY>(); // prod(1-smoothInd)
    T* vel_num = cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>();
    porosity = 1.-porosity; // 1-prod(1-smoothInd)
    for (int i=0; i<DESCRIPTOR::d; i++)  {
      u[i] += porosity * (vel_num[i] / vel_denom - u[i]);
    }
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  statistics.incrementStats(rho, uSqr);

  cell.template setField<descriptors::POROSITY>(_fieldTmp[0]);
  cell.template setField<descriptors::VELOCITY_NUMERATOR>({_fieldTmp[1], _fieldTmp[2]});
  cell.template setField<descriptors::VELOCITY_DENOMINATOR>(_fieldTmp[3]);
}

template<typename T, typename DESCRIPTOR>
T KrauseHBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void KrauseHBGKdynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  this->setOmega(omega);
  preFactor = computePreFactor(omega, smagoConst);
}

template<typename T, typename DESCRIPTOR>
T KrauseHBGKdynamics<T,DESCRIPTOR>::computePreFactor(T omega, T smagoConst)
{
  return (T)smagoConst*smagoConst*descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>()*2*sqrt(2);
}


template<typename T, typename DESCRIPTOR>
void KrauseHBGKdynamics<T,DESCRIPTOR>::computeOmega(T omega0, Cell<T,DESCRIPTOR>& cell, T preFactor, T rho,
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
    T tau_turb = 0.5*(sqrt(tau_mol*tau_mol+(preFactor*fNeq))-tau_mol);
    /// Effective realaxation time
    tau_eff = tau_mol + tau_turb;
    newOmega[iPop] = 1./tau_eff;
  }
}


//////////////////// Class ParticlePorousBGKdynamics ////////////////////
/*
template<typename T, typename DESCRIPTOR>
ParticlePorousBGKdynamics<T,DESCRIPTOR>::ParticlePorousBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_),
    omega(omega_)
{ }

template<typename T, typename DESCRIPTOR>
void ParticlePorousBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* porosity = cell.template getFieldPointer<descriptors::POROSITY>();
  T* localVelocity = cell.template getFieldPointer<descriptors::LOCAL_DRAG>();
  for (int i=0; i<DESCRIPTOR::d; i++)  {
    u[i] *= porosity[0];
    u[i] += localVelocity[i];
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  statistics.incrementStats(rho, uSqr);

//  *cell.template getFieldPointer<descriptors::POROSITY>() = 1;
//  *cell.template getFieldPointer<descriptors::LOCAL_DRAG>() = 0.;
//  *(cell.template getFieldPointer<descriptors::LOCAL_DRAG>()+1) = 0.;
}

template<typename T, typename DESCRIPTOR>
T ParticlePorousBGKdynamics<T,DESCRIPTOR>::getOmega() const {
  return omega;
}

template<typename T, typename DESCRIPTOR>
void ParticlePorousBGKdynamics<T,DESCRIPTOR>::setOmega(T omega_) {
  omega = omega_;
}
*/

//////////////////// Class SmallParticleBGKdynamics ////////////////////

template<typename T, typename DESCRIPTOR>
SmallParticleBGKdynamics<T,DESCRIPTOR>::SmallParticleBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_),
    omega(omega_)
{ }

template<typename T, typename DESCRIPTOR>
void SmallParticleBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* porosity = cell.template getFieldPointer<descriptors::POROSITY>();
  T* localVelocity = cell.template getFieldPointer<descriptors::LOCAL_DRAG>();

  //cout << porosity[0]  << " " <<   localVelocity[0]<< " " <<   localVelocity[1]<< " " <<   localVelocity[2]<<std::endl;
  for (int i=0; i<DESCRIPTOR::d; i++)  {
    u[i] *= porosity[0];
    u[i] += localVelocity[i];
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T SmallParticleBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void SmallParticleBGKdynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  omega = omega_;
}

////////////////////// Class PSMBGKdynamics //////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
PSMBGKdynamics<T,DESCRIPTOR>::PSMBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_, int mode_ )
  : BGKdynamics<T,DESCRIPTOR>(omega_, momenta_),
    omega(omega_), paramA(1. / omega_ - 0.5)
{
  mode = (Mode) mode_;
}

template<typename T, typename DESCRIPTOR>
void PSMBGKdynamics<T,DESCRIPTOR>::computeU (Cell<T,DESCRIPTOR> const& cell, T u[DESCRIPTOR::d] ) const
{
  T rho;
  this->_momenta.computeRhoU(cell, rho, u);
//  T epsilon = 1. - *(cell.template getFieldPointer<descriptors::POROSITY>());
//  // speed up paramB
//  T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
//  // speed up paramC
//  T paramC = (1. - paramB);
//  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
//    u[iVel] = paramC * u[iVel] +
//              paramB * cell.template getFieldPointer<descriptors::VELOCITY_SOLID>()[iVel];
//  }
}

template<typename T, typename DESCRIPTOR>
void PSMBGKdynamics<T,DESCRIPTOR>::computeRhoU (Cell<T,DESCRIPTOR> const& cell, T& rho, T u[DESCRIPTOR::d] ) const
{
    this->_momenta.computeRhoU(cell, rho, u);
//  T epsilon = 1. - *(cell.template getFieldPointer<descriptors::POROSITY>());
//  // speed up paramB
//  T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
//  // speed up paramC
//  T paramC = (1. - paramB);
//  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
//    u[iVel] = paramC * u[iVel] +
//              paramB * cell.template getFieldPointer<descriptors::VELOCITY_SOLID>()[iVel];
//  }
}

template<typename T, typename DESCRIPTOR>
void PSMBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d], uSqr;
  T epsilon = 1. - *(cell.template getFieldPointer<descriptors::POROSITY>());

  this->_momenta.computeRhoU(cell, rho, u);
  // velocity at the boundary
  T u_s[DESCRIPTOR::d];
  for (int i = 0; i < DESCRIPTOR::d; i++) {
    u_s[i] = (cell.template getFieldPointer<descriptors::VELOCITY_SOLID>())[i];
  }

  if (epsilon < 1e-5) {
    uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  } else {
    // speed up paramB
    T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
    // speed up paramC
    T paramC = (1. - paramB);

    T omega_s[DESCRIPTOR::q];
    T cell_tmp[DESCRIPTOR::q];

    const T uSqr_s = util::normSqr<T,DESCRIPTOR::d>(u_s);

    uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell_tmp[iPop] = cell[iPop];
      switch(mode){
        case M2: omega_s[iPop] = (lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u_s, uSqr_s ) - cell[iPop])
                         + (1 - omega) * (cell[iPop] - lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr )); break;
        case M3: omega_s[iPop] = (cell[descriptors::opposite<DESCRIPTOR>(iPop)] - lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(descriptors::opposite<DESCRIPTOR>(iPop), rho, u_s, uSqr_s ))
               - (cell[iPop] - lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u_s, uSqr_s ));
      }

    }

    uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = cell_tmp[iPop] + paramC * (cell[iPop] - cell_tmp[iPop]);
      cell[iPop] += paramB * omega_s[iPop];
    }
    for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
      u[iVel] = paramC * u[iVel] + paramB * u_s[iVel];
    }
    uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  }
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T PSMBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void PSMBGKdynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  paramA = (1. / omega_ - 0.5);
  omega = omega_;

}

////////////////////// Class ForcedPSMBGKdynamics //////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
ForcedPSMBGKdynamics<T,DESCRIPTOR>::ForcedPSMBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_, int mode_ )
  : ForcedBGKdynamics<T,DESCRIPTOR>(omega_, momenta_),
    omega(omega_), paramA(1. / omega_ - 0.5)
{
  mode = (Mode) mode_;
}


template<typename T, typename DESCRIPTOR>
void ForcedPSMBGKdynamics<T,DESCRIPTOR>::computeU (Cell<T,DESCRIPTOR> const& cell, T u[DESCRIPTOR::d] ) const
{
  T rho;
  this->_momenta.computeRhoU(cell, rho, u);
  T epsilon = 1. - *(cell.template getFieldPointer<descriptors::POROSITY>());
  // speed up paramB
  T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
  // speed up paramC
  T paramC = (1. - paramB);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] = paramC * (u[iVel] + cell.template getFieldPointer<descriptors::FORCE>()[iVel] / (T)2.) +
              paramB * cell.template getFieldPointer<descriptors::VELOCITY_SOLID>()[iVel];
  }
}

template<typename T, typename DESCRIPTOR>
void ForcedPSMBGKdynamics<T,DESCRIPTOR>::computeRhoU (Cell<T,DESCRIPTOR> const& cell, T& rho, T u[DESCRIPTOR::d] ) const
{
  this->_momenta.computeRhoU(cell, rho, u);
  T epsilon = 1. - *(cell.template getFieldPointer<descriptors::POROSITY>());
  // speed up paramB
  T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
  // speed up paramC
  T paramC = (1. - paramB);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] = paramC * (u[iVel] + cell.template getFieldPointer<descriptors::FORCE>()[iVel] / (T)2.) +
              paramB * cell.template getFieldPointer<descriptors::VELOCITY_SOLID>()[iVel];
  }
}

template<typename T, typename DESCRIPTOR>
void ForcedPSMBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d], uSqr;
  T epsilon = 1. - *(cell.template getFieldPointer<descriptors::POROSITY>());

  this->_momenta.computeRhoU(cell, rho, u);

  T* force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  // velocity at the boundary
  T u_s[DESCRIPTOR::d];
  for (int i = 0; i < DESCRIPTOR::d; i++) {
    u_s[i] = (cell.template getFieldPointer<descriptors::VELOCITY_SOLID>())[i];
  }

  if (epsilon < 1e-5) {
    uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
    lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, omega, rho);
  } else {
    // speed up paramB
    T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
    // speed up paramC
    T paramC = (1. - paramB);

    T omega_s[DESCRIPTOR::q];
    T cell_tmp[DESCRIPTOR::q];

    const T uSqr_s = util::normSqr<T,DESCRIPTOR::d>(u_s);

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell_tmp[iPop] = cell[iPop];
      switch(mode){
        case M2: omega_s[iPop] = (lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u_s, uSqr_s ) - cell[iPop])
                         + (1 - omega) * (cell[iPop] - lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr )); break;
        case M3: omega_s[iPop] = (cell[descriptors::opposite<DESCRIPTOR>(iPop)] - lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(descriptors::opposite<DESCRIPTOR>(iPop), rho, u_s, uSqr_s ))
               - (cell[iPop] - lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u_s, uSqr_s ));
      }
    }

    uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
    lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, omega, rho);

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = cell_tmp[iPop] + paramC * (cell[iPop] - cell_tmp[iPop]);
      cell[iPop] += paramB * omega_s[iPop];
    }
    for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
      u[iVel] = paramC * u[iVel] + paramB * u_s[iVel];
    }
    uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  }
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T ForcedPSMBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void ForcedPSMBGKdynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  paramA = (1. / omega_ - 0.5);
  omega = omega_;

}


} // olb

#endif
