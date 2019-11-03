/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Peter Weisbrod
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

#ifndef INTERACTION_POTENTIAL_HH
#define INTERACTION_POTENTIAL_HH


#include "dynamics/interactionPotential.h"


namespace olb {



template <typename T, typename S>
ShanChen93<T,S>::ShanChen93(T rhoZero) : AnalyticalF1D<T,S>(1), _rhoZero(rhoZero)
{
  this->getName() = "ShanChen93";
}

template <typename T, typename S>
bool ShanChen93<T,S>::operator()(T psi[], const S rho[])
{
  psi[0]=sqrt(_rhoZero)*(1-exp(-(rho[0]/_rhoZero)));
  return true;
}


template <typename T, typename S>
ShanChen94<T,S>::ShanChen94(T rhoZero, T psiZero) : AnalyticalF1D<T,S>(1), _rhoZero(rhoZero), _psiZero(psiZero)
{
  this->getName() = "ShanChen94";
}

template <typename T, typename S>
bool ShanChen94<T,S>::operator()(T psi[], const S rho[])
{
  psi[0]=_psiZero*exp(-_rhoZero/rho[0]);
  return true;
}


template <typename T, typename S>
PengRobinson<T,S>::PengRobinson(T G, T acentricFactor, T a, T b, T tr) : AnalyticalF1D<T,S>(1), _G(G), _acentricFactor(acentricFactor), _a(a), _b(b)
{
  _R = 1.;
  //a=0.45724*R*R*tc*tc/pc;
  //b=0.0778*R*tc/pc;
  _tc = 0.0778/0.45724*_a/_b/_R;
  //T pc = 0.0778*_R*tc/_b;
  //T rhoc = pc/0.307/_R/tc;
  _t = _tc*tr;
  //Zc=0.307 Tc=0.072922004 pc=0.059569985 rhoc=2.6609121
  _alpha = 1. + (0.37464+1.54226*_acentricFactor-0.26992*_acentricFactor*_acentricFactor)*(1.-sqrt(_t/_tc));
  _alpha = _alpha*_alpha;
  this->getName() = "PengRobinson";
}

template <typename T, typename S>
bool PengRobinson<T,S>::operator()(T psi[], const S rho[])
{
  T p = (rho[0]*_R*_t/(1.-_b*rho[0]))-(_a*_alpha*rho[0]*rho[0]/(1.+2.*_b*rho[0]-_b*_b*rho[0]*rho[0]));
  psi[0] = sqrt(6.*(p-rho[0]/3.)/_G);
  return true;
}

// second operator allows to incorporate temperature changes
template <typename T, typename S>
bool PengRobinson<T,S>::operator()(T psi[], const S rho[], const S t[])
{
  _t = t[0];
  _alpha = 1. + (0.37464+1.54226*_acentricFactor-0.26992*_acentricFactor*_acentricFactor)*(1.-sqrt(_t/_tc));
  _alpha = _alpha*_alpha;
  T p = (rho[0]*_R*_t/(1.-_b*rho[0]))-(_a*_alpha*rho[0]*rho[0]/(1.+2.*_b*rho[0]-_b*_b*rho[0]*rho[0]));
  psi[0] = sqrt(6.*(p-rho[0]/3.)/_G);
  return true;
}


template <typename T, typename S>
CarnahanStarling<T,S>::CarnahanStarling(T G, T a, T b, T tr) : AnalyticalF1D<T,S>(1), _G(G), _a(a), _b(b)
{
  _R = 1.;
  //a=0.4963*tc*tc*R*R/pc;
  //b=0.18727*R*tc/pc;
  T tc = 0.18727/0.4963*_a/_b/_R;
  //T pc = 0.18727*_R*tc/_b;
  //T rhoc = pc/0.35930763/_R/tc;
  _t = tc*tr;
  //Zc=0.35930763 Tc=0.094333065 pc=0.0044164383 rhoc=0.13029921
  this->getName() = "CarnahanStarling";
}

template <typename T, typename S>
bool CarnahanStarling<T,S>::operator()(T psi[], const S rho[])
{
  T c = _b*rho[0]/4.;
  T p = rho[0]*_R*_t*((1.+c+c*c-c*c*c)/(1.-c)/(1.-c)/(1.-c))-_a*rho[0]*rho[0];
  psi[0] = sqrt(6.*(p-rho[0]/3.)/_G);
  return true;
}

// second operator allows to incorporate temperature changes
template <typename T, typename S>
bool CarnahanStarling<T,S>::operator()(T psi[], const S rho[], const S t[])
{
  _t = t[0];
  T c = _b*rho[0]/4.;
  T p = rho[0]*_R*_t*((1.+c+c*c-c*c*c)/(1.-c)/(1.-c)/(1.-c))-_a*rho[0]*rho[0];
  psi[0] = sqrt(6.*(p-rho[0]/3.)/_G);
  return true;
}


template <typename T, typename S>
PsiEqualsRho<T,S>::PsiEqualsRho() : AnalyticalF1D<T,S>(1)
{
  this->getName() = "PsiEqualsRho";
}

template <typename T, typename S>
bool PsiEqualsRho<T,S>::operator()(T psi[], const S rho[])
{
  psi[0]=rho[0];
  return true;
}


template <typename T, typename S>
Krause<T,S>::Krause(T rhoZero, T psiZero) : AnalyticalF1D<T,S>(1), _rhoZero(rhoZero), _psiZero(psiZero)
{
  this->getName() = "Krause";
}

template <typename T, typename S>
bool Krause<T,S>::operator()(T psi[], const S rho[])
{
  psi[0]=_psiZero/1.77*1.414/_rhoZero*exp(-(_rhoZero-rho[0])*(_rhoZero-rho[0])/_rhoZero/_rhoZero);
  return true;
}


template <typename T, typename S>
WeisbrodKrause<T,S>::WeisbrodKrause(T rhoZero, T sigmu) : AnalyticalF1D<T,S>(1), _rhoZero(rhoZero), _sigmu(sigmu)
{
  _rhoZero=_rhoZero/1.005088;
  this->getName() = "WeisbrodKrause";
}

template <typename T, typename S>
bool WeisbrodKrause<T,S>::operator()(T psi[], const S rho[])
{
  psi[0]=sqrt(_rhoZero)*1.5179/_sigmu*exp(-(_sigmu-rho[0]/_rhoZero)*(_sigmu-rho[0]/_rhoZero)/_sigmu/_sigmu);
  return true;
}


template <typename T, typename S>
Normal<T,S>::Normal(T sigma, T mu) : AnalyticalF1D<T,S>(1), _sigma(sigma), _mu(mu)
{
  this->getName() = "Normal";
}

template <typename T, typename S>
bool Normal<T,S>::operator()(T psi[], const S rho[])
{
  psi[0]=1./2.507/_sigma*exp(-(rho[0]-_mu)*(rho[0]-_mu)/_sigma/_sigma/2.);
  return true;
}


} // end namespace olb

#endif
