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

#ifndef INTERACTION_POTENTIAL_H
#define INTERACTION_POTENTIAL_H


/**
 *  The functor dimensions are given by F: S^m -> T^n  (S=source, T=target)
 *  and are implemented via GenericF(n,m).
 *  Don't get confused by the flipped order of source and target.
 */

namespace olb {

// established -- original for both single- and multicomponent flow

template <typename T, typename S>
class ShanChen93 : public AnalyticalF1D<T,S> {
private:
  T _rhoZero;
public:
  ShanChen93(T rhoZero=1.);
  bool operator() (T psi[], const S rho[]);
};

// established -- only multicomponent flow

template <typename T, typename S>
class PsiEqualsRho : public AnalyticalF1D<T,S> {
private:
public:
  PsiEqualsRho();
  bool operator() (T psi[], const S rho[]) override;
};

// established -- only singlecomponent flow

template <typename T, typename S>
class ShanChen94 : public AnalyticalF1D<T,S> {
private:
  T _rhoZero;
  T _psiZero;
public:
  ShanChen94(T rhoZero=200., T psiZero=4.);
  bool operator() (T psi[], const S rho[]) override;
};

template <typename T, typename S>
class PengRobinson : public AnalyticalF1D<T,S> {
private:
  T _G;
  T _acentricFactor;
  T _a;
  T _b;
  T _R;
  T _alpha;
  T _t;
  T _tc;
public:
  PengRobinson(T G, T acentricFactor=0.334, T a=2./49., T b=2./21., T tr=.8);
  bool operator() (T psi[], const S rho[]);
  // second operator allows to incorporate temperature changes
  bool operator() (T psi[], const S rho[], const S t[]);
};

template <typename T, typename S>
class CarnahanStarling : public AnalyticalF1D<T,S> {
private:
  T _G;
  T _a;
  T _b;
  T _R;
  T _t;
public:
  CarnahanStarling(T G, T a=1., T b=4., T tr=.7);
  bool operator() (T psi[], const S rho[]) override;
  // second operator allows to incorporate temperature changes
  bool operator() (T psi[], const S rho[], const S t[]);
};

// under development -- for singlecomponent flow

// 0.5 -> psiZero=0.65
// 1 -> psiZero=1.9
// 1.5 -> psiZero=3.5
// 2. -> psiZero=5,45
template <typename T, typename S>
class Krause : public AnalyticalF1D<T,S> {
private:
  T _rhoZero;
  T _psiZero;
public:
  Krause(T rhoZero=1., T psiZero=1.9);
  bool operator() (T psi[], const S rho[]);
};

template <typename T, typename S> // density of liquid phase always equals rhoZero for G=-1
class WeisbrodKrause : public AnalyticalF1D<T,S> {
private:
  T _rhoZero;
  T _sigmu;
public:
  WeisbrodKrause(T rhoZero=1., T sigmu=1.);
  bool operator() (T psi[], const S rho[]);
};

template <typename T, typename S> // not very good
class Normal : public AnalyticalF1D<T,S> {
private:
  T _sigma;
  T _mu;
public:
  Normal(T sigma=1., T mu=1.);
  bool operator() (T psi[], const S rho[]);
};

} // end namespace olb

#endif
