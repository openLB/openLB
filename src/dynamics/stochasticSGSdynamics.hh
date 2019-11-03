/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Mathias J. Krause, Jonas Latt
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
 * MRT Dynamics with adjusted omega -- generic implementation.
 */
#ifndef STOCHASTIC_SGS_DYNAMICS_HH
#define STOCHASTIC_SGS_DYNAMICS_HH

#include <algorithm>
#include <limits>
#include "stochasticSGSdynamics.h"
#include "mrtDynamics.h"
#include "mrtHelpers.h"
#include "core/cell.h"
#include "core/util.h"
#include "math.h"

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>


using namespace std;
namespace olb {

/// Implementation of the Stochastic relaxation based on
/// " A stochastic subgrid model with application to turbulent flow and scalar mixing"; Phys. of Fluids 19; 2007

////////////////////// Class StochasticsSGSdynamics //////////////////////////

/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */



template<typename T, typename DESCRIPTOR>
StochasticSGSdynamics<T,DESCRIPTOR>::StochasticSGSdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_, T turbulenceInt_, T charU_, T smagoConst_, T dx_, T dt_)
  : MRTdynamics<T,DESCRIPTOR>(omega_, momenta_),
    turbulenceInt(turbulenceInt_),
    smagoConst(smagoConst_),
    charU(charU_),
    preFactor(computePreFactor(omega_,smagoConst_) )

{

  //    T invM_S_SGS[DESCRIPTOR::q][DESCRIPTOR::q];
  //    T rtSGS[DESCRIPTOR::q]; // relaxation times vector for SGS approach.
  //    for (int iPop  = 0; iPop < DESCRIPTOR::q; ++iPop)
  //    {
  //      rtSGS[iPop] = DESCRIPTOR::S[iPop];
  //    }
  //    for (int iPop  = 0; iPop < DESCRIPTOR::shearIndexes; ++iPop)
  //    {
  //      rtSGS[DESCRIPTOR::shearViscIndexes[iPop]] = omega;
  //    }
  //    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop)
  //    {
  //      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop)
  //      {
  //        invM_S_SGS[iPop][jPop] = T();
  //        for (int kPop = 0; kPop < DESCRIPTOR::q; ++kPop)
  //        {
  //          if (kPop == jPop)
  //          {
  //            invM_S_SGS[iPop][jPop] += DESCRIPTOR::invM[iPop][kPop] *
  //                                  rtSGS[kPop];
  //            cout << "wert"<<iPop <<jPop << "= "<<  invM_S_SGS[iPop][jPop]<< endl;
  //          }
  //        }
  //      }
  //    }
  //


}


template<typename T, typename DESCRIPTOR>
void StochasticSGSdynamics<T,DESCRIPTOR>::collide(
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];



  T drift  = computeTimeScale(preFactor, rho, pi, smagoConst, X_lang_n);
  T result = getRandBMTrans(cell, turbulenceInt, charU);

  //  cout << "vor neu setzen: "<<X_lang_n<< endl;
  X_lang_n = getRandomWalk(cell, drift, result);
  // cout << "nach neu setzen: "<<X_lang_n<< endl;
  //cout << X_lang_n<< endl;
  // cout << drift<< endl;
  // cout << result<< endl;


  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi, X_lang_n);


  T invM_S_SGS[DESCRIPTOR::q][DESCRIPTOR::q];
  T rtSGS[DESCRIPTOR::q]; // relaxation times vector for SGS approach.
  for (int iPop  = 0; iPop < DESCRIPTOR::q; ++iPop) {
    rtSGS[iPop] = DESCRIPTOR::S[iPop];
  }
  for (int iPop  = 0; iPop < DESCRIPTOR::shearIndexes; ++iPop) {
    rtSGS[DESCRIPTOR::shearViscIndexes[iPop]] = newOmega;
  }
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
      invM_S_SGS[iPop][jPop] = T();
      for (int kPop = 0; kPop < DESCRIPTOR::q; ++kPop) {
        if (kPop == jPop) {
          invM_S_SGS[iPop][jPop] += DESCRIPTOR::invM[iPop][kPop] *
                                    rtSGS[kPop];
          //cout << "wert"<<iPop <<jPop << "= "<<  invM_S_SGS[iPop][jPop]<< endl;
        }
      }
    }
  }

  T uSqr = mrtHelpers<T,DESCRIPTOR>::mrtSGSCollision(cell, rho, u, newOmega, invM_S_SGS);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
void StochasticSGSdynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  this->setOmega(omega);
  preFactor = computePreFactor(omega, smagoConst);
}

template<typename T, typename DESCRIPTOR>
T StochasticSGSdynamics<T,DESCRIPTOR>::getSmagorinskyOmega(Cell<T,DESCRIPTOR>& cell, T X_lang_n )
{
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi, X_lang_n);
  return newOmega;
}



template<typename T, typename DESCRIPTOR>
T StochasticSGSdynamics<T,DESCRIPTOR>::getRandBMTrans(
  Cell<T,DESCRIPTOR>& cell,
  T turbulenceInt, T CharU )
{
  /// Random number generator based on Box Müller transform to produuce random normal
  /// distributed numbers with zero mean and

  T mean = 0.;
  T TKE_ini = 1.5*turbulenceInt*turbulenceInt*charU*charU;
  T velStDev =  sqrt(2./3.*TKE_ini);
  static double n2 = 0.0;
  static int n2_cached = 0;
  if (!n2_cached) {
    double x, y, r;
    do {
      x = 2.0*rand()/RAND_MAX - 1;
      y = 2.0*rand()/RAND_MAX - 1;

      r = x*x + y*y;
    } while ( util::nearZero(r) || r > 1.0);
    {
      double d = sqrt(-2.0*log(r)/r);
      double n1 = x*d;
      n2 = y*d;
      double result = n1*velStDev + mean;
      n2_cached = 1;
      return result;
    }
  } else {
    n2_cached = 0;
    return n2*velStDev + mean;
  }
}


/// Create Random walk
template<typename T, typename DESCRIPTOR>
T StochasticSGSdynamics<T,DESCRIPTOR>::getRandomWalk(
  Cell<T,DESCRIPTOR>& cell,
  T drift, T result)
{
  /// initialisation of model standard variation, see Pope pp 484
  T sigma = 2.3;
  X_lang_n *=drift;

  X_lang_n +=  sigma*sqrt(drift*2)*result;

  return X_lang_n;


}
/// set random walk

// template<typename T, typename DESCRIPTOR>
// void StochasticSGSdynamics<T,DESCRIPTOR>::setRandomWalk(
//   Cell<T,DESCRIPTOR>& cell,
//   T CharU, T drift, T result )
// {
//   /// initialisation of model standard variation, see Pope pp 484
//   T X_lang_n = getRandomWalk(cell, CharU, drift, result);
// }


// /// get time sclae
template<typename T, typename DESCRIPTOR>
T StochasticSGSdynamics<T,DESCRIPTOR>::computeTimeScale(
  T preFactor, T rho, T pi[util::TensorVal<DESCRIPTOR >::n], T smagoConst, T X_lang_n  )
{
  T Const = 0.2;
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);

  /// SGS dissipation is calcualted directly withou any filter size due to effeiciency in tau
  // for post processing this has to be evaluated seperately with S_ij³
  T diss_corr = smagoConst*smagoConst*PiNeqNorm*PiNeqNorm*PiNeqNorm*(1+ X_lang_n);

  T tau= Const*pow(( 1. / diss_corr ), 1./3.);
  T drift = 1./tau;
  /// deterministic drift time scale T_L see Pope pp. 484
  return drift;

}


// // // /// set timescale
// template<typename T, typename DESCRIPTOR>
// void StochasticSGSdynamics<T,DESCRIPTOR>::setTimeScale(
//   T preFactor, T rho, T pi[util::TensorVal<DESCRIPTOR >::n], T smagoConst, T X_lang_n)
// {
//   T drift = computeTimeScale(preFactor, rho, pi, smagoConst, X_lang_n);
// }


template<typename T, typename DESCRIPTOR>
T StochasticSGSdynamics<T,DESCRIPTOR>::computePreFactor(T omega, T smagoConst)
{
  return (T)smagoConst*smagoConst*descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>()*2*sqrt(2);
}

template<typename T, typename DESCRIPTOR>
T StochasticSGSdynamics<T,DESCRIPTOR>::computeOmega(T omega0, T preFactor, T rho, T pi[util::TensorVal<DESCRIPTOR >::n], T X_lang_n)
{
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  /// Molecular realaxation time
  T tau_mol = 1. /omega0;
  /// Turbulent realaxation time
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol+(preFactor*tau_eff*PiNeqNorm*(1+X_lang_n)))-tau_mol);
  /// Effective realaxation time
  tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;
}

}

#endif
