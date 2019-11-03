/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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
 * Helper functions for the implementation of LB dynamics. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef LB_HELPERS_H
#define LB_HELPERS_H

#include "latticeDescriptors.h"
#include "core/cell.h"
#include "core/util.h"


namespace olb {


// Forward declarations
template<typename T, class DESCRIPTORBASE> struct lbDynamicsHelpers;
template<typename T, typename DESCRIPTOR> struct lbExternalHelpers;
template<typename T, typename DESCRIPTOR> struct lbLatticeHelpers;

/// This structure forwards the calls to the appropriate helper class
template<typename T, typename DESCRIPTOR>
struct lbHelpers {

  static T equilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], const T uSqr)
  {
    return lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>
           ::equilibrium(iPop, rho, u, uSqr);
  }

  static T equilibriumFirstOrder(int iPop, T rho, const T u[DESCRIPTOR::d])
  {
    return lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>
           ::equilibriumFirstOrder(iPop, rho, u);
  }

  static T incEquilibrium(int iPop, const T j[DESCRIPTOR::d], const T jSqr, const T pressure)
  {
    return lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>
           ::incEquilibrium(iPop, j, jSqr, pressure);
  }

  static void computeFneq ( Cell<T,DESCRIPTOR> const& cell,
                            T fNeq[DESCRIPTOR::q], T rho, const T u[DESCRIPTOR::d] )
  {
    lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>::computeFneq(cell, fNeq, rho, u);
  }

  static T bgkCollision(Cell<T,DESCRIPTOR>& cell, T const& rho, const T u[DESCRIPTOR::d], T const& omega)
  {
    return lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>
           ::bgkCollision(cell, rho, u, omega);
  }

  static T bgkCollision(Cell<T,DESCRIPTOR>& cell, T const& rho, const T u[DESCRIPTOR::d], const T omega[DESCRIPTOR::q])
  {
    const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] *= (T)1-omega[iPop];
      cell[iPop] += omega[iPop] * lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>::equilibrium (iPop, rho, u, uSqr );
    }
    return uSqr;
  }

  static T rlbCollision( Cell<T, DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T omega )
  {
    return lbDynamicsHelpers<T, typename DESCRIPTOR::BaseDescriptor>
           ::rlbCollision( cell, rho, u, omega );
  }

  static T mrtCollision( Cell<T, DESCRIPTOR>& cell, T const& rho, const T u[DESCRIPTOR::d], T invM_S[DESCRIPTOR::q][DESCRIPTOR::q] ) {

    return lbDynamicsHelpers<T, typename DESCRIPTOR::BaseDescriptor>
           ::mrtCollision(cell, rho, u, invM_S );
  }

  static T incBgkCollision(Cell<T,DESCRIPTOR>& cell, T pressure, const T j[DESCRIPTOR::d], T omega)
  {
    return lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>
           ::incBgkCollision(cell, pressure, j, omega);
  }

  static T constRhoBgkCollision(Cell<T,DESCRIPTOR>& cell,
                                T rho, const T u[DESCRIPTOR::d], T ratioRho, T omega)
  {
    return lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>
           ::constRhoBgkCollision(cell, rho, u, ratioRho, omega);
  }

  static T computeRho(Cell<T,DESCRIPTOR> const& cell)
  {
    return lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>
           ::computeRho(cell);
  }

  static void computeJ(Cell<T,DESCRIPTOR> const& cell, T j[DESCRIPTOR::d] )
  {
    lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>
    ::computeJ(cell, j);
  }

  static void computeRhoU(Cell<T,DESCRIPTOR> const& cell, T& rho, T u[DESCRIPTOR::d])
  {
    lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>
    ::computeRhoU(cell, rho, u);
  }

  static void computeFeq(Cell<T,DESCRIPTOR> const& cell, T fEq[DESCRIPTOR::q])
  {
    T rho{};
    T u[2] {};
    computeRhoU(cell, rho, u);
    const T uSqr = u[0]*u[0] + u[1]*u[1];
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      fEq[iPop] = equilibrium(iPop, rho, u, uSqr);
    }
  }

  static void computeFneq(Cell<T,DESCRIPTOR> const& cell, T fNeq[DESCRIPTOR::q])
  {
    T rho{};
    T u[2] {};
    computeRhoU(cell, rho, u);
    computeFneq(cell, fNeq, rho, u);
  }

  static void computeRhoJ(Cell<T,DESCRIPTOR> const& cell, T& rho, T j[DESCRIPTOR::d])
  {
    lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>
    ::computeRhoJ(cell, rho, j);
  }

  static void computeStress(Cell<T,DESCRIPTOR> const& cell, T rho, const T u[DESCRIPTOR::d],
                            T pi[util::TensorVal<DESCRIPTOR >::n] )
  {
    lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>
    ::computeStress(cell, rho, u, pi);
  }

  static void computeAllMomenta(Cell<T,DESCRIPTOR> const& cell, T& rho, T u[DESCRIPTOR::d],
                                T pi[util::TensorVal<DESCRIPTOR >::n] )
  {
    lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>
    ::computeAllMomenta(cell, rho, u, pi);
  }

  static T computePiNeqNormSqr(Cell<T,DESCRIPTOR> const& cell)
  {
    //return lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>
    //::computePiNeqNormSqr(cell);
    T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
    computeAllMomenta(cell, rho, u, pi);
    T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
    if (util::TensorVal<DESCRIPTOR >::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];
    }
    return PiNeqNormSqr;
  }

  static T computeForcedPiNeqNormSqr(Cell<T,DESCRIPTOR> const& cell)
  {
    //return lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>
    //::computeForcedPiNeqNormSqr(cell);
    T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
    computeAllMomenta(cell, rho, u, pi);
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
    T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
    if (util::TensorVal<DESCRIPTOR >::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];
    }
    return PiNeqNormSqr;
  }

  static void modifyVelocity(Cell<T,DESCRIPTOR>& cell, const T newU[DESCRIPTOR::d])
  {
    lbDynamicsHelpers<T,typename DESCRIPTOR::BaseDescriptor>
    ::modifyVelocity(cell, newU);
  }

  static void addExternalForce(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d], T omega, T amplitude=(T)1)
  {
    lbExternalHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, omega, amplitude);
  }

  static void swapAndStream2D(Cell<T,DESCRIPTOR> **grid, int iX, int iY)
  {
    lbLatticeHelpers<T,DESCRIPTOR>::swapAndStream2D(grid, iX, iY);
  }

  static void swapAndStream3D(Cell<T,DESCRIPTOR> ***grid, int iX, int iY, int iZ)
  {
    lbLatticeHelpers<T,DESCRIPTOR>::swapAndStream3D(grid, iX, iY, iZ);
  }

};  // struct lbHelpers


/// All helper functions are inside this structure
template<typename T, typename DESCRIPTORBASE>
struct lbDynamicsHelpers {
  /// Computation of equilibrium distribution, second order in u
  static T equilibrium(int iPop, T rho, const T u[DESCRIPTORBASE::d], const T uSqr)
  {
    T c_u = T();
    for (int iD=0; iD < DESCRIPTORBASE::d; ++iD) {
      c_u += descriptors::c<DESCRIPTORBASE>(iPop,iD)*u[iD];
    }
    return rho * descriptors::t<T,DESCRIPTORBASE>(iPop) * ( (T)1 + descriptors::invCs2<T,DESCRIPTORBASE>() * c_u
                                              + descriptors::invCs2<T,DESCRIPTORBASE>() * descriptors::invCs2<T,DESCRIPTORBASE>() * (T)0.5 * c_u *c_u
                                              - descriptors::invCs2<T,DESCRIPTORBASE>() * (T)0.5 * uSqr )
                                              - descriptors::t<T,DESCRIPTORBASE>(iPop);
  }

  /// Computation of equilibrium distribution, first order in u
  static T equilibriumFirstOrder(int iPop, T rho, const T u[DESCRIPTORBASE::d])
  {
    T c_u = T();
    for (int iD=0; iD < DESCRIPTORBASE::d; ++iD) {
      c_u += descriptors::c<DESCRIPTORBASE>(iPop,iD)*u[iD];
    }
    return rho * descriptors::t<T,DESCRIPTORBASE>(iPop) * ( (T)1 + c_u * descriptors::invCs2<T,DESCRIPTORBASE>() ) - descriptors::t<T,DESCRIPTORBASE>(iPop);
  }

  static T incEquilibrium( int iPop, const T j[DESCRIPTORBASE::d],
                           const T jSqr, const T pressure )
  {
    T c_j = T();
    for (int iD=0; iD < DESCRIPTORBASE::d; ++iD) {
      c_j += descriptors::c<DESCRIPTORBASE>(iPop,iD)*j[iD];
    }

    return descriptors::t<T,DESCRIPTORBASE>(iPop) * ( descriptors::invCs2<T,DESCRIPTORBASE>() * pressure +
                                   descriptors::invCs2<T,DESCRIPTORBASE>() * c_j +
                                   descriptors::invCs2<T,DESCRIPTORBASE>() * descriptors::invCs2<T,DESCRIPTORBASE>()/(T)2 * c_j*c_j -
                                   descriptors::invCs2<T,DESCRIPTORBASE>()/(T)2 * jSqr
                                 ) - descriptors::t<T,DESCRIPTORBASE>(iPop);
  }

  /// Computation of non-equilibrium distribution
  static void computeFneq(CellBase<T,DESCRIPTORBASE> const& cell, T fNeq[DESCRIPTORBASE::q], T rho, const T u[DESCRIPTORBASE::d])
  {
    const T uSqr = util::normSqr<T,DESCRIPTORBASE::d>(u);
    for (int iPop=0; iPop < DESCRIPTORBASE::q; ++iPop) {
      fNeq[iPop] = cell[iPop] - equilibrium(iPop, rho, u, uSqr);
    }
  }

  /// BGK collision step
  static T bgkCollision(CellBase<T,DESCRIPTORBASE>& cell, T const& rho, const T u[DESCRIPTORBASE::d], T const& omega)
  {
    const T uSqr = util::normSqr<T,DESCRIPTORBASE::d>(u);
    for (int iPop=0; iPop < DESCRIPTORBASE::q; ++iPop) {
      cell[iPop] *= (T)1-omega;
      cell[iPop] += omega * lbDynamicsHelpers<T,DESCRIPTORBASE>::equilibrium(iPop, rho, u, uSqr );
    }
    return uSqr;
  }

  /// Incompressible BGK collision step
  static T incBgkCollision(CellBase<T,DESCRIPTORBASE>& cell, T pressure, const T j[DESCRIPTORBASE::d], T omega)
  {
    const T jSqr = util::normSqr<T,DESCRIPTORBASE::d>(j);
    for (int iPop=0; iPop < DESCRIPTORBASE::q; ++iPop) {
      cell[iPop] *= (T)1-omega;
      cell[iPop] += omega * lbDynamicsHelpers<T,DESCRIPTORBASE>::incEquilibrium (
                      iPop, j, jSqr, pressure );
    }
    return jSqr;
  }

  /// BGK collision step with density correction
  static T constRhoBgkCollision(CellBase<T,DESCRIPTORBASE>& cell, T rho, const T u[DESCRIPTORBASE::d], T ratioRho, T omega)
  {
    const T uSqr = util::normSqr<T,DESCRIPTORBASE::d>(u);
    for (int iPop=0; iPop < DESCRIPTORBASE::q; ++iPop) {
      T feq = lbDynamicsHelpers<T,DESCRIPTORBASE>::equilibrium(iPop, rho, u, uSqr );
      cell[iPop] =
        ratioRho*(feq+descriptors::t<T,DESCRIPTORBASE>(iPop))-descriptors::t<T,DESCRIPTORBASE>(iPop) +
        ((T)1-omega)*(cell[iPop]-feq);
    }
    return uSqr;
  }

  /// RLB advection diffusion collision step
  static T rlbCollision(CellBase<T, DESCRIPTORBASE>& cell, T rho, const T u[DESCRIPTORBASE::d], T omega )
  {
    const T uSqr = util::normSqr<T, DESCRIPTORBASE::d>( u );
    // First-order moment for the regularization
    T j1[DESCRIPTORBASE::d];
    for ( int iD = 0; iD < DESCRIPTORBASE::d; ++iD ) {
      j1[iD] = T();
    }

    T fEq[DESCRIPTORBASE::q];
    for ( int iPop = 0; iPop < DESCRIPTORBASE::q; ++iPop ) {
      fEq[iPop] = lbDynamicsHelpers<T, DESCRIPTORBASE>::equilibriumFirstOrder( iPop, rho, u );
      for ( int iD = 0; iD < DESCRIPTORBASE::d; ++iD ) {
        j1[iD] += descriptors::c<DESCRIPTORBASE>(iPop,iD) * ( cell[iPop] - fEq[iPop] );
      }
    }

    // Collision step
    for ( int iPop = 0; iPop < DESCRIPTORBASE::q; ++iPop ) {
      T fNeq = T();
      for ( int iD = 0; iD < DESCRIPTORBASE::d; ++iD ) {
        fNeq += descriptors::c<DESCRIPTORBASE>(iPop,iD) * j1[iD];
      }
      fNeq *= descriptors::t<T,DESCRIPTORBASE>(iPop) * descriptors::invCs2<T,DESCRIPTORBASE>();
      cell[iPop] = fEq[iPop] + ( (T)1 - omega ) * fNeq;
    }
    return uSqr;
  }


///// Computation of all equilibrium distribution (in momenta space)
  static void computeMomentaEquilibrium( T momentaEq[DESCRIPTORBASE::q], T rho, const T u[DESCRIPTORBASE::d], T uSqr )
  {
    for (int iPop = 0; iPop < DESCRIPTORBASE::q; ++iPop) {
      momentaEq[iPop] = T();
      for (int jPop = 0; jPop < DESCRIPTORBASE::q; ++jPop) {
        momentaEq[iPop] += DESCRIPTORBASE::M[iPop][jPop] *
                (lbDynamicsHelpers<T, DESCRIPTORBASE>::equilibrium(jPop,rho,u,uSqr) + DESCRIPTORBASE::t[jPop]);
      }
    }
  }

  static void computeMomenta(T momenta[DESCRIPTORBASE::q], CellBase<T, DESCRIPTORBASE>& cell)
  {
    for (int iPop = 0; iPop < DESCRIPTORBASE::q; ++iPop) {
      momenta[iPop] = T();
      for (int jPop = 0; jPop < DESCRIPTORBASE::q; ++jPop) {
        momenta[iPop] += DESCRIPTORBASE::M[iPop][jPop] *
                         (cell[jPop] + DESCRIPTORBASE::t[jPop]);
      }
    }
  }

  static T mrtCollision( CellBase<T,DESCRIPTORBASE>& cell, T const& rho, const T u[DESCRIPTORBASE::d], T invM_S[DESCRIPTORBASE::q][DESCRIPTORBASE::q] )
  {
    //// Implemented in advectionDiffusionMRTlbHelpers2D.h and advectionDiffusionMRTlbHelpers3D.h
    T uSqr = util::normSqr<T, DESCRIPTORBASE::d>(u);
    T momenta[DESCRIPTORBASE::q];
    T momentaEq[DESCRIPTORBASE::q];

    computeMomenta(momenta, cell);
    computeMomentaEquilibrium(momentaEq, rho, u, uSqr);

//    std::cout << "momenta = ";
//    for (int i=0; i < DESCRIPTORBASE::q; ++i) {
//        std::cout << momenta[i] << ", ";
//    }
//    std::cout << std::endl;

//    std::cout << "momentaEq = ";
//    for (int i=0; i < DESCRIPTORBASE::q; ++i) {
//        std::cout << momentaEq[i] << ", ";
//    }
//    std::cout << std::endl;

    for (int iPop = 0; iPop < DESCRIPTORBASE::q; ++iPop) {
      T collisionTerm = T();
      for (int jPop = 0; jPop < DESCRIPTORBASE::q; ++jPop) {
        collisionTerm += invM_S[iPop][jPop] * (momenta[jPop] - momentaEq[jPop]);
      }
      cell[iPop] -= collisionTerm;
    }
    return uSqr;
  }

  /// Computation of density
  static T computeRho(CellBase<T,DESCRIPTORBASE> const& cell)
  {
    T rho = T();
    for (int iPop=0; iPop < DESCRIPTORBASE::q; ++iPop) {
      rho += cell[iPop];
    }
    rho += (T)1;
    return rho;
  }

  /// Computation of momentum
  static void computeJ(CellBase<T,DESCRIPTORBASE> const& cell, T j[DESCRIPTORBASE::d])
  {
    for (int iD=0; iD < DESCRIPTORBASE::d; ++iD) {
      j[iD] = T();
    }
    for (int iPop=0; iPop < DESCRIPTORBASE::q; ++iPop) {
      for (int iD=0; iD < DESCRIPTORBASE::d; ++iD) {
        j[iD] += cell[iPop]*descriptors::c<DESCRIPTORBASE>(iPop,iD);
      }
    }
  }

  /// Computation of hydrodynamic variables
  static void computeRhoU(CellBase<T,DESCRIPTORBASE> const& cell, T& rho, T u[DESCRIPTORBASE::d])
  {
    rho = T();
    for (int iD=0; iD < DESCRIPTORBASE::d; ++iD) {
      u[iD] = T();
    }
    for (int iPop=0; iPop < DESCRIPTORBASE::q; ++iPop) {
      rho += cell[iPop];
      for (int iD=0; iD < DESCRIPTORBASE::d; ++iD) {
        u[iD] += cell[iPop]*descriptors::c<DESCRIPTORBASE>(iPop,iD);
      }
    }
    rho += (T)1;
    for (int iD=0; iD < DESCRIPTORBASE::d; ++iD) {
      u[iD] /= rho;
    }
  }

  /// Computation of hydrodynamic variables
  static void computeRhoJ(CellBase<T,DESCRIPTORBASE> const& cell, T& rho, T j[DESCRIPTORBASE::d])
  {
    rho = T();
    for (int iD=0; iD < DESCRIPTORBASE::d; ++iD) {
      j[iD] = T();
    }
    for (int iPop=0; iPop < DESCRIPTORBASE::q; ++iPop) {
      rho += cell[iPop];
      for (int iD=0; iD < DESCRIPTORBASE::d; ++iD) {
        j[iD] += cell[iPop]*descriptors::c<DESCRIPTORBASE>(iPop,iD);
      }
    }
    rho += (T)1;
  }

  /// Computation of stress tensor
  static void computeStress(CellBase<T,DESCRIPTORBASE> const& cell, T rho, const T u[DESCRIPTORBASE::d],
                            T pi[util::TensorVal<DESCRIPTORBASE>::n] )
  {
    int iPi = 0;
    for (int iAlpha=0; iAlpha < DESCRIPTORBASE::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta < DESCRIPTORBASE::d; ++iBeta) {
        pi[iPi] = T();
        for (int iPop=0; iPop < DESCRIPTORBASE::q; ++iPop) {
          pi[iPi] += descriptors::c<DESCRIPTORBASE>(iPop,iAlpha)*
                     descriptors::c<DESCRIPTORBASE>(iPop,iBeta) * cell[iPop];
        }
        // stripe off equilibrium contribution
        pi[iPi] -= rho*u[iAlpha]*u[iBeta];
        if (iAlpha==iBeta) {
          pi[iPi] -= 1./descriptors::invCs2<T,DESCRIPTORBASE>()*(rho-(T)1);
        }
        ++iPi;
      }
    }
  }

  /// Computation of all hydrodynamic variables
  static void computeAllMomenta(CellBase<T,DESCRIPTORBASE> const& cell, T& rho, T u[DESCRIPTORBASE::d],
                                T pi[util::TensorVal<DESCRIPTORBASE>::n] )
  {
    computeRhoU(cell, rho, u);
    computeStress(cell, rho, u, pi);
  }

  /*
  /// Computes squared norm of non-equilibrium part of 2nd momentum
  static T computePiNeqNormSqr(CellBase<T,DESCRIPTORBASE> const& cell)
  {
    T rho, u[DESCRIPTORBASE::d], pi[util::TensorVal<DESCRIPTORBASE >::n];
    computeAllMomenta(cell, rho, u, pi);
    T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
    if (util::TensorVal<DESCRIPTORBASE >::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];
    }
    return PiNeqNormSqr;
  }

  /// Computes squared norm of forced non-equilibrium part of 2nd momentum
  static T computeForcedPiNeqNormSqr(CellBase<T,DESCRIPTORBASE> const& cell)
  {
    T rho, u[DESCRIPTORBASE::d], pi[util::TensorVal<DESCRIPTORBASE >::n];
    computeAllMomenta(cell, rho, u, pi);
    //Creation of body force tensor (rho/2.)*(G_alpha*U_beta + U_alpha*G_Beta)
    T ForceTensor[util::TensorVal<DESCRIPTORBASE >::n];
    int iPi = 0;
    for (int Alpha=0; Alpha<DESCRIPTORBASE::d; ++Alpha) {
      for (int Beta=Alpha; Beta<DESCRIPTORBASE::d; ++Beta) {
        ForceTensor[iPi] = rho/2.*(cell[DESCRIPTORBASE::template index<descriptors::FORCE>()][Alpha]*u[Beta] + u[Alpha]*cell[DESCRIPTORBASE::template index<descriptors::FORCE>()][Beta]);
        ++iPi;
      }
    }
    // Creation of second-order moment off-equilibrium tensor
    for (int iPi=0; iPi < util::TensorVal<DESCRIPTORBASE >::n; ++iPi){
      pi[iPi] += ForceTensor[iPi];
    }
    T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
    if (util::TensorVal<DESCRIPTORBASE >::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];
    }
    return PiNeqNormSqr;
  }
  */

  static void modifyVelocity(CellBase<T,DESCRIPTORBASE>& cell, const T newU[DESCRIPTORBASE::d])
  {
    T rho, oldU[DESCRIPTORBASE::d];
    computeRhoU(cell, rho, oldU);
    const T oldUSqr = util::normSqr<T,DESCRIPTORBASE::d>(oldU);
    const T newUSqr = util::normSqr<T,DESCRIPTORBASE::d>(newU);
    for (int iPop=0; iPop<DESCRIPTORBASE::q; ++iPop) {
      cell[iPop] = cell[iPop]
                   - equilibrium(iPop, rho, oldU, oldUSqr)
                   + equilibrium(iPop, rho, newU, newUSqr);
    }
  }

};  // struct lbDynamicsHelpers

/// Helper functions for dynamics that access external field
template<typename T, typename DESCRIPTOR>
struct lbExternalHelpers {
  /// Add a force term after BGK collision
  static void addExternalForce(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d], T omega, T amplitude)
  {
    const T* force = cell.template getFieldPointer<descriptors::FORCE>();
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      T c_u = T();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
      }
      c_u *= descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>();
      T forceTerm = T();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        forceTerm +=
          (   ((T)descriptors::c<DESCRIPTOR>(iPop,iD)-u[iD]) * descriptors::invCs2<T,DESCRIPTOR>()
              + c_u * descriptors::c<DESCRIPTOR>(iPop,iD)
          )
          * force[iD];
      }
      forceTerm *= descriptors::t<T,DESCRIPTOR>(iPop);
      forceTerm *= T(1) - omega/T(2);
      forceTerm *= amplitude;
      cell[iPop] += forceTerm;
    }
  }
};  // struct externalFieldHelpers

/// Helper functions with full-lattice access
template<typename T, typename DESCRIPTOR>
struct lbLatticeHelpers {
  /// Swap ("bounce-back") values of a cell (2D), and apply streaming step
  static void swapAndStream2D(Cell<T,DESCRIPTOR> **grid, int iX, int iY)
  {
    const int half = DESCRIPTOR::q/2;
    for (int iPop=1; iPop<=half; ++iPop) {
      int nextX = iX + descriptors::c<DESCRIPTOR>(iPop,0);
      int nextY = iY + descriptors::c<DESCRIPTOR>(iPop,1);
      T fTmp                   = grid[iX][iY][iPop];
      grid[iX][iY][iPop]       = grid[iX][iY][iPop+half];
      grid[iX][iY][iPop+half]  = grid[nextX][nextY][iPop];
      grid[nextX][nextY][iPop] = fTmp;
    }
  }

  /// Swap ("bounce-back") values of a cell (3D), and apply streaming step
  static void swapAndStream3D(Cell<T,DESCRIPTOR> ***grid,
                              int iX, int iY, int iZ)
  {
    const int half = DESCRIPTOR::q/2;
    for (int iPop=1; iPop<=half; ++iPop) {
      int nextX = iX + descriptors::c<DESCRIPTOR>(iPop,0);
      int nextY = iY + descriptors::c<DESCRIPTOR>(iPop,1);
      int nextZ = iZ + descriptors::c<DESCRIPTOR>(iPop,2);
      T fTmp                          = grid[iX][iY][iZ][iPop];
      grid[iX][iY][iZ][iPop]          = grid[iX][iY][iZ][iPop+half];
      grid[iX][iY][iZ][iPop+half]     = grid[nextX][nextY][nextZ][iPop];
      grid[nextX][nextY][nextZ][iPop] = fTmp;
    }
  }
};

/// All boundary helper functions are inside this structure
template<typename T, typename DESCRIPTOR, int direction, int orientation>
struct BoundaryHelpers {
  static void computeStress (
    Cell<T,DESCRIPTOR> const& cell, T rho, const T u[DESCRIPTOR::d],
    T pi[util::TensorVal<DESCRIPTOR >::n] )
  {
    typedef DESCRIPTOR L;
    const T uSqr = util::normSqr<T,L::d>(u);

    std::vector<int> const& onWallIndices = util::subIndex<L, direction, 0>();
    std::vector<int> const& normalIndices = util::subIndex<L, direction, orientation>();

    T fNeq[DESCRIPTOR::q];
    for (unsigned fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
      int iPop = onWallIndices[fIndex];
      if (iPop == 0) {
        fNeq[0] = T();  // fNeq[0] will not be used anyway
      } else {
        fNeq[iPop] =
          cell[iPop] -
          lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
      }
    }
    for (unsigned fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
      int iPop = normalIndices[fIndex];
      fNeq[iPop] =
        cell[iPop] -
        lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
    }

    int iPi = 0;
    for (int iAlpha=0; iAlpha<L::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta<L::d; ++iBeta) {
        pi[iPi] = T();
        for (unsigned fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
          const int iPop = onWallIndices[fIndex];
          pi[iPi] +=
            descriptors::c<L>(iPop,iAlpha)*descriptors::c<L>(iPop,iBeta)*fNeq[iPop];
        }
        for (unsigned fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
          const int iPop = normalIndices[fIndex];
          pi[iPi] += (T)2 * descriptors::c<L>(iPop,iAlpha)*descriptors::c<L>(iPop,iBeta)*
                     fNeq[iPop];
        }
        ++iPi;
      }
    }
  }

};  // struct boundaryHelpers

}  // namespace olb

// The specialized code is directly included. That is because we never want
// it to be precompiled so that in both the precompiled and the
// "include-everything" version, the compiler can apply all the
// optimizations it wants.
#include "lbHelpersD2Q5.h"
#include "lbHelpersD2Q9.h"
#include "lbHelpersD3Q7.h"
#include "lbHelpersD3Q15.h"
#include "lbHelpersD3Q19.h"
#include "lbHelpersD3Q27.h"

#endif
