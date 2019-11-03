/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Julius Woerner, Albert Mink
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
 * discretization of anisotrphic phasefunction
 * DOI 10.1016/j.ijheatmasstransfer.2011.11.009
 *
 */


#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>
#ifdef FEATURE_OPENBLAS
#include <openblas/lapacke.h>
#endif

#include "utilities/vectorHelpers.h"

namespace olb {

/** Function to compute Henyey Greenstein phase funtion
 * \param cosTheta cos(theta) of scattering event with net direction theta
 * \param g anisotropy factor
 */
double henyeyGreenstein(double cosTheta, double g)
{
  return (1-g*g) / std::pow(1+g*g-2*g*cosTheta ,1.5);
}

// check for energy conservation
template< int q, typename DESCRIPTOR>
std::array<double,q> testEnergyConservationColumn( const std::array<std::array<double,q>,q>&  phi )
{
  std::array<double,q> discIntegral;
  for( int iVec = 0; iVec < q; iVec++ ) {
    discIntegral[iVec] = 0;
    for( int iSum = 0; iSum < q; iSum++ ) {
      discIntegral[iVec] += descriptors::t<double,DESCRIPTOR>(iSum) * phi[iSum][iVec];
    }
  }
  return discIntegral;
}

// check for energy conservation
template<int q, typename DESCRIPTOR>
std::array<double,q> testEnergyConservationRow( const std::array<std::array<double,q>,q>& phi )
{
  std::array<double,q> discIntegral;
  for( int iVec = 0; iVec < q; iVec++ ) {
    discIntegral[iVec] = 0;
    for( int iSum = 0; iSum < q; iSum++ ) {
      discIntegral[iVec] += descriptors::t<double,DESCRIPTOR>(iSum) * phi[iVec][iSum];
    }
  }
  return discIntegral;
}

// check for anisotropy conservation
template< int q>
std::array<double,q> testAnisotropyConservationColumn( const std::array<std::array<double,q>,q>&  phi,
  const double weights[q], std::array<std::array<double,q>,q>& cosTheta)
{
  std::array<double,q> discIntegral;
  for( int iVec = 0; iVec < q; iVec++ ) {
    discIntegral[iVec] = 0;
    for( int iSum = 0; iSum < q; iSum++ ) {
      discIntegral[iVec] += weights[iSum] * cosTheta[iVec][iSum]* phi[iVec][iSum];
    }
  }
  return discIntegral;
}


/** Computes vector [a, a+stepsize, a+2*stepsize, ..., b]
 * \param stepsize
 * \param a
 * \param b
 */
std::vector<double> linespace( double const stepsize, double const a, double const b )
{
  std::vector<double> linspace{}; // initalize to empty
  if( util::nearZero( a-b ) ) {
    return linspace;
  }
  int const h = (int) ( (b-a) / stepsize);
  for (int n = 0; n <= h; n++) {
    linspace.push_back( a +double(n)/h * b );
  }
  return linspace;
}

/**
 * \tparam DESCRIPTOR a lattice stencil
 * \tparam q number of lattice stencils MINUS one, 0th direction not needed
 * \param stepsize choose precision
 * \param anisotropyFactor also known as g
 * \param solution for debugging
 * \param phi is anisotropy matrix(sym), element [i][j] probability of scattering from i into j (direction)
 * \param breakAfter to limit endless iteration for testing
 */
template<typename DESCRIPTOR>
void computeAnisotropyMatrix( double const stepsize, double const anisotropyFactor,
  double solution[(DESCRIPTOR::q-1)*((DESCRIPTOR::q-1)+1)/2],
  std::array<std::array<double,DESCRIPTOR::q-1>, DESCRIPTOR::q-1>& phi, int const breakAfter = -1)
{
  using namespace descriptors;

  OstreamManager clout( std::cout, "AnisotropyMatrix" );
  clout << "Call AnistorpiyMatrix()" << std::endl;
#ifdef FEATURE_OPENBLAS
  clout << "Compute anisotropy matrix ..." << std::endl;
  typedef DESCRIPTOR L;
  int const q = L::q-1;
  int const mm = 2*q;
  int const nn = (q+1)*q/2;

  // qxq matrix 0th row and 0th column are empty
  std::vector<std::vector<double>> angleProb;
  angleProb.resize(q);
  for ( int n = 0; n < q; n++ ) {
    angleProb[n].resize(q);
  }
  double dotProduct;
  double normI;
  double normJ;
  for (int iPop=0; iPop<q; iPop++) {
    for (int jPop=0; jPop<q; jPop++) {
      // shift by 1 due to notation in DESCRIPTOR/DESCRIPTOR
      // exclude 0th direction
      dotProduct = c<L>(iPop+1)*c<L>(jPop+1);
      normI = std::sqrt( c<L>(iPop+1)*c<L>(iPop+1) );
      normJ = std::sqrt( c<L>(jPop+1)*c<L>(jPop+1) );
      angleProb[iPop][jPop] = dotProduct / (normI*normJ);
    }
  }

  for (int i=0; i<q; i++) {
    for (int j=0; j<q; j++) {
      phi[i][j]=1.0;
    }
  }

  for( int i = 0; i < nn; i++ ) {
    solution[i] = 0;
  };

  int iter;
  double U[mm][nn];
  double U_array[nn*mm];
  std::vector<double> anisoIterVector = linespace( stepsize, 0, anisotropyFactor );

  // additional condition only for unit testing
  size_t index = 0;
  for ( ; index < anisoIterVector.size() && index != (std::size_t)(breakAfter); index++)
  {
    // wipe matrices and vectors
    for (int m = 0; m < mm; m++) {
      for (int n = 0; n < nn; n++) {
        U[m][n] = 0;
      }
    }

    // compute matrix U, iter handels symmetry
    iter = 0;
    for (int m=0; m<q; m++) {
      for (int n=m; n<q; n++) {
        U[m][iter]   = t<double,L>(m+1)*phi[n][m];
        U[n][iter]   = t<double,L>(n+1)*phi[m][n];
        U[m+q][iter] = t<double,L>(m+1)*phi[n][m]*angleProb[n][m];
        U[n+q][iter] = t<double,L>(n+1)*phi[m][n]*angleProb[m][n];
        iter++;
      }
    }

    double sum1;
    double sum2;
    // compute vector b
    for (int n=0; n<q; n++) {
      // get sum
      sum1 = 0;
      sum2 = 0;
      for (int k=0; k<q; k++) {
        sum1 += t<double,L>(k+1)*phi[k][n];
        sum2 += t<double,L>(k+1)*phi[k][n]*angleProb[k][n];
      }
      solution[n] = 1 - sum1;
      solution[q+n] = anisoIterVector[index] - sum2;
    }

    // transform 2d array to 1d array, column major
    for ( int n = 0; n < nn; ++n) {
      for ( int m = 0; m < mm; ++m ) {
        U_array[n*mm +m] = U[m][n];
      }
    }

    //util::print(b, nn, "b");

    // solve Ax = QRx = R'Q'x = b
    // equiv x = Q*R'\x
    LAPACKE_dgels( LAPACK_COL_MAJOR, 'N', mm,
                   nn, 1, U_array,
                   mm, solution, nn);
    //util::print(b, nn, "b after");

    iter = 0;
    for (int m=0; m<q; m++) {
      for (int n=m; n<q; n++) {
        phi[n][m] = phi[n][m]*(1+solution[iter]);
        phi[m][n] = phi[n][m];
        iter++;
      }
    }
    //util::print( testEnergyConservationColumn<q>(phi,L::t), "colum sum elementwise" );
    //util::print( testEnergyConservationRow<q>(phi,L::t), "row sum elementwise" );
    //util::print( testEnergyConservationSec<q>(phi,L::t,angleProb), "second Moment" );
  }
  clout << "Terminate after " << index << " iterations" << std::endl;
#endif
}


template<typename DESCRIPTOR>
void computeAnisotropyMatrixKimAndLee( double const anisotropyFactor,
  std::array<std::array<double,DESCRIPTOR::q>,DESCRIPTOR::q>& phi )
{
  OstreamManager clout( std::cout, "AnisotropyMatrix_KimAndLee" );
  clout << "Compute anisotropy matrix ..." << std::endl;
  typedef DESCRIPTOR L;
  int const q = L::q;

  // qxq matrix 0th row and 0th column are empty
  std::array<std::array<double,q>, q> cosTheta;
  double dotProduct;
  double normI;
  double normJ;
  for (int iPop=0; iPop<q; iPop++) {
    for (int jPop=0; jPop<q; jPop++) {
      dotProduct = descriptors::c<L>(iPop) * descriptors::c<L>(jPop);
      normI = std::sqrt( util::normSqr<int,3>(descriptors::c<L>(iPop)) );
      normJ = std::sqrt( util::normSqr<int,3>(descriptors::c<L>(jPop)) );
      cosTheta[iPop][jPop] = dotProduct /(normI*normJ);
      if( util::normSqr<int,3>(descriptors::c<L>(iPop)) == 0 || util::normSqr<int,3>(descriptors::c<L>(jPop)) == 0){
        cosTheta[iPop][jPop] = 0.0;
      }
    }
  }

  std::array<std::array<double,q>, q> phaseFunction;
  for (int m=0; m<q; m++) {
    for (int n=0; n<q; n++) {
      phaseFunction[m][n] = henyeyGreenstein(cosTheta[m][n], anisotropyFactor);
    }
  }

  for (int m=0; m<q; m++) {
    for (int n=0; n<q; n++) {
      dotProduct = 0;
      for (int i = 0; i < q; ++i) {
        dotProduct += phaseFunction[m][i]*descriptors::t<double,L>(i);
      }
      phi[m][n] = phaseFunction[m][n] / dotProduct;
    }
  }
  //util::print( testEnergyConservationColumn<q>(phi,L::t), "colum sum elementwise" );
  //util::print( testEnergyConservationRow<q>(phi,L::t), "row sum elementwise" );
  //util::print( testAnisotropyConservationColumn<q>(phi,L::t,cosTheta), "anisotropyConservation Moment" );
}


} // namespace olb
