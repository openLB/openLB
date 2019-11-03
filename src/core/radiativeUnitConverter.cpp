/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Albert Mink, Marc Haussmann
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
 * Function to extract refractive properties needed for boundary modeling.
 */


#include "core/radiativeUnitConverter.h"

// definition required only by cygwin
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/// All OpenLB code is contained in this namespace.
namespace olb {



/// Documentation of implemented functions found in 5.2.2 Biomedical Optics, Principles and Imaging; Wang 2007

double getThetaRefracted(double const& thetaIncident, double const& refractiveRelative)
{
  double thetaRefracted = M_PI/2.;
  if( refractiveRelative * sin(thetaIncident) < 1 ) {
    thetaRefracted = asin( refractiveRelative * sin(thetaIncident));  // eq.(5.118)
  }
  return thetaRefracted;
};

double getFresnelFunction(double const& theta, double const& refractiveRelative)
{
  double thetaRefracted = getThetaRefracted(theta, refractiveRelative);
  double rf_1 = 0.5 * pow((refractiveRelative * cos(thetaRefracted) - cos(theta)) /
                          (refractiveRelative * cos(thetaRefracted) + cos(theta)), 2.);
  double rf_2 = 0.5 * pow((refractiveRelative * cos(theta) - cos(thetaRefracted)) /
                          (refractiveRelative * cos(theta) + cos(thetaRefracted)), 2.);
  return rf_1 + rf_2;   // eq.(5.115)
};

double R_phi_diff (double const& theta, double const& refractiveRelative)
{
  return 2. * sin(theta) * cos(theta) * getFresnelFunction(theta,refractiveRelative);
};

double R_j_diff (double const& theta, double const& refractiveRelative)
{
  return 3. * sin(theta) * pow(cos(theta),2.) * getFresnelFunction(theta,refractiveRelative);
};

double getRefractionFunction(const double& refractiveRelative)
{
  int N = 10000.0;
  double h = (M_PI / 2.) /double(N);
  double R_phi = 0.0;
  double R_j = 0.0;
  for (int i = 0; i < N; i++) {
    R_phi += h*(R_phi_diff(0.5*h + h*i,refractiveRelative));
    R_j   += h*(R_j_diff  (0.5*h + h*i,refractiveRelative));
  }
  double R_eff = (R_phi + R_j) / (2 - R_phi + R_j);     // eq.(5.112)
  return (1 + R_eff) / (1 - R_eff);                     // eq.(5.111)    C_R = (1 + R_eff) / (1 - R_eff);
};

double getPartialBBCoefficient(double const& latticeDiffusionCoefficient, double const& relativeRefractiveIndex )
{
  double C_R = getRefractionFunction( relativeRefractiveIndex );
  return 2 - 2/(4*latticeDiffusionCoefficient*C_R +1);
};

}  // namespace olb
