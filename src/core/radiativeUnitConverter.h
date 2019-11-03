/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Max Gaedtke, Albert Mink
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
 * Unit conversion handling -- header file.
 */

#ifndef RADIATIVEUNITCONVERTER_H
#define RADIATIVEUNITCONVERTER_H


#include "io/ostreamManager.h"
#include "core/unitConverter.h"


/// All OpenLB code is contained in this namespace.
namespace olb {

double getThetaRefracted(double const& thetaIncident, double const& refractiveRelative);
double getFresnelFunction(double const& theta, double const& refractiveRelative);
double R_phi_diff (double const& theta, double const& refractiveRelative);
double R_j_diff (double const& theta, double const& refractiveRelative);
double getRefractionFunction(const double& refractiveRelative);
double getRefractionFunction(const double& refractiveRelative);
double getPartialBBCoefficient(double const& latticeDiffusionCoefficient, double const& relativeRefractiveIndex );

// forward declaration
template<typename T, typename DESCRIPTOR> class RadiativeUnitConverter;

// wrapper for above function
template <typename T, typename DESCRIPTOR>
double getPartialBBCoefficient(RadiativeUnitConverter<T,DESCRIPTOR> const& converter)
{
  return getPartialBBCoefficient( converter.getLatticeDiffusion(), converter.getRefractiveRelative() );
};

// wrapper for above function
template <typename T, typename DESCRIPTOR>
double getRefractionFunction(RadiativeUnitConverter<T,DESCRIPTOR> const& converter)
{
  return getRefractionFunction(converter.getRefractiveRelative());
};

/** Conversion between physical and lattice units, as well as discretization.
* Be aware of the nomenclature:
* We distingish between physical (dimensioned) and lattice (dimensionless) values.
* A specific conversion factor maps the two different scopes,
* e.g. __physLength = conversionLength * latticeLength__
*
*/
template <typename T, typename DESCRIPTOR>
class RadiativeUnitConverter : public UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> {
public:
  /** Documentation of constructor:
    *   \param resolution   is number of voxel per 1 meter
    *   \param latticeRelaxationTime    see class UnitConverterFromResolutionAndRelaxationTime
    *   \param physAbsorption   physical absorption in 1/meter
    *   \param physScattering   physical scattering in 1/meter
    */
  constexpr RadiativeUnitConverter( int resolution, T latticeRelaxationTime, T physAbsorption, T physScattering, T anisotropyFactor=0, T charPhysLength=1, T refractiveMedia=1, T refractiveAmbient=1 )
    : UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>( resolution, latticeRelaxationTime, charPhysLength, T(1), T(1), T(1) ),
      clout(std::cout, "RadiativeUnitConverter"),
      _physAbsorption(physAbsorption),
      _physScattering(physScattering),
      _anisotropyFactor(anisotropyFactor),
      _extinction( physAbsorption+physScattering ),
      _scatteringAlbedo( physScattering/(physAbsorption+physScattering) ),
      _physDiffusion( 1.0 / (3.0*(physAbsorption+physScattering)) ),
      _effectiveAttenuation( std::sqrt(3*physAbsorption*(physAbsorption+physScattering)) ),
      _refractiveRelative(refractiveMedia/refractiveAmbient),
      _latticeAbsorption( physAbsorption*this->getConversionFactorLength() ),
      _latticeScattering( physScattering*this->getConversionFactorLength() ),
      _latticeDiffusion(_physDiffusion/this->getConversionFactorLength())
  { };

  constexpr T getPhysAbsorption() const
  {
    return _physAbsorption;
  };

  constexpr T getPhysScattering() const
  {
    return _physScattering;
  };

  constexpr T getAnisotropyFactor() const
  {
    return _anisotropyFactor;
  };

  constexpr T getExtinction() const
  {
    return _extinction;
  };

  constexpr T getScatteringAlbedo() const
  {
    return _scatteringAlbedo;
  };

  constexpr T getPhysDiffusion() const
  {
    return _physDiffusion;
  };

  constexpr T getEffectiveAttenuation() const
  {
    return _effectiveAttenuation;
  };

  constexpr T getLatticeAbsorption() const
  {
    return _latticeAbsorption;
  };

  constexpr T getLatticeScattering() const
  {
    return _latticeScattering;
  };

  constexpr T getLatticeDiffusion() const
  {
    return _latticeDiffusion;
  };

  constexpr T getRefractiveRelative() const
  {
    return _refractiveRelative;
  };

  void print() const override;

private:
  mutable OstreamManager clout;

  double _physAbsorption;
  double _physScattering;
  double _anisotropyFactor;
  double _extinction;
  double _scatteringAlbedo;
  double _physDiffusion;
  double _effectiveAttenuation;

  double _refractiveRelative;

  double _latticeAbsorption;
  double _latticeScattering;
  double _latticeDiffusion;
};

template <typename T, typename DESCRIPTOR>
void RadiativeUnitConverter<T, DESCRIPTOR>::print() const
{
  clout << "----------------- UnitConverter information -----------------" << std::endl;
  clout << "-- Parameters:" << std::endl;
  clout << "Resolution:                       N=              " << this->getResolution() << std::endl;
  clout << "Lattice relaxation frequency:     omega=          " << this->getLatticeRelaxationFrequency() << std::endl;
  clout << "Lattice relaxation time:          tau=            " << this->getLatticeRelaxationTime() << std::endl;
  clout << "Characteristical length(m):       charL=          " << this->getCharPhysLength() << std::endl;
  clout << "Phys. density(kg/m^d):            charRho=        " << this->getPhysDensity() << std::endl;
  clout << "Phys. absorption(1/m):            mu_a=           " << getPhysAbsorption() << std::endl;
  clout << "Phys. scattering(1/m):            mu_s=           " << getPhysScattering() << std::endl;
  clout << "Extinction(1/m):                  mu_t=           " << getExtinction() << std::endl;
  clout << "Effective attenuation(1/m):       mu_eff=         " << getEffectiveAttenuation() << std::endl;
  clout << "Phys. diffusion(m):               D=              " << getPhysDiffusion() << std::endl;
  clout << "Single scattering albedo:         albedo=         " << getScatteringAlbedo() << std::endl;
  clout << "Anisotropy factor:                g=              " << getAnisotropyFactor() << std::endl;

  clout << std::endl;
  clout << "Lattice diffusion:                D^*=            " << getLatticeDiffusion() << std::endl;
  clout << "Lattice absorption:               absorption=     " << getLatticeAbsorption() << std::endl;
  clout << "Lattice scattering:               scattering=     " << getLatticeScattering() << std::endl;
  clout << "Lattice sink:                     sink=           " << 3./8.*getLatticeAbsorption()*(getLatticeScattering()+getLatticeAbsorption()) << std::endl;
  clout << "C_R: " << getRefractionFunction(getRefractiveRelative()) << std::endl;
  clout << "r_F: " << getPartialBBCoefficient(getLatticeDiffusion(),getRefractiveRelative()) << std::endl;

  clout << std::endl;
  clout << "-- Conversion factors:" << std::endl;
  clout << "Voxel length(m):                  physDeltaX=     " << this->getConversionFactorLength() << std::endl;
  clout << "Time step(s):                     physDeltaT=     " << this->getConversionFactorTime() << std::endl;
  clout << "Density factor(kg/m^3):           physDensity=    " << this->getConversionFactorDensity() <<  std::endl;
  clout << "-------------------------------------------------------------" << std::endl;

}


}  // namespace olb

#endif
