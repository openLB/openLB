/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017-2018 Max Gaedtke, Albert Mink, Davide Dapelo
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

#ifndef PL_UNITCONVERTER_H
#define PL_UNITCONVERTER_H


#include <math.h>
#include "io/ostreamManager.h"
#include "core/util.h"
#include "io/xmlReader.h"
#include "unitConverter.h"


/// All OpenLB code is contained in this namespace.
namespace olb {



/** Conversion between physical and lattice units, as well as discretization specialized for power-law rheology: \f$\nu=m^{n-1}\f$, with \f$m$ being the consistency coefficient and \f$n\f$ the power-law index.
 *  The Newtonian case is recovered for \f$n=1\f$.
 *  A characteristic (apparent) kinematic viscosity for problem similarity (e.g., to compute the Reynolds number) is obtained as __physViscosity = poweLawConsistencyCoeff * [physVelocity / (2*physLength)]^(n-1)__.
* Be aware of the nomenclature:
* We distingish between physical (dimensioned) and lattice (dimensionless) values.
* A specific conversion factor maps the two different scopes,
* e.g. __physLength = conversionLength * latticeLength__
*
* For pressure and temperature we first shift the physical values by a characteristic value to asure a lattice pressure and lattice temperature between 0 and 1, e.g. __physPressure - charPhysPressure = conversionPressure * latticePressure__
*
*  \param latticeRelaxationTime   relaxation time, have to be greater than 0.5!
*  - - -
*  \param physViscosity         physical kinematic viscosity in __m^2 / s__
*  \param poweLawIndex          power-law index
*  \param physConsistencyCoeff  physicsl consistency coefficient __m^2 s^(n-2)__
*  \param physDensity           physical density in __kg / m^3__
*  - - -
*  \param conversionLength      conversion factor for length __m__
*  \param conversionTime        conversion factor for time __s__
*  \param conversionMass        conversion factor for mass __kg__
*  - - -
*  \param conversionVelocity    conversion velocity __m / s__
*  \param conversionViscosity   conversion kinematic viscosity __m^2 / s__
*  \param conversionConsistencyCoeff conversion power-law consistency coefficient __m^2 s^(n-2)__
*  \param conversionDensity     conversion density __kg / m^3__
*  \param conversionForce       conversion force __kg m / s^2__
*  \param conversionPressure    conversion pressure __kg / m s^2__
*  - - -
*  \param resolution            number of grid points per charPhysLength
*  - - -
*  \param charLatticeVelocity
*/
template <typename T, typename DESCRIPTOR>
class PowerLawUnitConverter : public UnitConverter<T, DESCRIPTOR> {
public:
  /** Documentation of constructor:
    *  \param physDeltaX              spacing between two lattice cells in __m__
    *  \param physDeltaT              time step in __s__
    *  \param charPhysLength          reference/characteristic length of simulation geometry in __m__
    *  \param charPhysVelocity        maximal or highest expected velocity during simulation in __m / s__
    *  \param physConsistencyCoeff    physical power-law consistency coefficient in __m^2 s^(n-2)__
    *  \param powerLawIndex           Power-law index
    *  \param physDensity             physical density in __kg / m^3__
    *  \param charPhysPressure        reference/characteristic physical pressure in Pa = kg / m s^2
    */
  constexpr PowerLawUnitConverter( T physDeltaX, T physDeltaT, T charPhysLength,
                                   T charPhysVelocity, T physConsistencyCoeff, T powerLawIndex,
                                   T physDensity, T charPhysPressure = 0 )
    : UnitConverter<T, DESCRIPTOR>( physDeltaX, physDeltaT, charPhysLength, charPhysVelocity,
                                 physConsistencyCoeff * pow(charPhysVelocity / (2*charPhysLength), powerLawIndex-1),
                                 physDensity, charPhysPressure ),
      _conversionConsistencyCoeff(pow( physDeltaT,powerLawIndex-2 ) * pow( physDeltaX,2 ) ),
      _powerLawIndex(powerLawIndex),
      _physConsistencyCoeff(physConsistencyCoeff),
      clout(std::cout,"PowerLawUnitConverter")
  {
  }

  /// return consistency coefficient in physical units
  constexpr T getPhysConsistencyCoeff( ) const
  {
    return _physConsistencyCoeff;
  }
  /// conversion from lattice to  physical consistency coefficient
  constexpr T getPhysConsistencyCoeff( T latticeConsistencyCoeff ) const
  {
    return _conversionConsistencyCoeff * latticeConsistencyCoeff;
  }
  /// conversion from physical to lattice consistency coefficient
  constexpr T getLatticeConsistencyCoeff(  ) const
  {
    return _physConsistencyCoeff / _conversionConsistencyCoeff;
  }
  /// access (read-only) to private member variable
  constexpr T getConversionFactorConsistencyCoeff() const
  {
    return _conversionConsistencyCoeff;
  }

  /// access (read-only) to private member variable
  constexpr T getPowerLawIndex() const
  {
    return _powerLawIndex;
  }

  /// nice terminal output for conversion factors, characteristical and physical data
  virtual void print() const;

  void write(std::string const& fileName = "unitConverter") const;

protected:
  // conversion factors
  const T _conversionConsistencyCoeff;         // m^2 s^(n-2)

  // physical units, e.g characteristic or reference values
  const T _powerLawIndex;
  const T _physConsistencyCoeff;         // m^2 s^(n-2)

private:
  mutable OstreamManager clout;
};

/// creator function with data given by a XML file
template <typename T, typename DESCRIPTOR>
PowerLawUnitConverter<T, DESCRIPTOR>* createPowerLawUnitConverter(XMLreader const& params);

template <typename T, typename DESCRIPTOR>
class PowerLawUnitConverterFrom_Resolution_RelaxationTime_Reynolds_PLindex : public PowerLawUnitConverter<T, DESCRIPTOR> {
public:
  constexpr PowerLawUnitConverterFrom_Resolution_RelaxationTime_Reynolds_PLindex(
    int resolution,
    T latticeRelaxationTime,
    T charPhysLength,
    T charPhysVelocity,
    T Re,
    T powerLawIndex,
    T physDensity,
    T charPhysPressure = 0)
    : PowerLawUnitConverter<T, DESCRIPTOR>( (charPhysLength/resolution),
                                         (latticeRelaxationTime - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * pow((charPhysLength/resolution),2) / ( ( charPhysLength * charPhysVelocity * pow( charPhysVelocity / ( 2 * charPhysLength ), 1 - powerLawIndex ) / Re ) * pow( charPhysVelocity / (2 * charPhysLength ), powerLawIndex - 1 ) ),
                                         charPhysLength, charPhysVelocity,
                                         charPhysLength * charPhysVelocity * pow( charPhysVelocity / ( 2 * charPhysLength ), 1 - powerLawIndex ) / Re, powerLawIndex, physDensity, charPhysPressure )
  {
  }

};

}  // namespace olb

#endif
