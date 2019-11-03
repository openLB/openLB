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

#ifndef UNITCONVERTER_H
#define UNITCONVERTER_H


#include <math.h>
#include "io/ostreamManager.h"
#include "core/util.h"
#include "io/xmlReader.h"

// known design issues
//    1. How can we prevent abuse of constructur by mixing up parameters?
//    2. physical problems may have different names for viscosity, e.g. diffusity,  temperature conductivity
//    4. Feedback about stability or comment the chosen discretization
//    5. Explain why Desctiptor as template
//    6. Is it worth to introduce invConversionDensity to avoid division


/// All OpenLB code is contained in this namespace.
namespace olb {



/** Conversion between physical and lattice units, as well as discretization.
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
*  \param physDensity           physical density in __kg / m^3__
*  - - -
*  \param conversionLength      conversion factor for length __m__
*  \param conversionTime        conversion factor for time __s__
*  \param conversionMass        conversion factor for mass __kg__
*  - - -
*  \param conversionVelocity    conversion velocity __m / s__
*  \param conversionViscosity   conversion kinematic viscosity __m^2 / s__
*  \param conversionDensity     conversion density __kg / m^3__
*  \param conversionForce       conversion force __kg m / s^2__
*  \param conversionPressure    conversion pressure __kg / m s^2__
*  - - -
*  \param resolution            number of grid points per charPhysLength
*  - - -
*  \param charLatticeVelocity
*/
template <typename T, typename DESCRIPTOR>
class UnitConverter {
public:
  /** Documentation of constructor:
    *  \param physDeltaX              spacing between two lattice cells in __m__
    *  \param physDeltaT              time step in __s__
    *  \param charPhysLength          reference/characteristic length of simulation geometry in __m__
    *  \param charPhysVelocity        maximal or highest expected velocity during simulation in __m / s__
    *  \param physViscosity           physical kinematic viscosity in __m^2 / s__
    *  \param physDensity             physical density in __kg / m^3__
    *  \param charPhysPressure        reference/characteristic physical pressure in Pa = kg / m s^2
    */
  constexpr UnitConverter( T physDeltaX, T physDeltaT, T charPhysLength, T charPhysVelocity,
                           T physViscosity, T physDensity, T charPhysPressure = 0 )
    : _conversionLength(physDeltaX),
      _conversionTime(physDeltaT),
      _conversionVelocity(_conversionLength / _conversionTime),
      _conversionDensity(physDensity),
      _conversionMass( _conversionDensity * pow(_conversionLength, 3) ),
      _conversionViscosity(_conversionLength * _conversionLength / _conversionTime),
      _conversionForce( _conversionMass * _conversionLength / (_conversionTime * _conversionTime) ),
      _conversionPressure( _conversionForce / pow(_conversionLength, 2) ),
      _charPhysLength(charPhysLength),
      _charPhysVelocity(charPhysVelocity),
      _physViscosity(physViscosity),
      _physDensity(physDensity),
      _charPhysPressure(charPhysPressure),
      _resolution((int)(_charPhysLength / _conversionLength + 0.5)),
      _latticeRelaxationTime( (_physViscosity / _conversionViscosity * descriptors::invCs2<T,DESCRIPTOR>()) + 0.5 ),
      _charLatticeVelocity( _charPhysVelocity / _conversionVelocity ),
      clout(std::cout,"UnitConverter")
  {
  }

  virtual ~UnitConverter() = default;

  /// return resolution
  constexpr int getResolution(  ) const
  {
    return _resolution;
  }
  /// return relaxation time in lattice units
  constexpr T getLatticeRelaxationTime(  ) const
  {
    return _latticeRelaxationTime;
  }
  /// return relaxation frequency in lattice units
  constexpr T getLatticeRelaxationFrequency(  ) const
  {
    return 1./_latticeRelaxationTime;
  }
  /// return relaxation frequency in lattice units computed from given physical diffusivity in __m^2 / s__
  template <typename DESCRIPTOR_>
  constexpr T getLatticeRelaxationFrequencyFromDiffusivity( const T physDiffusivity ) const
  {
    T latticeDiffusivity = physDiffusivity / _conversionViscosity;
    return 1.0 / ( latticeDiffusivity * descriptors::invCs2<T,DESCRIPTOR_>() + 0.5 );
  }
  /// return characteristic length in physical units
  constexpr T getCharPhysLength(  ) const
  {
    return _charPhysLength;
  }
  /// return characteristic velocity in physical units
  constexpr T getCharPhysVelocity(  ) const
  {
    return _charPhysVelocity;
  }
  /// return characteristic velocity in lattice units
  constexpr T getCharLatticeVelocity(  ) const
  {
    return _charLatticeVelocity;
  }
  /// return viscosity in physical units
  constexpr T getPhysViscosity(  ) const
  {
    return _physViscosity;
  }
  /// return density in physical units
  constexpr T getPhysDensity(  ) const
  {
    return _physDensity;
  }
  /// return characteristic pressure in physical units
  constexpr T getCharPhysPressure(  ) const
  {
    return _charPhysPressure;
  }
  /// return Reynolds number
  constexpr T getReynoldsNumber(  ) const
  {
    return _charPhysVelocity * _charPhysLength / _physViscosity;
  }
  /// return Mach number
  constexpr T getMachNumber(  ) const
  {
    return getCharLatticeVelocity() * std::sqrt(descriptors::invCs2<T,DESCRIPTOR>());
  }
  /// return Knudsen number
  constexpr T getKnudsenNumber(  ) const
  {
    // This calculates the lattice Knudsen number.
    // See e.g. (7.22) in "The Lattice Boltzmann Method: Principles and Practice" [kruger2017lattice].
    return getMachNumber() / getReynoldsNumber();
  }
  /// conversion from lattice to  physical length
  constexpr T getPhysLength( int latticeLength ) const
  {
    return _conversionLength * latticeLength;
  }
  /// conversion from physical to lattice length, returns number of voxels for given physical length
  constexpr int getLatticeLength( T physLength ) const
  {
    return int( physLength / _conversionLength + 0.5 );
  }
  /// access (read-only) to private member variable
  constexpr T getConversionFactorLength() const
  {
    return _conversionLength;
  }
  /// returns grid spacing (voxel length) in __m__
  constexpr T getPhysDeltaX() const
  {
    return _conversionLength;
  }

  /// conversion from lattice to  physical time
  constexpr T getPhysTime( int latticeTime ) const
  {
    return _conversionTime * latticeTime;
  }
  /// conversion from physical to lattice time
  constexpr int getLatticeTime( T physTime ) const
  {
    return int(physTime / _conversionTime + 0.5);
  }
  /// access (read-only) to private member variable
  constexpr T getConversionFactorTime() const
  {
    return _conversionTime;
  }
  /// returns time spacing (timestep length) in __s__
  constexpr T getPhysDeltaT() const
  {
    return _conversionTime;
  }

  /// conversion from lattice to  physical velocity
  constexpr T getPhysVelocity( T latticeVelocity ) const
  {
    return _conversionVelocity * latticeVelocity;
  }
  /// conversion from physical to lattice velocity
  constexpr T getLatticeVelocity( T physVelocity ) const
  {
    return physVelocity / _conversionVelocity;
  }
  /// access (read-only) to private member variable
  constexpr T getConversionFactorVelocity() const
  {
    return _conversionVelocity;
  }

  /// conversion from lattice to  physical density
  constexpr T getPhysDensity( T latticeDensity ) const
  {
    return _conversionDensity * latticeDensity;
  }
  /// conversion from physical to lattice density
  constexpr T getLatticeDensity( T physDensity ) const
  {
    return physDensity / _conversionDensity;
  }
  constexpr T getLatticeDensityFromPhysPressure( T physPressure ) const
  {
    T latticePressure = getLatticePressure( physPressure );
    return util::densityFromPressure<T,DESCRIPTOR>( latticePressure);
  }
  /// access (read-only) to private member variable
  constexpr T getConversionFactorDensity() const
  {
    return _conversionDensity;
  }

  /// conversion from lattice to  physical mass
  constexpr T getPhysMass( T latticeMass ) const
  {
    return _conversionMass * latticeMass;
  }
  /// conversion from physical to lattice mass
  constexpr T getLatticeMass( T physMass ) const
  {
    return physMass / _conversionMass;
  }
  /// access (read-only) to private member variable
  constexpr T getConversionFactorMass() const
  {
    return _conversionMass;
  }

  /// conversion from lattice to  physical viscosity
  constexpr T getPhysViscosity( T latticeViscosity ) const
  {
    return _conversionViscosity * latticeViscosity;
  }
  /// conversion from physical to lattice viscosity
  constexpr T getLatticeViscosity(  ) const
  {
    return _physViscosity / _conversionViscosity;
  }
  /// access (read-only) to private member variable
  constexpr T getConversionFactorViscosity() const
  {
    return _conversionViscosity;
  }

  /// conversion from lattice to  physical force
  constexpr T getPhysForce( T latticeForce ) const
  {
    return _conversionForce * latticeForce;
  }
  /// conversion from physical to lattice force
  constexpr T getLatticeForce( T physForce ) const
  {
    return physForce / _conversionForce;
  }
  /// access (read-only) to private member variable
  constexpr T getConversionFactorForce() const
  {
    return _conversionForce;
  }

  /// conversion from lattice to  physical pressure
  constexpr T getPhysPressure( T latticePressure ) const
  {
    return _conversionPressure * latticePressure + _charPhysPressure;
  }
  /// conversion from physical to lattice pressure
  constexpr T getLatticePressure( T physPressure ) const
  {
    return ( physPressure - _charPhysPressure ) / _conversionPressure;
  }
  /// access (read-only) to private member variable
  constexpr T getConversionFactorPressure() const
  {
    return _conversionPressure;
  }
  /// nice terminal output for conversion factors, characteristical and physical data
  virtual void print() const;

  void write(std::string const& fileName = "unitConverter") const;

protected:
  // conversion factors
  const T _conversionLength;      // m
  const T _conversionTime;        // s
  const T _conversionVelocity;    // m / s
  const T _conversionDensity;     // kg / m^3
  const T _conversionMass;        // kg
  const T _conversionViscosity;   // m^2 / s
  const T _conversionForce;       // kg m / s^2
  const T _conversionPressure;    // kg / m s^2

  // physical units, e.g characteristic or reference values
  const T _charPhysLength;        // m
  const T _charPhysVelocity;      // m / s
  const T _physViscosity;         // m^2 / s
  const T _physDensity;           // kg / m^3
  const T _charPhysPressure;      // kg / m s^2

  // lattice units, discretization parameters
  const int _resolution;
  const T _latticeRelaxationTime;
  const T _charLatticeVelocity;   //

private:
  mutable OstreamManager clout;
};

/// creator function with data given by a XML file
template <typename T, typename DESCRIPTOR>
UnitConverter<T, DESCRIPTOR>* createUnitConverter(XMLreader const& params);

template <typename T, typename DESCRIPTOR>
class UnitConverterFromResolutionAndRelaxationTime : public UnitConverter<T, DESCRIPTOR> {
public:
  constexpr UnitConverterFromResolutionAndRelaxationTime(
    int resolution,
    T latticeRelaxationTime,
    T charPhysLength,
    T charPhysVelocity,
    T physViscosity,
    T physDensity,
    T charPhysPressure = 0) : UnitConverter<T, DESCRIPTOR>(
        (charPhysLength/resolution),
        (latticeRelaxationTime - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * pow((charPhysLength/resolution),2) / physViscosity,
        charPhysLength,
        charPhysVelocity,
        physViscosity,
        physDensity,
        charPhysPressure)
  {
  }
};

template <typename T, typename DESCRIPTOR>
class UnitConverterFromResolutionAndLatticeVelocity : public UnitConverter<T, DESCRIPTOR> {
public:
  constexpr UnitConverterFromResolutionAndLatticeVelocity(
    int resolution,
    T charLatticeVelocity,
    T charPhysLength,
    T charPhysVelocity,
    T physViscosity,
    T physDensity,
    T charPhysPressure = 0) : UnitConverter<T, DESCRIPTOR>(
        (charPhysLength/resolution),
        (charLatticeVelocity / charPhysVelocity * charPhysLength / resolution),
        charPhysLength,
        charPhysVelocity,
        physViscosity,
        physDensity,
        charPhysPressure)
  {
  }
};
template <typename T, typename DESCRIPTOR>
class UnitConverterFromRelaxationTimeAndLatticeVelocity : public UnitConverter<T, DESCRIPTOR> {
public:
  constexpr UnitConverterFromRelaxationTimeAndLatticeVelocity(
    T latticeRelaxationTime,
    T charLatticeVelocity,
    T charPhysLength,
    T charPhysVelocity,
    T physViscosity,
    T physDensity,
    T charPhysPressure = 0) : UnitConverter<T, DESCRIPTOR>(
        (physViscosity * charLatticeVelocity / charPhysVelocity * descriptors::invCs2<T,DESCRIPTOR>() / (latticeRelaxationTime - 0.5)),
        (physViscosity * charLatticeVelocity * charLatticeVelocity / charPhysVelocity / charPhysVelocity * descriptors::invCs2<T,DESCRIPTOR>() / (latticeRelaxationTime - 0.5)),
        charPhysLength,
        charPhysVelocity,
        physViscosity,
        physDensity,
        charPhysPressure)
  {
  }
};

}  // namespace olb

#include "thermalUnitConverter.h"
#endif
