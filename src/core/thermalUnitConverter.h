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

#ifndef THERMALUNITCONVERTER_H
#define THERMALUNITCONVERTER_H


#include <math.h>
#include "io/ostreamManager.h"
#include "core/util.h"
#include "io/xmlReader.h"
#include "core/unitConverter.h"

/// All OpenLB code is contained in this namespace.
namespace olb {



/** Conversion between physical and lattice units, as well as discretization specialized for thermal applications with boussinesq approximation.
* Be aware of the nomenclature:
* We distingish between physical (dimensioned) and lattice (dimensionless) values.
* A specific conversion factor maps the two different scopes,
* e.g. __physLength = conversionLength * latticeLength__
*
* For pressure and temperature we first shift the physical values by a characteristic value to asure a lattice pressure and between 0 and 1, e.g. __physPressure - charPhysPressure = conversionPressure * latticePressure__. For the temperature we set lattice values between 0.5 and 1.5 by __latticeTemperature = (physTemperature - charPhysLowTemperature) / conversionTemperature + 0.5 with conversionTemperature = charPhysHighTemperature - charPhysLowTemperature = charPhysTemperatureDifference
*
* TODO: Extend documentation for ThermalUnitConverter
*/
template <typename T, typename DESCRIPTOR, typename ThermalLattice>
class ThermalUnitConverter : public UnitConverter<T, DESCRIPTOR> {
public:
  /** Documentation of constructor:
    * TODO: Extend constructur documentation
    */
  constexpr ThermalUnitConverter(
    T physDeltaX,
    T physDeltaT,
    T charPhysLength,
    T charPhysVelocity,
    T physViscosity,
    T physDensity,
    T physThermalConductivity,
    T physSpecificHeatCapacity,
    T physThermalExpansionCoefficient,
    T charPhysLowTemperature,
    T charPhysHighTemperature,
    T charPhysPressure = 0 )
    : UnitConverter<T, DESCRIPTOR>(
        physDeltaX, physDeltaT, charPhysLength, charPhysVelocity,
        physViscosity, physDensity, charPhysPressure),
      _conversionTemperature(charPhysHighTemperature - charPhysLowTemperature),
      _conversionThermalDiffusivity(this->_conversionViscosity),
      _conversionSpecificHeatCapacity(this->_conversionVelocity * this->_conversionVelocity / _conversionTemperature),
      _conversionThermalConductivity(this->_conversionForce / this->_conversionTime / _conversionTemperature),
      _conversionHeatFlux(this->_conversionMass / pow(this->_conversionTime, 3)),
      _charPhysLowTemperature(charPhysLowTemperature),
      _charPhysHighTemperature(charPhysHighTemperature),
      _charPhysTemperatureDifference(charPhysHighTemperature - charPhysLowTemperature),
      _physThermalExpansionCoefficient(physThermalExpansionCoefficient),
      _physThermalDiffusivity(physThermalConductivity / (physDensity * physSpecificHeatCapacity)),
      _physSpecificHeatCapacity(physSpecificHeatCapacity),
      _physThermalConductivity(physThermalConductivity),
      _latticeThermalRelaxationTime( (_physThermalDiffusivity / _conversionThermalDiffusivity * descriptors::invCs2<T,ThermalLattice>()) + 0.5 ),
      clout(std::cout,"ThermalUnitConv")
  {
  };

  /// return thermal relaxation time in lattice units
  constexpr T getLatticeThermalRelaxationTime(  ) const
  {
    return _latticeThermalRelaxationTime;
  };
  /// return thermal relaxation frequency in lattice units
  constexpr T getLatticeThermalRelaxationFrequency(  ) const
  {
    return 1.0 / _latticeThermalRelaxationTime;
  };

  /// return characteristic low temperature in physical units
  constexpr T getCharPhysLowTemperature(  ) const
  {
    return _charPhysLowTemperature;
  };
  /// return characteristic high temperature in physical units
  constexpr T getCharPhysHighTemperature(  ) const
  {
    return _charPhysHighTemperature;
  };
  /// return characteristic temperature difference in physical units
  constexpr T getCharPhysTemperatureDifference(  ) const
  {
    return _charPhysTemperatureDifference;
  };
  /// return thermal expansion coefficient in physical units
  constexpr T getPhysThermalExpansionCoefficient(  ) const
  {
    return _physThermalExpansionCoefficient;
  };
  /// return thermal diffusivity in physical units
  constexpr T getPhysThermalDiffusivity(  ) const
  {
    return _physThermalDiffusivity;
  };
  /// return specific heat capacity in physical units
  constexpr T getPhysSpecificHeatCapacity(  ) const
  {
    return _physSpecificHeatCapacity;
  };
  /// return thermal conductivity in physical units
  constexpr T getThermalConductivity(  ) const
  {
    return _physThermalConductivity;
  };

  /// conversion from lattice to physical temperature
  constexpr T getPhysTemperature( T latticeTemperature ) const
  {
    return _conversionTemperature * (latticeTemperature - 0.5) + _charPhysLowTemperature;
  };
  /// conversion from physical to lattice temperature
  constexpr T getLatticeTemperature( T physTemperature ) const
  {
    return (physTemperature - _charPhysLowTemperature) / _conversionTemperature + 0.5;
  };
  /// access (read-only) to private member variable
  constexpr T getConversionFactorTemperature() const
  {
    return _conversionTemperature;
  };

  /// conversion from lattice to physical thermal diffusivity
  constexpr T getPhysThermalDiffusivity( T latticeThermalDiffusivity ) const
  {
    return _conversionThermalDiffusivity * latticeThermalDiffusivity;
  };
  /// conversion from physical to lattice thermal diffusivity
  constexpr T getLatticeThermalDiffusivity( T physThermalDiffusivity ) const
  {
    return physThermalDiffusivity / _conversionThermalDiffusivity;
  };
  /// access (read-only) to private member variable
  constexpr T getConversionFactorThermalDiffusivity() const
  {
    return _conversionThermalDiffusivity;
  };


  /// conversion from lattice to physical specific heat capacity
  constexpr T getPhysSpecificHeatCapacity( T latticeSpecificHeatCapacity ) const
  {
    return _conversionSpecificHeatCapacity * latticeSpecificHeatCapacity;
  };
  /// conversion from physical to lattice specific heat capacity
  constexpr T getLatticeSpecificHeatCapacity( T physSpecificHeatCapacity ) const
  {
    return physSpecificHeatCapacity / _conversionSpecificHeatCapacity;
  };
  /// access (read-only) to private member variable
  constexpr T getConversionFactorSpecificHeatCapacity() const
  {
    return _conversionSpecificHeatCapacity;
  };

  /// conversion from lattice to physical thermal  conductivity
  constexpr T getPhysThermalConductivity( T latticeThermalConductivity ) const
  {
    return _conversionThermalConductivity * latticeThermalConductivity;
  };
  /// conversion from physical to lattice thermal  conductivity
  constexpr T getLatticeThermalConductivity( T physThermalConductivity ) const
  {
    return physThermalConductivity / _conversionThermalConductivity;
  };
  /// access (read-only) to private member variable
  constexpr T getConversionFactorThermalConductivity() const
  {
    return _conversionThermalConductivity;
  };

  /// conversion from lattice to physical heat flux
  constexpr T getPhysHeatFlux( T latticeHeatFlux ) const
  {
    return _conversionHeatFlux * latticeHeatFlux;
  };
  /// conversion from physical to lattice heat flux
  constexpr T getLatticeHeatFlux( T physHeatFlux ) const
  {
    return physHeatFlux / _conversionHeatFlux;
  };
  /// access (read-only) to private member variable
  constexpr T getConversionFactorHeatFlux() const
  {
    return _conversionHeatFlux;
  };
  constexpr T getPrandtlNumber() const{
    return this->_physViscosity/_physThermalDiffusivity;
  };
  constexpr T getRayleighNumber() const{
    return 9.81 * _physThermalExpansionCoefficient/this->_physViscosity/_physThermalDiffusivity * (_charPhysHighTemperature - _charPhysLowTemperature) * pow(this->_charPhysLength,3);
  };
/// nice terminal output for conversion factors, characteristical and physical data
  virtual void print() const;

  void write(std::string const& fileName = "ThermalUnitConverter") const;



protected:
  // conversion factors
  const T _conversionTemperature; // K
  const T _conversionThermalDiffusivity; // m^2 / s
  const T _conversionSpecificHeatCapacity; // J / kg K = m^2 / s^2 K
  const T _conversionThermalConductivity; // W / m K = kg m / s^3 K
  const T _conversionHeatFlux; // W / m^2 = kg / s^3

  // physical units, e.g characteristic or reference values
  const T _charPhysLowTemperature; // K
  const T _charPhysHighTemperature; // K
  const T _charPhysTemperatureDifference; // K
  const T _physThermalExpansionCoefficient; // 1 / K
  const T _physThermalDiffusivity; // m^2 / s
  const T _physSpecificHeatCapacity; // J / kg K = m^2 / s^2 K
  const T _physThermalConductivity; // W / m K = kg m / s^3 K

  // lattice units, discretization parameters
  const T _latticeThermalRelaxationTime; // -

private:
  mutable OstreamManager clout;
};

}  // namespace olb

#endif
