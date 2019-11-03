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


#ifndef THERMALUNITCONVERTER_HH
#define THERMALUNITCONVERTER_HH

#include <fstream>
#include <iostream>
#include <unistd.h>
#include "core/singleton.h"
#include "io/fileName.h"

/// All OpenLB code is contained in this namespace.
namespace olb {

template <typename T, typename DESCRIPTOR, typename ThermalLattice>
void ThermalUnitConverter<T, DESCRIPTOR, ThermalLattice>::print() const
{
  clout << "----------------- UnitConverter information -----------------" << std::endl;
  clout << "-- Parameters:" << std::endl;
  clout << "Resolution:                                 N=                              " << this->getResolution() << std::endl;
  clout << "Lattice velocity:                           latticeU=                       " << this->getCharLatticeVelocity() << std::endl;
  clout << "Lattice relaxation frequency:               omega=                          " << this->getLatticeRelaxationFrequency() << std::endl;
  clout << "Lattice relaxation time:                    tau=                            " << this->getLatticeRelaxationTime() << std::endl;
  clout << "Thermal Lattice relaxation frequency:       omega_AD=                       " << this->getLatticeThermalRelaxationFrequency() << std::endl;
  clout << "Thermal Lattice relaxation time:            tau_AD=                         " << this->getLatticeThermalRelaxationTime() << std::endl;
  clout << "Characteristical length(m):                 charL=                          " << this->getCharPhysLength() << std::endl;
  clout << "Characteristical speed(m/s):                charU=                          " << this->getCharPhysVelocity() << std::endl;
  clout << "Phys. kinematic viscosity(m^2/s):           charNu=                         " << this->getPhysViscosity() << std::endl;
  clout << "Phys. density(kg/m^d):                      charRho=                        " << this->getPhysDensity() << std::endl;
  clout << "Characteristical pressure(N/m^2):           charPressure=                   " << this->getCharPhysPressure() << std::endl;
  clout << "Reynolds number:                            reynoldsNumber=                 " << this->getReynoldsNumber() << std::endl;

  clout << "-------------------------------------------------------------" << std::endl;

  clout << "----------------- ThermalUnitConverter information -----------------" << std::endl;
  clout << "-- Parameters:" << std::endl;
  clout << "Phys. Delta X(m):                           physDeltaX=                     " << this->getPhysDeltaX() << std::endl;
  clout << "Phys. Delta T(s):                           physDeltaT=                     " << this->getPhysDeltaT() << std::endl;
  clout << "Characteristical pressure(N/m^2):           charPressure=                   " << this->getCharPhysPressure() << std::endl;
  clout << "Phys. Thermal Conductivity(W/m/K):          physThermalCondcticity=         " << getThermalConductivity() << std::endl;
  clout << "Phys. specific Heat Capacity(J/kg/K):       physSpecificHeatCapacity=       " << getPhysSpecificHeatCapacity() << std::endl;
  clout << "Phys. Thermal Expasion Coefficent(K^-1):    physThermalExpansionCoefficent= " << getPhysThermalExpansionCoefficient() << std::endl;
  clout << "Characteristical Phys. low Temperature(K):  charPhysLowTemperature=         " << getCharPhysLowTemperature() << std::endl;
  clout << "Characteristical Phys. high Temperature(K): charPhysHighTemperature=        " << getCharPhysHighTemperature() << std::endl;
  clout << "Prandtl number:                             prandtlNumber=                  " << getPrandtlNumber() << std::endl;
  clout << "Rayleigh number:                            rayleighNumber=                 " << getRayleighNumber() << std::endl;


  clout << "-------------------------------------------------------------" << std::endl;

  clout << "----------------- Conversion factors:-----------------" << std::endl;
  clout << "Voxel length(m):                            physDeltaX=                     " << this->getConversionFactorLength() << std::endl;
  clout << "Time step(s):                               physDeltaT=                     " << this->getConversionFactorTime() << std::endl;
  clout << "Velocity factor(m/s):                       physVelocity=                   " << this->getConversionFactorVelocity() << std::endl;
  clout << "Density factor(kg/m^3):                     physDensity=                    " << this->getConversionFactorDensity() <<  std::endl;
  clout << "Mass factor(kg):                            physMass=                       " << this->getConversionFactorMass() << std::endl;
  clout << "Viscosity factor(m^2/s):                    physViscosity=                  " << this->getConversionFactorViscosity() << std::endl;
  clout << "Force factor(N):                            physForce=                      " << this->getConversionFactorForce() << std::endl;
  clout << "Pressure factor(N/m^2):                     physPressure=                   " << this->getConversionFactorPressure() << std::endl;

  clout << "-------------------------------------------------------------" << std::endl;

  clout << "----------------- ThermalConversion factors:-----------------" << std::endl;
  clout << "Temperature(K):                             temperature=                    " << getConversionFactorTemperature() << std::endl;
  clout << "Thermal Duffusity(m^2/s):                   physThermalDiffusity=           " << getConversionFactorThermalDiffusivity() << std::endl;
  clout << "specific Heat Capacity(J/kg):               physSpecificHeatCapacity=       " << getConversionFactorSpecificHeatCapacity() << std::endl;
  clout << "Thermal Coductivity(W/m/K):                 physThermalCondcticity=         " << getConversionFactorThermalConductivity() <<  std::endl;
  clout << "HeatFlux(W):                                physHeatFlux=                   " << getConversionFactorHeatFlux() << std::endl;

  clout << "-------------------------------------------------------------" << std::endl;

}
/*
template <typename T, typename DESCRIPTOR, typename ThermalLattice>
void ThermalUnitConverter<T, DESCRIPTOR, ThermalLattice>::writeDatFile(std::string const& title) const
{
  std::string dataFile = singleton::directories().getLogOutDir() + title + ".dat";

  if (singleton::mpi().isMainProcessor())
  {
    std::ofstream fout;
    fout.open(dataFile.c_str(), std::ios::trunc);

  fout << "----------------- ThermalUnitConverter information -----------------" << std::endl;
  fout << "-- Parameters:" << std::endl;
  fout << "Phys. Delta X:                    physDeltaX=              " << getPhysDeltaX() << std::endl;
  fout << "Phys. Delta T:                    physDeltaT=            " << getPhysDeltaT() << std::endl;
  fout << "Characteristical pressure(N/m^2): charPressure=   " << getCharPhysPressure() << std::endl;
  fout << "Phys. Thermal Conductivity:       physThermalCondcticity=   " << getPhysThermalCondcticity() << std::endl;
  fout << "Phys. specific Heat Capacity:     physSpecificHeatCapacity= " << getPhysSpecificHeatCapacity() << std::endl;
  fout << "Phys. Thermal Expasion Coefficent physThermalExpansionCoefficent= " << getPhysThermalExpasionCoefficent << std::endl;
  fout << "Characteristical Phys. low Temperature charPhysLowTemperature= " << getCharPhysLowTemperature << std::endl;
  fout << "Characteristical Phys. high Temperature charPhysHighTemperature= " << getCharPhysHighTemperature << std::endl;
  fout << std::endl;

  fout << "-- Conversion factors:" << std::endl;
  fout << "Temperature:                      temperature=       " << getConversionFactorTemperature() << std::endl;
  fout << "Thermal Duffusity:                physThermalDiffusity=       " << getConversionFactorThermalDiffusivity() << std::endl;
  fout << "specific Heat Capacity:           physSpecificHeatCapacity=   " << getConversionFactorSpecificHeatCapacity() << std::endl;
  fout << "Thermal Coductivity:              physThermalCondcticity=    " << getConversionFactorThermalConductivity() <<  std::endl;
  fout << "HeatFlux:                         physHeatFlux=       " << getConversionFactorHeatFlux() << std::endl;

  fout << "-------------------------------------------------------------" << std::endl;

    fout.close();
  }
}
*/
/*
template<typename T, typename DESCRIPTOR>
ThermalUnitConverter<T, DESCRIPTOR>* createThermalUnitConverter(XMLreader const& params)
{
  OstreamManager clout(std::cout,"createThermalUnitConverter");
  params.setWarningsOn(false);

  T physDeltaX;
  T physDeltaT;

  T charPhysHighTemperature;
  T charPhysLowTemperature;
  T physThermalCondcticity;
  T physDensity;
  T charPhysPressure = 0;

  int resolution;
  T latticeRelaxationTime;
  T charLatticeVelocity;

  // params[parameter].read(value) sets the value or returns false if the parameter can not be found
  params["Application"]["ThermalPhysParameters"]["CharPhysLowTemperature"].read(charPhysLowTemperature);
  params["Application"]["ThermalPhysParameters"]["charPhysHighTemperature"].read(charPhysHighTemperature);
  params["Application"]["ThermalPhysParameters"]["PhysThermalCondcticity"].read(physThermalCondcticity);
  params["Application"]["ThermalPhysParameters"]["PhysDensity"].read(physDensity);
  params["Application"]["ThermalPhysParameters"]["CharPhysPressure"].read(charPhysPressure);

  if (!params["Application"]["Discretization"]["PhysDeltaX"].read(physDeltaX,false)) {
    if (!params["Application"]["Discretization"]["Resolution"].read<int>(resolution,false)) {
      if (!params["Application"]["Discretization"]["CharLatticeVelocity"].read(charLatticeVelocity,false)) {
        // NOT found physDeltaX, resolution or charLatticeVelocity
        clout << "Error: Have not found PhysDeltaX, Resolution or CharLatticeVelocity in XML file."
              << std::endl;
        exit (1);
      }
      else {
        // found charLatticeVelocity
        if (params["Application"]["Discretization"]["PhysDeltaT"].read(physDeltaT,false)) {
          physDeltaX = charPhysVelocity / charLatticeVelocity * physDeltaT;
        }
        else if (params["Application"]["Discretization"]["LatticeRelaxationTime"].read(latticeRelaxationTime,false)) {
          physDeltaX = physViscosity * charLatticeVelocity / charPhysVelocity * descriptors::invCs2<T,DESCRIPTOR>() / (latticeRelaxationTime - 0.5);
        }
      }
    }
    else {
      // found resolution
      physDeltaX = charPhysLength / resolution;
    }
  }
  // found physDeltaX
  if (!params["Application"]["Discretization"]["PhysDeltaT"].read(physDeltaT,false)) {
    if (!params["Application"]["Discretization"]["LatticeRelaxationTime"].read(latticeRelaxationTime,false)) {
      if (!params["Application"]["Discretization"]["CharLatticeVelocity"].read(charLatticeVelocity,false)) {
        // NOT found physDeltaT, latticeRelaxationTime and charLatticeVelocity
        clout << "Error: Have not found PhysDeltaT, LatticeRelaxationTime or CharLatticeVelocity in XML file."
              << std::endl;
        exit (1);
      }
      else {
        // found charLatticeVelocity
        physDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX;
      }
    }
    else {
      // found latticeRelaxationTime
      physDeltaT = (latticeRelaxationTime - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * physDeltaX * physDeltaX / physViscosity;
    }
  }

  return new ThermalUnitConverter<T, DESCRIPTOR, ThermalLattice>(physDeltaX, physDeltaT, charPhysLength, charPhysVelocity, physViscosity, physDensity,T physThermalConductivity,T physSpecificHeatCapacity,T physThermalExpansionCoefficient,T charPhysLowTemperature,T charPhysHighTemperature, charPhysPressure);
}*/

}  // namespace olb

#endif
