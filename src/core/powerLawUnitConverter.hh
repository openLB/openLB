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


#ifndef PL_UNITCONVERTER_HH
#define PL_UNITCONVERTER_HH

#include <fstream>
#include <iostream>
#include <unistd.h>
#include "core/singleton.h"
#include "io/fileName.h"

/// All OpenLB code is contained in this namespace.
namespace olb {

template <typename T, typename DESCRIPTOR>
void PowerLawUnitConverter<T, DESCRIPTOR>::print() const
{
  clout << "----------------- UnitConverter information ------------------" << std::endl;
  clout << "-- Parameters:" << std::endl;
  clout << "Resolution:                           N=              " << this->getResolution() << std::endl;
  clout << "DESCRIPTOR velocity:                     latticeU=       " << this->getCharLatticeVelocity() << std::endl;
  clout << "DESCRIPTOR relaxation frequency:         omega=          " << this->getLatticeRelaxationFrequency(  ) << std::endl;
  clout << "DESCRIPTOR relaxation time:              tau=            " << this->getLatticeRelaxationTime() << std::endl;
  clout << "Characteristical length(m):           charL=          " << this->getCharPhysLength() << std::endl;
  clout << "Characteristical speed(m/s):          charU=          " << this->getCharPhysVelocity() << std::endl;
  clout << "Phys. char kinematic visco(m^2/s):    charNu=         " << this->getPhysViscosity() << std::endl;
  clout << "Phys. consistency coeff(m^2 s^(n-2)): charM=          " << this->getPhysConsistencyCoeff() << std::endl;
  clout << "Power-law index:                      n=              " << this->getPowerLawIndex() << std::endl;
  clout << "Phys. density(kg/m^d):                charRho=        " << this->getPhysDensity() << std::endl;
  clout << "Characteristical pressure(N/m^2):     charPressure=   " << this->getCharPhysPressure() << std::endl;
  clout << "Reynolds number:                      reynoldsNumber= " << this->getReynoldsNumber() << std::endl;

  clout << std::endl;
  clout << "-- Conversion factors:" << std::endl;
  clout << "Voxel length(m):                      physDeltaX=     " << this->getConversionFactorLength() << std::endl;
  clout << "Time step(s):                         physDeltaT=     " << this->getConversionFactorTime() << std::endl;
  clout << "Velocity factor(m/s):                 physVelocity=   " << this->getConversionFactorVelocity() << std::endl;
  clout << "Density factor(kg/m^3):               physDensity=    " << this->getConversionFactorDensity() <<  std::endl;
  clout << "Mass factor(kg):                      physMass=       " << this->getConversionFactorMass() << std::endl;
  clout << "Viscosity factor(m^2/s):              physViscosity=  " << this->getConversionFactorViscosity() << std::endl;
  clout << "Force factor(N):                      physForce=      " << this->getConversionFactorForce() << std::endl;
  clout << "Pressure factor(N/m^2):               physPressure=   " << this->getConversionFactorPressure() << std::endl;

  clout << "--------------------------------------------------------------" << std::endl;

}

template <typename T, typename DESCRIPTOR>
void PowerLawUnitConverter<T, DESCRIPTOR>::write(std::string const& fileName) const
{
  std::string dataFile = singleton::directories().getLogOutDir() + fileName + ".dat";

  if (singleton::mpi().isMainProcessor())
  {
    std::ofstream fout;
    fout.open(dataFile.c_str(), std::ios::trunc);

    fout << "UnitConverter information\n\n";
    fout << "----------------- UnitConverter information ------------------\n";
    fout << "-- Parameters:" << std::endl;
    fout << "Resolution:                           N=              " << this->getResolution()                << "\n";
    fout << "DESCRIPTOR velocity:                     latticeU=       " << this->getCharLatticeVelocity()       << "\n";
    fout << "DESCRIPTOR relaxation frequency:         omega=          " << this->getLatticeRelaxationFrequency(  ) << std::endl;
    fout << "DESCRIPTOR relaxation time:              tau=            " << this->getLatticeRelaxationTime()     << "\n";
    fout << "Characteristical length(m):           charL=          " << this->getCharPhysLength()            << "\n";
    fout << "Characteristical speed(m/s):          charU=          " << this->getCharPhysVelocity()          << "\n";
    fout << "Phys. char kinematic visco(m^2/s):    charNu=         " << this->getPhysViscosity() << std::endl;
    fout << "Phys. consistency coeff(m^2 s^(n-2)): charM=          " << this->getPhysConsistencyCoeff() << std::endl;
    fout << "Power-law index:                      n=              " << this->getPowerLawIndex() << std::endl;
    fout << "Phys. density(kg/m^d):                charRho=        " << this->getPhysDensity()               << "\n";
    fout << "Characteristical pressure(N/m^2):     charPressure=   " << this->getCharPhysPressure()          << "\n";
    fout << "Reynolds number:                      reynoldsNumber= " << this->getReynoldsNumber() << std::endl;
    fout << "\n";
    fout << "-- Conversion factors:"                                                               << "\n";
    fout << "Voxel length(m):                      physDeltaX=     " << this->getConversionFactorLength() << std::endl;
    fout << "Time step(s):                         physDeltaT=     " << this->getConversionFactorTime()      << "\n";
    fout << "Velocity factor(m/s):                 physVelocity=   " << this->getConversionFactorVelocity()  << "\n";
    fout << "Density factor(kg/m^3):               physDensity=    " << this->getConversionFactorDensity()   << "\n";
    fout << "Mass factor(kg):                      physMass=       " << this->getConversionFactorMass()      << "\n";
    fout << "Viscosity factor(m^2/s):              physViscosity=  " << this->getConversionFactorViscosity() << "\n";
    fout << "Force factor(N):                      physForce=      " << this->getConversionFactorForce()     << "\n";
    fout << "Pressure factor(N/m^2):               physPressure=   " << this->getConversionFactorPressure()  << "\n";

    fout << "--------------------------------------------------------------" << "\n";

    fout.close();
  }
}

template<typename T, typename DESCRIPTOR>
PowerLawUnitConverter<T, DESCRIPTOR>* createPowerLawUnitConverter(XMLreader const& params)
{
  OstreamManager clout(std::cout,"createUnitConverter");
  params.setWarningsOn(false);

  T physDeltaX;
  T physDeltaT;

  T charPhysLength;
  T charPhysVelocity;
  T physViscosity;
  T physConsistencyCoeff;
  T powerLawIndex;
  T physDensity;
  T charPhysPressure = 0;

  int resolution;
  T latticeRelaxationTime;
  T charLatticeVelocity;

  // params[parameter].read(value) sets the value or returns false if the parameter can not be found
  params["Application"]["PhysParameters"]["CharPhysLength"].read(charPhysLength);
  params["Application"]["PhysParameters"]["CharPhysVelocity"].read(charPhysVelocity);
  params["Application"]["PhysParameters"]["PhysConsistencyCoeff"].read(physConsistencyCoeff);
  params["Application"]["PhysParameters"]["powerLawIndex"].read(powerLawIndex);
  params["Application"]["PhysParameters"]["PhysDensity"].read(physDensity);
  params["Application"]["PhysParameters"]["CharPhysPressure"].read(charPhysPressure);

  physViscosity = physConsistencyCoeff * pow(charPhysVelocity / (2*charPhysLength), powerLawIndex-1);

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

  return new PowerLawUnitConverter<T, DESCRIPTOR>(physDeltaX, physDeltaT, charPhysLength, charPhysVelocity, physConsistencyCoeff, powerLawIndex, physDensity, charPhysPressure);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/*
template <typename T, typename DESCRIPTOR>
constexpr PowerLawUnitConverterFrom_Resolution_RelaxationTime_Reynolds_PLindex<T, DESCRIPTOR>::
PowerLawUnitConverterFrom_Resolution_RelaxationTime_Reynolds_PLindex(
    int resolution,
    T latticeRelaxationTime,
    T charPhysLength,
    T charPhysVelocity,
    T Re,
   T powerLawIndex,
    T physDensity,
    T charPhysPressure)
{
  T physDeltaX = (charPhysLength/resolution);
  T physConsistencyCoeff = charPhysLength * charPhysVelocity * pow( charPhysVelocity / ( 2 * charPhysLength ), 1 - powerLawIndex ) / Re;
  T physViscosity = physConsistencyCoeff * pow( charPhysVelocity / (2 * charPhysLength ), powerLawIndex - 1 );
  T physDeltaT = (latticeRelaxationTime - 0.5) / descriptors::invCs2<T,DESCRIPTOR>() * pow((charPhysLength/resolution),2) / physViscosity;

PowerLawUnitConverter<T, DESCRIPTOR>( physDeltaX, physDeltaT, charPhysLength, charPhysVelocity,
              physConsistencyCoeff, powerLawIndex, physDensity, charPhysPressure );
}
*/




}  // namespace olb

#endif
