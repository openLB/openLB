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


#ifndef UNITCONVERTER_HH
#define UNITCONVERTER_HH

#include <fstream>
#include <iostream>
#include <unistd.h>
#include "core/singleton.h"
#include "io/fileName.h"
#include "unitConverter.h"

/// All OpenLB code is contained in this namespace.
namespace olb {

template <typename T, typename DESCRIPTOR>
void UnitConverter<T, DESCRIPTOR>::print() const
{
  clout << "----------------- UnitConverter information -----------------" << std::endl;
  clout << "-- Parameters:" << std::endl;
  clout << "Resolution:                       N=              " << getResolution() << std::endl;
  clout << "Lattice velocity:                 latticeU=       " << getCharLatticeVelocity() << std::endl;
  clout << "Lattice relaxation frequency:     omega=          " << getLatticeRelaxationFrequency(  ) << std::endl;
  clout << "Lattice relaxation time:          tau=            " << getLatticeRelaxationTime() << std::endl;
  clout << "Characteristical length(m):       charL=          " << getCharPhysLength() << std::endl;
  clout << "Characteristical speed(m/s):      charU=          " << getCharPhysVelocity() << std::endl;
  clout << "Phys. kinematic viscosity(m^2/s): charNu=         " << getPhysViscosity() << std::endl;
  clout << "Phys. density(kg/m^d):            charRho=        " << getPhysDensity() << std::endl;
  clout << "Characteristical pressure(N/m^2): charPressure=   " << getCharPhysPressure() << std::endl;
  clout << "Mach number:                      machNumber=     " << getMachNumber() << std::endl;
  clout << "Reynolds number:                  reynoldsNumber= " << getReynoldsNumber() << std::endl;
  clout << "Knudsen number:                   knudsenNumber=  " << getKnudsenNumber() << std::endl;

  clout << std::endl;
  clout << "-- Conversion factors:" << std::endl;
  clout << "Voxel length(m):                  physDeltaX=     " << getConversionFactorLength() << std::endl;
  clout << "Time step(s):                     physDeltaT=     " << getConversionFactorTime() << std::endl;
  clout << "Velocity factor(m/s):             physVelocity=   " << getConversionFactorVelocity() << std::endl;
  clout << "Density factor(kg/m^3):           physDensity=    " << getConversionFactorDensity() <<  std::endl;
  clout << "Mass factor(kg):                  physMass=       " << getConversionFactorMass() << std::endl;
  clout << "Viscosity factor(m^2/s):          physViscosity=  " << getConversionFactorViscosity() << std::endl;
  clout << "Force factor(N):                  physForce=      " << getConversionFactorForce() << std::endl;
  clout << "Pressure factor(N/m^2):           physPressure=   " << getConversionFactorPressure() << std::endl;

  clout << "-------------------------------------------------------------" << std::endl;

}

template <typename T, typename DESCRIPTOR>
void UnitConverter<T, DESCRIPTOR>::write(std::string const& fileName) const
{
  std::string dataFile = singleton::directories().getLogOutDir() + fileName + ".dat";

  if (singleton::mpi().isMainProcessor())
  {
    std::ofstream fout;
    fout.open(dataFile.c_str(), std::ios::trunc);

    fout << "UnitConverter information\n\n";
    fout << "----------------- UnitConverter information -----------------\n";
    fout << "-- Parameters:" << std::endl;
    fout << "Resolution:                       N=              " << getResolution()                << "\n";
    fout << "Lattice velocity:                 latticeU=       " << getCharLatticeVelocity()       << "\n";
    fout << "Lattice relaxation frequency:     omega=          " << getLatticeRelaxationFrequency(  ) << std::endl;
    fout << "Lattice relaxation time:          tau=            " << getLatticeRelaxationTime()     << "\n";
    fout << "Characteristical length(m):       charL=          " << getCharPhysLength()            << "\n";
    fout << "Characteristical speed(m/s):      charU=          " << getCharPhysVelocity()          << "\n";
    fout << "Phys. kinematic viscosity(m^2/s): charNu=         " << getPhysViscosity()             << "\n";
    fout << "Phys. density(kg/m^d):            charRho=        " << getPhysDensity()               << "\n";
    fout << "Characteristical pressure(N/m^2): charPressure=   " << getCharPhysPressure()          << "\n";
    fout << "Mach number:                      machNumber=     " << getMachNumber()                << "\n";
    fout << "Reynolds number:                  reynoldsNumber= " << getReynoldsNumber()            << "\n";
    fout << "Knudsen number:                   knudsenNumber=  " << getKnudsenNumber()             << std::endl;
    fout << "\n";
    fout << "-- Conversion factors:"                                                               << "\n";
    fout << "Voxel length(m):                  physDeltaX=     " << getConversionFactorLength() << std::endl;
    fout << "Time step(s):                     physDeltaT=     " << getConversionFactorTime()      << "\n";
    fout << "Velocity factor(m/s):             physVelocity=   " << getConversionFactorVelocity()  << "\n";
    fout << "Density factor(kg/m^3):           physDensity=    " << getConversionFactorDensity()   << "\n";
    fout << "Mass factor(kg):                  physMass=       " << getConversionFactorMass()      << "\n";
    fout << "Viscosity factor(m^2/s):          physViscosity=  " << getConversionFactorViscosity() << "\n";
    fout << "Force factor(N):                  physForce=      " << getConversionFactorForce()     << "\n";
    fout << "Pressure factor(N/m^2):           physPressure=   " << getConversionFactorPressure()  << "\n";

    fout << "-------------------------------------------------------------" << "\n";

    fout.close();
  }
}

template<typename T, typename DESCRIPTOR>
UnitConverter<T, DESCRIPTOR>* createUnitConverter(XMLreader const& params)
{
  OstreamManager clout(std::cout,"createUnitConverter");
  params.setWarningsOn(false);

  T physDeltaX;
  T physDeltaT;

  T charPhysLength;
  T charPhysVelocity;
  T physViscosity;
  T physDensity;
  T charPhysPressure = 0;

  int resolution;
  T latticeRelaxationTime;
  T charLatticeVelocity;

  // params[parameter].read(value) sets the value or returns false if the parameter can not be found
  params["Application"]["PhysParameters"]["CharPhysLength"].read(charPhysLength);
  params["Application"]["PhysParameters"]["CharPhysVelocity"].read(charPhysVelocity);
  params["Application"]["PhysParameters"]["PhysViscosity"].read(physViscosity);
  params["Application"]["PhysParameters"]["PhysDensity"].read(physDensity);
  params["Application"]["PhysParameters"]["CharPhysPressure"].read(charPhysPressure);

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

      if (params["Application"]["Discretization"]["CharLatticeVelocity"].read(charLatticeVelocity,false)) {
        latticeRelaxationTime = physViscosity * charLatticeVelocity * descriptors::invCs2<T,DESCRIPTOR>() * resolution + 0.5;
      }
      else {
        if (!params["Application"]["Discretization"]["LatticeRelaxationTime"].read(latticeRelaxationTime,false)) {
          clout << "Error: Have not found LatticeRelaxationTime and was not able to derive it using CharLatticeVelocity"
                << std::endl;
          exit (1);
        }
      }
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

  return new UnitConverter<T, DESCRIPTOR>(physDeltaX, physDeltaT, charPhysLength, charPhysVelocity, physViscosity, physDensity, charPhysPressure);
}

}  // namespace olb

#endif
