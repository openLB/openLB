/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Robin Trunk
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

#ifndef ADVECTION_DIFFUSION_FORCES_H
#define ADVECTION_DIFFUSION_FORCES_H

#include "core/unitConverter.h"

namespace olb {

template<typename T, typename DESCRIPTOR,
typename ADLattice=descriptors::D3Q7<descriptors::VELOCITY,descriptors::VELOCITY2>>
class AdvectionDiffusionForce3D {
public:
  AdvectionDiffusionForce3D()
  {
    initArg = 0;
  };
  virtual ~AdvectionDiffusionForce3D() {};
  virtual void applyForce(T force[], Cell<T,DESCRIPTOR> *nsCell, Cell<T,ADLattice> *adCell, T vel[], int latticeR[])=0;
  int getInitArg()
  {
    return initArg;
  }
private:
  int initArg;
};

template<typename T, typename DESCRIPTOR,
typename ADLattice=descriptors::D3Q7<descriptors::VELOCITY,descriptors::VELOCITY2>>
class AdvDiffDragForce3D : public AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice> {
public:
  AdvDiffDragForce3D(UnitConverter<T,DESCRIPTOR> const& converter_, T St_);
  AdvDiffDragForce3D(UnitConverter<T,DESCRIPTOR> const& converter_, T pRadius_, T pRho_);
  ~AdvDiffDragForce3D() override {};
  void applyForce(T force[], Cell<T,DESCRIPTOR> *nsCell, Cell<T,ADLattice> *adCell, T vel[], int latticeR[]) override;

private:
  int initArg;
  T dragCoeff;
};

template<typename T, typename DESCRIPTOR,
typename ADLattice=descriptors::D3Q7<descriptors::VELOCITY,descriptors::VELOCITY2>>
class AdvDiffRotatingForce3D : public AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice> {
public:
  AdvDiffRotatingForce3D(SuperGeometry3D<T>& superGeometry_,
                         const UnitConverter<T,DESCRIPTOR>& converter_,
                         std::vector<T> axisPoint_,
                         std::vector<T> axisDirection_,
                         T w_, T* frac_,
                         bool centrifugeForceOn_ = true,
                         bool coriolisForceOn_ = true);
  AdvDiffRotatingForce3D(UnitConverter<T,DESCRIPTOR> const& converter_, T pRadius_, T pRho_);
  virtual ~AdvDiffRotatingForce3D() {};
  void applyForce(T force[], Cell<T,DESCRIPTOR> *nsCell, Cell<T,ADLattice> *adCell, T vel[], int latticeR[]);

protected:
  SuperGeometry3D<T>& sg;
  std::vector<T> axisPoint;
  std::vector<T> axisDirection;
  T invMassLessForce;
  T w;
  T* frac;
  bool centrifugeForceOn;
  bool coriolisForceOn;

};

template<typename T, typename DESCRIPTOR,
typename ADLattice=descriptors::D3Q7<descriptors::VELOCITY,descriptors::VELOCITY2>>
class AdvDiffMagneticWireForce3D : public AdvectionDiffusionForce3D<T,DESCRIPTOR,ADLattice> {
public:
  AdvDiffMagneticWireForce3D(SuperGeometry3D<T>& superGeometry_, UnitConverter<T,DESCRIPTOR> const& converter_, T pMass, AnalyticalF3D<T, T>& getMagForce);
  ~AdvDiffMagneticWireForce3D() override {};
  void applyForce(T force[], Cell<T,DESCRIPTOR> *nsCell, Cell<T,ADLattice> *adCell, T vel[], int latticeR[]) override;

private:
  SuperGeometry3D<T>& sg;
  int initArg;
  T _pMass;
  T _conversionVelocity;
  AnalyticalF3D<T, T>& _getMagForce;
};

}

#endif
