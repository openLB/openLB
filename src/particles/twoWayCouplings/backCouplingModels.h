/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2019 Davide Dapelo
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

/* Models for Lagrangian back-coupling methods -- header file.
 */

#ifndef LB_BACK_COUPLING_MODELS_H
#define LB_BACK_COUPLING_MODELS_H

#include "twoWayHelperFunctionals.h"

namespace olb {

/** Abstact base class for BaseBackCouplingModel.
  * Its raison d'etre consists of not being templetized in Lattice.
  */
template<typename T, template<typename V> class Particle>
class BackCouplingModel {
public:
  /// Class operator to apply the coupling, for overload.
  virtual bool operator() (Particle<T>* p, int globic, int material, int subCycles=1)=0;
  /// Resets external field
  virtual void resetExternalField(int material)=0;
};

/** Abstact class for all the back-coupling models,
  * viz., momentum coupling from particle to fluid.
  * Input parameters in attice units.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class BaseBackCouplingModel : public BackCouplingModel<T,Particle> {
public:
  /// Resets external field
  virtual void resetExternalField(int material) override;
protected:
  /// Constructor
  BaseBackCouplingModel ( UnitConverter<T, Lattice>& converter,
                         SuperLattice3D<T, Lattice>& sLattice,
                         SuperGeometry3D<T>& sGeometry );
  UnitConverter<T, Lattice>& _converter; // reference to a UnitConverter
  SuperGeometry3D<T>& _sGeometry;
  SuperLattice3D<T, Lattice>& _sLattice; // reference to a lattice
private:
  // Pointers to functions to reset fluid force
  std::shared_ptr<AnalyticalConst3D<T, T> > _zeroAnalytical;
  std::shared_ptr<AnalyticalComposed3D<T, T> > _zeroField;
};

/** Back-coupling is performed on the cell containing the particle
  * and its neighbours within a cube of one lattice spacing, as per in Maier et al. (2017).
  * Input parameters in attice units.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class CubicDeltaBackCouplingModel : public BaseBackCouplingModel<T,Lattice,Particle> {
public:
  /// Constructor
  CubicDeltaBackCouplingModel ( UnitConverter<T, Lattice>& converter,
                         SuperLattice3D<T, Lattice>& sLattice,
                         SuperGeometry3D<T>& sGeometry );
  /// Class operator to apply the coupling.
  virtual bool operator() (Particle<T>* p, int globic, int material, int subCycles=1) override;
protected:
  int _range = 1;
  T _delta[4][4][4] = { T() };
  std::shared_ptr<SuperLatticeSmoothDiracDelta3D<T, Lattice> > _cubicDeltaFunctional;
};

/** Back-coupling is performed only on the cell containing the particle.
  * Input parameters in attice units.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class LocalBackCouplingModel : public BaseBackCouplingModel<T,Lattice,Particle> {
public:
  /// Constructor
  LocalBackCouplingModel ( UnitConverter<T, Lattice>& converter,
                         SuperLattice3D<T, Lattice>& sLattice,
                         SuperGeometry3D<T>& sGeometry );
  /// Class operator to apply the coupling.
  virtual bool operator() (Particle<T>* p, int globic, int material, int subCycles=1) override;
};

/** Class for a generic non-local back-coupling model,
  * viz., momentum coupling from particle to fluid.
  * It reproduces the characteristics (viz., smoothing) of an input forward coupling model.
  * Input parameters in attice units.
  */
template<typename T, typename Lattice, template<typename V> class Particle>
class NonLocalBaseBackCouplingModel : public BaseBackCouplingModel<T,Lattice,Particle> {
public:
  /// Constructor
  NonLocalBaseBackCouplingModel ( UnitConverter<T, Lattice>& converter,
                         SuperLattice3D<T, Lattice>& sLattice,
                         SuperGeometry3D<T>& sGeometry,
                         SmoothingFunctional<T, Lattice>& smoothingFunctional );
  /// Class operator to apply the coupling.
  virtual bool operator() (Particle<T>* p, int globic, int material, int subCycles=1) override;
protected:
  SmoothingFunctional<T, Lattice>& _smoothingFunctional; // Functional to treat non-local smoothing
};

}

#endif
