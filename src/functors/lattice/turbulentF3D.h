/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013-2015 Patrick Nathen, Mathias J. Krause
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

#ifndef TURBULENT_F_3D_H
#define TURBULENT_F_3D_H

#include <list>

#include "blockBaseF3D.h"
#include "superBaseF3D.h"
#include "blockLatticeLocalF3D.h"
#include "superLatticeLocalF3D.h"
#include "core/unitConverter.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"


/** These are functors used for turbulent flows. Some like AMD have an execute member
 *  function which writes the data into the external field of a lattice descriptor.
 */

namespace olb {

/// functor to get pointwise yPlus from rho, shear stress and local density on local lattices
template <typename T, typename DESCRIPTOR>
class SuperLatticeYplus3D : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  IndicatorF3D<T>&    _indicator;
  const int           _material;
public:
  SuperLatticeYplus3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter,
                      SuperGeometry3D<T>& superGeometry, IndicatorF3D<T>& indicator,
                      const int material );
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise explicit filtering on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
/*template <typename T, typename DESCRIPTOR>
class BlockLatticeADM3D : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  T _sigma;
  int _order;
  bool _adaptive;
  const UnitConverter<T,DESCRIPTOR>& _converter;

public:
  BlockLatticeADM3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, T sigma, int order, bool adaptive, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
  void execute(const int input[]);
  void execute();

};

/// functor to get pointwise ecplicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperLatticeADM3D : public SuperLatticeF3D<T,DESCRIPTOR> {
protected:
  T _sigma;
  int _order;
  bool _adaptive;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticeADM3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, T sigma, int order, bool adaptive, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
  void execute(SuperGeometry3D<T>& superGeometry, const int material);
};
*/

/// functor to get pointwise finite difference Dissipation on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T>
class BlockFiniteDifference3D : public BlockF3D<T> {
private:
  BlockGeometryStructure3D<T>& _blockGeometry;
  BlockF3D<T>& _blockFunctor;
  std::list<int>& _matNumber;
  int _targetDim;
  int _n[3];

public:
  BlockFiniteDifference3D(BlockGeometryStructure3D<T>& blockGeometry, BlockF3D<T>& blockFunctor, std::list<int>& matNumber);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T>
class SuperFiniteDifference3D : public SuperF3D<T> {
private:
  SuperGeometry3D<T>& _sGeometry;
  SuperF3D<T>& _sFunctor;
  std::list<int>& _matNumber;
public:
  SuperFiniteDifference3D(SuperGeometry3D<T>& sGeometry, SuperF3D<T>& sFunctor, std::list<int>& matNumber);
  bool operator() (T output[], const int input[]) override;
};

template <typename T, typename DESCRIPTOR>
class BlockPhysFiniteDifference3D : public BlockF3D<T> {
private:
  BlockF3D<T>& _blockFinDiff;
  int _targetDim;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  BlockPhysFiniteDifference3D(BlockF3D<T>& blockFunctor, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperPhysFiniteDifference3D : public SuperF3D<T> {
private:
  SuperFiniteDifference3D<T> _sFinDiff;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperPhysFiniteDifference3D(SuperGeometry3D<T>& sGeometry, SuperF3D<T>& sFunctor, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

template <typename T, typename DESCRIPTOR>
class BlockLatticeVelocityGradientFD3D : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockF3D<T>& _blockFinDiff;

public:
  BlockLatticeVelocityGradientFD3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFunctor);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class BlockLatticeExternalVelocityGradientFD3D : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockF3D<T>& _blockFinDiff;

public:
  BlockLatticeExternalVelocityGradientFD3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFunctor);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperLatticeVelocityGradientFD3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticeVelocity3D<T,DESCRIPTOR> _sVelocity;
  SuperFiniteDifference3D<T> _sFinDiff;
public:
  SuperLatticeVelocityGradientFD3D(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperLatticeExternalVelocityGradientFD3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticeExternalVelocity3D<T,DESCRIPTOR> _sVelocity;
  SuperFiniteDifference3D<T> _sFinDiff;
public:
  SuperLatticeExternalVelocityGradientFD3D(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber);
  bool operator() (T output[], const int input[]);
};


template <typename T, typename DESCRIPTOR>
class BlockLatticePhysVelocityGradientFD3D : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockF3D<T>& _blockFinDiff;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  BlockLatticePhysVelocityGradientFD3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFunctor, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysVelocityGradientFD3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> _sVelocity;
  SuperPhysFiniteDifference3D<T,DESCRIPTOR> _sFinDiff;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticePhysVelocityGradientFD3D(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

template <typename T, typename DESCRIPTOR>
class BlockLatticeStrainRateFD3D : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockF3D<T>& _blockVeloGrad;
public:
  BlockLatticeStrainRateFD3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFunctor);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperLatticeStrainRateFD3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticeVelocityGradientFD3D<T,DESCRIPTOR> _sVeloGrad;
public:
  SuperLatticeStrainRateFD3D(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class BlockLatticePhysStrainRateFD3D : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockF3D<T>& _blockVeloGrad;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  BlockLatticePhysStrainRateFD3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFunctor, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysStrainRateFD3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticePhysVelocityGradientFD3D<T,DESCRIPTOR> _sVeloGrad;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticePhysStrainRateFD3D(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class BlockLatticeDissipationFD3D : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockF3D<T>& _blockVeloGrad;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  BlockLatticeDissipationFD3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFunctor, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperLatticeDissipationFD3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticeVelocityGradientFD3D<T,DESCRIPTOR> _sVeloGrad;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticeDissipationFD3D(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class BlockLatticePhysDissipationFD3D : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockF3D<T>& _blockVeloGrad;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  BlockLatticePhysDissipationFD3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFunctor, const UnitConverter<T,DESCRIPTOR>& _converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysDissipationFD3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticePhysVelocityGradientFD3D<T,DESCRIPTOR> _sVeloGrad;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticePhysDissipationFD3D(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

template <typename T, typename DESCRIPTOR>
class BlockLatticeEffectiveDissipationFD3D : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockF3D<T>& _blockVeloGrad;
  const UnitConverter<T,DESCRIPTOR>& _converter;
  LESDynamics<T, DESCRIPTOR>& _LESdynamics;
public:
  BlockLatticeEffectiveDissipationFD3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFunctor,
                                       const UnitConverter<T,DESCRIPTOR>& converter, LESDynamics<T, DESCRIPTOR>& LESdynamics);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperLatticeEffectiveDissipationFD3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticeVelocityGradientFD3D<T,DESCRIPTOR> _sVeloGrad;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticeEffectiveDissipationFD3D(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber,
                                       const UnitConverter<T,DESCRIPTOR>& converter, LESDynamics<T, DESCRIPTOR>& LESdynamics);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class BlockLatticePhysEffectiveDissipationFD3D : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockF3D<T>& _blockVeloGrad;
  const UnitConverter<T,DESCRIPTOR>& _converter;
  LESDynamics<T, DESCRIPTOR>& _LESdynamics;
public:
  BlockLatticePhysEffectiveDissipationFD3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFunctor,
      const UnitConverter<T,DESCRIPTOR>& converter, LESDynamics<T, DESCRIPTOR>& LESdynamics);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysEffectiveDissipationFD3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticePhysVelocityGradientFD3D<T,DESCRIPTOR> _sVeloGrad;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticePhysEffectiveDissipationFD3D(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber,
      const UnitConverter<T,DESCRIPTOR>& converter, LESDynamics<T, DESCRIPTOR>& LESdynamics);
  bool operator() (T output[], const int input[]) override;
};

template <typename T, typename DESCRIPTOR>
class BlockLatticeVorticityFD3D : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockF3D<T>& _blockVeloGrad;
public:
  BlockLatticeVorticityFD3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFunctor);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperLatticeVorticityFD3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticeVelocityGradientFD3D<T,DESCRIPTOR> _sVeloGrad;
public:
  SuperLatticeVorticityFD3D(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class BlockLatticePhysVorticityFD3D : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockF3D<T>& _blockVeloGrad;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  BlockLatticePhysVorticityFD3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFunctor, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysVorticityFD3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticePhysVelocityGradientFD3D<T,DESCRIPTOR> _sVeloGrad;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticePhysVorticityFD3D(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor that returns pointwise the enstrophy
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysEnstrophyFD3D : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockF3D<T>& _blockVeloGrad;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  BlockLatticePhysEnstrophyFD3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockVeloGrad, const  UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class SuperLatticePhysEnstrophyFD3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticePhysVelocityGradientFD3D<T,DESCRIPTOR> _sVeloGrad;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticePhysEnstrophyFD3D(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const  UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class BlockLatticePhysStressFD3D : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockF3D<T>& _blockStrainRate;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  BlockLatticePhysStressFD3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& blockFunctor, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysStressFD3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticePhysStrainRateFD3D<T,DESCRIPTOR> _sStrainRate;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticePhysStressFD3D(SuperGeometry3D<T>& sGeometry, SuperLattice3D<T,DESCRIPTOR>& sLattice, std::list<int>& matNumber, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor that returns pointwise the turbulent, kinetic energy
template <typename T, typename DESCRIPTOR>
class BlockIsotropicHomogeneousTKE3D : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockF3D<T>& _blockVelocity;

public:
  BlockIsotropicHomogeneousTKE3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, BlockF3D<T>& f);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class SuperIsotropicHomogeneousTKE3D : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> _sVelocity;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperIsotropicHomogeneousTKE3D( SuperLattice3D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

/*
template <typename T, typename DESCRIPTOR>
class BlockLatticeSigmaADM3D : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
public:
  BlockLatticeSigmaADM3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise explicit filter on local lattice, if globIC is not on
/// the local processor, the returned vector is empty
template <typename T, typename DESCRIPTOR>
class SuperLatticeSigmaADM3D : public SuperLatticeF3D<T,DESCRIPTOR> {
protected:
public:
  SuperLatticeSigmaADM3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]);
};
*/

} // end namespace olb

#endif


