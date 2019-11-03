/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#ifndef BLOCK_LATTICE_LOCAL_F_3D_H
#define BLOCK_LATTICE_LOCAL_F_3D_H

#include "blockBaseF3D.h"
#include "geometry/blockGeometry3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/blockIndicatorBaseF3D.h"
#include "dynamics/smagorinskyBGKdynamics.h"
#include "dynamics/porousBGKdynamics.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

template<typename T, typename DESCRIPTOR> class blockLatticeStructure3D;
template<typename T, typename S> class ParticleIndicatorF3D;

////////////////////////////////////////////////////////////////////////////////
//// if globIC is not on the local processor, the returned vector is empty//////
////////////////////////////////////////////////////////////////////////////////

/// functor returns pointwise f population on local lattices
template <typename T, typename DESCRIPTOR>
class BlockLatticeFpop3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticeFpop3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise dissipation density on local lattices
template <typename T, typename DESCRIPTOR>
class BlockLatticeDissipation3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  BlockLatticeDissipation3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise dissipation density on local lattices
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysDissipation3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  const int             _overlap;
  const UnitConverter<T, DESCRIPTOR>& _converter;
public:
  BlockLatticePhysDissipation3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                int overlap,
                                const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise turbulent dissipation density on local lattices
template <typename T, typename DESCRIPTOR>
class BlockLatticeEffevtiveDissipation3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  const UnitConverter<T,DESCRIPTOR>& _converter;
  T _smagoConst;
  LESDynamics<T, DESCRIPTOR>& _LESdynamics;
public:
  BlockLatticeEffevtiveDissipation3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst,
                                     LESDynamics<T, DESCRIPTOR>& LESdynamics);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise turbulent dissipation density on local lattices
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysEffevtiveDissipation3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  const UnitConverter<T,DESCRIPTOR>& _converter;
  T _smagoConst;
  LESDynamics<T, DESCRIPTOR>& _LESdynamics;
public:
  BlockLatticePhysEffevtiveDissipation3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                         const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst,
                                         LESDynamics<T, DESCRIPTOR>& LESdynamics);
  bool operator() (T output[], const int input[]) override;
};


/// functor returns pointwise density rho on local lattices
template <typename T, typename DESCRIPTOR>
class BlockLatticeDensity3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticeDensity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]) override;
};


/// functor returns pointwise velocity on local lattice
template <typename T, typename DESCRIPTOR>
class BlockLatticeVelocity3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticeVelocity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise external velocity (external field) on local lattice
template <typename T, typename DESCRIPTOR>
class BlockLatticeExternalVelocity3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticeExternalVelocity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]);
};

/// functor returns pointwise the material no. presenting the geometry on local lattice
template <typename T, typename DESCRIPTOR>
class BlockLatticeGeometry3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
  BlockGeometryStructure3D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticeGeometry3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                         BlockGeometryStructure3D<T>& blockGeometry, int material = -1);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise the rank no. + 1 on local lattice
template <typename T, typename DESCRIPTOR>
class BlockLatticeRank3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticeRank3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise the cuboid no. + 1 on local lattice
template <typename T, typename DESCRIPTOR>
class BlockLatticeCuboid3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  // holds cuboid nmb of current block
  int _iC;
public:
  BlockLatticeCuboid3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, int iC);
  bool operator() (T output[], const int input[]) override;
};


/// functor returns pointwise phys pressure from rho on local lattices
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysPressure3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  const int _overlap;
public:
  BlockLatticePhysPressure3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                             int overlap,
                             const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// functor returns pointwise phys velocity on local lattice
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysVelocity3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  const int  _overlap;
  const bool _print;
public:
  BlockLatticePhysVelocity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                             int overlap,
                             const UnitConverter<T,DESCRIPTOR>& converter,
                             bool print=false);
  bool operator() (T output[], const int input[]) override;
};

template <typename T, typename DESCRIPTOR>
class BlockLatticePhysExternalVelocity3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
public:
  BlockLatticePhysExternalVelocity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class BlockLatticePhysExternalPorosity3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  const int _overlap;
public:
  BlockLatticePhysExternalPorosity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                     int overlap,
                                     const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class BlockLatticePhysExternalParticleVelocity3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
public:
  BlockLatticePhysExternalParticleVelocity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class BlockLatticePhysExternal3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticePhysExternal3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                             T convFactorToPhysUnits, int offset, int size);
  bool operator() (T output[], const int input[]) override;
private:
  T _convFactorToPhysUnits;
  int _offset, _size;
};

/// functor returns pointwise lattice flux on local lattice
template <typename T, typename DESCRIPTOR>
class BlockLatticeFlux3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
public:
  BlockLatticeFlux3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[3], const int input[3]) override;
};

/// functor returns pointwise strain rate on local lattice, s_ij = 1/2*(du_idr_j + du_jdr_i)
template <typename T, typename DESCRIPTOR>
class BlockLatticeStrainRate3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
public:
  BlockLatticeStrainRate3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                           const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise phys strain rate on local lattice, s_ij = 1/2*(du_idr_j + du_jdr_i)
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysStrainRate3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  const int _overlap;
public:
  BlockLatticePhysStrainRate3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                               int overlap,
                               const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor returns pointwise phys force acting on a boundary with a given material on local lattice
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysBoundaryForce3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  BlockIndicatorF3D<T>&        _indicatorF;
  BlockGeometryStructure3D<T>& _blockGeometry;
public:
  BlockLatticePhysBoundaryForce3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                  BlockIndicatorF3D<T>& indicatorF,
                                  const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise phys force for PSM dynamics
template <typename T, typename DESCRIPTOR>
class BlockLatticePSMPhysForce3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
public:
  BlockLatticePSMPhysForce3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                             const UnitConverter<T,DESCRIPTOR>& converter, int mode_);
  bool operator() (T output[], const int input[]) override;
private:
  Mode mode;
};

/// functor returns pointwise phys wall shear stress acting on a boundary with a given material on local lattice
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysWallShearStress3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure3D<T>& _blockGeometry;
  const int _overlap;
  const int _material;
  std::vector<std::vector<std::vector<std::vector<int>>>> _discreteNormal;
  std::vector<std::vector<std::vector<std::vector<T>>>> _normal;
  T _physFactor;
public:
  BlockLatticePhysWallShearStress3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                    BlockGeometryStructure3D<T>& blockGeometry,
                                    int overlap,
                                    int material,
                                    const UnitConverter<T,DESCRIPTOR>& converter,
                                    IndicatorF3D<T>& indicator);
  bool operator() (T output[], const int input[]) override;
};

/**
 *  functor returns pointwise phys force acting on a indicated boundary on local lattice
 *  see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
 */
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysCorrBoundaryForce3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  BlockIndicatorF3D<T>&        _indicatorF;
  BlockGeometryStructure3D<T>& _blockGeometry;
public:
  BlockLatticePhysCorrBoundaryForce3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                      BlockIndicatorF3D<T>& indicatorF,
                                      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise, lattice-dependent external field
template <typename T, typename DESCRIPTOR, typename FIELD>
class BlockLatticeField3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticeField3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]) override;
};

/**
 *  functor returns pointwise, lattice-dependent porosity values in [0,1]
 *  in combination with (Extended)PorousBGKdynamics: 0->solid, 1->fluid
 */
template <typename T, typename DESCRIPTOR>
class BlockLatticePorosity3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticePorosity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]);
};

/**
 *  functor returns pointwise an approximation for the volume fraction
 */
template <typename T, typename DESCRIPTOR>
class BlockLatticeVolumeFractionApproximation3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure3D<T>& _blockGeometry;
  IndicatorF3D<T>& _indicator;
  int _refinementLevel;
  const UnitConverter<T,DESCRIPTOR>& _converter;
  bool _insideOut;
  T _physSubGridMinPhysRshift;
  T _physSubGridDeltaX;
  T _latticeSubGridVolume;
public:
  BlockLatticeVolumeFractionApproximation3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
      BlockGeometryStructure3D<T>& blockGeometry,
      IndicatorF3D<T>& indicator,
      int refinementLevel,
      const UnitConverter<T,DESCRIPTOR>& converter, bool insideOut);
  bool operator() (T output[], const int input[]);
};

/**
 *  functor to get pointwise mesh-independent permeability values in (0,inf)
 *  in combination with (Extended)PorousBGKdynamics
 *  note: result is cropped to 999999
 */
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysPermeability3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
public:
  BlockLatticePhysPermeability3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};


/**
 *  functor to get pointwise mesh-independent permeability values in (0,inf)
 *  in combination with (Extended)PorousBGKdynamics
 *  note: result is cropped to 1
 */
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysCroppedPermeability3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
public:
  BlockLatticePhysCroppedPermeability3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

//TODO: consistency with 2D (181219)
/*
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysPermeability3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  BlockGeometry3D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticePhysPermeability3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                 BlockGeometry3D<T>& blockGeometry,
                                 int material, const UnitConverter<T>& converter);
  bool operator() (T output[], const int input[]);
};*/


/// functor returns pointwise -nu/K*u on the lattice, can be used with BlockSum3D as objective
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysDarcyForce3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  BlockGeometry3D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticePhysDarcyForce3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                               BlockGeometry3D<T>& blockGeometry,
                               int material, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};


/// functor returns pointwise the l2-norm, e.g. of a velocity
template <typename T, typename DESCRIPTOR>
class BlockEuklidNorm3D final : public BlockF3D<T> {
protected:
  BlockF3D<T>& _f;
public:
  BlockEuklidNorm3D(BlockF3D<T>& f);
  bool operator() (T output[], const int input[]) override;
};


template <typename T, typename DESCRIPTOR>
class BlockLatticeInterpPhysVelocity3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
protected:
  Cuboid3D<T>* _cuboid;
  int _overlap;
public:
  BlockLatticeInterpPhysVelocity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                   const UnitConverter<T,DESCRIPTOR>& conv, Cuboid3D<T>* c, int overlap);
  BlockLatticeInterpPhysVelocity3D(const BlockLatticeInterpPhysVelocity3D<T,DESCRIPTOR>& rhs);
  bool operator() (T output[3], const int input[3]) override
  {
    return false;
  }
  void operator() (T output[3], const T input[3]);
};

/** Functor that returns forces acting on a particle surface, returns data in output for every particle in a row(described are return values for the first particle).
 * \return output[0]-output[2] translational force - physical units
 * \return output[3]-output[5] torque - physical units
 * \return output[7] number of voxels
 */

template <typename T, typename DESCRIPTOR>
class BlockLatticePorousMomentumLossForce3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure3D<T>& _blockGeometry;
  std::vector<SmoothIndicatorF3D<T,T,true>* >& _vectorOfIndicator;  
public:
  BlockLatticePorousMomentumLossForce3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                      BlockGeometryStructure3D<T>& blockGeometry,
                                      std::vector<SmoothIndicatorF3D<T,T,true>* >& indicator,
                                      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class BlockLatticePhysTemperature3D final : public BlockLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR> {
public:
  BlockLatticePhysTemperature3D(BlockLatticeStructure3D<T,TDESCRIPTOR>& blockLattice,
                                ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter);
  bool operator() (T output[], const int input[]);
};

/// BlockLatticePhysHeatFlux3D returns pointwise phys heat flux on local lattice.
template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class BlockLatticePhysHeatFlux3D final : public BlockLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR> {
public:
  BlockLatticePhysHeatFlux3D(BlockLatticeStructure3D<T,TDESCRIPTOR>& blockLattice,
                             const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
private:
  T _temp; // contains latticeSpecificHeatCapacity * (tau - 0.5) / tau
};

/// functor returns pointwise minimum distance to boundary given by indicators
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysBoundaryDistance3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure3D<T>& _blockGeometry;
  std::shared_ptr<IndicatorF3D<T>> _tmpIndicator = nullptr;
  std::vector<std::shared_ptr<IndicatorF3D<T>>> _indicatorList;
public:
  BlockLatticePhysBoundaryDistance3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                     BlockGeometryStructure3D<T>& blockGeometry,
                                     XMLreader const& xmlReader);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise pore radius for packings of spheres given by indicators
/// returns NAN for non-pore voxels
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysPoreSizeDistribution3D final  : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure3D<T>& _blockGeometry;
  int _material;
  std::shared_ptr<IndicatorF3D<T>> _tmpIndicator = nullptr;
  std::vector<std::shared_ptr<IndicatorF3D<T>>> _indicatorList;
  BlockLatticePhysBoundaryDistance3D<T,DESCRIPTOR> _distanceFunctor;
  BlockData3D<T,T> _distanceCache;
public:
  BlockLatticePhysPoreSizeDistribution3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                         BlockGeometryStructure3D<T>& blockGeometry, int material,
                                         XMLreader const& xmlReader);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise pore radius for packings of spheres given by indicators
/// returns NAN for non-pore voxels
template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class BlockLatticePhysTauFromBoundaryDistance3D final  : public BlockLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR> {
private:
  BlockGeometryStructure3D<T>& _blockGeometry;
  BlockLatticePhysBoundaryDistance3D<T,DESCRIPTOR> _distanceFunctor;
  const T _tmp1, _tmp2;
public:
  BlockLatticePhysTauFromBoundaryDistance3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                         BlockGeometryStructure3D<T>& blockGeometry,
                                         XMLreader const& xmlReader,
                                         ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter,
                                         const T p, const T T_avg, const T c_p, const T beta, const T lambda_0, const T sigma, const T p_0, const T n_0);
  bool operator() (T output[], const int input[]) override;
};

/// functor that returns 1 if SmoothIndicatorF A intersects IndicatorF B; otherwise, 0
template <typename T, typename DESCRIPTOR, bool HLBM>
class BlockLatticeIndicatorSmoothIndicatorIntersection3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure3D<T>& _blockGeometry;
  IndicatorF3D<T>& _normalInd;
  SmoothIndicatorF3D<T,T,HLBM>& _smoothInd;
public:
  BlockLatticeIndicatorSmoothIndicatorIntersection3D( 
                      BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                      BlockGeometryStructure3D<T>& blockGeometry,
                      IndicatorF3D<T>& normalInd, 
                      SmoothIndicatorF3D<T,T,HLBM>& smoothInd );
  bool operator() (T output[], const int input[]) override;
};

/// Returns pointwise porosity on local lattices for Guo & Zhao (2002)'s model.
template <typename T, typename DESCRIPTOR>
class BlockLatticeGuoZhaoEpsilon3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticeGuoZhaoEpsilon3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]) override;
};

/// Returns pointwise porous conductivity on local lattices for Guo & Zhao (2002)'s model.
template <typename T, typename DESCRIPTOR>
class BlockLatticeGuoZhaoPhysK3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
public:
  BlockLatticeGuoZhaoPhysK3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                             const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// Returns pointwise body force on local lattices for Guo & Zhao (2002)'s model.
template <typename T, typename DESCRIPTOR>
class BlockLatticeGuoZhaoPhysBodyForce3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
public:
  BlockLatticeGuoZhaoPhysBodyForce3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to scale particle distributions to a time step
template <typename T, typename DESCRIPTOR>
class BlockLatticeTimeStepScale3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  T _tau_old;
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  BlockLatticeTimeStepScale3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                              T oldTau, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

} // end namespace olb

#endif
