/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#ifndef SUPER_LATTICE_LOCAL_F_3D_H
#define SUPER_LATTICE_LOCAL_F_3D_H

#include<vector>

#include "superBaseF3D.h"
#include "superCalcF3D.h"
#include "blockLatticeLocalF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "core/superLattice3D.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


////////////////////////////////////////////////////////////////////////////////
//////if globIC is not on the local processor, the returned vector is empty/////
////////////////////////////////////////////////////////////////////////////////

/// functor to get pointwise f population on local lattices
template <typename T, typename DESCRIPTOR>
class SuperLatticeFpop3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeFpop3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise dissipation density on local lattices
template <typename T, typename DESCRIPTOR>
class SuperLatticeDissipation3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticeDissipation3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise dissipation density on local lattices
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysDissipation3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysDissipation3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise turbulent dissipation density on local lattices
template <typename T, typename DESCRIPTOR>
class SuperLatticeEffevtiveDissipation3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticeEffevtiveDissipation3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst,
                                     LESDynamics<T, DESCRIPTOR>& LESdynamics);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise turbulent dissipation density on local lattices
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysEffevtiveDissipation3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysEffevtiveDissipation3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                         const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst,
                                         LESDynamics<T, DESCRIPTOR>& LESdynamics);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise density rho on local lattices
template <typename T, typename DESCRIPTOR>
class SuperLatticeDensity3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeDensity3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise velocity on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticeVelocity3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeVelocity3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise velocity on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticeExternalVelocity3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeExternalVelocity3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise strain rate on local lattice
/// s_ij = 1/2*(du_idr_j + du_jdr_i)
template <typename T, typename DESCRIPTOR>
class SuperLatticeStrainRate3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticeStrainRate3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                           const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise flux on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticeFlux3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeFlux3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys strain rate on local lattice
/// s_ij = 1/2*(du_idr_j + du_jdr_i)
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysStrainRate3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysStrainRate3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                               const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise the material no. presenting the geometry on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticeGeometry3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticeGeometry3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                         SuperGeometry3D<T>& superGeometry, const int material = -1);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise the rank no. + 1 on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticeRank3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeRank3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise the cuboid no. + 1 on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticeCuboid3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeCuboid3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise phys pressure from rho on local lattices
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysPressure3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysPressure3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                             const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise phys velocity on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysVelocity3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysVelocity3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                             const UnitConverter<T,DESCRIPTOR>& converter, bool print=false);
  bool operator() (T output[], const int input[]) override;
private:
  bool _print;
};

template <typename T, typename DESCRIPTOR>
class SuperLatticePhysExternalPorosity3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysExternalPorosity3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class SuperLatticePhysExternalVelocity3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysExternalVelocity3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class SuperLatticePhysExternalParticleVelocity3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysExternalParticleVelocity3D(SuperLattice3D<T,DESCRIPTOR>& blockLattice,
      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class SuperLatticePhysExternal3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysExternal3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                             T convFactorToPhysUnits,
                             int offset, int size);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys force acting on a boundary with a given material on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysBoundaryForce3D : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
public:
  SuperLatticePhysBoundaryForce3D(SuperLattice3D<T,DESCRIPTOR>&      sLattice,
                                  FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
                                  const UnitConverter<T,DESCRIPTOR>& converter);
  SuperLatticePhysBoundaryForce3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                  SuperGeometry3D<T>& superGeometry, const int material,
                                  const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys force for the PSM dynamics
template <typename T, typename DESCRIPTOR>
class SuperLatticePSMPhysForce3D : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePSMPhysForce3D(SuperLattice3D<T,DESCRIPTOR>&      sLattice,
                             const UnitConverter<T,DESCRIPTOR>& converter,
                             int mode_=0);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys wall shear stress with a given material on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysWallShearStress3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysWallShearStress3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                    SuperGeometry3D<T>& superGeometry, const int material,
                                    const UnitConverter<T,DESCRIPTOR>& converter,
                                    IndicatorF3D<T>& indicator);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys force acting on a boundary with a given material on local lattice
/// see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysCorrBoundaryForce3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
public:
  SuperLatticePhysCorrBoundaryForce3D(SuperLattice3D<T,DESCRIPTOR>&      sLattice,
                                      FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF,
                                      const UnitConverter<T,DESCRIPTOR>& converter);
  SuperLatticePhysCorrBoundaryForce3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                      SuperGeometry3D<T>& superGeometry, const int material,
                                      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise, lattice-dependent external field
template <typename T, typename DESCRIPTOR, typename FIELD>
class SuperLatticeField3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeField3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise, lattice-dependent porosity values in [0,1]
/// in combination with (Extended)PorousBGKdynamics: 0->solid, 1->fluid
template <typename T, typename DESCRIPTOR>
class SuperLatticePorosity3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
public:
  SuperLatticePorosity3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise an approx. for the volume fraction
template <typename T, typename DESCRIPTOR>
class SuperLatticeVolumeFractionApproximation3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeVolumeFractionApproximation3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, SuperGeometry3D<T>& superGeometry,
      IndicatorF3D<T>& indicator, int refinementLevel, const UnitConverter<T,DESCRIPTOR>& converter, bool insideOut = false);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise mesh-independent permeability values in (0,inf) in combination with (Extended)PorousBGKdynamics
/// note: result is cropped to 999999
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysPermeability3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysPermeability3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise mesh-independent permeability values in (0,inf) in combination with (Extended)PorousBGKdynamics
/// note: result is cropped to 1
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysCroppedPermeability3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysCroppedPermeability3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};


/// computes pointwise -nu/K*u on the lattice, can be used with SuperSum3D as objective
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysDarcyForce3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysDarcyForce3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                               SuperGeometry3D<T>& superGeometry,
                               const int material, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};


/// functor that returns pointwise the l2-norm, e.g. of a velocity
template <typename T, typename DESCRIPTOR>
class SuperEuklidNorm3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticeF3D<T,DESCRIPTOR>& _f;
public:
  SuperEuklidNorm3D(SuperLatticeF3D<T,DESCRIPTOR>& f);
  bool operator() (T output[], const int input[]) override;
};


template <typename T, typename DESCRIPTOR>
class SuperLatticeInterpPhysVelocity3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticeInterpPhysVelocity3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, UnitConverter<T,DESCRIPTOR> const& converter);
  bool operator()(T output[], const int input[]) override;
  void operator()(T output[], const T input[], const int iC);
};

/// functor to get pointwise phys temperature from rho on local lattices
template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class SuperLatticePhysTemperature3D final : public SuperLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR> { // templatename before <
public:
  SuperLatticePhysTemperature3D(SuperLattice3D<T,TDESCRIPTOR>& sLattice,
                                ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise heat flux on local lattice
template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class SuperLatticePhysHeatFlux3D final : public SuperLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR> {
public:
  SuperLatticePhysHeatFlux3D(SuperLattice3D<T,TDESCRIPTOR>& sLattice,
                             const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/** Functor that returns forces acting on a particle surface, returns data in output for every particle in a row(described are return values for the first particle).
 * \return output[0]-output[2] translational force - physical units
 * \return output[3]-output[5] torque - physical units
 * \return output[7] number of voxels
 */
template <typename T, typename DESCRIPTOR>
class SuperLatticePorousMomentumLossForce3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePorousMomentumLossForce3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                      SuperGeometry3D<T>& superGeometry,
                                      std::vector<SmoothIndicatorF3D<T,T,true>* >& indicator,
                                      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// functor that returns the minimum distance (in m) to a set of indicators given by an xmlReader
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysBoundaryDistance3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
public:
  SuperLatticePhysBoundaryDistance3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                     SuperGeometry3D<T>& superGeometry,
                                     XMLreader const& xmlReader);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise pore radius (in m) for packings of spheres given by an xmlReader
/// returns NAN for non-pore voxels
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysPoreSizeDistribution3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
public:
  SuperLatticePhysPoreSizeDistribution3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                         SuperGeometry3D<T>& superGeometry,
                                         int material,
                                         XMLreader const& xmlReader);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise pore radius (in m) for packings of spheres given by an xmlReader
/// returns NAN for non-pore voxels
template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class SuperLatticePhysTauFromBoundaryDistance3D final : public SuperLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR> {
public:
  SuperLatticePhysTauFromBoundaryDistance3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                         SuperGeometry3D<T>& sGeometry,
                                         XMLreader const& xmlReader,
                                         ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter,
                                         const T p, const T T_avg, const T c_p, const T beta, const T lambda_0, const T sigma, const T p_0, const T n_0);
  bool operator() (T output[], const int input[]) override;
};

/// functor that returns 1 if SmoothIndicatorF A intersects IndicatorF B; otherwise, 0
template <typename T, typename DESCRIPTOR, bool HLBM>
class SuperLatticeIndicatorSmoothIndicatorIntersection3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeIndicatorSmoothIndicatorIntersection3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                             SuperGeometry3D<T>& superGeometry,
                             IndicatorF3D<T>& normalInd, SmoothIndicatorF3D<T,T,HLBM>& smoothInd );
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise porosity on local lattice for Guo & Zhao (2002)'s model
template <typename T, typename DESCRIPTOR>
class SuperLatticeGuoZhaoEpsilon3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeGuoZhaoEpsilon3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise porous conductivity on local lattice for Guo & Zhao (2002)'s model
template <typename T, typename DESCRIPTOR>
class SuperLatticeGuoZhaoPhysK3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticeGuoZhaoPhysK3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                             const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise body force on local lattice for Guo & Zhao (2002)'s model
template <typename T, typename DESCRIPTOR>
class SuperLatticeGuoZhaoPhysBodyForce3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticeGuoZhaoPhysBodyForce3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to scale particle distributions to a time step
template <typename T, typename DESCRIPTOR>
class SuperLatticeTimeStepScale3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeTimeStepScale3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                             T oldTau, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

} // end namespace olb

#endif
