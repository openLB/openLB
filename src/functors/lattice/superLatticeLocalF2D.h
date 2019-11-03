/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause
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

#ifndef SUPER_LATTICE_LOCAL_F_2D_H
#define SUPER_LATTICE_LOCAL_F_2D_H

#include <vector>

#include "superBaseF2D.h"
#include "core/superLattice2D.h"
#include "indicator/superIndicatorBaseF2D.h"
#include "utilities/functorPtr.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


template<typename T> class SuperGeometry2D;

////////////////////////////////////////////////////////////////////////////////
//////if globIC is not on the local processor, the returned vector is empty/////
////////////////////////////////////////////////////////////////////////////////


/// functor to get pointwise dissipation density on local lattices
template <typename T, typename DESCRIPTOR>
class SuperLatticeDissipation2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
private:
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticeDissipation2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise dissipation density on local lattices
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysDissipation2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
public:
  SuperLatticePhysDissipation2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                                const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise density rho on local lattices
template <typename T, typename DESCRIPTOR>
class SuperLatticeDensity2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
public:
  SuperLatticeDensity2D(SuperLattice2D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise velocity on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticeVelocity2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
public:
  SuperLatticeVelocity2D(SuperLattice2D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys strain rate on local lattice
/// s_ij = 1/2*(du_idr_j + du_jdr_i)
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysStrainRate2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
public:
  SuperLatticePhysStrainRate2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                               const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys wall shear stress with a given material on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysWallShearStress2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysWallShearStress2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                                    SuperGeometry2D<T>& superGeometry, const int material,
                                    const UnitConverter<T,DESCRIPTOR>& converter,
                                    IndicatorF2D<T>& indicator);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise the material no. presenting the geometry on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticeGeometry2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
private:
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticeGeometry2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                         SuperGeometry2D<T>& superGeometry, const int material = -1);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise the rank no. + 1 on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticeRank2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
public:
  SuperLatticeRank2D(SuperLattice2D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise the cuboid no. + 1 on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticeCuboid2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
public:
  SuperLatticeCuboid2D(SuperLattice2D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise phys pressure from rho on local lattices
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysPressure2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
public:
  SuperLatticePhysPressure2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                             const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise phys velocity on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysVelocity2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
public:
  SuperLatticePhysVelocity2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                             const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

template <typename T, typename DESCRIPTOR>
class SuperLatticePhysExternalPorosity2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
public:
  SuperLatticePhysExternalPorosity2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR, typename FIELD>
class SuperLatticeField2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
public:
  SuperLatticeField2D(SuperLattice2D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class SuperLatticePhysExternalVelocity2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
public:
  SuperLatticePhysExternalVelocity2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, typename DESCRIPTOR>
class SuperLatticePhysExternalParticleVelocity2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
public:
  SuperLatticePhysExternalParticleVelocity2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise phys force acting on a boundary with a given material on local lattice
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysBoundaryForce2D : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  FunctorPtr<SuperIndicatorF2D<T>> _indicatorF;
public:
  SuperLatticePhysBoundaryForce2D(SuperLattice2D<T,DESCRIPTOR>&      sLattice,
                                  FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF,
                                  const UnitConverter<T,DESCRIPTOR>& converter);
  SuperLatticePhysBoundaryForce2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                                  SuperGeometry2D<T>& superGeometry, const int material,
                                  const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys force for the PSM dynamics
template <typename T, typename DESCRIPTOR>
class SuperLatticePSMPhysForce2D : public SuperLatticePhysF2D<T,DESCRIPTOR> {
public:
  SuperLatticePSMPhysForce2D(SuperLattice2D<T,DESCRIPTOR>&      sLattice,
                             const UnitConverter<T,DESCRIPTOR>& converter,
                             int mode_);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys force acting on a boundary with a given material on local lattice
/// see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysCorrBoundaryForce2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysCorrBoundaryForce2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                                      SuperGeometry2D<T>& superGeometry, const int material,
                                      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise, lattice-dependent porosity values in [0,1]
/// in combination with (Extended)PorousBGKdynamics: 0->solid, 1->fluid
template <typename T, typename DESCRIPTOR>
class SuperLatticePorosity2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePorosity2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                         SuperGeometry2D<T>& superGeometry, const int material,
                         const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise an approx. for the volume fraction
template <typename T, typename DESCRIPTOR>
class SuperLatticeVolumeFractionApproximation2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
public:
  SuperLatticeVolumeFractionApproximation2D(SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
      IndicatorF2D<T>& indicator, int refinementLevel, const UnitConverter<T,DESCRIPTOR>& converter, bool insideOut = false);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise an approx. for the volume fraction
template <typename T, typename DESCRIPTOR>
class SuperLatticeVolumeFractionPolygonApproximation2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
public:
  SuperLatticeVolumeFractionPolygonApproximation2D(SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
      IndicatorF2D<T>& indicator, const UnitConverter<T,DESCRIPTOR>& converter, bool insideOut = false);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise mesh-independent permeability values in (0,inf)  in combination with (Extended)PorousBGKdynamics
/// note: result is cropped to 999999
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysPermeability2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysPermeability2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                                 SuperGeometry2D<T>& superGeometry,
                                 const int material, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// computes pointwise -nu/K*u on the lattice, can be used with SuperSum2D as objective
template <typename T, typename DESCRIPTOR>
class SuperLatticePhysDarcyForce2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
private:
  SuperGeometry2D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysDarcyForce2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                               SuperGeometry2D<T>& superGeometry, const int material,
                               const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// functor that returns pointwise the l2-norm, e.g. of a velocity
template <typename T, typename DESCRIPTOR>
class SuperEuklidNorm2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
private:
  SuperLatticeF2D<T,DESCRIPTOR>& _f;
public:
  SuperEuklidNorm2D(SuperLatticeF2D<T,DESCRIPTOR>& f);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise phys temperature from rho on local lattices
template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class SuperLatticePhysTemperature2D final : public SuperLatticeThermalPhysF2D<T,DESCRIPTOR,TDESCRIPTOR> {
public:
  SuperLatticePhysTemperature2D(SuperLattice2D<T,TDESCRIPTOR>& sLattice,
                                ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter);
  bool operator() (T output[], const int input[]);
};


/** Functor that returns forces acting on a particle surface, returns data in output for every particle in a row (described are return values for the first particle).
 * \return output[0]-output[1] translational force - physical units
 * \return output[3] torque - physical units
 * \return output[4] number of voxels
 */
template <typename T, typename DESCRIPTOR>
class SuperLatticePorousMomentumLossForce2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
public:
  SuperLatticePorousMomentumLossForce2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                                      SuperGeometry2D<T>& superGeometry,
                                      std::vector<SmoothIndicatorF2D<T,T,true>* >& indicator,
                                      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise heat flux on local lattice
template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class SuperLatticePhysHeatFlux2D final : public SuperLatticeThermalPhysF2D<T,DESCRIPTOR,TDESCRIPTOR> {
public:
  SuperLatticePhysHeatFlux2D(SuperLattice2D<T,TDESCRIPTOR>& sLattice,
                             const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor that returns 1 if SmoothIndicatorF A intersects IndicatorF B; otherwise, 0
template <typename T, typename DESCRIPTOR, bool HLBM>
class SuperLatticeIndicatorSmoothIndicatorIntersection2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
public:
  SuperLatticeIndicatorSmoothIndicatorIntersection2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                             SuperGeometry2D<T>& superGeometry,
                             IndicatorF2D<T>& normalInd, SmoothIndicatorF2D<T,T,HLBM>& smoothInd );
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise porosity on local lattice for Guo & Zhao (2002)'s model
template <typename T, typename DESCRIPTOR>
class SuperLatticeGuoZhaoEpsilon2D final : public SuperLatticeF2D<T,DESCRIPTOR> {
public:
  SuperLatticeGuoZhaoEpsilon2D(SuperLattice2D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise porous conductivity on local lattice for Guo & Zhao (2002)'s model
template <typename T, typename DESCRIPTOR>
class SuperLatticeGuoZhaoPhysK2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
public:
  SuperLatticeGuoZhaoPhysK2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                             const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise body force on local lattice for Guo & Zhao (2002)'s model
template <typename T, typename DESCRIPTOR>
class SuperLatticeGuoZhaoPhysBodyForce2D final : public SuperLatticePhysF2D<T,DESCRIPTOR> {
public:
  SuperLatticeGuoZhaoPhysBodyForce2D(SuperLattice2D<T,DESCRIPTOR>& sLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

} // end namespace olb
#endif
