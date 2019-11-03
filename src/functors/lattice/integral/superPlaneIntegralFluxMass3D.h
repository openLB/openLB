/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerlaender
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

#ifndef SUPER_PLANE_INTEGRAL_FLUX_MASS_3D_H
#define SUPER_PLANE_INTEGRAL_FLUX_MASS_3D_H

#include "superPlaneIntegralF3D.h"

namespace olb {


/// Mass flux plane integral
/**
 * Calculates the flux integral of a 3-dimensional velocity functor multiplied
 * by a 1-dimensional density functor i.e. mass flux.
 *
 * Flux calculation is performed by SuperPlaneIntegralF3D<T>.
 * This class adds a print method and mass flux specific constructor wrappers.
 *
 * See SuperPlaneIntegralF3D for further documentation.
 **/
template<typename T>
class SuperPlaneIntegralFluxMass3D final : public SuperPlaneIntegralF3D<T> {
private:
  /// Velocity functor
  FunctorPtr<SuperF3D<T>> _velocityF;
  /// Density functor
  FunctorPtr<SuperF3D<T>> _densityF;

  /// Mass conversation factor
  const T _conversationFactorMass;
  /// Mass conversation factor
  const T _conversationFactorTime;

public:
  /// Primary constructor
  /**
   * \param velocityF
   *        (non-)owning pointer or reference to velocity functor.
   * \param densityF
   *        (non-)owning pointer or reference to density functor.
   * \param conversationFactorMass
   *        Mass conversation factor
   * \param conversationFactorTime
   *        Time conversation factor e.g. `converter.getConversionFactorTime()`
   * \param hyperplaneLattice
   *        Parametrization of the hyperplane lattice to be interpolated.
   * \param integrationIndicator
   *        (non-)owning pointer or reference to SuperIndicatorF3D<T>.
   *        Describes the set of lattice points relevant for integration.
   * \param subplaneIndicator
   *        (non-)owning pointer or reference to IndicatorF2D<T>.
   *        Describes the relevant subplane of the interpolated hyperplane.
   * \param mode
   *        Defines how the values of the discrete hyperplane are determined.
   *        i.e. if they are interpolated or read directly from lattice points.
   *        Note: BlockDataReductionMode::Analytical imposes restrictions on
   *        hyperplane definition and discretization. If you are not sure
   *        consider providing only a hyperplane defintion instead of both a
   *        definition and a discretization.
   **/
  SuperPlaneIntegralFluxMass3D(FunctorPtr<SuperF3D<T>>&& velocityF,
                               FunctorPtr<SuperF3D<T>>&& densityF,
                               SuperGeometry3D<T>&       geometry,
                               T conversationFactorMass,
                               T conversationFactorTime,
                               const HyperplaneLattice3D<T>& hyperplaneLattice,
                               FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
                               FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
                               BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  /// Constructor providing automatic lattice generation
  /**
   * \param velocityF
   *        (non-)owning pointer or reference to velocity functor.
   * \param densityF
   *        (non-)owning pointer or reference to density functor.
   * \param conversationFactorMass
   *        Mass conversation factor
   * \param conversationFactorTime
   *        Time conversation factor e.g. `converter.getConversionFactorTime()`
   * \param hyperplane
   *        Parametrization of the hyperplane to be integrated.
   *        The lattice resolution is set to CuboidGeometry3D<T>::getMinDeltaR.
   * \param integrationIndicator
   *        (non-)owning pointer or reference to SuperIndicatorF3D<T>.
   *        Describes the set of lattice points relevant for integration.
   * \param subplaneIndicator
   *        (non-)owning pointer or reference to IndicatorF2D<T>.
   *        Describes the relevant subplane of the interpolated hyperplane.
   * \param mode
   *        Defines how the values of the discrete hyperplane are determined.
   *        i.e. if they are interpolated or read directly from lattice points.
   **/
  SuperPlaneIntegralFluxMass3D(FunctorPtr<SuperF3D<T>>&& velocityF,
                               FunctorPtr<SuperF3D<T>>&& densityF,
                               SuperGeometry3D<T>&       geometry,
                               T conversationFactorMass,
                               T conversationFactorTime,
                               const Hyperplane3D<T>& hyperplane,
                               FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
                               FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
                               BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  /// Constructor providing automatic lattice generation and omitting subplane restriction
  /**
   * i.e. the intersection between geometry and hyperplane is integrated wherever _integrationIndicatorF allows.
   *
   * \param velocityF
   *        (non-)owning pointer or reference to velocity functor.
   * \param densityF
   *        (non-)owning pointer or reference to density functor.
   * \param conversationFactorMass
   *        Mass conversation factor
   * \param conversationFactorTime
   *        Time conversation factor e.g. `converter.getConversionFactorTime()`
   * \param hyperplane
   *        Parametrization of the hyperplane to be integrated.
   *        The lattice resolution is set to the cuboid geometry's minDeltaR.
   * \param integrationIndicator
   *        (non-)owning pointer or reference to SuperIndicatorF3D<T>.
   *        Describes the set of lattice points relevant for integration.
   * \param mode
   *        Defines how the values of the discrete hyperplane are determined.
   *        i.e. if they are interpolated or read directly from lattice points.
   **/
  SuperPlaneIntegralFluxMass3D(FunctorPtr<SuperF3D<T>>&& velocityF,
                               FunctorPtr<SuperF3D<T>>&& densityF,
                               SuperGeometry3D<T>&       geometry,
                               T conversationFactorMass,
                               T conversationFactorTime,
                               const Hyperplane3D<T>&  hyperplane,
                               FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
                               BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  /// Constructor providing automatic lattice and material indicator instantiation
  /**
   * \param velocityF (non-)owning pointer or reference to velocity functor.
   * \param densityF  (non-)owning pointer or reference to density functor.
   * \param conversationFactorMass Mass conversation factor
   * \param conversationFactorTime Time conversation factor
   * \param origin    hyperplane origin
   * \param u         hyperplane span vector
   * \param v         hyperplane span vector
   * \param materials material numbers relevant for hyperplane integration
   * \param mode      defines how the values of the discrete hyperplane are determined
   **/
  SuperPlaneIntegralFluxMass3D(FunctorPtr<SuperF3D<T>>&& velocityF,
                               FunctorPtr<SuperF3D<T>>&& densityF,
                               SuperGeometry3D<T>&       geometry,
                               T conversationFactorMass,
                               T conversationFactorTime,
                               const Vector<T,3>& origin,
                               const Vector<T,3>& u, const Vector<T,3>& v,
                               std::vector<int> materials,
                               BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  /// Constructor providing automatic lattice parametrization, only interpolating material 1
  /**
   * \param velocityF (non-)owning pointer or reference to velocity functor.
   * \param densityF  (non-)owning pointer or reference to density functor.
   * \param conversationFactorMass Mass conversation factor
   * \param conversationFactorTime Time conversation factor
   * \param origin hyperplane origin
   * \param u      hyperplane span vector
   * \param v      hyperplane span vector
   * \param mode   defines how the values of the discrete hyperplane are determined
   **/
  SuperPlaneIntegralFluxMass3D(FunctorPtr<SuperF3D<T>>&& velocityF,
                               FunctorPtr<SuperF3D<T>>&& densityF,
                               SuperGeometry3D<T>&       geometry,
                               T conversationFactorMass,
                               T conversationFactorTime,
                               const Vector<T,3>& origin,
                               const Vector<T,3>& u, const Vector<T,3>& v,
                               BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  /// Constructor providing automatic lattice parametrization to fit a given circle
  /**
   * \param velocityF (non-)owning pointer or reference to velocity functor.
   * \param densityF  (non-)owning pointer or reference to density functor.
   * \param conversationFactorMass Mass conversation factor
   * \param conversationFactorTime Time conversation factor
   * \param circle    circle indicator to be used for hyperplane subset parametrization
   * \param materials material numbers relevant for hyperplane interpolation
   * \param mode      defines how the values of the discrete hyperplane are determined
   **/
  SuperPlaneIntegralFluxMass3D(FunctorPtr<SuperF3D<T>>&& velocityF,
                               FunctorPtr<SuperF3D<T>>&& densityF,
                               SuperGeometry3D<T>&       geometry,
                               T conversationFactorMass,
                               T conversationFactorTime,
                               IndicatorCircle3D<T>& circle,
                               std::vector<int> materials,
                               BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  /// Constructor providing automatic lattice parametrization to fit a given circle and only interpolating material 1
  /**
   * \param velocityF (non-)owning pointer or reference to velocity functor.
   * \param densityF  (non-)owning pointer or reference to density functor.
   * \param conversationFactorMass Mass conversation factor
   * \param conversationFactorTime Time conversation factor
   * \param circle    circle indicator to be used for hyperplane subset parametrization
   * \param materials material numbers relevant for hyperplane interpolation
   * \param mode      defines how the values of the discrete hyperplane are determined
   **/
  SuperPlaneIntegralFluxMass3D(FunctorPtr<SuperF3D<T>>&& velocityF,
                               FunctorPtr<SuperF3D<T>>&& densityF,
                               SuperGeometry3D<T>&       geometry,
                               T conversationFactorMass,
                               T conversationFactorTime,
                               IndicatorCircle3D<T>& circle,
                               BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  bool operator() (T output[], const int input[]) override;

  void print(std::string regionName, std::string massFluxSiScaleName);

};


}

#endif
