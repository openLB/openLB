/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#ifndef SUPER_PLANE_INTEGRAL_F_2D_H
#define SUPER_PLANE_INTEGRAL_F_2D_H

#include "core/vector.h"
#include "io/ostreamManager.h"
#include "geometry/superGeometry2D.h"
#include "functors/lattice/superBaseF2D.h"
#include "functors/analytical/indicator/indicatorBaseF2D.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"
#include "functors/lattice/indicator/superIndicatorBaseF2D.h"
#include "functors/lattice/blockReduction2D1D.h"
#include "utilities/functorPtr.h"
#include "utilities/hyperplane2D.h"
#include "utilities/hyperplaneLattice2D.h"

namespace olb {


/// Surface integral of a subset of a interpolated hyperplane
template<typename T>
class SuperPlaneIntegralF2D : public SuperF2D<T> {
protected:
  SuperGeometry2D<T>& _geometry;

  /// Functor to be integrated on the line
  FunctorPtr<SuperF2D<T>> _f;
  /// Indicator describing relevant discrete integration points
  FunctorPtr<SuperIndicatorF2D<T>> _integrationIndicatorF;
  /// Indicator describing the relevant subset of the interpolated hyperplane
  FunctorPtr<IndicatorF2D<T>> _subplaneIndicatorF;
  /// Functor describing line to be interpolated and integrated
  BlockReduction2D1D<T> _reductionF;

  /// Origin vector as given by hyperplane definition, (0,0) in respect to the
  /// subplane indicator _subplaneIndicatorF
  /**
   * Note: The reduced plane _reductionF calculates its own origin based on _origin
   * as its spans the maximum size possible given the direction vector, origin and
   * geometry size.
   **/
  Vector<T,2> _origin;
  /// Direction vector u as given by hyperplane definition, normalized to h
  Vector<T,2> _u;
  /// Orthogonal vector to _u
  Vector<T,2> _normal;

  /// Subset of the discrete line points given by _reductionF as indicated
  /// by _integrationIndicatorF. i.e. the points used to interpolate the hyperplane
  std::vector<int> _rankLocalSubplane;

  /// \returns true iff the given physical position is to be integrated
  /**
   * This is determined using the _integrationIndicatorF indicated subset of the 2d
   * hyperplane reduced by _reductionF.
   **/
  bool isToBeIntegrated(const Vector<T,2>& physR, int iC);

public:
  /// Primary constructor
  /**
   * All other constructors defer the actual construction to this constructor.
   *
   * \param f
   *        (non-)owning pointer or reference to SuperF2D<T>.
   * \param hyperplaneLattice
   *        Parametrization of the hyperplane lattice to be interpolated.
   * \param integrationIndicator
   *        (non-)owning pointer or reference to SuperIndicatorF2D<T>.
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
  SuperPlaneIntegralF2D(FunctorPtr<SuperF2D<T>>&& f,
                        SuperGeometry2D<T>&       geometry,
                        const HyperplaneLattice2D<T>& hyperplaneLattice,
                        FunctorPtr<SuperIndicatorF2D<T>>&& integrationIndicator,
                        FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
                        BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  /// Constructor providing automatic lattice generation
  /**
   * \param f
   *        (non-)owning pointer or reference to SuperF2D<T>.
   * \param hyperplane
   *        Parametrization of the hyperplane to be integrated.
   *        The lattice resolution is set to CuboidGeometry2D<T>::getMinDeltaR.
   * \param integrationIndicator
   *        (non-)owning pointer or reference to SuperIndicatorF2D<T>.
   *        Describes the set of lattice points relevant for integration.
   * \param subplaneIndicator
   *        (non-)owning pointer or reference to IndicatorF2D<T>.
   *        Describes the relevant subplane of the interpolated hyperplane.
   * \param mode
   *        Defines how the values of the discrete hyperplane are determined.
   *        i.e. if they are interpolated or read directly from lattice points.
   **/
  SuperPlaneIntegralF2D(FunctorPtr<SuperF2D<T>>&& f,
                        SuperGeometry2D<T>&       geometry,
                        const Hyperplane2D<T>&    hyperplane,
                        FunctorPtr<SuperIndicatorF2D<T>>&& integrationIndicator,
                        FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
                        BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  /// Constructor providing automatic lattice generation and omitting subplane restriction
  /**
   * i.e. the intersection between geometry and hyperplane is integrated wherever _integrationIndicatorF allows.
   *
   * \param f
   *        (non-)owning pointer or reference to SuperF2D<T>.
   * \param hyperplane
   *        Parametrization of the hyperplane to be integrated.
   *        The lattice resolution is set to the cuboid geometry's minDeltaR.
   * \param integrationIndicator
   *        (non-)owning pointer or reference to SuperIndicatorF2D<T>.
   *        Describes the set of lattice points relevant for integration.
   * \param mode
   *        Defines how the values of the discrete hyperplane are determined.
   *        i.e. if they are interpolated or read directly from lattice points.
   **/
  SuperPlaneIntegralF2D(FunctorPtr<SuperF2D<T>>&& f,
                        SuperGeometry2D<T>&       geometry,
                        const Hyperplane2D<T>&    hyperplane,
                        FunctorPtr<SuperIndicatorF2D<T>>&& integrationIndicator,
                        BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  /// Constructor providing automatic lattice and material indicator instantiation
  /**
   * \param f         (non-)owning pointer or reference to SuperF2D<T>.
   * \param origin    hyperplane origin origin
   * \param u         hyperplane direction vector
   * \param materials material numbers relevant for hyperplane integration
   * \param mode      defines how the values of the discrete hyperplane are determined
   **/
  SuperPlaneIntegralF2D(FunctorPtr<SuperF2D<T>>&& f,
                        SuperGeometry2D<T>& geometry,
                        const Vector<T,2>& origin,
                        const Vector<T,2>& u,
                        std::vector<int> materials,
                        BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  /// Constructor providing automatic lattice parametrization, only interpolating material 1
  /**
   * \param f      (non-)owning pointer or reference to SuperF2D<T>.
   * \param origin hyperplane origin
   * \param u      hyperplane direction vector
   * \param mode   defines how the values of the discrete hyperplane are determined
   **/
  SuperPlaneIntegralF2D(FunctorPtr<SuperF2D<T>>&& f,
                        SuperGeometry2D<T>& geometry,
                        const Vector<T,2>& origin,
                        const Vector<T,2>& u,
                        BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  /**
   * Returns the line integral in the following structure:
   *
   * \code
   * output[0] = integral, e.g. flow[0] * h for 1-dimensional target sizes
   * output[1] = #voxels * h i.e. length of line
   * output[2..2+f.getTargetSize()] = flow
   * \endcode
   *
   * Note: output[0] contains the flux value if applicable
   *
   * \param input irrelevant
   **/
  bool operator() (T output[], const int input[]) override;

};


}

#endif
