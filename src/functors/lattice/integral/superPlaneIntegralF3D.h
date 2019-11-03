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

#ifndef SUPER_PLANE_INTEGRAL_F_3D_H
#define SUPER_PLANE_INTEGRAL_F_3D_H

#include "core/superLattice3D.h"
#include "core/vector.h"
#include "geometry/superGeometry3D.h"
#include "functors/lattice/superBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "functors/lattice/indicator/superIndicatorF3D.h"
#include "functors/lattice/indicator/superIndicatorBaseF3D.h"
#include "functors/lattice/blockReduction3D2D.h"
#include "functors/lattice/indicator/indicator2D.h"
#include "utilities/functorPtr.h"
#include "utilities/hyperplane3D.h"
#include "utilities/hyperplaneLattice3D.h"

namespace olb {


/// Surface integral of a subset of a interpolated hyperplane
template<typename T>
class SuperPlaneIntegralF3D : public SuperF3D<T> {
protected:
  SuperGeometry3D<T>& _geometry;

  /// Functor to be integrated on the plane
  FunctorPtr<SuperF3D<T>> _f;
  /// Indicator describing relevant discrete integration points
  FunctorPtr<SuperIndicatorF3D<T>> _integrationIndicatorF;
  /// Indicator describing the relevant subset of the interpolated hyperplane
  FunctorPtr<IndicatorF2D<T>> _subplaneIndicatorF;
  /// Functor describing plane to be interpolated and integrated
  BlockReduction3D2D<T> _reductionF;

  /// Origin vector as given by hyperplane definition, (0,0) in respect to the
  /// subplane indicator _subplaneIndicatorF
  /**
   * Note: The reduced plane _reductionF calculates its own origin based on _origin
   * as its spans the maximum size possible given the span vectors, origin and geometry
   * size.
   **/
  Vector<T,3> _origin;
  /// Span vector u as given by hyperplane definition, normalized to h
  Vector<T,3> _u;
  /// Span vector v as given by hyperplane definition, normalized to h
  Vector<T,3> _v;
  /// Orthogonal vector to _u and _v
  Vector<T,3> _normal;

  /// Subset of the discrete plane points given by _reductionF as indicated
  /// by _integrationIndicatorF. i.e. the points used to interpolate the hyperplane
  std::vector<std::tuple<int, int>> _rankLocalSubplane;

  /// \returns true iff the given physical position is to be integrated
  /**
   * This is determined using the _integrationIndicatorF indicated subset of the 2d
   * plane reduced by _reductionF.
   **/
  bool isToBeIntegrated(const Vector<T,3>& physR, int iC);

public:
  /// Primary constructor
  /**
   * All other constructors defer the actual construction to this constructor.
   *
   * \param f
   *        (non-)owning pointer or reference to SuperF3D<T>.
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
  SuperPlaneIntegralF3D(FunctorPtr<SuperF3D<T>>&& f,
                        SuperGeometry3D<T>&       geometry,
                        const HyperplaneLattice3D<T>& hyperplaneLattice,
                        FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
                        FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
                        BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  /// Constructor providing automatic lattice generation
  /**
   * \param f
   *        (non-)owning pointer or reference to SuperF3D<T>.
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
  SuperPlaneIntegralF3D(FunctorPtr<SuperF3D<T>>&& f,
                        SuperGeometry3D<T>&       geometry,
                        const Hyperplane3D<T>&    hyperplane,
                        FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
                        FunctorPtr<IndicatorF2D<T>>&&      subplaneIndicator,
                        BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  /// Constructor providing automatic lattice generation and omitting subplane restriction
  /**
   * i.e. the intersection between geometry and hyperplane is integrated wherever _integrationIndicatorF allows.
   *
   * \param f
   *        (non-)owning pointer or reference to SuperF3D<T>.
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
  SuperPlaneIntegralF3D(FunctorPtr<SuperF3D<T>>&& f,
                        SuperGeometry3D<T>&       geometry,
                        const Hyperplane3D<T>&    hyperplane,
                        FunctorPtr<SuperIndicatorF3D<T>>&& integrationIndicator,
                        BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  /// Constructor providing automatic lattice and material indicator instantiation
  /**
   * \param f         (non-)owning pointer or reference to SuperF3D<T>.
   * \param origin    hyperplane origin
   * \param u         hyperplane span vector
   * \param v         hyperplane span vector
   * \param materials material numbers relevant for hyperplane integration
   * \param mode      defines how the values of the discrete hyperplane are determined
   **/
  SuperPlaneIntegralF3D(FunctorPtr<SuperF3D<T>>&& f,
                        SuperGeometry3D<T>& geometry,
                        const Vector<T,3>& origin,
                        const Vector<T,3>& u, const Vector<T,3>& v,
                        std::vector<int> materials,
                        BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  /// Constructor providing automatic lattice parametrization, only interpolating material 1
  /**
   * \param f      (non-)owning pointer or reference to SuperF3D<T>.
   * \param origin hyperplane origin
   * \param u      hyperplane span vector
   * \param v      hyperplane span vector
   * \param mode   Defines how the values of the discrete hyperplane are determined.
   **/
  SuperPlaneIntegralF3D(FunctorPtr<SuperF3D<T>>&& f,
                        SuperGeometry3D<T>& geometry,
                        const Vector<T,3>& origin,
                        const Vector<T,3>& u, const Vector<T,3>& v,
                        BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  /// Constructor providing automatic lattice and material indicator instantiation
  /**
   * \param f         (non-)owning pointer or reference to SuperF3D<T>.
   * \param origin    hyperplane origin
   * \param normal    hyperplane normal
   * \param materials material numbers relevant for hyperplane integration
   * \param mode      defines how the values of the discrete hyperplane are determined
   **/
  SuperPlaneIntegralF3D(FunctorPtr<SuperF3D<T>>&& f,
                        SuperGeometry3D<T>& geometry,
                        const Vector<T,3>& origin,
                        const Vector<T,3>& normal,
                        std::vector<int> materials,
                        BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  /// Constructor providing automatic lattice parametrization, only interpolating material 1
  /**
   * \param f      (non-)owning pointer or reference to SuperF3D<T>.
   * \param origin hyperplane origin
   * \param normal hyperplane normal
   * \param mode   defines how the values of the discrete hyperplane are determined
   **/
  SuperPlaneIntegralF3D(FunctorPtr<SuperF3D<T>>&& f,
                        SuperGeometry3D<T>& geometry,
                        const Vector<T,3>& origin,
                        const Vector<T,3>& normal,
                        BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  /// Constructor providing automatic lattice and material indicator instantiation
  /**
   * \param f         (non-)owning pointer or reference to SuperF3D<T>.
   * \param normal    hyperplane normal (centered in mother cuboid)
   * \param materials material numbers relevant for hyperplane integration
   * \param mode      defines how the values of the discrete hyperplane are determined
   **/
  SuperPlaneIntegralF3D(FunctorPtr<SuperF3D<T>>&& f,
                        SuperGeometry3D<T>& geometry,
                        const Vector<T,3>& normal,
                        std::vector<int> materials,
                        BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  /// Constructor providing automatic lattice parametrization, only interpolating material 1
  /**
   * \param f      (non-)owning pointer or reference to SuperF3D<T>.
   * \param normal hyperplane normal (centered in mother cuboid)
   * \param mode   defines how the values of the discrete hyperplane are determined
   **/
  SuperPlaneIntegralF3D(FunctorPtr<SuperF3D<T>>&& f,
                        SuperGeometry3D<T>& geometry,
                        const Vector<T,3>& normal,
                        BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  /// Constructor providing automatic lattice parametrization to fit a given circle
  /**
   * \param f         (non-)owning pointer or reference to SuperF3D<T>.
   * \param circle    circle indicator to be used for hyperplane subset parametrization
   * \param materials material numbers relevant for hyperplane interpolation
   * \param mode      defines how the values of the discrete hyperplane are determined
   **/
  SuperPlaneIntegralF3D(FunctorPtr<SuperF3D<T>>&& f,
                        SuperGeometry3D<T>& geometry,
                        const IndicatorCircle3D<T>& circle,
                        std::vector<int> materials,
                        BlockDataReductionMode mode=BlockDataReductionMode::Analytical);
  /// Constructor providing automatic lattice parametrization to fit a given circle and only interpolating material 1
  /**
   * \param f      (non-)owning pointer or reference to SuperF3D<T>.
   * \param circle circle indicator to be used for hyperplane subset parametrization
   * \param mode   defines how the values of the discrete hyperplane are determined
   **/
  SuperPlaneIntegralF3D(FunctorPtr<SuperF3D<T>>&& f,
                        SuperGeometry3D<T>& geometry,
                        const IndicatorCircle3D<T>& circle,
                        BlockDataReductionMode mode=BlockDataReductionMode::Analytical);

  /**
   * Returns the plane integral in the following structure:
   *
   * \code
   * output[0] = integral, e.g. flow[0] * h^2 for 1-dimensional target sizes
   * output[1] = #voxels * h^2 i.e. area
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
