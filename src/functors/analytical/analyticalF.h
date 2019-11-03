/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause, Albert Mink
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

#ifndef ANALYTICAL_F_H
#define ANALYTICAL_F_H

#include <vector>

#include "analyticalBaseF.h"
#include "indicator/smoothIndicatorF2D.h"
#include "indicator/smoothIndicatorF3D.h"


/**
 *  The functor dimensions are given by F: S^m -> T^n  (S=source, T=target)
 *  and are implemented via GenericF(n,m).
 *  Don't get confused by the flipped order of source and target.
 */

namespace olb {

template<typename T, typename S, bool> class SmoothIndicatorSphere3D;
template<typename T,  typename DESCRIPTOR> class RadiativeUnitConverter;

////////////////////////////////////////////////////////////////////////////////
////////implementation of several 1d,2d,3d functors (analyticalFXD)/////////////
////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////1D////////////////////////////////////////////

/// AnalyticalConst1D: 1D -> XD, where XD is defined by value.size()
template <typename T, typename S>
class AnalyticalConst1D : public AnalyticalF1D<T,S> {
private:
  // is constant return value of operator()
  std::vector<T> _c;
public:
  AnalyticalConst1D(T value);
  AnalyticalConst1D(const std::vector<T>& value);
  bool operator() (T output[], const S x[]) override;
};

/// AnalyticalLinear1D: 1D -> 1D troughout given points (x0,v0) and (x1,v1)
//  Punktsteigungsform
template <typename T, typename S>
class AnalyticalLinear1D : public AnalyticalF1D<T,S> {
private:
  T _a;
  T _b;
public:
  AnalyticalLinear1D(T a, T b);
  AnalyticalLinear1D(S x0, T v0, S x1, T v1);
  bool operator() (T output[], const S x[]) override; ///< returns line _a*x + _b
};

/// AnalyticalRandom1D: 1D -> 1D with random image in (0,1)
template <typename T, typename S>
class AnalyticalRandom1D : public AnalyticalF1D<T,S> {
public:
  AnalyticalRandom1D();
  bool operator() (T output[], const S x[]) override;
};


/// represents an inverse parabola profile like it is used in Poiseuille inflow
/// note: output depends only on first parameter, maps 1D,2D,3D->1D
template <typename T, typename S>
class AnalyticalSquare1D : public AnalyticalF1D<T,S> {
private:
  S _cp;
  S _r;
  T _maxi;
public:
  AnalyticalSquare1D(S cp, S r, T maxi);
  bool operator() (T output[], const S x[]) override;
};


/// SinusStartScale: 1D -> 1D a start curve based on sinus for a continuous transition at 0 and 1
template <typename T, typename S>
class SinusStartScale : public AnalyticalF1D<T,S> {
protected:
  S _numTimeSteps;
  T _maxValue;
  T _pi;
public:
  SinusStartScale(int numTimeSteps=1, T maxValue=1);
  bool operator() (T output[], const S x[]) override;
};


/// PolynomialStartScale: 1D -> 1D a start curve based on a polynomial fifth order for a continuous transition at 0 and 1: maxValue*(6*y^5-15*y^4+10*y^3)
template <typename T, typename S>
class PolynomialStartScale : public AnalyticalF1D<T,S> {
protected:
  S _numTimeSteps;
  T _maxValue;
public:
  PolynomialStartScale(S numTimeSteps=S(1), T maxValue=T(1));
  bool operator() (T output[], const S x[]) override;
};

/// Derivative of a given 1D functor computed with a finite difference
template <typename T>
class AnalyticalDiffFD1D : public AnalyticalF1D<T,T> {
protected:
  AnalyticalF1D<T,T>& _f;
  T _eps;
public:
  AnalyticalDiffFD1D(AnalyticalF1D<T,T>& f, T eps = 1.e-10);
  bool operator() (T output[], const T input[]) override;
};

//////////////////////////////////2D////////////////////////////////////////////

template <typename T, typename S>
class AnalyticalComposed2D final : public AnalyticalF2D<T,S> {
private:
  AnalyticalF2D<T,S>& _f0;
  AnalyticalF2D<T,S>& _f1;
public:
  AnalyticalComposed2D(AnalyticalF2D<T,S>& f0, AnalyticalF2D<T,S>& f1);
  bool operator() (T output[], const S x[]) override;
};

/// AnalyticalConst2D: 2D -> XD, where XD is defined by value.size()
template <typename T, typename S>
class AnalyticalConst2D final : public AnalyticalF2D<T,S> {
private:
  // is constant return value of operator()
  std::vector<T> _c;
public:
  AnalyticalConst2D(T value);
  AnalyticalConst2D(T value0, T value1);
  AnalyticalConst2D(T value0, T value1, T value2);
  AnalyticalConst2D(const std::vector<T>& value);
  bool operator() (T output[], const S x[]) override;
};

/// AnalyticalLinear2D: 2D -> 1D troughout given points (x0,y0,v0), (x1,y1,v1), (x2,y2,v2)
template <typename T, typename S>
class AnalyticalLinear2D final : public AnalyticalF2D<T,S> {
protected:
  T _a;
  T _b;
  T _c;
public:
  AnalyticalLinear2D(T a, T b, T c);
  AnalyticalLinear2D(S x0, S y0, T v0, S x1, S y1, T v1, S x2, S y2, T v2);
  bool operator() (T output[], const S x[]) override;
};

/// AnalyticalRandom2D: 2D -> 1D with random image in (0,1)
template <typename T, typename S>
class AnalyticalRandom2D final : public AnalyticalF2D<T,S> {
public:
  AnalyticalRandom2D();
  bool operator() (T output[], const S x[]) override;
};

/// AnalyticalRandom2D: 2D -> 1D with maxValue in the center decreasing linearly with the distrance to the center to zero at the radius and zero outside
template <typename T, typename S>
class AnalyticalParticleAdsorptionLinear2D final : public AnalyticalF2D<T,S> {
protected:
  T _center[2];
  T _radius;
  T _maxValue;
public:
  AnalyticalParticleAdsorptionLinear2D(T center[], T radius, T maxValue);
  bool operator() (T output[], const S x[]);
};

/** Computes resulting velocity of an object from translational and rotational velocity.
 * \param indicator Class defining the object (needs to be a SmoothIndicatorF2D<T,T,true>)
 * \param u translational velocity of the object - expected in lattice units
 * \param omega rotational velocity of the object - expected in lattice units
 */
template <typename T, typename S, typename DESCRIPTOR>
class ParticleU2D : public AnalyticalF2D<T,S> {
protected:
  SmoothIndicatorF2D<T,T,true>& _indicator;
  UnitConverter<T,DESCRIPTOR> const& _converter;
public:
  ParticleU2D(SmoothIndicatorF2D<T,T,true>& indicator, UnitConverter<T,DESCRIPTOR> const& converter);
  bool operator()(T output[], const S input[]);
};


//////////////////////////////////3D////////////////////////////////////////////
template <typename T, typename S>
class AnalyticalComposed3D final : public AnalyticalF3D<T,S> {
private:
  AnalyticalF3D<T,S>& _f0;
  AnalyticalF3D<T,S>& _f1;
  AnalyticalF3D<T,S>& _f2;
public:
  AnalyticalComposed3D(AnalyticalF3D<T,S>& f0, AnalyticalF3D<T,S>& f1, AnalyticalF3D<T,S>& f2);
  bool operator() (T output[], const S x[]) override;
};

/// AnalyticalConst3D: 3D -> XD, where XD is defined by value.size()
template <typename T, typename S>
class AnalyticalConst3D final : public AnalyticalF3D<T,S> {
private:
  // is constant return value of operator()
  std::vector<T> _c;
public:
  AnalyticalConst3D(T value);
  AnalyticalConst3D(T value0, T value1);
  AnalyticalConst3D(T value0, T value1, T value2);
  AnalyticalConst3D(const std::vector<T>& value);
  AnalyticalConst3D(const Vector<T,3>& value);
  bool operator() (T output[], const S x[]) override;
};

/// AnalyticalLinear3D: 3D -> 1D troughout given points (x0,y0,z0,v0), (x1,y1,z1,v1), (x2,y2,z2,v2), (x3,y3,z3,v3)
template <typename T, typename S>
class AnalyticalLinear3D final : public AnalyticalF3D<T,S> {
protected:
  T _a;
  T _b;
  T _c;
  T _d;
public:
  AnalyticalLinear3D(T a, T b, T c, T d);
  AnalyticalLinear3D(S x0, S y0, S z0, T v0, S x1, S y1, S z1, T v1, S x2, S y2,
                     S z2, T v2, S x3, S y3, S z3, T v3);
  bool operator() (T output[], const S x[]) override;
};

/// AnalyticalRandom3D: 3D -> 1D with random image in (0,1)
template <typename T, typename S>
class AnalyticalRandom3D final : public AnalyticalF3D<T,S> {
public:
  AnalyticalRandom3D();
  bool operator() (T output[], const S x[]) override;
};

/// AnalyticalScaled3D: 3D -> Image(AnalyticalF) scales AnalyticalF by _scale
template <typename T, typename S>
class AnalyticalScaled3D final : public AnalyticalF3D<T,S> {
private:
  AnalyticalF3D<T,S>& _f;
  T _scale;
public:
  AnalyticalScaled3D(AnalyticalF3D<T,S>& f, T scale);
  bool operator() (T output[], const S x[]) override;
};

/// see Mink et al. 2016 in Sec.3.1.
template <typename T, typename S, typename DESCRIPTOR>
class PLSsolution3D : public AnalyticalF3D<T,S> {
private:
  T _physSigmaEff;
  T _physDiffusionCoefficient;
public:
  PLSsolution3D(RadiativeUnitConverter<T,DESCRIPTOR> const& converter);
  bool operator()(T output[1], const S x[3]) override;
};

/// light source as a cylinder along z-axis
template <typename T, typename S, typename DESCRIPTOR>
class LightSourceCylindrical3D : public AnalyticalF3D<T,S> {
private:
  T _physSigmaEff;
  T _physDiffusionCoefficient;
  Vector<T,3> _center;
public:
  LightSourceCylindrical3D(RadiativeUnitConverter<T,DESCRIPTOR> const& converter, Vector<T,3> center = {T(0), T(0), T(0)});
  bool operator()(T output[1], const S x[3]) override;
};

/**
* \param _position of light source
* \param _orientation direction of light source (normalized)
* \param _falloff is power of the cosine
*/
template <typename T, typename S>
class Spotlight : public AnalyticalF3D<T,S> {
private:
  Vector<T,3> const _position;
  Vector<T,3> const _orientation;
  T const _falloff;
public:
  Spotlight(Vector<T,3> position, Vector<T,3> direction, T falloff);
  bool operator()(T output[1], const S x[3]) override;
};

/// 8.6.1 Gauss Hill inital values
template <typename T, typename S>
class GaussianHill2D : public AnalyticalF2D<T,S> {
private:
  T _sigma;
  Vector<T,2> _x0;
  T _c0;
public:
  GaussianHill2D(T sigma, Vector<T,2> x0, T c0);
  bool operator()(T output[1], const S x[2]) override;
};

/// 8.6.1 Gauss Hill time evolution
template <typename T, typename S>
class GaussianHillTimeEvolution2D : public AnalyticalF2D<T,S> {
private:
  T _sigma02;
  T _D;
  T _t;
  Vector<T,2> _x0;
  Vector<T,2> _u;
  T _c0;
public:
  GaussianHillTimeEvolution2D(T sigma0, T D, T t, Vector<T,2> x0, Vector<T,2> u, T c0);
  bool operator()(T output[1], const S x[2]) override;
};

/** Computes resulting velocity of an object from translational and rotational velocity.
 * \param indicator Class defining the object (needs to be a SmoothIndicatorF3D)
 * \param u translational velocity of the object - expected in lattice units
 * \param omega rotational velocity of the object - expected in lattice units
 */
template <typename T, typename S, typename DESCRIPTOR>
class ParticleU3D : public AnalyticalF3D<T,S> {
protected:
  SmoothIndicatorF3D<T, T, true>& _indicator;
  UnitConverter<T,DESCRIPTOR> const& _converter;
public:
  ParticleU3D(SmoothIndicatorF3D<T, T, true>& indicator, UnitConverter<T,DESCRIPTOR> const& converter);
  bool operator()(T output[], const S input[]);
};

} // end namespace olb
#endif
