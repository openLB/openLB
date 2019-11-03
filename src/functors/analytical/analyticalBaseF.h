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

#ifndef ANALYTICAL_BASE_F_H
#define ANALYTICAL_BASE_F_H

#include "functors/genericF.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"

/**
 *  The functor dimensions are given by F: S^m -> T^n  (S=source, T=target)
 *  and are implemented via GenericF(n,m).
 *  Don't get confused by the flipped order of source and target.
 */

namespace olb {

////////////////////////////////////////////////////////////////////////////////
// 2nd level classes
// note: for LatticeFunctions the number indicates the SOURCE dimension,
//       target dim depends on return variable type, so std::vector<T> is used

template<typename T, typename S> class AnalyticalIdentity2D;
template<typename T, typename S> class AnalyticalIdentity3D;


/// AnalyticalF1D are applications from 1D to XD, where X is set by the constructor.
template <typename T, typename S>
class AnalyticalF1D : public GenericF<T,S> {
protected:
  // n denotes the target dimension
  AnalyticalF1D(int n);
public:
  AnalyticalF1D<T,S>& operator-(AnalyticalF1D<T,S>& rhs);
  AnalyticalF1D<T,S>& operator+(AnalyticalF1D<T,S>& rhs);
  AnalyticalF1D<T,S>& operator*(AnalyticalF1D<T,S>& rhs);
  AnalyticalF1D<T,S>& operator/(AnalyticalF1D<T,S>& rhs);
};

/// AnalyticalF2D are applications from 2D to XD, where X is set by the constructor.
template <typename T, typename S>
class AnalyticalF2D : public GenericF<T,S> {
protected:
  // n denotes the target dimension
  AnalyticalF2D(int n);
public:
  using identity_functor_type = AnalyticalIdentity2D<T,S>;

  AnalyticalF2D<T,S>& operator-(AnalyticalF2D<T,S>& rhs);
  AnalyticalF2D<T,S>& operator+(AnalyticalF2D<T,S>& rhs);
  AnalyticalF2D<T,S>& operator*(AnalyticalF2D<T,S>& rhs);
  AnalyticalF2D<T,S>& operator/(AnalyticalF2D<T,S>& rhs);
};

/// AnalyticalF3D are applications from 3D to XD, where X is set by the constructor.
template <typename T, typename S>
class AnalyticalF3D : public GenericF<T,S> {
protected:
  // n denotes the target dimension
  AnalyticalF3D(int n);
public:
  using identity_functor_type = AnalyticalIdentity3D<T,S>;

  AnalyticalF3D<T,S>& operator-(AnalyticalF3D<T,S>& rhs);
  AnalyticalF3D<T,S>& operator+(AnalyticalF3D<T,S>& rhs);
  AnalyticalF3D<T,S>& operator*(AnalyticalF3D<T,S>& rhs);
  AnalyticalF3D<T,S>& operator/(AnalyticalF3D<T,S>& rhs);
};

/// AnalyticalIdentity1D stores vectors, result of addition,multiplication, ...
template <typename T, typename S>
class AnalyticalIdentity1D final : public AnalyticalF1D<T,S> {
protected:
  AnalyticalF1D<T,S>& _f;
public:
  AnalyticalIdentity1D(AnalyticalF1D<T,S>& f);
  bool operator() (T output[], const S input[]) override;
};

/// AnalyticalIdentity2D stores vectors, result of addition,multiplication, ...
template <typename T, typename S>
class AnalyticalIdentity2D final : public AnalyticalF2D<T,S> {
protected:
  AnalyticalF2D<T,S>& _f;
public:
  AnalyticalIdentity2D(AnalyticalF2D<T,S>& f);
  bool operator() (T output[], const S input[]) override;
};

/// AnalyticalIdentity3D stores vectors, result of addition,multiplication, ...
template <typename T, typename S>
class AnalyticalIdentity3D final : public AnalyticalF3D<T,S> {
protected:
  AnalyticalF3D<T,S>& _f;
public:
  AnalyticalIdentity3D(AnalyticalF3D<T,S>& f);
  bool operator() (T output[], const S input[]) override;
};


/// Converts IndicatorF to AnalyticalF (used for Analytical operands for Identity)
template <typename T, typename S>
class AnalyticalFfromIndicatorF3D : public AnalyticalF3D<T,S> {
protected:
  IndicatorF3D<T>& _indicatorF;
public:
  AnalyticalFfromIndicatorF3D(IndicatorF3D<T>& indicatorF);
  bool operator() (T output[], const S input[]) override;
};



} // end namespace olb

#endif
