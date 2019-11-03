/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Mathias J. Krause, Benjamin Förster
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

#ifndef INDICATOR_BASE_F_2D_H
#define INDICATOR_BASE_F_2D_H

#include <vector>

#include "core/vector.h"
#include "functors/genericF.h"

namespace olb {


/** IndicatorF1D is an application from \f$ \Omega \subset R^3 \to {0,1} \f$.
  * \param _myMin   holds minimal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _myMax   holds maximal(component wise) vector of the domain \f$ \Omega \f$.
  */
template <typename S>
class IndicatorF1D : public GenericF<bool,S> {
protected:
  IndicatorF1D();
  Vector<S,1> _myMin;
  Vector<S,1> _myMax;
public:
  virtual Vector<S,1>& getMin();
  virtual Vector<S,1>& getMax();
  IndicatorF1D<S>& operator+(IndicatorF1D<S>& rhs);
  IndicatorF1D<S>& operator-(IndicatorF1D<S>& rhs);
  IndicatorF1D<S>& operator*(IndicatorF1D<S>& rhs);
};


/** IndicatorF2D is an application from \f$ \Omega \subset R^3 \to {0,1} \f$.
  * \param _myMin   holds minimal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _myMax   holds maximal(component wise) vector of the domain \f$ \Omega \f$.
  */
template <typename S>
class IndicatorF2D : public GenericF<bool,S> {
protected:
  IndicatorF2D();
  Vector<S,2> _myMin;
  Vector<S,2> _myMax;
public:
  using GenericF<bool,S>::operator();

  virtual Vector<S,2>& getMin();
  virtual Vector<S,2>& getMax();

  /// returns false or true and pos. distance if there was one found for an given origin and direction
  /**
   * (mind that the default computation is done by a numerical approximation which searches .. [TODO])
   */
  virtual bool distance(S& distance, const Vector<S,2>& origin, const Vector<S,2>& direction, int iC=-1);

  /// returns true and the normal if there was one found for an given origin and direction
  /**
   * (mind that the default computation is done by a numerical approximation which searches .. [TODO])
   */
  virtual bool normal(Vector<S,2>& normal, const Vector<S,2>& origin, const Vector<S,2>& direction, int iC=-1);

  /// Returns true if `point` is inside a cube with corners `_myMin` and `_myMax`
  bool isInsideBox(Vector<S,2> point);

  /// Indicator specific function operator overload
  /**
   * \return Domain indicator i.e. `true` iff the input lies within the described domain.
   **/
  virtual bool operator() (const S input[]);
};


/////////////////////////////////////IdentityF//////////////////////////////////
template <typename S>
class IndicatorIdentity2D : public IndicatorF2D<S> {
public:
  std::shared_ptr<IndicatorF2D<S>> _f;

  IndicatorIdentity2D(std::shared_ptr<IndicatorF2D<S>> f);
  bool operator() (bool output[1], const S input[2]) override;
};


}

#endif
