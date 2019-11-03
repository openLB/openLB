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

#ifndef SMOOTH_INDICATOR_BASE_F_2D_H
#define SMOOTH_INDICATOR_BASE_F_2D_H

#include <vector>

#include "core/vector.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalBaseF.h"

namespace olb {

template <typename T, typename S, bool HLBM=false>
class SmoothIndicatorF2D;

/** SmoothIndicatorF2D is an application from \f$ \Omega \subset R^3 \to [0,1] \f$.
  * \param _myMin   holds minimal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _myMax   holds maximal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _center
  * \param _diam
  */
template <typename T, typename S>
class SmoothIndicatorF2D<T,S,false> : public AnalyticalF2D<T,S> {
protected:
  SmoothIndicatorF2D();
  Vector<S,2> _myMin;
  Vector<S,2> _myMax;
  Vector<S,2> _pos;
  Vector<S,4> _rotMat;  //saved values of rotation matrix
  S _circumRadius;
  S _theta;
  S _epsilon;
  std::string _name = "smoothIndicator2D";
public:
  void init(T theta, Vector<S,2> vel, T mass, T mofi);
  const Vector<S,2>& getMin() const;
  const Vector<S,2>& getMax() const;
  const Vector<S,2>& getPos() const;
  const Vector<S,4>& getRotationMatrix() const;
  const S& getCircumRadius() const;
  const S& getTheta() const;
  const S& getEpsilon() const;
  std::string name();
  void setRotationMatrix(Vector<S,4> rotMat);
  void setTheta(S theta);
  void setEpsilon(S epsilon);

  SmoothIndicatorF2D<T,S,false>& operator+(SmoothIndicatorF2D<T,S,false>& rhs);
};


template <typename T, typename S>
class SmoothIndicatorIdentity2D : public SmoothIndicatorF2D<T,S,false> {
protected:
  SmoothIndicatorF2D<T,S,false>& _f;
public:
  SmoothIndicatorIdentity2D(SmoothIndicatorF2D<T,S,false>& f);
  bool operator() (T output[], const S input[]) override;
};

/** SmoothIndicatorF2D is an application from \f$ \Omega \subset R^3 \to [0,1] \f$.
  * Base class for specific SmoothIndicator implementation providing common data.
  */
template <typename T, typename S>
class SmoothIndicatorF2D<T,S,true> : public AnalyticalF2D<T,S> {
protected:
  SmoothIndicatorF2D();
  Vector<S,2> _myMin;
  Vector<S,2> _myMax;
  Vector<S,2> _pos;
  Vector<S,2> _vel;
  Vector<S,2> _acc;
  Vector<S,2> _acc2;
  Vector<S,2> _force;
  Vector<S,4> _rotMat;  //saved values of rotation matrix
  S _circumRadius;
  S _theta;
  S _omega;
  S _alpha;
  S _alpha2;
  S _mass;
  S _mofi; //Moment of Inertia
  S _epsilon;
  std::string _name = "HLBMobject2D";

public:
  void init(T theta, Vector<S,2> vel, T mass, T mofi);
  const Vector<S,2>& getMin() const;
  const Vector<S,2>& getMax() const;
  const Vector<S,2>& getPos() const;
  const Vector<S,2>& getVel() const;
  const Vector<S,2>& getAcc() const;
  const Vector<S,2>& getAcc2() const;
  const Vector<S,2>& getForce() const;
  const Vector<S,4>& getRotationMatrix() const;
  const S& getCircumRadius() const;
  const S& getTheta() const;
  const S& getOmega() const;
  const S& getAlpha() const;
  const S& getAlpha2() const;
  const S& getMass() const;
  const S& getMofi() const;
  const S& getEpsilon() const;
  std::string name();
  void setPos(Vector<S,2> pos);
  void setVel(Vector<S,2> vel);
  void setAcc(Vector<S,2> acc);
  void setAcc2(Vector<S,2> acc2);
  void setForce(Vector<S,2> force);
  void setRotationMatrix(Vector<S,4> rotMat);
  void setTheta(S theta);
  void setOmega(S omega);
  void setAlpha(S alpha);
  void setAlpha2(S alpha2);
  void setMass(S mass);
  void setMofi(S mofi);
  void setEpsilon(S epsilon);
};

}

#endif
