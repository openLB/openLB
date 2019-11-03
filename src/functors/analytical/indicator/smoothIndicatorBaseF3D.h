/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Albert Mink, Mathias J. Krause, Benjamin Förster
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

#ifndef SMOOTH_INDICATOR_BASE_F_3D_H
#define SMOOTH_INDICATOR_BASE_F_3D_H

#include <vector>

#include "core/vector.h"
#include "functors/analytical/analyticalBaseF.h"
#include "functors/genericF.h"

namespace olb {


template <typename T, typename S, bool HLBM=false>
class SmoothIndicatorF3D;

/** SmoothIndicatorF3D is an application from \f$ \Omega \subset R^3 \to [0,1] \f$.
  * \param _myMin   holds minimal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _myMax   holds maximal(component wise) vector of the domain \f$ \Omega \f$.
  * \param _center
  * \param _diam_
  */
template <typename T, typename S>
class SmoothIndicatorF3D<T,S,false> : public AnalyticalF3D<T,S> {
protected:
  SmoothIndicatorF3D();
  //~SmoothIndicatorF3D() override {};
  Vector<S,3> _myMin;
  Vector<S,3> _myMax;
  Vector<S,3> _pos;
  Vector<S,9> _rotMat;  //saved values of rotation matrix
  S _circumRadius;
  Vector<S,3>  _theta;
  S _epsilon;
  std::string _name = "smoothIndicator3D";
public:
  void init(Vector<S,3> theta, Vector<S,3> vel, T mass, Vector<S,3> mofi);
  const Vector<S,3>& getMin() const;
  const Vector<S,3>& getMax() const;
  const Vector<S,3>& getPos() const;
  const Vector<S,9>& getRotationMatrix() const;
  const Vector<S,3>& getTheta() const;
  const S& getCircumRadius() const;
  const S& getEpsilon() const;
  std::string name();
  void setRotationMatrix(Vector<S,9> rotMat);
  void setTheta(Vector<S,3> theta);  
  void setEpsilon(S epsilon);

  SmoothIndicatorF3D<T,S,false>& operator+(SmoothIndicatorF3D<T,S,false>& rhs);
};


template <typename T, typename S>
class SmoothIndicatorIdentity3D : public SmoothIndicatorF3D<T, S, false> {
protected:
  SmoothIndicatorF3D<T, S, false>& _f;
public:
  SmoothIndicatorIdentity3D(SmoothIndicatorF3D<T, S, false>& f);
  bool operator() (T output[], const S input[]) override;
};

/** SmoothIndicatorF3D is an application from \f$ \Omega \subset R^3 \to [0,1] \f$.
  * Base class for specific SmoothIndicator implementation providing common data..
  */
template <typename T, typename S>
class SmoothIndicatorF3D<T, S, true> : public AnalyticalF3D<T,S> {
protected:
  SmoothIndicatorF3D();
  Vector<S, 3> _myMin;
  Vector<S, 3> _myMax;
  Vector<S,3> _pos;
  Vector<S,3> _vel;
  Vector<S,3> _acc;
  Vector<S,3> _acc2;
  Vector<S,3> _force;
  Vector<S,9> _rotMat;  //saved values of rotation matrix
  S _circumRadius;
  Vector<S,3> _theta;
  Vector<S,3> _omega;
  Vector<S,3> _alpha;
  Vector<S,3> _alpha2;
  S _mass;
  Vector<S,3> _mofi;  //moment of inertia
  S _epsilon;
  std::string _name="HLBMobject3D";

public:
  void init(Vector<S,3> theta, Vector<S,3> vel, T mass, Vector<S,3> mofi);
  const Vector<S,3>& getMin() const;    
  const Vector<S,3>& getMax() const;    
  const Vector<S,3>& getPos() const;
  const Vector<S,3>& getVel() const;
  const Vector<S,3>& getAcc() const;
  const Vector<S,3>& getAcc2() const;
  const Vector<S,3>& getForce() const;
  const Vector<S,9>& getRotationMatrix() const;
  const Vector<S, 3>& getTheta() const;
  const Vector<S, 3>& getOmega() const;
  const Vector<S, 3>& getAlpha() const;
  const Vector<S, 3>& getAlpha2() const;
  const Vector<S,3>& getMofi() const;
  const S& getCircumRadius() const;
  const S& getMass() const;
  const S& getEpsilon() const;
  std::string name();
  void setPos(Vector<S,3> pos);
  void setVel(Vector<S,3> vel);
  void setAcc(Vector<S,3> acc);
  void setAcc2(Vector<S,3> acc2);
  void setForce(Vector<S,3> force);
  void setRotationMatrix(Vector<S,9> rotMat);
  void setTheta(Vector<S, 3> theta);
  void setOmega(Vector<S, 3> omega);
  void setAlpha(Vector<S, 3> alpha);
  void setAlpha2(Vector<S,3> alpha2);
  void setMofi(Vector<S, 3> mofi);
  void setMass(S mass);
  void setEpsilon(S epsilon);
};

}

#endif
