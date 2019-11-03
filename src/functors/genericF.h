/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2013 Lukas Baron, Mathias J. Krause, Albert Mink
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

/** \file genericF.h
 * The description of a generic interface for all functor classes
 * -- header file.
 */


#ifndef GENERIC_F_H
#define GENERIC_F_H

#include<string>
#include<memory>

namespace olb {

/**
 *  GenericF is a base class, that can represent continuous as well as discrete
 *  functions.
 *  Please take care about the source and target dimensions in the constructor.
 *                F: S^m -> T^n (S=source, T=target)
 *
 *  \param _m     source dimension
 *  \param _n     target dimension
 *  \param _name  is functor name e.g. velocity, pressure
 *  \param _ptrCalcC  stores arithmeticClasses
 */
template <typename T, typename S>
class GenericF {
protected:
  // constructor
  GenericF(int targetDim, int sourceDim);

private:
  std::string _name;
  int _n;
  int _m;
  // private copy constructor
  GenericF(const GenericF&) = delete;
  // copy assignment operator
  GenericF& operator=(const GenericF&) = delete;

public:
  // virtual destructor
  virtual ~GenericF();
  /// memory management, frees resouces (calcClass)
  std::shared_ptr< GenericF<T,S> > _ptrCalcC;
  /// read only access to member variable _m
  int getSourceDim() const;
  /// read only access to member variable _n
  int getTargetDim() const;
  /// read and write access to name
  std::string& getName();
  /// read only access to name
  std::string const& getName() const;
  /// has to be implemented for 'every' derived class
  virtual bool operator() (T output[], const S input[])=0;
  /// wrapper that call the pure virtual operator() (T output[], const S input[]) from above
  // it is aimed that it even calls the implemented pure virtual operator() of derived classes
  // how to use: derived class handles the overloaded operators by including:
  // using GenericF<T,S>::operator();
  bool operator() (T output[]);
  bool operator() (T output[], S input0);
  bool operator() (T output[], S input0, S input1);
  bool operator() (T output[], S input0, S input1, S input2);
  bool operator() (T output[], S input0, S input1, S input2, S input3);
};

} // end namespace olb

#endif
