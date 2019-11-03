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

#ifndef GENERIC_F_HH
#define GENERIC_F_HH


/** \file
 * The description of a generic interface for all functor classes
 * -- generic implementation.
 */

#include"genericF.h"


namespace olb {

template <typename T, typename S>
GenericF<T,S>::GenericF (int targetDim, int sourceDim) : _n(targetDim), _m(sourceDim)
{
}

template <typename T, typename S>
GenericF<T,S>::~GenericF () {}

template <typename T, typename S>
int GenericF<T,S>::getSourceDim() const
{
  return _m;
}

template <typename T, typename S>
int GenericF<T,S>::getTargetDim() const
{
  return _n;
}

template <typename T, typename S>
std::string& GenericF<T,S>::getName()
{
  return _name;
}

template <typename T, typename S>
std::string const& GenericF<T,S>::getName() const
{
  return _name;
}

template <typename T, typename S>
bool GenericF<T,S>::operator() (T output[])
{
  S tmp[0];
  return operator()(output,tmp);
}

template <typename T, typename S>
bool GenericF<T,S>::operator() (T output[], S input0)
{
  S tmp[1] = {input0};
  return operator()(output,tmp);
}

template <typename T, typename S>
bool GenericF<T,S>::operator() (T output[], S input0, S input1)
{
  S tmp[2] = {input0,input1};
  return operator()(output,tmp);
}

template <typename T, typename S>
bool GenericF<T,S>::operator() (T output[], S input0, S input1, S input2)
{
  S tmp[3] = {input0,input1,input2};
  return operator()(output,tmp);
}

template <typename T, typename S>
bool GenericF<T,S>::operator() (T output[], S input0, S input1, S input2, S input3)
{
  S tmp[4] = {input0,input1,input2,input3};
  return operator()(output,tmp);
}

} // end namespace olb

#endif
