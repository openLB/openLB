/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
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

#ifndef DESCRIPTOR_FUNCTION_H
#define DESCRIPTOR_FUNCTION_H

#include <type_traits>
#include <stdexcept>

#include "descriptorTag.h"
#include "utilities/fraction.h"

namespace olb {

template <typename, unsigned> class Vector;

namespace descriptors {

/// \defgroup descriptor_interface Descriptor functions
/// \ingroup descriptor
//@{

/// \defgroup descriptor_interface_details Descriptor data
/// \ingroup descriptor_interface
//@{

namespace data {

using utilities::Fraction;

template <unsigned D, unsigned Q>
constexpr int vicinity = {};

template <unsigned D, unsigned Q>
constexpr int c[Q][D] = {};

template <unsigned D, unsigned Q>
constexpr int opposite[Q] = {};

template <unsigned D, unsigned Q>
constexpr Fraction t[Q] = {};

template <unsigned D, unsigned Q>
constexpr Fraction cs2 = {};

template <unsigned D, unsigned Q>
constexpr Fraction lambda_e = {};

template <unsigned D, unsigned Q>
constexpr Fraction lambda_h = {};

}

template <unsigned D, unsigned Q>
constexpr int vicinity()
{
  return data::vicinity<D,Q>;
}

template <unsigned D, unsigned Q>
constexpr int c(unsigned iPop, unsigned iDim)
{
  return data::c<D,Q>[iPop][iDim];
}

template <unsigned D, unsigned Q>
constexpr Vector<int,D> c(unsigned iPop)
{
  return Vector<int,D>(data::c<D,Q>[iPop]);
}

template <unsigned D, unsigned Q>
constexpr int opposite(unsigned iPop)
{
  return data::opposite<D,Q>[iPop];
}

template <typename T, unsigned D, unsigned Q>
constexpr T t(unsigned iPop, tag::DEFAULT)
{
  return data::t<D,Q>[iPop].template as<T>();
}

template <typename T, unsigned D, unsigned Q>
constexpr T invCs2()
{
  return data::cs2<D,Q>.template inverseAs<T>();
}

template <typename T, unsigned D, unsigned Q>
constexpr T lambda_e()
{
  return data::lambda_e<D,Q>.template inverseAs<T>();
}

template <typename T, unsigned D, unsigned Q>
constexpr T lambda_h()
{
  return data::lambda_h<D,Q>.template inverseAs<T>();
}

//@}

template <typename DESCRIPTOR>
constexpr int d()
{
  return DESCRIPTOR::d;
}

template <typename DESCRIPTOR>
constexpr int q()
{
  return DESCRIPTOR::q;
}

template <typename DESCRIPTOR>
constexpr int vicinity()
{
  return vicinity<DESCRIPTOR::d, DESCRIPTOR::q>();
}

template <typename DESCRIPTOR>
constexpr int c(unsigned iPop, unsigned iDim)
{
  return c<DESCRIPTOR::d, DESCRIPTOR::q>(iPop, iDim);
}

template <typename DESCRIPTOR>
constexpr Vector<int,DESCRIPTOR::d> c(unsigned iPop)
{
  return c<DESCRIPTOR::d, DESCRIPTOR::q>(iPop);
}

template <typename DESCRIPTOR>
constexpr int opposite(unsigned iPop)
{
  return opposite<DESCRIPTOR::d, DESCRIPTOR::q>(iPop);
}

template <typename T, typename DESCRIPTOR>
constexpr T t(unsigned iPop)
{
  return t<T, DESCRIPTOR::d, DESCRIPTOR::q>(iPop, typename DESCRIPTOR::category_tag());
}

template <typename T, typename DESCRIPTOR>
constexpr T invCs2()
{
  return invCs2<T, DESCRIPTOR::d, DESCRIPTOR::q>();
}

template <typename T, typename DESCRIPTOR>
constexpr T lambda_e()
{
  return lambda_e<T, DESCRIPTOR::d, DESCRIPTOR::q>();
}

template <typename T, typename DESCRIPTOR>
constexpr T lambda_h()
{
  return lambda_h<T, DESCRIPTOR::d, DESCRIPTOR::q>();
}

//@}

}

}

#endif
