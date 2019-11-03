/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Jonas Latt
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

#include "base64.h"
#include "base64.hh"

namespace olb {

// All of the following is a workaround to the following problem: on a
// 32-bit machine where unsigned int is the same as size_t, Base64Encoder
// needs to be instantiated on one integer type only. On some 64-bit
// platforms however, size_t is not equal to unsigned int. In that case,
// Base64Encoder needs to be instantiated twice. It is however not
// possible to instantiate this class first on unsigned int and then
// on size_t, because this yields a double instantiation, and thus
// an error, where these types are the same. To avoid this problem,
// the chosen instantiation types are unsigned int and size_t where
// these types are different, and unsigned char and unsigned int where
// they are similar. A template-based if-then-else construct is used
// to distinguish the two cases.

template<bool areEqual> struct DistinctUint;

template<>
struct DistinctUint<true> {
  typedef unsigned char T1;
  typedef size_t T2;
};

template<>
struct DistinctUint<false> {
  typedef unsigned int T1;
  typedef size_t T2;
};

typedef DistinctUint<sizeof(unsigned int)==sizeof(size_t)>::T1 T1;
typedef DistinctUint<sizeof(unsigned int)==sizeof(size_t)>::T2 T2;

template class Base64Encoder<bool>;
template class Base64Encoder<float>;
template class Base64Encoder<double>;
template class Base64Encoder<T1>;
template class Base64Encoder<T2>;
template class Base64Encoder<unsigned char>;

template class Base64Decoder<bool>;
template class Base64Decoder<float>;
template class Base64Decoder<double>;
template class Base64Decoder<T1>;
template class Base64Decoder<T2>;

}
