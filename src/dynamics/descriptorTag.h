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

#ifndef DESCRIPTOR_TAG_H
#define DESCRIPTOR_TAG_H

#include "utilities/meta.h"

namespace olb {

namespace descriptors {

/// \defgroup descriptor
//@{

/// Base of a descriptor tag
struct DESCRIPTOR_TAG {
  template <unsigned, unsigned>
  static constexpr unsigned size()
  {
    return 0; // a tag doesn't have a size
  }
};

namespace tag {

/// Base of all tags describing the _category_ of a descriptor
/**
 * e.g. RTLBM, MRT
 **/
struct CATEGORY { };

/// Implicit default category of _normal_ descriptors
struct DEFAULT : public CATEGORY, public DESCRIPTOR_TAG { };

/// Returns first item of FIELDS type list that is derived from BASE.
/**
 * If such a type list item doesn't exist, FALLBACK is _returned_.
 **/
template <typename BASE, typename FALLBACK, typename... FIELDS>
using field_with_base = typename std::conditional<
  std::is_void<typename utilities::meta::list_item_with_base<BASE, FIELDS...>::type>::value,
  FALLBACK,
  typename utilities::meta::list_item_with_base<BASE, FIELDS...>::type
>::type;

}

//@}

}

}

#endif
