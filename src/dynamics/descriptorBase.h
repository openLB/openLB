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

#ifndef DESCRIPTOR_BASE_H
#define DESCRIPTOR_BASE_H

#include "descriptorTag.h"
#include "descriptorField.h"

namespace olb {

namespace descriptors {

/// \defgroup descriptor
//@{

/// Base descriptor of a D-dimensional lattice with Q directions and a list of additional fields
template <unsigned D, unsigned Q, typename... FIELDS>
struct DESCRIPTOR_BASE {

  /// Deleted constructor to enforce pure usage as type and prevent implicit narrowing conversions
  DESCRIPTOR_BASE() = delete;

  /// Tag that describes the category of the descriptor
  /**
   * e.g. tag::RTLBM for RTLBM descriptors or tag::DEFAULT if no category tag-field is found.
   *
   * This is needed to enable tag-dispatching of descriptor function overloads.
   **/
  using category_tag = tag::field_with_base<tag::CATEGORY, tag::DEFAULT, FIELDS...>;

  /// Number of dimensions
  static constexpr int d = D;
  /// Number of velocities
  static constexpr int q = Q;

  /// Returns index of WANTED_FIELD
  /**
   * Fails compilation if WANTED_FIELD is not contained in FIELDS.
   * Branching that depends on this information can be realized using `provides`.
   **/
  template <typename WANTED_FIELD>
  static constexpr int index(const unsigned localIndex=0)
  {
    return getIndexFromFieldList<D,Q,WANTED_FIELD,FIELDS...>()
           + WANTED_FIELD::template getLocalIndex<D,Q>(localIndex);
  }

  /// Returns size of FIELD
  /**
   * Works for all fields, even if they are not provided by this instantiation of DESCRIPTOR_BASE.
   **/
  template <
    typename FIELD = void,
    utilities::meta::enable_if_t<!std::is_same<FIELD,void>::value, int> = 0
  >
  static constexpr int size()
  {
    return FIELD::template size<D,Q>();
  }

  /// Returns total size of descriptor
  /**
   * This implements DESCRIPTOR_BASE::size()
   **/
  template <
    typename FIELD = void,
    utilities::meta::enable_if_t<std::is_same<FIELD,void>::value, int> = 0
  >
  static constexpr int size()
  {
    return getFieldListSize<D,Q,FIELDS...>();
  }

  /// Returns whether WANTED_FIELD is contained in FIELDS
  template <typename WANTED_FIELD>
  static constexpr bool provides()
  {
    return utilities::meta::list_contains_item<WANTED_FIELD,FIELDS...>::type::value;
  }

};

//@}

}

}

#endif
