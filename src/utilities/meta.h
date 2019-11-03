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

#ifndef UTILITIES_META_H
#define UTILITIES_META_H

#include <type_traits>

namespace olb {

namespace utilities {

namespace meta {

/// Provide implementation of std::enable_if_t
/**
 * Weirdly version 18.0.3 of the Intel C++ compiler doesn't provide this
 * useful template alias despite otherwise supporting C++14.
 **/
template<bool B, class T = void>
using enable_if_t = typename std::enable_if<B,T>::type;

/// Check whether a given type list contains WANTED
template <
  typename WANTED,
  typename HEAD = void, // Default argument in case the list is empty
  typename... TAIL
>
struct list_contains_item {
  using type = typename std::conditional<
    std::is_same<WANTED, HEAD>::value,
    std::true_type,
    typename list_contains_item<WANTED, TAIL...>::type
  >::type;
};

template <typename WANTED, typename HEAD>
struct list_contains_item<WANTED, HEAD> {
  using type = typename std::is_same<WANTED, HEAD>::type;
};


/// Get first type based on BASE contained in a given type list
/**
 * If no such list item exists, type is void.
 **/
template <
  typename BASE,
  typename HEAD = void, // Default argument in case the list is empty
  typename... TAIL
>
struct list_item_with_base {
  using type = typename std::conditional<
    std::is_base_of<BASE, HEAD>::value,
    HEAD,
    typename list_item_with_base<BASE, TAIL...>::type
  >::type;
};

template <typename BASE, typename HEAD>
struct list_item_with_base<BASE, HEAD> {
  using type = typename std::conditional<
    std::is_base_of<BASE, HEAD>::value,
    HEAD,
    void
  >::type;
};

}

}

}

#endif
