/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerländer
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

#ifndef UTILITIES_NAMED_TYPE_H
#define UTILITIES_NAMED_TYPE_H

namespace olb {

namespace utilities {

template <typename T, typename Identificator>
class NamedType {
private:
  T _value;

public:
  explicit NamedType(const T& value):
    _value(value) {}
  explicit NamedType(T&& value):
    _value(std::move(value)) {}

  operator T() const
  {
    return _value;
  }
};

}

}

#endif
