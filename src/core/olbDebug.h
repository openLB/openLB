/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

#ifndef OLB_DEBUG_H
#define OLB_DEBUG_H

#include <cassert>
#include <iostream>
#include <ostream>

#ifdef OLB_DEBUG

#define OLB_ASSERT( COND, MESSAGE )             \
  if ( !(COND) ) {                            \
    std::cerr << (MESSAGE) << std::endl;    \
  }                                           \
  assert( COND );

#define OLB_PRECONDITION( COND )  assert( COND );
#define OLB_POSTCONDITION( COND ) assert( COND );
#define OLB_STATECHECK( A,B )     assert( (A) == (B) );

#else

#define OLB_ASSERT( COND, MESSAGE )
#define OLB_PRECONDITION( COND )
#define OLB_POSTCONDITION( COND )
#define OLB_STATECHECK( A,B )

#endif  // OLB_DEBUG

#endif  // OLB_DEBUG_H
