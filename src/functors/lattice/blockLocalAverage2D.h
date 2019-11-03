/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#ifndef BLOCK_LOCAL_AVERAGE_2D_H
#define BLOCK_LOCAL_AVERAGE_2D_H

#include "blockBaseF2D.h"
#include "indicator/blockIndicatorBaseF2D.h"

namespace olb {


/// Averages given functor inside the local sphere
/**
 * Note that the average is calculated over the intersection between local sphere
 * and block cuboid. i.e. results will differ compared to  SuperLocalAverage2D if
 * the local sphere intersects the neighboring blocks.
 **/
template <typename T, typename W = T>
class BlockLocalAverage2D final : public BlockF2D<W> {
private:
  BlockF2D<W>&          _f;
  BlockIndicatorF2D<T>& _indicatorF;
  const T _radius;
public:
  /// Primary constructor
  /**
   * \param f          Functor to be locally averaged
   * \param indicatorF Indicator describing relevant cells
   * \param radius     Radius of the locality sphere
   **/
  BlockLocalAverage2D(BlockF2D<W>&          f,
                      BlockIndicatorF2D<T>& indicatorF,
                      T radius);

  /**
   * Returns average of functor \e _f evaluated on all cells both inside a sphere
   * of \e _radius around \e input and indicated by \e _indicatorF.
   **/
  bool operator() (W output[], const int input[]) override;
};


}

#endif
