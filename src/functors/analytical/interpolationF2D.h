/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2018 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Adrian Kummerlaender
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

#ifndef INTERPOLATION_F_2D_H
#define INTERPOLATION_F_2D_H

#include "analyticalF.h"
#include "functors/lattice/blockBaseF2D.h"
#include "functors/lattice/superBaseF2D.h"
#include "geometry/cuboidGeometry2D.h"
#include "geometry/blockGeometry2D.h"
#include "geometry/superGeometry2D.h"

namespace olb {


/// Converts block functors to analytical functors
template <typename T, typename W = T>
class AnalyticalFfromBlockF2D final : public AnalyticalF2D<T,W> {
protected:
  BlockF2D<W>& _f;
  Cuboid2D<T>& _cuboid;
  const int    _overlap;
public:
  AnalyticalFfromBlockF2D(BlockF2D<W>& f, Cuboid2D<T>& cuboid, const int overlap);
  bool operator() (W output[], const T physC[]) override;
};

/// Converts super functions to analytical functions
template <typename T, typename W = T>
class AnalyticalFfromSuperF2D final : public AnalyticalF2D<T,W> {
protected:
  const bool _communicateToAll;
  const bool _communicateOverlap;

  SuperF2D<T>&         _f;
  CuboidGeometry2D<T>& _cuboidGeometry;
  int                  _overlap;

  std::vector<std::unique_ptr<AnalyticalFfromBlockF2D<T,W>>> _blockF;
public:
  AnalyticalFfromSuperF2D(SuperF2D<T>& f,
                          bool communicateToAll=false,
                          int overlap=-1,
                          bool communicateOverlap=true);
  bool operator() (T output[], const T physC[]) override;

  /// \return Size of _blockF vector
  int getBlockFSize() const;
  /// \return _blockF[iCloc]
  AnalyticalFfromBlockF2D<T,W>& getBlockF(int iCloc);
};


} // end namespace olb

#endif
