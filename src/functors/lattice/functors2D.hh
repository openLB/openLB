/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause
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

/** \file
 * Groups all include files for the directory genericFunctions.
 */

#include "blockBaseF2D.hh"
#include "blockCalcF2D.hh"
#include "blockLatticeIntegralF2D.hh"
#include "blockLatticeLocalF2D.hh"
#include "blockLatticeRefinementMetricF2D.hh"
#include "blockLocalAverage2D.hh"
#include "blockReduction2D2D.hh"
#include "blockReduction2D1D.hh"
#include "latticeFrameChangeF3D.hh"
#include "reductionF2D.hh"
#include "superBaseF2D.hh"
#include "superCalcF2D.hh"
#include "superConst2D.hh"
#include "superLatticeIntegralF2D.hh"
#include "superLatticeLocalF2D.hh"
#include "superLatticeRefinementMetricF2D.hh"
#include "superLocalAverage2D.hh"
#include "superErrorNorm2D.hh"
#include "indicator/indicator2D.hh"
#include "integral/integral2D.hh"
#include "timeAveraged/superLatticeTimeAveraged2D.hh"
