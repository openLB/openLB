/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink0
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
#include "blockBaseF3D.hh"
#include "blockCalcF3D.hh"
#include "blockLatticeIntegralF3D.hh"
#include "blockMin3D.hh"
#include "blockMax3D.hh"
#include "blockAverage3D.hh"
#include "blockGeometryFaces3D.hh"
#include "blockLatticeLocalF3D.hh"
#include "blockLocalAverage3D.hh"
#include "latticeFrameChangeF3D.hh"
#include "reductionF3D.hh"
#include "blockReduction3D2D.hh"
#include "superBaseF3D.hh"
#include "superCalcF3D.hh"
#include "superConst3D.hh"
#include "superLatticeIntegralF3D.hh"
#include "superLatticeLocalF3D.hh"
#include "superLocalAverage3D.hh"
#include "superMin3D.hh"
#include "superMax3D.hh"
#include "superAverage3D.hh"
#include "superGeometryFaces3D.hh"
#include "turbulentF3D.hh"
#include "superErrorNorm3D.hh"
#include "indicator/indicator3D.hh"
#include "integral/integral3D.hh"
#include "timeAveraged/superLatticeTimeAveraged3D.hh"
